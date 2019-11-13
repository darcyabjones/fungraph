#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/


def helpMessage() {
    log.info"""
    # fungraph

    ## Exit codes

    - 0: All ok.
    - 1: Incomplete parameter inputs.

    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}


params.genomes = false
params.complete_genomes = false
params.minigraph_preset = "ggs"
params.max_ns = 50
params.min_alnlen_one = 50000
params.min_alnlen_two = 5000
params.min_contig = params.min_alnlen_two
params.min_coverage = 0.9
params.min_bubble = 30


if ( !(params.genomes || params.complete_genomes) ) {
    log.error "Please provide some genomes."
    exit 1
}

if ( params.genomes ) {
    Channel
        .fromPath(params.genomes, checkIfExists: true, type: "file")
        .map { g -> [g.simpleName, g] }
        .set { genomes }
} else {
    genomes = Channel.empty()
}

if ( params.complete_genomes ) {
    Channel
        .fromPath(params.complete_genomes, checkIfExists: true, type: "file")
        .map { g -> [g.simpleName, g] }
        .set { completeGenomes }
} else {
    completeGenomes = Channel.empty()
}

completeGenomes.into {
    completeGenomes4PreprocessGenomes;
    completeGenomes4CallMinimapVariants;
}


genomes.into {
    genomes4PreprocessGenomes;
    genomes4CallMinimapVariants;
}


/*
 * Add the genome name to the beginning of each sequence.
 */
process preprocessGenomes {

    label "posix"
    label "small_task"
    time "1h"

    tag "${name}"

    input:
    set val(name), file("in.fasta") from genomes4PreprocessGenomes

    output:
    set val(name), file("out.fasta") into preprocessedGenomes

    script:
    """
    awk -v name="${name}" '
        /^>/ { print ">" name "." substr(\$1, 2) }
        \$0 !~ />/ { print }
    ' < "in.fasta" \
    > "out.fasta"
    """
}


preprocessedGenomes.set { preprocessedGenomes4UnscaffoldGenomes }


/*
 */
process unscaffoldGenomes {

    label "python3"
    label "small_task"
    time "1h"

    publishDir "${params.outdir}/processed"

    tag "${name}"

    input:
    set val(name),
        file("in.fasta") from preprocessedGenomes4UnscaffoldGenomes

    output:
    set val(name),
        file("${name}.fasta") into unscaffoldedGenomes

    script:
    """
    split_at_n_stretch.py \
      --nsize "${params.max_ns}" \
      --min-length "${params.min_contig}" \
      -o "${name}.fasta" \
      in.fasta
    """
}


/*
 * Add the genome name to the beginning of each sequence.
 */
process preprocessCompleteGenomes {

    label "posix"
    label "small_task"
    time "1h"

    publishDir "${params.outdir}/processed"

    tag "${name}"

    input:
    set val(name), file("in.fasta") from completeGenomes4PreprocessGenomes

    output:
    set val(name), file("${name}.fasta") into preprocessedCompleteGenomes

    script:
    """
    awk -v name="${name}" '
        /^>/ { print ">" name "." substr(\$1, 2) }
        \$0 !~ />/ { print }
    ' < "in.fasta" \
    > "${name}.fasta"
    """
}


unscaffoldedGenomes.set { genomes4CombineGenomes }
preprocessedCompleteGenomes.set { completeGenomes4CombineGenomes }


/*
 * All genomes need to be in the same file for minimap
 * all-vs-all alignment.
 */
process combineGenomes {

    label "python3"
    label "small_task"

    time "3h"

    input:
    file "in/*.fasta" from genomes4CombineGenomes
        .map { n, f -> f }
        .collect()

    output:
    file "soft_masked.fasta" into softMasked
    file "soft_masked.bed" into softMaskedBed
    file "hard_masked.fasta" into hardMasked

    script:
    """
    cat in/*.fasta > soft_masked.fasta
    softmask_to_bed.py soft_masked.fasta > soft_masked.bed
    sed -E '/^>/!s/[atgcryswkmbdhvn]/N/g' soft_masked.fasta > hard_masked.fasta
    """
}


softMasked.into {
    softMasked4AlignScaffoldsAllVAll;
    softMasked4RealignScaffoldsToComponents;
    softMasked4SelectScaffoldsToComponents;
    softMasked4SquishAlignmentsFine;
}


/*
 * All genomes need to be in the same file for minimap
 * all-vs-all alignment.
 */
process combineCompleteGenomes {

    label "python3"
    label "small_task"

    time "3h"

    input:
    file "in/*.fasta" from completeGenomes4CombineGenomes
        .map { n, f -> f }
        .collect()

    output:
    file "soft_masked.fasta" into softMaskedComplete
    file "soft_masked.bed" into softMaskedBedComplete

    script:
    """
    cat in/*.fasta > soft_masked.fasta
    softmask_to_bed.py soft_masked.fasta > soft_masked.bed
    sed -E '/^>/!s/[atgcryswkmbdhvn]/N/g' soft_masked.fasta > hard_masked.fasta
    """
}


softMaskedComplete.into {
    softMaskedComplete4AlignCompleteScaffoldsAllVAll;
    softMaskedComplete4SquishAlignmentsCoarse;
    softMaskedComplete4RealignScaffoldsToComponents;
    softMaskedComplete4SelectScaffoldsToComponents;
    softMaskedComplete4AlignScaffoldsAllVAll;
    softMaskedComplete4SquishAlignmentsFine;
}


// Minimap2
// url: https://github.com/lh3/minimap2
// doi: 10.1093/bioinformatics/bty191
//
// Does pairwise alignments of genomes.
process alignCompleteScaffoldsAllVAll {

    label "minimap2"
    label "big_task"
    time "12h"

    input:
    file "in.fasta" from softMasked4AlignCompleteScaffoldsAllVAll

    output:
    file "aligned.paf" into alignedCompleteScaffoldsAllVAll

    script:
    """
    minimap2 \
      -cx asm20 \
      -DP \
      --dual=no \
      -t "${task.cpus}" \
      in.fasta in.fasta \
    > aligned.paf
    """
}


// fpa
// url: https://github.com/natir/fpa
//
// Filter minimap output to only include long alignments.
// Useful to reduce complexity of seqwish graph construction.
process filterAlignedScaffoldsAllvAllCoarse {

    label "python3"
    label "small_task"

    input:
    file "aligned.paf" from alignedScaffoldsAllVAll
    file "repeats.bed" from softMaskedBed

    output:
    file "filtered.paf" into filteredAlignedScaffoldsAllVAll

    script:
    """
    filter_paf.py repeats --sep "." aligned.paf > duplications.bed
    cat repeats.bed duplications.bed > combined_repeats.bed

    filter_paf.py filter \
      --min-length "${params.min_alnlen_one}" \
      --prop-overlap 0.8 \
      combined_repeats.bed \
      aligned.paf \
    > filtered.paf

    #       --sep "."
    """
}


/*
 * seqwish
 * url: https://github.com/ekg/seqwish
 *
 * This "squishes" alignments into a graph.
 */
process squishAlignmentsCoarse {

    label "seqwish"
    label "big_task"
    time "1d"

    publishDir "${params.outdir}"

    input:
    set file("genomes.fasta"),
        file("alignments.paf") from softMasked4SquishAlignmentsCoarse
            .combine(filteredAlignedScaffoldsAllVAll)
            .first()

    output:
    file "pan.gfa" into squishedAlignments

    script:
    """
    seqwish \
      -s "genomes.fasta" \
      -p "alignments.paf" \
      -g pan.gfa \
      --threads "${task.cpus}"
    """
}


// Splits the graph up into connected components.
process gfa2InitialComponentGFAs {

    label "vg"
    label "medium_task"
    time "5h"

    publishDir "${params.outdir}/initial/component_gfas"

    input:
    file "pan.gfa" from squishedAlignments

    output:
    file "initial/*.gfa" into initialComponentGFAs mode flatten

    script:
    // Parallelisation is limited for this step.
    """
    vg view \
      --vg \
      --gfa-in pan.gfa \
      --threads "${task.cpus}" \
    > pan.vg

    vg explode \
      --threads "${task.cpus}" \
      pan.vg \
      initial

    for v in initial/*.vg
    do
      vg view --gfa --vg-in "\${v}" > "\${v%.vg}.gfa"
    done

    rm pan.vg
    rm initial/*.vg
    """
}


// Realign the contigs to the GFA to find ones to re-assemble.
process findComponentContigsAlign {

    label "minigraph"
    label "medium_task"
    time "5h"

    tag "${component}"

    publishDir "${params.outdir}/initial/component_gafs"

    input:
    set val(component),
        file("component.gfa"),
        file("complete_genomes.fasta") from initialComponentGFAs
        file("genomes.fasta") from initialComponentGFAs
            .map { g -> [g.baseName, g] }
            .combine(softMaskedComplete4RealignScaffoldsToComponents)
            .combine(softMasked4RealignScaffoldsToComponents)

    output:
    set val(component),
        file("${component}.gaf") into alignedComponentContigs

    script:
    """
    minigraph \
      -x lr \
      -N 0 \
      -t "${task.cpus}" \
      -o "${component}.gaf" \
      component.gfa \
      complete_genomes.fasta \
      genomes.fasta
    """
}


// Finds which components have the best coverage for each contig.
// We'll use this to re-assemble the contigs individually.
process selectBestContigsForComponents {

    label "python3"
    label "small_task"
    time "3h"

    publishDir "${params.outdir}/initial"

    input:
    set file("complete_genomes.fasta"),
        file("genomes.fasta")
        file("gafs/*") from softMaskedComplete4SelectScaffoldsToComponents
            .combine(softMasked4SelectScaffoldsToComponents)
            .combine(alignedComponentContigs.map { c, g -> g }.collect().toList())
            .first()

    output:
    file "fasta_components/*.fasta" into bestContigsForComponents mode flatten
    file "unplaced_contigs.fasta" into unplacedContigs optional true
    file "contigs_to_components.tsv"

    script:
    """
    select_best_contigs_for_components.py \
      --min-coverage "${params.min_coverage}" \
      --outfile "contigs_to_components.tsv" \
      gafs/*


    mkdir -p "fasta_components"

    select_sequences.py \
      --outdir "fasta_components" \
      "contigs_to_components.tsv" \
      in.fasta

    if [ -s "fasta_components/unplaced.fasta" ]
    then
      mv "fasta_components/unplaced.fasta" "unplaced_contigs.fasta"
    fi
    """
}


bestContigsForComponents
    .map { s -> [s.baseName, s] }
    .into { bestContigsForComponents4Align; bestContigsForComponents4Squish }


// Minimap2
// url: https://github.com/lh3/minimap2
// doi: 10.1093/bioinformatics/bty191
//
// Does pairwise alignments of genomes.
process alignComponentScaffoldsAllVAll {

    label "minimap2"
    label "big_task"
    time "12h"

    tag "${component}"

    input:
    set val(component), file("in.fasta") from bestContigsForComponents4Align

    output:
    set val(component),
        file("aligned.paf") into alignedComponentScaffoldsAllVAll

    script:
    """
    minimap2 \
      -cx asm20 \
      -DP \
      --dual=no \
      -t "${task.cpus}" \
      in.fasta in.fasta \
    > aligned.paf
    """
}


// fpa
// url: https://github.com/natir/fpa
//
// Filter minimap output to only include long alignments.
// Useful to reduce complexity of seqwish graph construction.
process filterAlignedScaffoldsAllvAllFine {

    label "python3"
    label "small_task"
    time "2h"

    tag "${component}"

    input:
    set val(component),
        file("aligned.paf"),
        file("complete_repeats.bed"),
        file("repeats.bed") from alignedComponentScaffoldsAllVAll
            .combine(softMaskedBedComplete)
            .combine(softMaskedBed)

    output:
    set val(component),
        file("filtered.paf") into filteredAlignedComponentScaffoldsAllVAll

    script:
    """
    filter_paf.py repeats --sep "." aligned.paf > duplications.bed
    cat complete_repeats.bed repeats.bed duplications.bed > combined_repeats.bed

    filter_paf.py filter \
      --min-length "${params.min_alnlen_two}" \
      --prop-overlap 0.8 \
      combined_repeats.bed \
      aligned.paf \
    > filtered.paf

    #       --sep "."
    """
}


/*
 * seqwish
 * url: https://github.com/ekg/seqwish
 *
 * This "squishes" alignments into a graph.
 */
process squishAlignmentsFine {

    label "seqwish"
    label "big_task"
    time "1d"

    publishDir "${params.outdir}/fine"

    tag "${component}"

    input:
    set val(component),
        file("genomes.fasta"),
        file("alignments.paf") from bestContigsForComponents4Squish
            .combine(filteredAlignedComponentScaffoldsAllVAll, by: 0)

    output:
    set val(component), file("${component}.gfa") into squishedComponentAlignments

    script:
    """
    seqwish \
      -s "genomes.fasta" \
      -p "alignments.paf" \
      -g "${component}.gfa" \
      --threads "${task.cpus}"
    """
}


/*
 * ODGI
 * url: https://github.com/vgteam/odgi
 */
process gfa2ODGI {

    label "odgi"
    label "small_task"
    time "5h"

    publishDir "${params.outdir}/fine"
    tag "${component}"

    input:
    set val(component),
        file("pan.gfa") from squishedComponentAlignments

    output:
    set val(component),
        file("${component}.dg") into odgiGraph

    script:
    """
    odgi build -g pan.gfa -o - -p | odgi sort -i - -o "${component}.dg"
    """
}


odgiGraph.into {
    odgiGraph4VisualiseGraph;
    odgiGraph4GetBins;
}

/*
 * ODGI
 * url: https://github.com/vgteam/odgi
 *
 * Makes this interesting plot.
 */
process visualiseGraph {

    label "odgi"
    label "small_task"
    time "2h"

    publishDir "${params.outdir}/fine"

    tag "${component}"

    input:
    set val(component),
        file("pan.dg") from odgiGraph4VisualiseGraph

    output:
    file "${component}.png"

    script:
    """
    odgi viz \
      --threads "${task.cpus}" \
      -i pan.dg \
      -x 4000 \
      -y 800 \
      -L 0 \
      -X 1 \
      -P 10 \
      -R \
      -o "${component}.png"
    """
}


process getBins {

    label "odgi"
    label "small_task"
    time "2h"

    publishDir "${params.outdir}/fine/bins"
    tag "${component}"

    input:
    set val(component),
        file("pan.dg") from odgiGraph4GetBins

    output:
    set val(component),
        file("${component}.json")

    script:
    """
    odgi bin \
        -i pan.dg \
        --json \
        --bin-width 10000 \
        --path-delim "." \
    > "${component}.json"
    """
}
