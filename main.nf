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
params.min_alnlen = 10000
params.min_contig = 2 * params.min_alnlen
params.min_cluster_size = 3
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


/*
 * Add the genome name to the beginning of each sequence.
 */
process preprocessGenomes {

    label "python3"
    label "small_task"
    time "1h"

    publishDir "${params.outdir}/processed"

    tag "${name}"

    input:
    set val(name), file("in.fasta") from genomes

    output:
    set val(name), file("${name}.fasta") into unscaffoldedGenomes

    script:
    """
    awk -v name="${name}" '
        /^>/ { print ">" name "." substr(\$1, 2) }
        \$0 !~ />/ { print }
    ' < "in.fasta" \
    > "out.fasta"

    split_at_n_stretch.py \
      --nsize "${params.max_ns}" \
      --min-length "${params.min_contig}" \
      -o "${name}.fasta" \
      in.fasta

    rm out.fasta
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
    set val(name), file("in.fasta") from completeGenomes

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

    label "ppg"
    label "small_task"
    time "3h"

    input:
    file "in/*.fasta" from genomes4CombineGenomes
        .mix(completeGenomes4CombineGenomes)
        .map { n, f -> f }
        .collect()

    output:
    file "genomes.fasta" into softMasked
    file "soft_masked.bed" into softMaskedBed
    file "hard_masked.fasta" into hardMasked

    script:
    """
    cat in/*.fasta > soft_masked.fasta

    ppg unsoftmask \
      --outbed soft_masked.bed \
      --outfasta genomes.fasta \
      soft_masked.fasta

    sed -E '/^>/!s/[atgcryswkmbdhvn]/N/g' soft_masked.fasta > hard_masked.fasta
    """
}


softMasked.into {
    softMasked4AlignScaffolds;
    softMasked4FindClusters;
}


// Minimap2
// url: https://github.com/lh3/minimap2
// doi: 10.1093/bioinformatics/bty191
//
// Does pairwise alignments of genomes.
process alignScaffolds {

    label "minimap2"
    label "big_task"
    time "12h"

    publishDir "${params.outdir}/alignment1"

    input:
    file "in.fasta" from softMasked4AlignScaffolds

    output:
    file "aligned.paf" into alignedScaffolds

    script:
    """
    minimap2 \
      -cx asm10 \
      -DP \
      --dual=no \
      -t "${task.cpus}" \
      in.fasta in.fasta \
    > aligned.paf
    """
}


// Filter out alignments that are only in repeat regions or are shorter than the min length.
process filterAlignedScaffolds {

    label "ppg"
    label "small_task"
    time "6h"

    publishDir "${params.outdir}/alignment1"

    input:
    file "aligned.paf" from alignedScaffolds
    file "repeats.bed" from softMaskedBed

    output:
    file "filtered.paf" into filteredAlignedScaffolds

    script:
    """
    ppg repeats --sep "." aligned.paf > duplications.bed
    cat repeats.bed duplications.bed > combined_repeats.bed

    ppg filter \
      --min-length "${params.min_alnlen}" \
      --prop-overlap 0.8 \
      combined_repeats.bed \
      aligned.paf \
    > filtered.paf
    """
}


// Cluster the alignments to find components.
process findClusters {

    label "ppg"
    label "small_task"
    time "12h"

    publishDir "${params.outdir}/alignment1"

    input:
    file "in.paf" from filteredAlignedScaffolds
    file "in.fasta" from softMasked4FindClusters

    output:
    file "clusters/*.fasta" into clusters mode flatten
    file "unplaced.fasta" into unplacedScaffolds optional true
    file "clusters.png"
    file "clusters.tsv"

    script:
    """
    ppg cluster \
      --inflation 1.4 \
      --expansion 4 \
      --plot clusters.png \
      --outfile clusters.tsv \
      in.paf

    ppg selectseqs \
      --prefix "component" \
      --outdir clusters \
      --unplaced "unplaced" \
      --min-size "${params.min_cluster_size}" \
      clusters.tsv \
      in.fasta

    if [ -e clusters/unplaced.fasta ]
    then
      mv clusters/unplaced.fasta unplaced.fasta
    fi
    """
}


clusters
    .map { f -> [f.simpleName, f] }
    .into {
        clusters4AlignComponents;
        clusters4SquishComponentAlignments;
    }


// Minimap2
// url: https://github.com/lh3/minimap2
// doi: 10.1093/bioinformatics/bty191
//
// Does pairwise alignments of genomes within components.
process alignComponents {

    label "minimap2"
    label "big_task"
    time "12h"

    publishDir "${params.outdir}/alignment2"

    tag "${component}"

    input:
    set val(component),
        file("in.fasta") from clusters4AlignComponents

    output:
    set val(component),
        file("${component}.paf") into alignedComponents

    script:
    """
    minimap2 \
      -cx asm20 \
      -DP \
      --dual=no \
      -t "${task.cpus}" \
      in.fasta in.fasta \
    > "${component}.paf"
    """
}


// Filter out alignments that are only in repeat regions or are shorter than the min length.
process filterAlignedComponents {

    label "ppg"
    label "small_task"
    time "6h"

    publishDir "${params.outdir}/alignment2"

    tag "${component}"

    input:
    set val(component),
        file("aligned.paf"),
    	file("repeats.bed") from alignedComponents.combine(softMaskedBed)

    output:
    set val(component),
        file("${component}_filtered.paf") into filteredAlignedComponents

    script:
    """
    ppg repeats --sep "." aligned.paf > duplications.bed
    cat repeats.bed duplications.bed > combined_repeats.bed

    ppg filter \
      --min-length "${params.min_alnlen}" \
      --prop-overlap 1.0 \
      combined_repeats.bed \
      aligned.paf \
    > "${component}_filtered.paf"
    """
}


/*
 * seqwish
 * url: https://github.com/ekg/seqwish
 *
 * This "squishes" alignments into a graph.
 */
process squishAlignments {

    label "seqwish"
    label "big_task"
    time "1d"

    publishDir "${params.outdir}/alignment2"

    tag "${component}"

    input:
    set val(component),
        file("in.fasta"),
        file("alignments.paf") from clusters4SquishComponentAlignments
            .combine(filteredAlignedComponents, by: 0)

    output:
    set val(component),
        file("${component}.gfa") into squishedAlignments

    script:
    """
    seqwish \
      -s "in.fasta" \
      -p "alignments.paf" \
      -g pan.gfa \
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

    publishDir "${params.outdir}/alignment2"
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

    publishDir "${params.outdir}/alignment2"

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

    publishDir "${params.outdir}/alignment2"
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
