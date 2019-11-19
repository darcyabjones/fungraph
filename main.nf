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
      out.fasta

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
        clusters4RealignAgainstGraph;
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
    label "biggish_task"
    time "1d"


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
      -g "${component}.gfa" \
      --threads "${task.cpus}"
    """
}


squishedAlignments.set {
    squishedAlignments4ODGI;
}


/*
 * The squished graph has ~2x the sequence content and is too complex for vg.
 * This removes simple bubbles, then bubbles shorter than 50bp.
 * I do seem to get different results with this, rather than just using -B.
 * `-B` is used to avoid tip-trimming.
 */
process popBubbles {

    label "gfatools"
    label "small_task"
    time "1h"

    tag "${component}"
    input:
    set val(component),
        file("in.gfa") from squishedAlignments

    output:
    set val(component),
        file("${component}.gfa") into poppedGraph

    script:
    """
    gfatools asm -s -l 50 -B in.gfa > "${component}.gfa"
    """
}


process simplifyGraph {

    label "vg"
    label "small_task"
    time "2h"

    tag "${component}"

    input:
    set val(component),
        file("in.gfa") from poppedGraph

    output:
    set val(component),
        file("${component}.vg") into simplifiedGraph

    script:
    """
    vg view --vg --gfa-in in.gfa > in.vg
    vg mod --unfold 1 in.vg > unfolded.vg
    vg mod --dagify-step 1 unfolded.vg > dagified.vg
    vg mod --until-normal 10 dagified.vg > "${component}.vg"

    rm in.vg unfolded.vg dagified.vg
    """
}


process combineGraph {

    label "vg"
    label "small_task"
    time "2h"

    input:
    file "components/*.vg" from simplifiedGraph.map { c, g -> g }.collect()

    output:
    file "combined.vg" into combinedGraph

    script:
    """
    # This is to avoid mutating and messing up checkpoints.
    cp -rL components new_components

    vg ids --join --compact --sort new_components/*.vg
    cat new_components/*.vg > combined.vg
    """
}


combinedGraph.into {
    combinedGraph4RealignAgainstGraph;
    combinedGraph4AugmentGraph;
}


process realignAgainstGraph {

    label "graphaligner"
    label "medium_task"
    time "12h"

    tag "${component}"

    input:
    set val(component),
        file("contigs.fasta"),
        file("combined.vg") from clusters4RealignAgainstGraph
            .mix( unplacedScaffolds.map { f -> [f.simpleName, f] } )
            .combine(combinedGraph)

    output:
    set val(component),
        file("out.gam") into realignedAgainstGraph

    script:
    """
    GraphAligner \
      --graph combined.vg \
      --reads contigs.fasta \
      --alignments-out out.gam \
      --threads "${task.cpus}" \
      --seeds-mem-count -1 \
      --seeds-mxm-length 20 \
      --tangle-effort 100000
    """
}


process augmentGraph {

    label "vg"
    label "medium_task"
    time "12h"

    publishDir "${params.outdir}/alignment2"

    input:
    file "in.vg" from combinedGraph4AugmentGraph
    file "alignments/*.gam" from realignAgainstGraph
        .map { c, g -> g }
        .collect()

    output:
    file "pan.gfa"
    file "components/*.gfa" into augmentedComponents mode flatten
    file "pan.vg"

    script:
    """
    cat alignments/*.gam > alignments.gam

    vg augment \
      --label-paths \
      --progress \
      --threads "${task.cpus}" \
      in.vg \
      alignments.gam \
    > pan.vg

    vg view --gfa --vg-in pan.vg > pan.gfa

    vg explode pan.vg augmented

    mkdir -p components
    for f in augmented/*.vg
    do
      FNAME="\$(basename \${f%.*}}"
      vg view --gfa --vg-in "\${f}" > "components/\${FNAME}.gfa"
    done
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
        file("pan.gfa") from augmentedComponents
            .map { f -> [f.simpleName, f] }

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
