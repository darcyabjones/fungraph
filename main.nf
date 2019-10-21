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

params.minalign = 10000


if ( params.genomes ) {
    Channel
        .fromPath(params.genomes, checkIfExists: true, type: "file")
        .map { g -> [g.simpleName, g] }
        .set { genomes }
} else {
    log.error "Please provide some genomes."
    exit 1
}


genomes.set {
    genomes4PreprocessGenomes
}


/*
 * Add the genome name to the beginning of each sequence.
 */
process preprocessGenomes {

    label "posix"
    label "small_task"

    tag "${name}"

    input:
    set val(name), file(fasta) from genomes4PreprocessGenomes

    output:
    set val(name), file("${name}.fasta.gz") into preprocessedGenomes

    script:
    """
    awk -v name="${name}" '
        /^>/ { print ">" name "." substr(\$1, 2) }
        \$0 !~ />/ { print toupper(\$0) }
    ' "${fasta}" \
    | gzip \
    > "${name}.fasta.gz"
    """
}


/*
 * All genomes need to be in the same file for minimap
 * all-vs-all alignment.
 */
process joinGenomes {

    label "posix"
    label "small_task"

    input:
    file "*.fasta.gz" from preprocessedGenomes
        .map { n, f -> f }
        .collect()

    output:
    file "combined.fasta.gz" into combinedGenomes

    script:
    """
    cat *.fasta.gz > combined.fasta.gz
    """
}

combinedGenomes.into {
    combinedGenomes4Align;
    combinedGenomes4Squish;
}


/*
 * Minimap2
 * url: https://github.com/lh3/minimap2
 * doi: 10.1093/bioinformatics/bty191
 *
 * Does pairwise alignments of genomes.
 */
process alignSelf {

    label "minimap2"
    label "big_task"

    input:
    file fasta from combinedGenomes4Align

    output:
    file "aligned.paf.gz" into aligned

    script:
    """
    minimap2 -cx asm20 -X -t "${task.cpus}" "${fasta}" "${fasta}" \
    | gzip \
    > "aligned.paf.gz"
    """
}

aligned.set {aligned4Filter}


/*
 * fpa
 * url: https://github.com/natir/fpa
 *
 * Filter minimap output to only include long alignments.
 * Useful to reduce complexity of seqwish graph construction.
 */
process filterAlignment {

    label "fpa"
    label "small_task"

    input:
    file "aligned.paf.gz" from aligned4Filter

    output:
    file "filtered.paf.gz" into filtered

    script:
    """
    zcat aligned.paf.gz | fpa drop -l "${params.minalign}" | gzip > filtered.paf.gz
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

    publishDir "${params.outdir}"

    input:
    set file("genomes.fasta.gz"),
        file("alignments.paf.gz") from combinedGenomes4Squish
            .combine(filtered)
            .first()

    output:
    file "pan.gfa" into squishedAlignments
    file "pan.vgp" optional true

    script:
    """
    seqwish \
        -s "genomes.fasta.gz" \
        -p "alignments.paf.gz" \
        -b pan.graph \
        -o pan.vgp \
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

    input:
    file "pan.gfa" from squishedAlignments

    output:
    file "pan.dg" into odgiGraph

    script:
    """
    odgi build -g pan.gfa -o - -p | odgi sort -i - -o pan.dg
    """
}


/*
 * ODGI
 * url: https://github.com/vgteam/odgi
 *
 * Makes this interesting plot.
 */
process visualiseGraph {

    label "odgi"
    label "medium_task"

    publishDir "${params.outdir}"

    input:
    file "pan.dg" from odgiGraph

    output:
    file "pan.png"

    script:
    """
    odgi viz \
      --threads "${task.cpus}" \
      -i pan.dg \
      -x 4000 \
      -y 800 \
      -L 0 \
      -A=SN15 \
      --show-strand \
      -X 1 \
      -P 10 \
      -R \
      -o pan.png
    """
}


/*
 */
process gfa2VG {

    label "vg"
    label "small_task"

    publishDir "${params.outdir}"

    input:
    file "pan.gfa" from squishedAlignments

    output:
    file "pan.vg" into vg

    script:
    // Parallelisation is limited for this step.
    """
    vg view \
      --vg \
      --gfa-in pan.gfa \
      --threads "${task.cpus}" \
    > pan.vg
    """
}


/*
 * Drop the path information and "chops" nodes to be no longer than 32 bp.
 */
process simplifyVG {

    label "vg"
    label "small_task"

    publishDir "${params.outdir}"

    input:
    file "pan.vg" from vg

    output:
    file "simplified.vg" into simplifiedVG

    script:
    """
    vg mod -X 32 pan.vg | vg ids -s - > simplified.vg 
    """
}


/*
 */
process getVGXG {

    label "vg"
    label "big_task"

    publishDir "${params.outdir}"

    input:
    file "simplified.vg" from simplifiedVG

    output:
    set file("simplified.vg"),
        file("simplified.xg") into vgWithXG

    script:
    """
    mkdir tmp
    TMPDIR="\${PWD}/tmp"
    vg index \
      -x simplified.xg \
      --temp-dir ./tmp \
      --threads "${task.cpus}" \
      simplified.vg

    rm -rf -- tmp
    """
}


process pruneGraph {

    label "vg"
    label "big_task"

    publishDir "${params.outdir}"

    input:
    set file("simplified.vg"),
        file("simplified.xg") from vgWithXG

    output:
    set file("simplified.vg"),
        file("simplified.xg"),
        file("simplified.gcsa") into prunedGraph

    script:
    """
    mkdir tmp
    TMPDIR="\${PWD}/tmp"

    vg mod -M 32 simplified.vg > resimplified.vg

    vg prune \
      -u resimplified.vg \
      -m node_mapping \
      --threads "${task.cpus}" \
    > pruned.vg

    vg index \
      -g simplified.gcsa \
      -f node_mapping \
      --temp-dir ./tmp \
      --threads "${task.cpus}" \
      --progress \
      pruned.vg
    """
}
