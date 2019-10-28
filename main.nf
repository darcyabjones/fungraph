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
params.min_contig = 1000

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
        \$0 !~ />/ { print toupper(\$0) }
    ' in.fasta \
    > "out.fasta"
    """
}


preprocessedGenomes.set {
    preprocessedGenomes4UnscaffoldGenomes;
}


/*
 */
process unscaffoldGenomes {

    label "python3"
    label "small_task"
    time "1h"

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

    tag "${name}"

    input:
    set val(name), file("in.fasta") from completeGenomes4PreprocessGenomes

    output:
    set val(name), file("${name}.fasta") into preprocessedCompleteGenomes

    script:
    """
    awk -v name="${name}" '
        /^>/ { print ">" name "." substr(\$1, 2) }
        \$0 !~ />/ { print toupper(\$0) }
    ' in.fasta \
    > "${name}.fasta"
    """
}


preprocessedCompleteGenomes.into {
    completeGenomes4CombineChannel;
    completeGenomes4AssembleMinigraph;
}


unscaffoldedGenomes.into {
    genomes4CombineChannel;
    genomes4AssembleMinigraph;
}


completeGenomes4CombineChannel
    .concat(genomes4CombineChannel)
    .into {
        genomes4RealignScaffoldsToRGFAAssembly;
        genomes4RealignScaffoldsToLinearisedRGFAAssembly;
        genomes4RealignScaffoldsToGraph;
}


/*
 */
process assembleMinigraph {

    label "minigraph"
    label "big_task"
    time "1d"

    publishDir "${params.outdir}/initial"

    input:
    set file("genomes/*"),
        file("complete_genomes/*") from genomes4AssembleMinigraph
            .map { n, f -> f }
            .collect()
            .toList()
            .combine(
                completeGenomes4AssembleMinigraph
                    .map { n, f -> f }
                    .collect()
                    .toList()
            )
            .collect()

    output:
    file "pan_minigraph.gfa" into minigraphAssembly
    file "pan_minigraph.rgfa" into minigraphRGFAAssembly

    script:
    """
    COMPLETE_ORDER=\$(order_genome_by_average_length.sh complete_genomes/*)
    ORDER=\$(order_genome_by_average_length.sh genomes/*)

    minigraph \
      -x "${params.minigraph_preset}" \
      -t "${task.cpus}" \
      -o pan_minigraph.rgfa \
      \${COMPLETE_ORDER} \
      \${ORDER}

    # This returns a canonical gfa that can be put into vg and odgi.
    echo -e "H\\tVN:Z:1.0" > pan_minigraph.gfa
    awk '
      BEGIN {OFS="\\t"}
      \$0 ~ /^S/ {
        \$2=gensub(/^s/, "", "g", \$2);
        \$0=\$1 "\\t" \$2 "\\t" \$3
      }
      \$0 ~ /^L/ {
        \$2=gensub(/^s/, "", "g", \$2);
        \$4=gensub(/^s/, "", "g", \$4);
        \$0=\$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t" \$5 "\\t" \$6
      }
      { print }
    ' pan_minigraph.rgfa \
    >> pan_minigraph.gfa
    """
}


minigraphRGFAAssembly.into {
    minigraphRGFAAssembly4RealignScaffoldsToRGFAAssembly;
    minigraphRGFAAssembly4LineariseRGFAAssembly;
}


/*
 */
process realignScaffoldsToRGFAAssembly {

    label "minigraph"
    label "medium_task"
    time "5h"

    publishDir "${params.outdir}/initial/gafs"
    tag "${name}"

    input:
    set val(name),
        file("genome.fasta"),
        file("pan.rgfa") from genomes4RealignScaffoldsToRGFAAssembly
            .combine(minigraphRGFAAssembly4RealignScaffoldsToRGFAAssembly)

    output:
    set val(name),
        file("${name}.gaf") into realignedScaffoldsToRGFA

    script:
    """
    minigraph \
      -x lr \
      -N 0 \
      -t "${task.cpus}" \
      -o "${name}.gaf" \
      pan.rgfa \
      genome.fasta
    """
}


process lineariseRGFAAssembly {

    label "gfatools"
    label "small_task"

    publishDir "${params.outdir}/initial"

    input:
    file "pan_minigraph.rgfa" from minigraphRGFAAssembly4LineariseRGFAAssembly

    output:
    file "pan_minigraph.fasta"

    script:
    """
    gfatools gfa2fa -s pan_minigraph.rgfa > pan_minigraph.fasta
    """
}


/*
 */
process gfa2InitialVG {

    label "vg"
    label "medium_task"
    time "5h"

    publishDir "${params.outdir}/initial"

    input:
    file "pan.gfa" from minigraphAssembly

    output:
    file "pan.vg" into initialVg

    script:
    // Parallelisation is limited for this step.
    """
    vg view \
      --vg \
      --gfa-in pan.gfa \
      --threads "${task.cpus}" \
    > pan_first.vg

    vg mod -X 32 pan_first.vg > pan.vg
    vg ids -s -j pan.vg

    rm -f pan_first.vg
    """
}


/*
 * We are aligning each genome to the graph and updating
 * the graph before aligning the next one.
 * To do this I add a cycle.
 */
realignmentAccumulator = Channel.create()

/*
 * Compute the indices needed to align against the genome.
 */
process indexRealignedVG {

    label "vg"
    label "biggish_task"
    time "12h"

    tag "${component}"

    // This will be overwritten at each iteration
    publishDir "${params.outdir}/realigned"

    input:
    file "realigned.vg" initialVg.mix(realignmentAccumulator)

    output:
    set file("realigned.vg"),
        file("realigned.xg"),
        file("realigned.gcsa") into indexedRealignmentAccumulator

    // This is where the output will go.
    set file("realigned.vg"),
        file("realigned.xg"),
        file("realigned.gcsa") into finalRealignedGraph

    script:
    """
    mkdir -p tmp
    TMPDIR="\${PWD}/tmp"

    vg index \
      -x "output.xg" \
      --temp-dir ./tmp \
      --threads "${task.cpus}" \
      realigned.vg

    vg mod -M 32 realigned.vg > simplified.vg

    vg prune \
      -u simplified.vg \
      -m node_mapping \
      --threads "${task.cpus}" \
    > extra_simplified.vg

    vg index \
      -g realigned.gcsa \
      -f node_mapping \
      --temp-dir ./tmp \
      --threads "${task.cpus}" \
      --progress \
      extra_simplified.vg

    rm -rf -- tmp
    rm -f simplified.vg extra_simplified.vg node_mapping
    """
}


/*
 * Take the updated indexed graph and align the genome to it, keeping
 * only the best alignments.
 */
process realignScaffoldsToGraph {

    label "vg"
    label "big_task"
    time "12h"

    input:
    set val(name),
        file("genome.fasta"),
        file("input.vg"),
        file("input.xg"),
        file("input.gcsa") from genomes4RealignScaffoldsToGraph
            .merge(indexedRealignmentAccumulator)

    output:
    // This feeds back into the accumulator, to be re-indexed.
    file "output.vg" into realignmentAccumulator

    script:
    """
    vg align \
      -x in.xg \
      -g in.gcsa \
      --fasta genome.fasta \
      --threads "${task.cpus}" \
      --hit-max 500 \
      --band-multi 128 \
      --xdrop-alignment \
      --alignment-model long \
      --max-multimaps 1 \
      --debug \
    > out.gam

    vg augment \
      --include-paths \
      --cut-softclips \
      input.vg \
      out.gam \
    > augmented.vg

    vg mod -X 32 augmented.vg > output.vg

    rm -f out.gam augmented.vg
    """
}


/*
 * This is just to split connected components.
 */
process explodeRealignedGraph {

    label "vg"
    label "small_task"

    publishDir "${params.outdir}/realigned"

    input:
    set file("realigned.vg"),
        file("realigned.xg"),
        file("realigned.gcsa") from finalRealignedGraph.last()

    output:
    file "exploded"

    script:
    """
    vg explode realigned.vg exploded
    """
}


/*
process callMinimapVariants {

    label "minimap2"
    label "medium_task"
    time "4h"

    publishDir "${params.outdir}/paf_call"

    input:
    file "pan.fasta" from panFasta
    set val(name), file("genome.fasta") from genomes4CallMinimapVariants

    output:
    file "${name}_minimap2.paf"
    file "${name}_minimap2.var"

    script:
    """
    mkdir -p tmp
    minimap2 -cx asm20 -t${task.cpus} --cs pan.fasta genome.fasta \
    | sort -k6,6 -k8,8n --temporary-directory=tmp \
    > "${name}_minimap2.paf"

    paftools.js call \
      -l 1000 \
      -L 10000 \
      -f pan.fasta \
      -s "${name}" \
      "${name}_minimap2.paf" \
    > "${name}_minimap2.var"

    rm -rf -- tmp
    """
}
*/




// /*
//  * All genomes need to be in the same file for minimap
//  * all-vs-all alignment.
//  */
// process joinGenomes {
//
//     label "posix"
//     label "small_task"
//     time "1h"
//
//     input:
//     file "*.fasta.gz" from preprocessedGenomes
//     .map { n, f -> f }
//     .collect()
//
//     output:
//     file "combined.fasta.gz" into combinedGenomes
//
//     script:
//     """
//     cat *.fasta.gz > combined.fasta.gz
//     """
// }
//
// combinedGenomes.into {
//     combinedGenomes4Align;
//     combinedGenomes4Squish;
// }
//
//
// /*
//  * Minimap2
//  * url: https://github.com/lh3/minimap2
//  * doi: 10.1093/bioinformatics/bty191
//  *
//  * Does pairwise alignments of genomes.
//  */
// process alignSelf {
//
//     label "minimap2"
//     label "big_task"
//     time "1d"
//
//     input:
//     file fasta from combinedGenomes4Align
//
//     output:
//     file "aligned.paf.gz" into aligned
//
//     script:
//     """
//     minimap2 -cx asm20 -X -t "${task.cpus}" "${fasta}" "${fasta}" \
//     | gzip \
//     > "aligned.paf.gz"
//     """
// }
//
// aligned.set {aligned4Filter}
//
//
// /*
//  * fpa
//  * url: https://github.com/natir/fpa
//  *
//  * Filter minimap output to only include long alignments.
//  * Useful to reduce complexity of seqwish graph construction.
//  */
// process filterAlignment {
//
//     label "fpa"
//     label "small_task"
//     time "5h"
//
//     input:
//     file "aligned.paf.gz" from aligned4Filter
//
//     output:
//     file "filtered.paf.gz" into filtered
//
//     script:
//     """
//     zcat aligned.paf.gz | fpa drop -l "${params.minalign}" | gzip > filtered.paf.gz
//     """
// }
//
//
// /*
//  * seqwish
//  * url: https://github.com/ekg/seqwish
//  *
//  * This "squishes" alignments into a graph.
//  */
// process squishAlignments {
//
//     label "seqwish"
//     label "big_task"
//     time "1d"
//
//     publishDir "${params.outdir}"
//
//     input:
//     set file("genomes.fasta.gz"),
//     file("alignments.paf.gz") from combinedGenomes4Squish
//         .combine(filtered)
//         .first()
//
//     output:
//     file "pan.gfa" into squishedAlignments
//     file "pan.vgp" optional true
//
//     script:
//     """
//     seqwish \
//     -s "genomes.fasta.gz" \
//     -p "alignments.paf.gz" \
//     -b pan.graph \
//     -o pan.vgp \
//     -g pan.gfa \
//     --threads "${task.cpus}"
//     """
// }
//
//
// /*
//  * ODGI
//  * url: https://github.com/vgteam/odgi
//  */
// process gfa2ODGI {
//
//     label "odgi"
//     label "small_task"
//     time "5h"
//
//     publishDir "${params.outdir}"
//
//     when:
//     params.assembler == "seqwish"
//
//     input:
//     file "pan.gfa" from squishedAlignments
//
//     output:
//     file "pan.dg" into odgiGraph
//
//     script:
//     """
//     odgi build -g pan.gfa -o - -p | odgi sort -i - -o pan.dg
//     """
// }
//
//
// odgiGraph.into {
//     odgiGraph4VisualiseGraph;
//     odgiGraph4SimplifyGraph;
// }
//
// /*
//  * ODGI
//  * url: https://github.com/vgteam/odgi
//  *
//  * Makes this interesting plot.
//  */
// process visualiseGraph {
//
//     label "odgi"
//     label "big_task"
//     time "2h"
//
//     publishDir "${params.outdir}"
//
//     when:
//     params.assembler == "seqwish"
//
//     input:
//     file "pan.dg" from odgiGraph4VisualiseGraph
//
//     output:
//     file "pan.png"
//
//     script:
//     """
//     odgi viz \
//       --threads "${task.cpus}" \
//       -i pan.dg \
//       -x 4000 \
//       -y 800 \
//       -L 0 \
//       -A=SN15 \
//       --show-strand \
//       -X 1 \
//       -P 10 \
//       -R \
//       -o pan.png
//     """
// }
