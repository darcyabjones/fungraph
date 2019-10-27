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
params.minigraph_preset = "ggs"
params.max_ns = 50
params.min_contig = 1000

if ( params.genomes ) {
    Channel
        .fromPath(params.genomes, checkIfExists: true, type: "file")
        .map { g -> [g.simpleName, g] }
        .set { genomes }
} else {
    log.error "Please provide some genomes."
    exit 1
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
        file("out.fasta") into unscaffoldedGenomes

    script:
    """
    split_at_n_stretch.py \
      --nsize "${params.max_ns}" \
      --min-length "${params.min_contig}" \
      -o out.fasta \
      in.fasta
    """
}


unscaffoldedGenomes.into {
    unscaffoldGenomes4AssembleMinigraph;
    unscaffoldedGenomes4FindComponentContigs;
    unscaffoldedGenomes4SelectBestContigsForComponents;
    unscaffoldedGenomes4RealignScaffoldsToGraph;
}


process assembleMinigraph {

    label "minigraph"
    label "big_task"
    time "1d"

    publishDir "${params.outdir}/initial"

    input:
    file "genomes/*.fasta" from unscaffoldGenomes4AssembleMinigraph
        .map { n, f -> f }
        .collect()

    output:
    file "pan_minigraph.gfa" into minigraphAssembly
    file "pan_minigraph.rgfa"

    script:
    """
    minigraph \
      -x "${params.minigraph_preset}" \
      -t "${task.cpus}" \
      -o pan_minigraph.rgfa \
      genomes/*

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


/*
 */
process gfa2InitialVG {

    label "vg"
    label "medium_task"
    time "5h"

    publishDir "${params.outdir}/initial/gfas"

    input:
    file "pan.gfa" from minigraphAssembly

    output:
    file "initial" into initialVg
    file "gfas/*" into initialGfas mode flatten

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

    mkdir -p gfas
    find initial -name "*.vg" -printf '%f\\0' \
    | xargs -0 -P "${task.cpus}" -I {} -- \
        bash -eu -c '
            F="{}";
            B="\${F%.*}";
            vg view --vg-in --gfa "initial/{}" > "gfas/\${B}.gfa"
        '
    """
}


/*
 */
process chopInitialVG {

    label "vg"
    label "big_task"
    time "5h"

    publishDir "${params.outdir}/initial/vg"

    input:
    file "in" from initialVg

    output:
    file "initial/*.vg" into choppedInitialVg mode flatten

    script:
    """
    mkdir -p initial

    find in/ -name "*.vg" -printf '%f\\0' \
    | xargs -0 -P 1 -I {} -- \
        bash -eu -c 'vg mod -X 32 in/{} > initial/{}'

    vg ids -s -j initial/*
    """
}


process findComponentContigs {

    label "minigraph"
    label "medium_task"
    time "5h"

    publishDir "${params.outdir}/initial/gafs"
    tag "${component}"

    input:
    set val(component),
        file("component.gfa"),
        file("genomes/*.fasta") from initialGfas
            .map { f -> [f.baseName, f] }
            .combine(
                unscaffoldedGenomes4FindComponentContigs
                    .map { n, f -> f }
                    .collect()
                    .toList()
            )

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
      genomes/*.fasta
    """
}


process selectBestContigsForComponents {

    label "python3"
    label "small_task"
    time "3h"

    publishDir "${params.outdir}/initial/selected_contigs"

    input:
    file "gafs/*" from alignedComponentContigs
        .map { c, g -> g }
        .collect()

    file "fastas/*.fasta" from unscaffoldedGenomes4SelectBestContigsForComponents
        .map { n, f -> f }
        .collect()

    output:
    file "components/*.fasta" into bestContigsForComponents mode flatten
    file "contigs_to_components.tsv"

    script:
    """
    select_best_contigs_for_components.py \
      --min-coverage 0.5 \
      --outfile contigs_to_components.tsv \
      gafs/*

    cat fastas/* > combined.fasta

    mkdir -p components
    select_sequences.py \
        --outdir components \
        contigs_to_components.tsv \
        combined.fasta

    rm combined.fasta
    """
}


/*
 * It might be better to do this as an MGSA?
 */
process realignScaffoldsToInitialGraph {

    label "vg"
    label "small_task"
    time "1d"

    tag "${component}"

    publishDir "${params.outdir}/initial/realigned"

    input:
    set val(component),
        file("in.vg"),
        file("in.fasta") from choppedInitialVg
            .map { v -> [v.baseName, v] }
            .join(
                bestContigsForComponents
                    .map { f -> [f.baseName, f] },
                by: 0
            )

    output:
    set val(component),
        file("${component}.vg") into realignedScaffoldsToInitialGraph

    script:
    """
    # MGSA appears to modify this, so checkpointing gets
    # screwed.
    cp -L in.vg in.vg.tmp

    mkdir tmp
    TMPDIR="\${PWD}/tmp"

    vg msga \
      --graph in.vg.tmp \
      --from in.fasta \
      --threads "${task.cpus}" \
      --band-multi 128 \
      --try-at-least 1 \
      --hit-max 500 \
      --max-multimaps 0 \
      --idx-kmer-size 16 \
      --idx-edge-max 3 \
      --idx-doublings 2 \
      --node-max 32 \
      --hit-max 5 \
      --debug \
    > "${component}.vg"

    # --xdrop-alignment

    rm -rf -- tmp
    """
}


/*
process addPathsToInitialGraph {

    label "vg"
    label "small_task"
    time "1d"

    publishDir "${params.outdir}/with_paths"

    tag "${component}"

    input:
    set val(component),
        file("component.vg"),
        file("component.xg"),
        file("component.gcsa"),
        file("gams/*.gam") from realignedScaffoldsToInitialGraph
            .map { c, n, g -> [c, g] }
            .groupTuple(by: 0)

    output:

    script:
    """
    cat gams/*.gam > single.gam

    vg augment \
      --include-paths \
      --cut-softclips \
      component.vg \
      single.gam \
    > tmp.vg

    vg mod \
      --remove-non-path \
      --X 32 \
      --threads "${task.cpus}" \
      tmp.vg \
    > "${component}_with_paths.vg"

    rm -f single.gam tmp.vg
    """
}
 */


/*
process getInitialXG {

    label "vg"
    label "biggish_task"
    time "12h"

    publishDir "${params.outdir}/initial"

    tag "${component}"

    input:
    set val(component),
        file("input.vg") from explodedInitialVg4GetInitialVGXG

    output:
    set val(component),
        file("${component}.xg") into explodedInitialXG

    script:
    """
    mkdir tmp
    TMPDIR="\${PWD}/tmp"
    vg index \
      -x "${component}.xg" \
      --temp-dir ./tmp \
      --threads "${task.cpus}" \
      input.vg

    rm -rf -- tmp
    """
}
 */


/*
process pruneInitialVG {

    label "vg"
    label "big_task"
    time "1d"

    tag "${component}"

    input:
    set val(component),
        file("input.vg") from explodedInitialVg4PruneInitialVG

    output:
    set val(component),
        file("pruned.vg"),
        file("node_mapping") into explodedInitialPrunedVg

    script:
    """
    vg mod -M 32 input.vg > simplified.vg

    vg prune \
      -u simplified.vg \
      -m node_mapping \
      --threads "${task.cpus}" \
    > pruned.vg

    rm -f simplified.vg
    """
}
 */


/*
process getInitialVGGCSA {

    label "vg"
    label "big_task"
    time "1d"

    publishDir "${params.outdir}/initial"

    tag "${component}"

    input:
    set val(component),
        file("pruned.vg"),
        file("node_mapping") from explodedInitialPrunedVg

    output:
    set val(component),
        val("${component}.gcsa") into explodedInitialGCSA

    script:
    """
    mkdir tmp
    TMPDIR="\${PWD}/tmp"

    vg index \
      -g initial.gcsa \
      -f node_mapping \
      --temp-dir ./tmp \
      --threads "${task.cpus}" \
      --progress \
      pruned.vg

    rm -rf -- tmp
    """
}
 */



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
