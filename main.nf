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
        genomes4FindComponentContigsAlign;
        genomes4SelectBestContigsForComponents;
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
    // This garbage is all needed to so that the two channels are merged
    // even if one is empty.
    set file("genomes/*"),
        file("complete_genomes/*") from genomes4AssembleMinigraph
            .map { n, f -> f }
            .collect()
            .ifEmpty([])
            .toList()
            .merge(
                completeGenomes4AssembleMinigraph
                    .map { n, f -> f }
                    .collect()
                    .ifEmpty([])
                    .toList()
            )
            .collect()



    output:
    file "pan_minigraph.gfa" into minigraphAssembly
    file "pan_minigraph.rgfa" into minigraphRGFAAssembly
    file "genome_order.tsv" into genomeOrder

    script:
    """
    if [ -d complete_genomes ]
    then
      COMPLETE_ORDER=\$(order_genome_by_average_length.sh complete_genomes/*)
    else
      COMPLETE_ORDER=""
    fi

    if [ -d genomes ]
    then
      ORDER=\$(order_genome_by_average_length.sh genomes/*)
    else
      ORDER=""
    fi

    ORDER_ARR_FILES=\${COMPLETE_ORDER} \${ORDER}

    i=0
    touch genome_order.tsv
    for n in \${ORDER_ARR_FILES}
    do
      NAME=\$(basename \${n%.*})
      echo -e "\${NAME}\\t\${i}" >> genome_order.tsv
      i=\$(( i+1 ))
    done

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
process gfa2InitialComponentVG {

    label "vg"
    label "medium_task"
    time "5h"

    input:
    file "pan.gfa" from minigraphAssembly

    output:
    file "initial/*" into initialComponentVgs mode flatten

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
    """
}


initialComponentVgs
    .map { v -> [v.baseName, v] }
    .into {
        initialComponentVgs4GetInitialComponentGFAs;
        initialComponentVgs4RealignScaffoldsToInitialGraph;
    }


process getInitialComponentGFAs {

    label "vg"
    label "small_task"
    time "2h"

    tag "${component}"

    publishDir "${params.outdir}/initial/component_gfas"

    input:
    set val(component),
        file("component.vg") from initialComponentVgs4GetInitialComponentGFAs

    output:
    set val(component),
        file("${component}.gfa") into initialComponentGFAs

    script:
    """
    vg view --vg-in --gfa "component.vg" > "${component}.gfa"
    """
}


process findComponentContigsAlign {

    label "minigraph"
    label "medium_task"
    time "5h"

    tag "${component}"

    publishDir "${params.outdir}/initial/component_gafs"

    input:
    set val(component),
        file("component.gfa"),
        file("genomes/*.fasta") from initialComponentGFAs
            .combine(
                genomes4FindComponentContigsAlign
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

    tag "${name}"

    publishDir "${params.outdir}/initial"

    input:
    set val(name),
        file("in.fasta"),
        file("gafs/*") from genomes4SelectBestContigsForComponents
            .combine(alignedComponentContigs.map { c, g -> g }.collect().toList())

    output:
    set val(name),
        file("fasta_components_${name}/*.fasta") into bestContigsForComponents mode flatten

    set val(name),
        file("unplaced_${name}.fasta") into unplacedContigs optional true

    file "contigs_to_components_${name}.tsv"

    script:
    """
    select_best_contigs_for_components.py \
      --min-coverage 0.7 \
      --outfile "contigs_to_components_${name}.tsv" \
      gafs/*


    mkdir -p "fasta_components_${name}"

    select_sequences.py \
      --outdir "fasta_components_${name}" \
      "contigs_to_components_${name}.tsv" \
      in.fasta

    if [ -s "fasta_components_${name}/unplaced.fasta" ]
    then
      mv "fasta_components_${name}/unplaced.fasta" "unplaced_${name}.fasta"
    fi
    """
}


// This sorts the channel so that contigs are in the same order that the
// original minigraph assembly was in.
// This does three things, 1 we get the major contributors to the graph first.
// and 2 we align the more contiguous assemblies first, which should help with
// the alignments via paths. and 3 maintains the order between runs so checkpointing
// is preserved.

bestContigsForComponents
    .map { n, c -> [n, c.baseName, c] }
    .combine(genomeOrder.splitCsv(sep: "\t"), by: 0)
    .toSortedList( { a, b -> a[3] <=> b[3] } )
    .flatMap { li -> li.collect { n, c, f, o -> [c, n, f] } }
    .into {
        sortedBestContigsForComponents;
        sortedBestContigsForComponents4FindLength;
        sortedBestContigsForComponents4FindLast;
    }

sortedBestContigsForComponents4FindLength
    .count()
    .set { numRealignments }

sortedBestContigsForComponents4FindLast
    .map { c, n, f -> [c, n] }
    .groupTuple(by: 0)
    .map { c, ns -> [c, ns[-1]] }
    .set { lastComponentNames }


// The component filename will match the component label.

/*
 * We are aligning each genome to the graph and updating
 * the graph before aligning the next one.
 * To do this I add a cycle.
 */
realignmentAccumulator = Channel.create()

process realignScaffoldsToInitialGraph {

    label "vg"
    label "small_task"
    time "1d"


    tag "${component} - ${name}"

    input:
    // The merge channel here is just to provide a stop condition for the recursion.
    set val(component),
        file("in.vg"),
        val(name),
        file("in.fasta") from initialComponentVgs4RealignScaffoldsToInitialGraph
            .mix(realignmentAccumulator)
            .join(sortedBestContigsForComponents, by: 0)
            .merge(Channel.from( 1..numRealignments.val )
            .map { c, v, n, f, i -> [c, v, n, f] }

    output:
    set val(component),
        file("out.vg") into realignmentAccumulator

    set val(component),
        val(name),
        file("out.vg") into realignedScaffoldsToInitialGraph

    script:
    """
    mkdir -p tmp
    TMPDIR="\${PWD}/tmp"

    vg mod -X 512 in.vg > chopped.vg

    vg msga \
      --graph chopped.vg \
      --from in.fasta \
      --threads "${task.cpus}" \
      --band-multi 128 \
      --try-at-least 1 \
      --hit-max 500 \
      --max-multimaps 1 \
      --idx-kmer-size 16 \
      --idx-edge-max 3 \
      --idx-doublings 3 \
      --node-max 512 \
      --debug \
    > msga.vg

    vg mod --unchop msga.vg > unchopped.vg
    vg mod --until-normal 10 unchopped.vg > normalised.vg
    vg mod --compact-ids normalised.vg > out.vg

    rm -f msga.vg unchopped.vg normalised.vg
    rm -rf -- tmp
    """
}


process normaliseRealignedScaffolds {

    label "vg"
    label "small_task"
    time "4h"

    tag "${component}"

    input:
    set val(component),
        file("in.vg") from realignedScaffoldsToInitialGraph
            .combine(lastComponentNames, by: [0, 1])
            .filter { c, n, v, l -> n == l }
            .map { c, n, v, l -> [c, v] }

    output:
    set val(component),
        file("out.vg") into realignedComponents

    script:
    """
    vg mod --remove-non-path in.vg > path_only.vg
    vg mod --unchop path_only.vg > unchopped.vg
    vg mod --until-normal 10 unchopped.vg > normalised.vg
    vg mod --compact-ids normalised.vg > out.vg

    rm -f path_only.vg unchopped.vg normalised.vg
    """
}


process combineRealignedComponents {

    label "vg"
    label "small_task"

    input:
    file "in/*" from realignedComponents
        .map { c, v -> v }
        .collect()

    output:
    file "combined.vg" into combinedRealignedComponents

    script:
    """
    cp -rL in out

    vg ids -j out/*.vg
    cat out/*.vg > combined.vg
    """
}


unplacedContigs
    .tap { unplacedContigs4RealignScaffolds }
    .count()
    .set { numCombinedRealignments }

if ( numCombinedRealignments.val > 0 ) {

    combinedRealignmentAccumulator = Channel.create()

    process realignScaffoldsToCombinedComponents {

        label "vg"
        label "small_task"
        time "1d"

        tag "${name}"

        input:
        // The merge channel here is just to provide a stop condition for the recursion.
        set file("in.vg"),
            val(name),
            file("in.fasta") from combineRealignedComponents
                .mix(combineRealignedAccumulator)
                .combine(unplacedContigsForComponents)
                .merge(Channel.from( 1..numCombinedRealignments.val )
                .map { v, n, f, i -> [v, n, f] }

        output:
        file "out.vg"  into combinedRealignmentAccumulator
        file "out.vg"  into realignedScaffoldsToCombinedGraph

        script:
        """
        mkdir -p tmp
        TMPDIR="\${PWD}/tmp"

        vg mod -X 512 in.vg > chopped.vg

        vg msga \
          --graph chopped.vg \
          --from in.fasta \
          --threads "${task.cpus}" \
          --band-multi 128 \
          --try-at-least 1 \
          --hit-max 500 \
          --max-multimaps 1 \
          --idx-kmer-size 16 \
          --idx-edge-max 3 \
          --idx-doublings 3 \
          --node-max 512 \
          --debug \
        > msga.vg

        vg mod --unchop msga.vg > unchopped.vg
        vg mod --until-normal 10 unchopped.vg > normalised.vg
        vg mod --compact-ids normalised.vg > out.vg

        rm -f msga.vg unchopped.vg normalised.vg
        rm -rf -- tmp
        """
    }

    realignedScaffoldsToCombinedGraph.last().set {finalisedRealignedGraph}

} else {

    combineRealignedComponent.set {finalisedRealignedGraph}

}

finalisedRealignedGraph.set {
    finalisedRealignedGraph4Split;
    finalisedRealignedGraph4Index;
}


process splitFinalisedGraph {

    label "vg"
    label "medium_task"
    time "6h"

    input:
    file "realigned.vg" from finalisedRealignedGraph4Split

    output:
    file "components/*.vg" into finalisedRealignedComponents mode flatten

    script:
    """
    vg view --gfa --vg-in realigned.vg > realigned.gfa
    vg explode --threads "${task.cpus}" realigned.vg components

    for f in components/*
    do
      NOEXT="\${f%.*}"
      vg view --gfa --vg-in "\${f}" > "\${NOEXT}.gfa"
    done
    """
}


/*
 * Compute the indices needed to align against the genome.
 */
process indexRealignedVG {

    label "vg"
    label "biggish_task"
    time "12h"

    // This will be overwritten at each iteration
    publishDir "${params.outdir}/realigned"

    input:
    file "realigned.vg" from finalisedRealignedGraph4Index

    output:
    set file("realigned_chopped.vg"),
        file("realigned_chopped.xg"),
        file("realigned_chopped.gcsa"),
        file("realigned_chopped.gcsa.lcp") into indexedRealignedGraph

    file "realigned.vg"

    script:
    """
    mkdir -p tmp
    TMPDIR="\${PWD}/tmp"

    vg mod -X 32 realigned.vg > realigned_chopped.vg

    vg index \
      -x "realigned_chopped.xg" \
      --temp-dir ./tmp \
      --threads "${task.cpus}" \
      realigned_chopped.vg

    vg mod -M 32 realigned_chopped.vg > simplified.vg

    vg prune \
      -u simplified.vg \
      -m node_mapping \
      --threads "${task.cpus}" \
    > extra_simplified.vg

    vg index \
      -g realigned_chopped.gcsa \
      -f node_mapping \
      --temp-dir ./tmp \
      --threads "${task.cpus}" \
      --progress \
      --kmer-size 16 \
      --doubling-steps 3 \
      extra_simplified.vg

    rm -rf -- tmp
    rm -f simplified.vg extra_simplified.vg node_mapping
    """
}



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
//       -L 0.5 \
//       -A=SN15 \
//       --show-strand \
//       -X 1 \
//       -P 10 \
//       -R \
//       -o pan.png
//     """
// }
