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


process filterAlignment {

    label "fpa"
    label "small_task"

    input:
    file "aligned.paf.gz" from aligned4Filter

    output:
    file "filtered.paf.gz" into filtered

    script:
    """
    zcat aligned.paf.gz | fpa drop -l 10000 | gzip > filtered.paf.gz
    """
}


process squishAlignments {

    label "seqwish"
    label "big_task"

    input:
    set file("genomes.fasta.gz"),
        file("alignments.paf.gz") from combinedGenomes4Squish
            .combine(filtered)

    output:
    file "x.gfa"

    script:
    """
    seqwish \
        -s "genomes.fasta.gz" \
        -p "alignments.paf.gz" \
        -b x.graph \
        -o x.vgp \
        -g x.gfa \
        --threads "${task.cpus}"
    """
}
