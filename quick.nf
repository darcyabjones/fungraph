params.gfa = false
params.genomes = false

Channel
  .fromPath(params.gfa, checkIfExists: true, type: "file")
  .first()
  .set { gfa }

Channel
  .fromPath(params.genomes, checkIfExists: true, type: "file")
  .map { g -> [g.simpleName, g] }
  .set { genomes }


process gfa2InitialComponentVG {

    label "vg"
    label "medium_task"
    time "5h"

    input:
    file "pan.gfa" from gfa

    output:
    set file("pan.vg"),
        file("pan.xg"),
        file("pan.gcsa"),
        file("pan.gcsa.lcp") into initialVg

    script:
    // Parallelisation is limited for this step.
    """
    vg view \
      --vg \
      --gfa-in pan.gfa \
      --threads "${task.cpus}" \
    > in.vg

    mkdir -p tmp
    TMPDIR="\${PWD}/tmp"

    vg mod -X 32 in.vg > pan.vg

    vg index \
      -x "pan.xg" \
      --temp-dir ./tmp \
      --threads "${task.cpus}" \
      pan.vg

    vg mod -M 32 pan.vg > simplified.vg

    vg prune \
      -u simplified.vg \
      -m node_mapping \
      --threads "${task.cpus}" \
    > extra_simplified.vg

    vg index \
      -g pan.gcsa \
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

initialVg.into {
   initialVg4Align;
   initialVg4Augment;
}

process alignGenomes {

    label "vg"
    label "medium_task"
    time "12h"

    tag "${name}"
    publishDir "${params.outdir}/aligned"

    input:
    set val(name),
        file("in.fasta"),
        file("pan.vg"),
        file("pan.xg"),
        file("pan.gcsa"),
        file("pan.gcsa.lcp") from genomes.combine(initialVg4Align)

    output:
    set val(name),
        file("${name}.gam") into alignedGams

    script:
    """
    vg map \
      -x pan.xg \
      -g pan.gcsa \
      -F in.fasta \
      --threads "${task.cpus}" \
      --band-width 512 \
      -u 2 \
      -L 63 \
      -q 1 \
      -z 2 \
      -o 2 \
      -y 1 \
      -M 1 \
    > "${name}.gam"
    """
}


alignedGams.into {
    alignedGams4Augment;
    alignedGams4Call;
}

process augmentGraph {

    label "vg"
    label "small_task"
    time "6h"

    input:
    set file("pan.vg"),
        file("pan.xg"),
        file("pan.gcsa"),
        file("pan.gcsa.lcp"),
        file("gams/*.gam") from initialVg4Augment
            .merge(alignedGams4Augment.collect().toList())

    output:
    file "augmented_labeled.vg" 
    file "augmented.vg" 

    script:
    """
    vg augment -C --label-paths pan.vg gams/*.gam > tmp1.vg
    vg mod --remove-non-path tmp1.vg > augmented_labeled.vg

    vg augment -C --include-paths pan.vg gams/*.gam > tmp2.vg
    vg mod --remove-non-path tmp2.vg > augmented.vg
    """
}
