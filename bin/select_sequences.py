#!/usr/bin/env python3

import sys
import argparse
from os.path import join as pjoin
from collections import defaultdict

from Bio import SeqIO


def cli(prog, args):
    parser = argparse.ArgumentParser(
        prog=prog,
        description="""
        Select seqs.
        """
    )

    parser.add_argument(
        "table",
        type=argparse.FileType("r"),
        help="Input table",
    )

    parser.add_argument(
        "infiles",
        nargs="+",
        type=argparse.FileType('r'),
        help="Input fasta files.",
    )

    parser.add_argument(
        "-o", "--outdir",
        default="",
        type=str,
        help="Output text file path. Default stdout.",
    )

    return parser.parse_args(args)


def parse_tsv(handle):
    out = dict()

    for line in handle:
        sline = line.strip().split("\t")
        out[sline[1]] = sline[0]

    return out


def main():
    args = cli(sys.argv[0], sys.argv[1:])

    seqid_to_component = parse_tsv(args.table)

    touched_components = {c: False for s, c in seqid_to_component.items()}

    for infile in args.infiles:
        components = defaultdict(list)
        seqs = SeqIO.parse(infile, format="fasta")

        for seq in seqs:
            component = seqid_to_component.get(seq.id, "unplaced")
            components[component].append(seq)

        for component in components:
            if touched_components[component]:
                mode = "a"
            else:
                mode = "w"

            touched_components[component] = True

            filename = pjoin(args.outdir, f"{component}.fasta")
            with open(filename, mode) as handle:
                SeqIO.write(components[component], handle, format="fasta")

    return


if __name__ == "__main__":
    main()
