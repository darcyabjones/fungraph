#!/usr/bin/env python3

import sys
import argparse
from collections import namedtuple

from Bio import SeqIO

BED = namedtuple("BED", ["seqid", "start", "end"])


def print_bed(handle, row):
    print(f"{row.seqid}\t{row.start}\t{row.end}", file=handle)
    return


def cli(prog, args):
    parser = argparse.ArgumentParser(
        prog=prog,
        description="""
        Finds regions of a fasta file with lowercase characters.
        """
    )

    parser.add_argument(
        "infile",
        default=sys.stdin,
        type=argparse.FileType('r'),
        help="Input fasta file. Use '-' for stdin.",
    )

    parser.add_argument(
        "-o", "--outfile",
        default=sys.stdout,
        type=argparse.FileType('w'),
        help="Output fasta file path. Default stdout.",
    )

    return parser.parse_args(args)


def find_lowercase_stretches(sr):

    start = None
    pos = 0

    for base in sr.seq:
        if base.islower() and start is None:
            start = pos
        elif base.isupper() and start is not None:
            yield BED(sr.id, start, pos)
            start = None

        pos += 1

    if start is not None:
        yield BED(sr.id, start, pos)

    return


def main():
    args = cli(sys.argv[0], sys.argv[1:])

    seqs = SeqIO.parse(args.infile, format="fasta")
    for seq in seqs:
        for bed_row in find_lowercase_stretches(seq):
            print_bed(args.outfile, bed_row)

    return


if __name__ == "__main__":
    main()
