#!/usr/bin/env python3

import sys
import argparse

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def cli(prog, args):
    parser = argparse.ArgumentParser(
        prog=prog,
        description="""
        Splits a fasta file at long stretches of Ns.
        """
    )

    parser.add_argument(
        "infile",
        default=sys.stdin,
        type=argparse.FileType('r'),
        help="Input fasta files. Use '-' for stdin.",
    )

    parser.add_argument(
        "-o", "--outfile",
        default=sys.stdout,
        type=argparse.FileType('w'),
        help="Output fasta file path. Default stdout.",
    )

    parser.add_argument(
        "-n", "--nsize",
        default=50,
        type=int,
        help=("The maximum number of Ns allowable to keep a sequence as a "
              "single contig."),
    )

    parser.add_argument(
        "-m", "--min-length",
        default=1,
        type=int,
        help=("The minimum length that a contig is allowed to be after "
              "splitting."),
    )

    return parser.parse_args(args)


def new_sequence(seq, start, end):
    assert start != end
    assert seq[start].upper() != "N"
    assert seq[end - 1].upper() != "N"

    if start > 0:
        assert seq[start - 1].upper() == "N"

    try:
        assert seq[end].upper() == "N"
    except IndexError:
        pass

    new_id = f"{seq.id}.{start}_{end}"
    new_seq = SeqRecord(
        id=new_id,
        name=new_id,
        description=seq.description if seq.description != seq.id else new_id,
        seq=seq.seq[start:end]
    )
    return new_seq


def split_at_n(seq, length):
    size = 0

    start = None
    end = 0
    for i, base in enumerate(seq):
        if base.upper() == "N":
            size += 1

        # This handles the case where a sequence begins with Ns
        elif start is None:
            start = i
            end = i + 1
            size = 0

        # If we've just finished a stretch of Ns
        elif size >= length:
            yield new_sequence(seq, start, end)

            size = 0
            start = i
            end = i + 1

        # If it's just a regular base reset the N counter and update the end.
        else:
            end = i + 1
            size = 0

    # At the end of the sequence we print the last block.
    # If start is none handles case where sequence is all Ns.
    if start is not None:
        yield new_sequence(seq, start, end)

    return


def main():
    args = cli(sys.argv[0], sys.argv[1:])

    seqs = SeqIO.parse(args.infile, format="fasta")
    for seq in seqs:
        for split_seq in split_at_n(seq, args.nsize):
            if len(split_seq) >= args.min_length:
                SeqIO.write(split_seq, args.outfile, format="fasta")

    return


if __name__ == "__main__":
    main()
