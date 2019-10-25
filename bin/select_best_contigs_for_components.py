#!/usr/bin/env python3

import sys
import argparse
from os.path import split as psplit
from os.path import splitext

from collections import defaultdict
from intervaltree import Interval, IntervalTree


def cli(prog, args):
    parser = argparse.ArgumentParser(
        prog=prog,
        description="""
        Selects the best matches for a GAF by coverage.
        """
    )

    parser.add_argument(
        "infiles",
        nargs="+",
        type=str,
        help="Input fasta files.",
    )

    parser.add_argument(
        "-o", "--outfile",
        default=sys.stdout,
        type=argparse.FileType('w'),
        help="Output text file path. Default stdout.",
    )

    parser.add_argument(
        "-m", "--min-coverage",
        default=0.5,
        type=float,
        help=("The minimum coverage that a contig should have for it to "
              "be included. Excluded contig names will be printed to stderr."),
    )

    return parser.parse_args(args)


def parse_gaf(handle):
    columns = ["query", "qlen", "qstart", "qend", "strand", "path",
               "plen", "pstart", "pend", "nmatch", "alilen", "mq"]

    for line in handle:
        sline = line.strip().split("\t")
        dline = dict(zip(columns, sline))

        dline["qlen"] = int(dline["qlen"])
        dline["qstart"] = int(dline["qstart"])
        dline["qend"] = int(dline["qend"])

        yield dline["query"], dline["qlen"], dline["qstart"], dline["qend"]

    return


def main():
    args = cli(sys.argv[0], sys.argv[1:])
    qlens = dict()
    qintervals = defaultdict(list)

    for infile in args.infiles:
        component = psplit(splitext(infile)[0])[-1]
        with open(infile, "r") as handle:
            for query, qlen, qstart, qend in parse_gaf(handle):
                qlens[query] = qlen
                qintervals[(component, query)].append(Interval(qstart, qend))

    coverages = defaultdict(list)

    for (component, query), intervals in qintervals.items():
        itree = IntervalTree(intervals)
        itree.merge_overlaps()

        count = sum((i.end - i.begin for i in itree))
        coverages[query].append((component, count / qlens[query]))

    for query, components in coverages.items():
        max_cov = max((cv for cp, cv in components))

        if max_cov < args.min_coverage:
            print(f"unplaced\t{query}\t{max_cov}", file=args.outfile)
        else:
            max_comp = [cp for cp, cv in components if cv == max_cov][0]
            print(f"{max_comp}\t{query}\t{max_cov}", file=args.outfile)

    return


if __name__ == "__main__":
    main()
