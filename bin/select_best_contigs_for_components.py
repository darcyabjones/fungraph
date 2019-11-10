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

    parser.add_argument(
        "-n", "--min-scaffolds",
        default=2,
        type=int,
        help="The minimum number of scaffolds assigned to a contig."
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

    best_components = defaultdict(list)

    for query, components in coverages.items():
        matches = [cv for cp, cv in components]

        if len(matches) > 0:
            max_cov = max(matches)
        else:
            max_cov = 0

        if max_cov < args.min_coverage:
            best_components["unplaced"].append((query, max_cov))
        else:
            max_comp = [cp for cp, cv in components if cv == max_cov][0]
            best_components[max_comp].append((query, max_cov))

    to_reassign = []
    to_drop = set()
    for component, assigned in best_components.items():
        if component == "unplaced":
            continue

        if len(assigned) < args.min_scaffolds:
            to_reassign.extend([s for s, c in assigned])
            to_drop.add(component)

    for component in to_drop:
        del best_components[component]

    for query in to_reassign:
        components = coverages[query]
        matches = [cv for cp, cv in components if cp not in to_drop]
        if len(matches) > 0:
            max_cov = max(matches)
        else:
            max_cov = 0

        if max_cov < args.min_coverage:
            best_components["unplaced"].append((query, max_cov))
        else:
            max_comp = [cp for cp, cv in components if cv == max_cov][0]
            best_components[max_comp].append((query, max_cov))

    for component, coverages in best_components.items():
        for query, max_cov in coverages:
            print(f"{component}\t{query}\t{max_cov}", file=args.outfile)

    return


if __name__ == "__main__":
    main()
