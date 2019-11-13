#!/usr/bin/env python3

import sys
import argparse

from collections import defaultdict

from intervaltree import Interval, IntervalTree


class PAF(object):

    columns = [
        "query",
        "qlen",
        "qstart",
        "qend",
        "strand",
        "target",
        "tlen",
        "tstart",
        "tend",
        "nmatch",
        "alilen",
        "mq",
    ]

    def __init__(
            self,
            query,
            qlen,
            qstart,
            qend,
            strand,
            target,
            tlen,
            tstart,
            tend,
            nmatch,
            alilen,
            mq,
            attrs=[],
    ):
        self.query = query
        self.qlen = qlen
        self.qstart = qstart
        self.qend = qend
        self.strand = strand
        self.target = target
        self.tlen = tlen
        self.tstart = tstart
        self.tend = tend
        self.nmatch = nmatch
        self.alilen = alilen
        self.mq = mq
        self.attrs = attrs
        return

    @classmethod
    def parse(cls, line):
        sline = line.strip().split("\t")
        if len(sline) < len(cls.columns):
            raise AssertionError(f"Wrong number of columns {line}")

        dline = dict(zip(cls.columns, sline))
        qlen = int(dline["qlen"])
        qstart = int(dline["qstart"])
        qend = int(dline["qend"])
        tlen = int(dline["tlen"])
        tstart = int(dline["tstart"])
        tend = int(dline["tend"])
        nmatch = int(dline["nmatch"])
        alilen = int(dline["alilen"])
        mq = int(dline["mq"])
        attrs = sline[len(cls.columns):]
        return cls(dline["query"], qlen, qstart, qend, dline["strand"],
                   dline["target"], tlen, tstart, tend, nmatch, alilen,
                   mq, attrs)

    @classmethod
    def from_file(cls, handle):
        for line in handle:
            yield cls.parse(line)
        return

    def __str__(self):
        line = [str(getattr(self, c)) for c in self.columns]
        line.extend(self.attrs)
        return "\t".join(line)

    def __repr__(self):
        args = ", ".join(str(getattr(self, c)) for c in self.columns)
        return f"PAF({args}, attrs={self.attrs})"

    def query_as_interval(self):
        start = min([self.qstart, self.qend])
        end = max([self.qstart, self.qend])
        return self.query, Interval(start, end)

    def target_as_interval(self):
        start = min([self.tstart, self.tend])
        end = max([self.tstart, self.tend])
        return self.target, Interval(start, end)


class BED(object):

    def __init__(self, seqid, start, end):
        self.seqid = seqid
        self.start = start
        self.end = end
        return

    def __str__(self):
        return f"{self.seqid}\t{self.start}\t{self.end}"

    def __repr__(self):
        return f"BED({self.seqid}, {self.start}, {self.end})"

    @classmethod
    def parse(cls, line):
        sline = line.strip().split("\t")
        return cls(sline[0], int(sline[1]), int(sline[2]))

    @classmethod
    def from_file(cls, handle):
        for line in handle:
            yield cls.parse(line)
        return

    def as_interval(self):
        return self.seqid, Interval(self.start, self.end)


def bed_to_itree(beds):
    intervals = defaultdict(list)

    for record in beds:
        seqid, interval = record.as_interval()
        intervals[seqid].append(interval)

    itree = {s: IntervalTree(i) for s, i in intervals.items()}

    for it in itree.values():
        it.merge_overlaps()
    return itree


def intersect(left, right):
    lstart = min([left.begin, left.end])
    lend = max([left.begin, left.end])
    rstart = min([right.begin, right.end])
    rend = max([right.begin, right.end])
    lbound = max([lstart, rstart])
    rbound = min([lend, rend])
    return Interval(lbound, rbound)


def union(left, right):
    lstart = min([left.begin, left.end])
    lend = max([left.begin, left.end])
    rstart = min([right.begin, right.end])
    rend = max([right.begin, right.end])

    if not ((lstart < rend) and (lend > rstart)):
        raise ValueError("Left and right must be overlapping.")

    lbound = min([lstart, rstart])
    rbound = max([lend, rend])
    return Interval(lbound, rbound)


def diff(left, right):
    if not left.overlaps(right):
        return [left, right]

    lstart = min([left.begin, left.end])
    lend = max([left.begin, left.end])
    rstart = min([right.begin, right.end])
    rend = max([right.begin, right.end])

    out = []
    if rstart > lstart:
        out.append(Interval(lstart, min([rstart, lend])))

    if rend < lend:
        out.append(Interval(max([lstart, rend]), lend))

    return out


def sym_diff(left, right):
    inter = intersect(left, right)
    uni = union(left, right)
    return diff(uni, inter)


def total_intersection(itree, interval):
    if interval.length() <= 0:
        return 0

    total = 0
    ovlps = IntervalTree(itree.overlap(interval))
    ovlps.merge_overlaps()
    for ovlp in ovlps:
        inter = intersect(interval, ovlp)
        total += inter.length()

    return total


def get_genome_name(string, sep):
    return string.split(sep, maxsplit=1)[0]


def filter_by_interval(itree, interval, min_length, prop_coverage):
    len_intersect = total_intersection(itree, interval)

    if interval.length() <= 0:
        prop_intersect = 0
    else:
        prop_intersect = len_intersect / interval.length()

    return ((interval.length() - len_intersect) < min_length) or (prop_intersect >= prop_coverage)


def cli(prog, args):

    parser = argparse.ArgumentParser(
        prog=prog,
        description=""
    )

    subparsers = parser.add_subparsers(dest='subparser_name')

    filter_subparser = subparsers.add_parser(
        "filter",
        help="Filter matches."
    )

    filter_cli(filter_subparser)

    repeats_subparser = subparsers.add_parser(
        "repeats",
        help="Find segmental duplications."
    )

    repeats_cli(repeats_subparser)

    parsed = parser.parse_args(args)

    if parsed.subparser_name is None:
        parser.print_help()
        sys.exit(0)

    return parsed


def repeats_cli(parser):
    parser.add_argument(
        "inpaf",
        default=sys.stdin,
        type=argparse.FileType('r'),
        help="Input paf file. Use '-' for stdin.",
    )

    parser.add_argument(
        "-o", "--outfile",
        default=sys.stdout,
        type=argparse.FileType('w'),
        help="Output bed file path. Default stdout.",
    )

    parser.add_argument(
        "-s", "--sep",
        default=".",
        type=str,
        help=(
            "Split sequence ids by this separator and exclude matches where "
            "both members have the same prefix."
        )
    )
    return


def filter_cli(parser):
    parser.add_argument(
        "inbed",
        default=sys.stdin,
        type=argparse.FileType('r'),
        help="Input bed file. Use '-' for stdin.",
    )
    parser.add_argument(
        "inpaf",
        default=sys.stdin,
        type=argparse.FileType('r'),
        help="Input paf file. Use '-' for stdin.",
    )

    parser.add_argument(
        "-o", "--outfile",
        default=sys.stdout,
        type=argparse.FileType('w'),
        help="Output paf file path. Default stdout.",
    )

    parser.add_argument(
        "-m", "--min-length",
        default=1,
        type=int,
        help="The minimum allowed match length excluding repeats.",
    )

    parser.add_argument(
        "-s", "--sep",
        default=None,
        type=str,
        help=(
            "Split sequence ids by this separator and exclude matches where "
            "both members have the same prefix."
        )
    )

    parser.add_argument(
        "-p", "--prop-overlap",
        default=0.5,
        type=float,
        help="The maximum proportion of bases allowed to be in repeat regions."
    )

    return


def filter_main(args):

    bed = BED.from_file(args.inbed)
    bed_tree = bed_to_itree(bed)

    paf = PAF.from_file(args.inpaf)
    for p in paf:
        if args.sep is not None:
            qgenome = get_genome_name(p.query, args.sep)
            tgenome = get_genome_name(p.target, args.sep)

            if qgenome == tgenome:
                continue

        if p.alilen < args.min_length:
            continue

        query, qinterval = p.query_as_interval()
        if filter_by_interval(bed_tree[query], qinterval,
                              args.min_length, args.prop_overlap):
            continue

        target, tinterval = p.target_as_interval()
        if filter_by_interval(bed_tree[target], tinterval,
                              args.min_length, args.prop_overlap):
            continue

        print(p, file=args.outfile)

    return


def repeats_main(args):
    paf = PAF.from_file(args.inpaf)

    repeats = defaultdict(list)

    for p in paf:
        qgenome = get_genome_name(p.query, args.sep)
        tgenome = get_genome_name(p.target, args.sep)

        if qgenome != tgenome:
            continue

        query, qinterval = p.query_as_interval()
        target, tinterval = p.target_as_interval()

        if (query == target) and qinterval.overlaps(tinterval):
            filtered = sym_diff(qinterval, tinterval)
            repeats[query].extend(filtered)
        else:
            repeats[query].append(qinterval)
            repeats[target].append(tinterval)

    for seqid, intervals in repeats.items():
        itree = IntervalTree(intervals)
        itree.merge_overlaps()
        for interval in itree:
            bed = BED(seqid, interval.begin, interval.end)
            print(bed, file=args.outfile)
    return


def main():
    args = cli(prog=sys.argv[0], args=sys.argv[1:])
    try:
        if args.subparser_name == "filter":
            filter_main(args)
        elif args.subparser_name == "repeats":
            repeats_main(args)
        else:
            raise ValueError("I shouldn't reach this point ever")

    except BrokenPipeError:
        # Pipes get closed and that's normal
        sys.exit(0)

    except KeyboardInterrupt:
        print("Received keyboard interrupt. Exiting.", file=sys.stderr)
        sys.exit(1)

    except EnvironmentError as e:
        print((
            "Encountered a system error.\n"
            "We can't control these, and they're usually related to your OS.\n"
            "Try running again."
        ), file=sys.stderr)
        raise e

    except Exception as e:
        print((
            "I'm so sorry, but we've encountered an unexpected error.\n"
            "This shouldn't happen, so please file a bug report with the "
            "authors.\nWe will be extremely grateful!\n\n"
        ), file=sys.stderr)
        raise e

    return


if __name__ == "__main__":
    main()
