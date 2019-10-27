#!/usr/bin/env bash

set -euo pipefail


fasta_to_tsv() {
  awk '
    /^>/ {
      b=gensub(/^>\s*(\S+).*$/, "\\1", "g", $0);
      printf("%s%s\t", (N>0?"\n":""), b);
      N++;
      next;
    }
    {
      printf("%s", $0)
    }
    END {
      printf("\n");
    }
  '
}


tsv_to_lengths() {
  awk -F '\t' -v fname="${1}" '
    BEGIN {total=0; nseqs=0}
    {
      total = total + length($2);
      nseqs++;
    }
    END { printf("%s\t%f\n", fname, total / nseqs) }
  '
}


for f in $@
do
  cat "${f}" \
  | fasta_to_tsv \
  | tsv_to_lengths "${f}"
done \
| sort -k2,2nr \
| awk '{print $1}'
