#!/usr/bin/env python3
"""Remap the MolBreeding GATK tables from AGPv3 to B73 NAM v5 positions.

Joins gatk_table_<set>_v3.tsv (CONTIG/POSITION = v3) to the liftover map
sites_v5_<set>.tsv (chr_v3,pos_v3 -> chr_v5,pos_v5) produced by
molbreeding_liftover.R. Markers that did not lift uniquely to chr1-10 are dropped.

Output: data/molbreeding_45k/gatk_table_<set>_v5.tsv
  - CONTIG rewritten to "chr1".."chr10" (project v5 convention; skim/wideseq use it)
  - sorted by SAMPLE, then chromosome, then position  (matches allelic_counts50K.tsv,
    so it feeds the same per-sample RTIGER splitter)

Usage: python3 molbreeding_gatk_table_to_v5.py [SNP mSNP ...]   (default: both)
"""
import sys
from pathlib import Path

DIR = Path("data/molbreeding_45k")


def load_lift(setname):
    """(chr_v3, pos_v3) -> (chr_v5_int, pos_v5_int)."""
    f = DIR / f"sites_v5_{setname}.tsv"
    lift = {}
    with f.open() as fh:
        head = fh.readline().rstrip("\n").split("\t")
        ix = {name: i for i, name in enumerate(head)}
        for line in fh:
            c = line.rstrip("\n").split("\t")
            lift[(c[ix["chr_v3"]], c[ix["pos_v3"]])] = (
                int(c[ix["chr_v5"]]), int(c[ix["pos_v5"]]))
    return lift


def remap(setname):
    lift = load_lift(setname)
    src = DIR / f"gatk_table_{setname}_v3.tsv"
    out = DIR / f"gatk_table_{setname}_v5.tsv"

    rows = []
    sites_in = set()
    sites_kept = set()
    with src.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            c = line.rstrip("\n").split("\t")
            sample, contig, pos, ref, alt, rc, ac, tc = c
            sites_in.add((contig, pos))
            v5 = lift.get((contig, pos))
            if v5 is None:
                continue
            chr5, pos5 = v5
            sites_kept.add((contig, pos))
            # sort key (sample, chr, pos); store formatted row
            rows.append((sample, chr5, pos5, ref, alt, rc, ac, tc))

    rows.sort(key=lambda r: (r[0], r[1], r[2]))
    with out.open("w") as w:
        w.write("\t".join(header) + "\n")
        for sample, chr5, pos5, ref, alt, rc, ac, tc in rows:
            w.write(f"{sample}\tchr{chr5}\t{pos5}\t{ref}\t{alt}\t{rc}\t{ac}\t{tc}\n")

    print(f"{setname}: wrote {out}")
    print(f"  v3 sites in table : {len(sites_in):,}")
    print(f"  lifted to v5      : {len(sites_kept):,}  ({100*len(sites_kept)/len(sites_in):.1f}%)")
    print(f"  dropped (no lift) : {len(sites_in)-len(sites_kept):,}")
    print(f"  rows written      : {len(rows):,}")


def main():
    for s in (sys.argv[1:] or ["SNP", "mSNP"]):
        remap(s)


if __name__ == "__main__":
    sys.exit(main())