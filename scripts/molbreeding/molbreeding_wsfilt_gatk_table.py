#!/usr/bin/env python3
"""Filter the MolBreeding v5 GATK table to wideseq-panel positions.

The wideseq panel (Schnable 2023, teosinte-vs-B73 segregating sites, B73 v5) is
the project's canonical informative-marker set with consistent REF=B73/ALT=teo
polarity. Restricting the MolBreeding markers to it gives a cross-source-
consistent ancestry-marker subset.

Input : data/molbreeding_45k/gatk_table_SNP_v5.tsv         (real-count GATK table, v5)
        wideseq panel positions  (chr<TAB>pos, "chr1".."chr10"; default on rsstu mount)
Output: data/molbreeding_45k/gatk_table_SNP_wsfilt_v5.tsv  (same schema, kept rows only)
        data/molbreeding_45k/wideseq_keep_v5.tsv           (chr,pos kept; provenance)

Usage:  python3 molbreeding_wsfilt_gatk_table.py [gatk_table_v5.tsv] [panel.pos]
"""
import os
import sys
from pathlib import Path

GT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("data/molbreeding_45k/gatk_table_SNP_v5.tsv")
PANEL = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(os.environ.get(
    "WIDESEQ_ALL_POS",
    "/Volumes/rsstu/users/r/rrellan/BZea/bzeaseq/wideseq_ref/wideseq_all.pos"))
OUT = Path("data/molbreeding_45k/gatk_table_SNP_wsfilt_v5.tsv")
KEEP = Path("data/molbreeding_45k/wideseq_keep_v5.tsv")


def main():
    if not PANEL.exists():
        sys.exit(f"ERROR: wideseq panel not found: {PANEL} (mount not available?)")

    # 1. unique MolBreeding (CONTIG, POSITION) keys
    mb = set()
    with GT.open() as fh:
        fh.readline()
        for line in fh:
            c = line.split("\t", 3)
            mb.add((c[1], c[2]))
    print(f"MolBreeding v5 sites: {len(mb):,}")

    # 2. intersect with the wideseq panel (single scan of the big .pos file)
    keep = set()
    scanned = 0
    with PANEL.open() as fh:
        for line in fh:
            scanned += 1
            chrom, _, pos = line.rstrip("\n").partition("\t")
            if (chrom, pos) in mb:
                keep.add((chrom, pos))
    print(f"wideseq panel positions scanned: {scanned:,}")
    print(f"overlap (MolBreeding sites that are wideseq SNPs): "
          f"{len(keep):,} ({100*len(keep)/len(mb):.1f}%)")

    # 3. write kept positions (sorted) + filtered GATK table
    def sort_key(k):                          # numeric chr if "chrN", else stable string order
        chrom, pos = k
        n = chrom[3:] if chrom.startswith("chr") else chrom
        return (int(n) if n.isdigit() else float("inf"), chrom, int(pos))
    with KEEP.open("w") as w:
        for chrom, pos in sorted(keep, key=sort_key):
            w.write(f"{chrom}\t{pos}\n")

    rows = 0
    with GT.open() as fh, OUT.open("w") as w:
        w.write(fh.readline())  # header
        for line in fh:
            c = line.split("\t", 3)
            if (c[1], c[2]) in keep:
                w.write(line)
                rows += 1
    print(f"wrote {OUT} ({rows:,} rows) and {KEEP}")


if __name__ == "__main__":
    sys.exit(main())