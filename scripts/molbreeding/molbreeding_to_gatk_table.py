#!/usr/bin/env python3
"""Melt the MolBreeding GBTS delivery (wide, .xls=TSV) into the canonical
GATK table (tidy/long readcount format) used across the BZea analyses.

Canonical "GATK table" columns (CAPS, from GATK CollectAllelicCounts; the
input format for the ancestry analysis at rpubs.com/faustovrz/1337797):

    SAMPLE  CONTIG  POSITION  REF_NUCLEOTIDE  ALT_NUCLEOTIDE  REF_COUNT  ALT_COUNT  TOTAL_COUNT

Source wide file: <set>/All.Genotype.xls
    ID  chrom  position  ref  [genotype(S) ref_depth(S) alt_depth(S) mutation_frequence(S)]xN

TRUST THE COUNT DATA. Everything is read from All.Genotype.xls itself:
    REF_NUCLEOTIDE <- ref column
    ALT_NUCLEOTIDE <- the non-ref base observed across the genotype calls
                      (a per-SITE property; biallelic SNPs. A hom-ref sample at
                      a variant site is REF_COUNT=total, ALT_COUNT=0 with ALT
                      still defined.)
    REF_COUNT <- ref_depth   ALT_COUNT <- alt_depth   TOTAL_COUNT <- ref+alt
    CONTIG <- chrom (AGPv3; liftover to v5 chr# downstream)   POSITION <- position

snp.annot.xls is NOT used here — it was only a cross-check
(agent/validate_molbreeding_alt.py: REF/ALT and depth polarization agree).

Only clean BIALLELIC SNP markers are kept. A site is DROPPED if:
  - REF is not a single A/C/G/T base;
  - nonvariant: no sample's genotype carries a non-ref base;
  - multiallelic: >1 distinct non-ref base across the genotypes;
  - any genotype base is not A/C/G/T (indel/MNP).

Output: data/molbreeding_45k/gatk_table_<set>_v3.tsv  (build tagged in name)
Usage:  python3 molbreeding_to_gatk_table.py [SNP|mSNP ...]   (default: both)
"""
import sys
from pathlib import Path

BASE = Path("data/GSER2026030032P01")
OUTDIR = Path("data/molbreeding_45k")
HEADER = ["SAMPLE", "CONTIG", "POSITION", "REF_NUCLEOTIDE",
          "ALT_NUCLEOTIDE", "REF_COUNT", "ALT_COUNT", "TOTAL_COUNT"]
ACGT = {"A", "C", "G", "T"}
MISSING_GT = {"", "NN", "NA", "./.", ".", "--"}


def melt(setname):
    src = BASE / setname / "All.Genotype.xls"
    out = OUTDIR / f"gatk_table_{setname}_v3.tsv"
    OUTDIR.mkdir(parents=True, exist_ok=True)

    d = dict(sites=0, kept=0, ref_non_acgt=0, nonvariant=0,
             multiallelic=0, gt_non_acgt=0)
    with src.open() as fh, out.open("w") as w:
        header = fh.readline().rstrip("\n").split("\t")
        samples = [h[len("genotype("):-1] for h in header[4::4]
                   if h.startswith("genotype(")]
        n = len(samples)
        w.write("\t".join(HEADER) + "\n")
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 4 + 4 * n:
                continue
            d["sites"] += 1
            contig, pos, ref = c[1], c[2], c[3]
            if ref not in ACGT:
                d["ref_non_acgt"] += 1
                continue
            # ALT = non-ref base(s) observed across all genotype calls
            alt_bases = set()
            bad_base = False
            depths = []  # (ref_depth, alt_depth) per sample
            for i in range(n):
                b = 4 + i * 4
                gt, rd, ad = c[b], c[b + 1], c[b + 2]
                depths.append((rd, ad))
                if not gt or gt.upper() in MISSING_GT:
                    continue
                for base in gt:
                    if base not in ACGT:
                        bad_base = True
                    elif base != ref:
                        alt_bases.add(base)
            if bad_base:
                d["gt_non_acgt"] += 1
                continue
            if not alt_bases:
                d["nonvariant"] += 1
                continue
            if len(alt_bases) > 1:
                d["multiallelic"] += 1
                continue
            alt = next(iter(alt_bases))
            for i in range(n):
                rd, ad = depths[i]
                try:
                    tot = int(rd) + int(ad)
                except ValueError:
                    rd, ad, tot = "0", "0", 0
                w.write(f"{samples[i]}\t{contig}\t{pos}\t{ref}\t{alt}\t{rd}\t{ad}\t{tot}\n")
            d["kept"] += 1

    print(f"{setname}: wrote {out}")
    print(f"  sites read       : {d['sites']:,}")
    print(f"  kept (biallelic) : {d['kept']:,}  -> {d['kept'] * n:,} rows ({n} samples)")
    print(f"  dropped nonvariant  : {d['nonvariant']:,}")
    print(f"  dropped multiallelic: {d['multiallelic']:,}")
    print(f"  dropped REF non-ACGT: {d['ref_non_acgt']:,}")
    print(f"  dropped GT  non-ACGT: {d['gt_non_acgt']:,}")


def main():
    for s in (sys.argv[1:] or ["SNP", "mSNP"]):
        melt(s)


if __name__ == "__main__":
    sys.exit(main())