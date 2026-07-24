#!/usr/bin/env bash
# Stage the MolBreeding target-seq + SNP50K-skim calibration inputs from the
# EXPLORATORY zealtiger repo into the PUBLICATION repo (zealhmm), so the
# molbreeding-truth vs skim nnil calibration reproduces with no cross-repo
# dependency. Provenance is recorded in DATA.md. Idempotent (cp -f).
#
# The full skim count set is ~1.5 GB in zealtiger; only the 16 calibration
# samples (the molbreeding-shared lines, from the 3-way correspondence) are copied.
#
# Override the source repo with ZEALTIGER=/path. Usage: bash scripts/stage_molbreeding_from_zealtiger.sh
set -euo pipefail
ZT="${ZEALTIGER:-$HOME/repos/zealtiger}"
ZH="$HOME/repos/zealhmm"
MB="$ZH/data/zeal/molbreeding"
SK="$ZH/data/zeal/skim"
CO="$ZH/data/zeal/correspondence"
SC="$ZH/scripts/molbreeding"
mkdir -p "$MB/source" "$MB/counts_targetseq" "$SK/counts" "$CO" "$SC"

# --- raw vendor hard-call source (SNP set; the genotype() columns = GBTS hard calls)
cp -f "$ZT/data/GSER2026030032P01/SNP/All.Genotype.xls" "$MB/source/All.Genotype.xls"

# --- v5 wsfilt count table (also the site REF=B73/ALT=donor polarity authority) + well->PN->pedigree map
cp -f "$ZT/data/molbreeding_45k/gatk_table_SNP_wsfilt_v5.tsv" "$MB/gatk_table_SNP_wsfilt_v5.tsv"
cp -f "$ZT/data/molbreeding_45k_sample_map.tsv" "$MB/molbreeding_sample_map.tsv"

# --- MolBreeding target-seq TRUTH counts (~110x, 16 samples, wsfilt grid)
cp -f "$ZT"/results/sim_calibration/molb_calls/SNP_wsfilt/counts/*.tsv "$MB/counts_targetseq/"

# --- cross-source correspondence + metadata (join key = canonical pedigree; PN_SID is per-experiment)
cp -f "$ZT/results/sample_correspondence/molbreeding_3way_correspondence.csv" "$CO/"
cp -f "$ZT/results/sample_correspondence/skim_brbseq_correspondence.csv" "$CO/"
cp -f "$ZT/data/skim_sample_pedigree.csv" "$CO/"
cp -f "$ZT/data/sample_metadata_master.csv" "$CO/"
cp -f "$ZT/data/brbseq_metadata_master.csv" "$CO/"

# --- skim: seqlengths + RTIGER reference calls + ONLY the 16 calibration count files
cp -f "$ZT/data/rtiger_50K/seqlengths.csv" "$SK/seqlengths.csv"
cp -f "$ZT/data/rtiger_50K/calls_taxa_r5.csv" "$SK/calls_taxa_r5.csv"
# skim ids of the molbreeding-shared samples = skim_prefix (col 4) of the 3-way correspondence
skim_ids=$(tail -n +2 "$CO/molbreeding_3way_correspondence.csv" | cut -d, -f4)
for id in $skim_ids; do
  src=$(find "$ZT/data/rtiger_50K/counts" -name "${id}.tsv" | head -1)
  if [ -n "$src" ]; then cp -f "$src" "$SK/counts/${id}.tsv"; else echo "WARN: skim counts not found for $id"; fi
done

# --- provenance scripts (how the v3 melt / v5 lift / wsfilt / truth counts were produced)
cp -f "$ZT"/molbreeding_to_gatk_table.py "$ZT"/molbreeding_gatk_table_to_v5.py \
      "$ZT"/molbreeding_wsfilt_gatk_table.py "$ZT"/molbreeding_liftover.R \
      "$ZT"/fit_rtiger_molbreeding.R "$SC/"

echo "=== staged manifest ==="
find "$MB/source" "$MB/counts_targetseq" "$SK" "$CO" "$SC" -type f | sort
echo "molb truth counts: $(ls "$MB/counts_targetseq" | wc -l | tr -d ' ')  |  skim counts: $(ls "$SK/counts" | wc -l | tr -d ' ')"
