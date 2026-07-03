# Data provenance

`data/`, `results/`, and `sim/` are **gitignored** (bulk inputs / generated
outputs — not tracked). Large inputs live on the mount and are referenced here;
resolve the mount roots with `R/staging.R::bzea_mounts()` (`pick()` fallback).
This file is the tracked record of where everything comes from.

## Mount roots

| Alias | Resolved by `bzea_mounts()` | Holds |
|-------|-----------------------------|-------|
| `CAS` (cassini) | `/Volumes/rsstu/.../tlaloc/cassini` → `/Volumes/tlaloc/cassini` → `/Volumes/cassini` | pangene counts, gene→pangene map, B73 gene coords (ATLAS) |
| `MNT` (ancestry) | `/Volumes/rsstu/.../BZea/bzeaseq/ancestry` → `/Volumes/BZea/bzeaseq/ancestry` | per-bin `ALT_FREQ` `*_bin_genotypes.tsv` (dense skim, binhmm) |

## Data sources (naming `<Source>-<Caller>`)

| Source | Modality | Marker set | Notes |
|--------|----------|------------|-------|
| Skim | low-cov WGS (~0.4×) | 50K panel (thinned) + dense ~27.6 M wideseq | TP teosinte–B73 wideseq filter; MAF ≥ 0.05 |
| BrB | 3′ RNA-seq (BRB-seq) | wideseq-thinned counts + cassini pangene | expression-driven depth; ASE |
| Target-seq | MolBreeding 45K, high cov | 45K array | completeness / Holland-style cross-check only |

The ~27.6 M wideseq set = Schnable et al. 2023 teosinte/*Tripsacum* VCF →
biallelic → MAF ≥ 0.05 (`wideseq_ref`); GATK `CollectAllelicCounts` tallies
skim/BRB reads at those positions. (memory `snp50k-cohort-provenance`,
`gatk-table-readcount-standard`.)

## Staged subset: skim sweep (`data/`, gitignored, ~33 MB)

`stage_sanity_paint.R` copies the 11-sample vary-skim cohort (B73 control + 10
NILs) into `data/`, organized by **(source × input type)** and consumed by
`analysis/nilhmm_sanity_paint.qmd`. Verified: the note reproduces the original
zealtiger output **row-for-row** (1709 segments, all 6 tracks identical).

```
data/
  ref/                                genome-level, dataset-independent
    gene_to_pangene.tsv               ATLAS: pangene <-> B73 gene map
    b73_gene_coords.tsv               ATLAS: B73 gene coordinates
  skimsweep/                          the vary-skim cohort
    coverage_sweep_members.csv        full sweep table (note filters to vary_skim)
    skim/counts_50k/<skim>.tsv        -> Skim-nNIL, Skim-rtiger
    skim/bins/<skim>_bin_genotypes.tsv-> Skim-binhmm  (PER-SAMPLE only, never the aggregate)
    brb/counts_wideseq/<brb>.tsv      -> BrB-nNIL, BrB-rtiger
    brb/pangene/<species>/<brb>.pangene_counts.tsv  -> BrB-atlas
```

| Staged path | Source |
|-------------|--------|
| `skimsweep/coverage_sweep_members.csv`, `skim/counts_50k/`, `brb/counts_wideseq/` | zealtiger repo |
| `skim/bins/<skim>_bin_genotypes.tsv` | rsstu `BZea/bzeaseq/ancestry/<skim>_bin_genotypes.tsv` (per-sample) |
| `brb/pangene/`, `ref/` | rsstu `tlaloc/cassini/` (`results/<species>/pangene/`, `data/pangene/`, `data/meta/`) |

Re-stage with `Rscript stage_sanity_paint.R` **run with the sandbox disabled** —
both shares live under the `rsstu` automount, which the default sandbox cannot
see (it reports the files as missing).

## Still to freeze (plan B5 — blocks B2.3 / B2.4)

- [ ] **~400 paired BrB + skim cohort manifest** — IDs, donor species, per-source
      depth. Trace from `sample_metatada.csv` / `bzea-sample-donor-metadata` +
      the BRB and skim sample maps. Blocks `03-source-method-comparison.qmd`.
- [ ] Which MolBreeding target-seq samples for the calibration cross-check.
- [ ] B1 anthocyanin phenotype table + NIL panel. Blocks
      `04-b1-mapping-benchmark.qmd`.

## Migration from zealtiger (plan B3)

The sanity-paint sweep inputs currently read from the zealtiger working repo
(`results/sim_calibration/coverage_sweep_members.csv`, `data/rtiger_50K/`,
`results/sim_calibration/brbseq_ks_wideseq/counts/`). Only the **validated,
paper-bound** subsets migrate onto the mount / local `data/`; the `agent/`
scratch tree is not copied.
