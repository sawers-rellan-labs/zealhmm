# Data manifest — provenance for every input

Large inputs stay on the mount and are referenced here (not committed); only
small inputs + this manifest live under `data/`. Mount roots drift across
machines — resolve them with `R/staging.R::bzea_mounts()` (`pick()` fallback).

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

## Still to freeze (plan B5 — blocks B2.3 / B2.4)

- [ ] **~400 paired BrB + skim cohort manifest** — IDs, donor species, per-source
      depth. Trace from `sample_metatada.csv` / `bzea-sample-donor-metadata` +
      the BRB and skim sample maps. Blocks `03-source-method-comparison.qmd`.
- [ ] Which MolBreeding target-seq samples for the calibration cross-check.
- [ ] B1 anthocyanin phenotype table + NIL panel. Blocks
      `04-b1-mapping-benchmark.qmd`.

## Migration from zealtiger (plan B3)

The sanity-paint sweep inputs currently read from the zealtiger working repo
(`results/sim_calibration/coverage_sweep_members.csv`,
`data/rtiger_50K/`, `results/sim_calibration/brbseq_ks_wideseq/counts/`). Only
the **validated, paper-bound** subsets migrate here as fixtures; the `agent/`
scratch tree is not copied.
