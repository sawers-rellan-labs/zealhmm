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

`stage_sanity_check_paint.R` copies the 11-sample vary-skim cohort (B73 control + 10
NILs) into `data/`, organized by **(source × input type)** and consumed by
`analysis/nilhmm_sanity_check_paint.qmd`. Verified: the note reproduces the original
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

Re-stage with `Rscript stage_sanity_check_paint.R` **run with the sandbox disabled** —
both shares live under the `rsstu` automount, which the default sandbox cannot
see (it reports the files as missing).

## TeoNAM genotypes (`data/teonam/`, gitignored, ~236 MB)

Recoded TeoNAM SNP genotypes (Chen et al. 2019, *Genetics* 213:1065; `chen2019teonam`).
The genotypes cited in the paper at CyVerse `/iplant/home/shared/panzea/genotypes/GBS/TeosinteNAM`
are **no longer reachable** — the whole `panzea` collection returns "not found" for our
account and for anonymous CyVerse access (panzea.org still routes there and requires a
CyVerse login, so it is likely permission-gated to a group we're not in, not necessarily
deleted). Figshare (9820178 / 9250682) has phenotypes/supplement only; panzea.org's Cornell
archive has only the older BC2S3 population (Yang et al. 2019 PNAS, `yang2019genetic`).

**Source used instead:** the EasiGP repo (Tomura et al. 2025, bioRxiv `10.1101/2025.07.15.664852`,
`tomura2025ensemble`) vendors a recoded copy.

| | |
|-|-|
| Paper | Tomura et al. 2025, bioRxiv `10.1101/2025.07.15.664852` (`tomura2025ensemble`) |
| Repo | https://github.com/ShunichiroT/EasiGP |
| Commit | `259c242006bce415a8995be1833ef44399655718` (2026-06-16) |
| Path in repo | `Data/TeoNAM/` |
| `TeoNAM_genotype.zip` sha256 | `3a8b5adcdff6ed162032347ecaabe25fed5384ea8ed387f0b22fe492b2c6f918` |

Staged files (copied verbatim; `TeoNAM_genotype.csv` is the unzipped `TeoNAM_genotype.zip`):

```
data/teonam/
  TeoNAM_genotype.csv        2,434 lines x ~51,483 markers; 0/1/2 coding
                             (0 = W22 hom, 1 = het, 2 = teosinte hom — paper's coding)
  TeoNAM_phenotype.csv       2,434 lines; traits DTA, ASI only (2 of the paper's 22)
  marker_info.csv            51,544 markers; chromosome,name,start,end
                             names like S1_10045 = chr1 / AGPv2 position 10045 (start/end are 0)
  chrom.csv                  per-family chromosome spans
  gene_info.csv              flowering-gene markers (Dong 2012, Wisser 2019)
  W22TIL{01,03,11,14,25}_genotype.csv / _phenotype.csv   the 5 teosinte-parent subpopulations
  W22TIL01_subset_*          small subset used by EasiGP as a quick example
  TeoNAM_genotype.zip        original artifact (kept for integrity)
```

**Caveats (verify before analysis):**
- Numeric **dosage matrix only** — no raw HapMap / ref-alt alleles. Fine for ancestry-HMM
  (0/1/2 *is* the W22↔teosinte ancestry dosage); insufficient if raw allele calls or an
  AGPv2 VCF are needed (positions survive only in the marker names).
- **Line count 2,434 exceeds the paper's 1,257 RILs — resolved (2026-07-04):** the file
  is EasiGP's per-family files concatenated, and each family file already doubles every
  line. The two copies of a line differ **only in `2` vs `2.0` integer/decimal
  formatting** — numerically identical (0/1,197 disagree). Cleaned to `TeoNAM_genotype_clean.csv`
  (see characterization below): dedup → **1,237** unique lines, all cells decimal.
- "Adjusted example data" per EasiGP's `Data/README` — treat as a re-derived convenience
  copy, not the authoritative release. For the authoritative genotypes, request access to
  the CyVerse `panzea` collection from the authors (Qiuyue Chen / John Doebley).

**Genotype-matrix characterization (measured 2026-07-04, full-scan of the raw and
cleaned files).** Raw `TeoNAM_genotype.csv`: 2,434 rows × 51,482 marker columns
(cols 4–51485; first 3 = `ID, population, factor`), CRLF endings. Cleaned artifact
`TeoNAM_genotype_clean.csv`: **1,237 rows × 51,482 markers, all decimals, LF**
(dedup + `2`→`2.0` normalization; original preserved).

- **The 2,434→1,237 doubling is a formatting artifact, not replicates.** Each line
  appears twice (byte-identical ID); the two copies are numerically identical and
  differ only in integer vs decimal spelling (`2` vs `2.0`). Raw matrix had 328,590
  integer-formatted cells (0.26%) vs 33.4M decimal; the clean file has **0** integer
  cells. Dedup keeps one copy → 1,237 unique lines.
- **Imputed / complete within family — proven, not inferred.** For every family,
  each marker column is present in either *all* of the family's lines or *none*
  (0 partially-present columns) → within-family missingness **0%**; hard `0/1/2`
  calls only. This is the paper's per-subpopulation FSFHap-imputed output — not raw
  GBS, and not the cross-imputed/thinned 4,578-marker JLM matrix.
- **Block-sparse union (73% blank, 100% structural):** each line is blank at other
  families' markers. Per-family line/marker counts (cleaned, one row per line):

  | Family | nLines (clean) | markers (present per line) |
  |---|---|---|
  | W22TIL01 | 220 | 13,042 |
  | W22TIL03 | 267 | 16,076 |
  | W22TIL11 | 218 | 13,152 |
  | W22TIL14 | 222 | 11,375 |
  | W22TIL25 | 310 | 14,857 |
  | total / mean | **1,237** | 13,700 |

  Mean 13,700 markers ≈ Chen 2019's "average 13,733 high-quality SNPs per subpopulation."
  Clean grid check: `1,237 × 51,482 = 63,683,234 = 17,159,588 present + 46,523,646 empty`.
- **Genotype frequencies** among present calls (dedup-invariant): **0 = 76.65%,
  1 (het) = 8.15%, 2 (teo hom) = 15.20%** — matches the paper's "15% homozygous
  teosinte, 8% heterozygous."
- **Reconciliation with phenotype (1,257 canonical RILs):** inner-join on
  `(family, line)` = **1,237** — i.e. **20 phenotyped RILs have no genotype in this
  EasiGP copy** (scattered, not truncation; per family missing: TIL01 3, TIL03 3,
  TIL11 1, **TIL14 13**, TIL25 0). Genuine incompleteness of the convenience copy.
- **Implication:** use each family as a dense `nLines × family-markers` truth block
  (simulate reads per family at that family's markers); no missingness to model on
  the truth side. Downstream densification to a common grid = `interpolate_genotype()`
  (nilHMM); see `agent/teonam-qtl-recovery-plan.md` §4.

## TeoNAM marker sets (the thinning/densification hierarchy)

All derive from the same cleaned genotypes; each row is a distinct marker count used
downstream. Counts verified 2026-07-05.

**Naming system — `teonam_<panel>_<assembly>_<role>`.** `panel` = `map` (all our data
derives from the composite MAP panel); `assembly` = `v2` (marker-name coords) | `v5`
(lifted); `role` = `family_imputed` (the base FSFHap within-family imputed calls) |
`markers` (annotation) | `gwas` (GWAS base) | `jlm` (JLM base — thinned directly from `gwas`).
(`gwas_nr` is a **retired** label, not a pipeline role: it does **not** feed JLM; kept only
as the duplicate-cM record — see the table note below.)
**All map sets share the same base calls, which are FSFHap family-imputed** (Chen's map
construction — complete *within family*; verified 0% partial-missing). So `family_imputed` is
the base state, not a distinguishing step; the later suffixes describe downstream *selection*.
The `family_imputed` label is deliberate: FSFHap imputes **within each family**, leaving
cross-family gaps NA — distinct from the `interpolated` matrix, where our cross-family
step-interpolation (Tian 2011 / Chen densifier) fills those gaps to supply the complete
genotypes for `jlm` (and the interpolated GWAS matrix). Interpolation is a step, not a named panel.

| name | markers | assembly | role / provenance |
|---|---|---|---|
| `teonam_map_v2_family_imputed` | 51,482 | v2 | **FSFHap family-imputed** recoded dosages (Chen map construction; complete *within family* — verified 0% partial-missing); EasiGP `TeoNAM_genotype_clean.csv`; coords in marker names |
| `teonam_map_v5_markers` | 51,065 | v5 | annotation (lift v2→v5 + consensus cM); `marker_info_v5_cm.tsv` (`R/teonam_liftover.R`) |
| `teonam_map_v5_gwas` | 51,004 | v5 | **GWAS base** — per-marker, typed lines only; `stam_gwas_scan_family_imputed.csv` |
| `teonam_map_v5_gwas_nr` | 47,750 | v5 | non-redundant (one marker per unique cM). **Retained only as the record of the duplicate-cM markers (`agent/notes-redundant-markers.md`); NOT in any active path.** |
| `teonam_map_v5_jlm` | 6,049 | v5 | **JLM base** — FastIndep (**cM-distance** @0.1) directly on `_gwas` (51,004); genotypes step-interpolated; `data/teonam/tassel/geno.hmp.txt` (`agent/teonam_jlm_build.R`) |

Notes:
- **`teonam_map_v5_jlm` thinning is a cM-distance thin, not an r² LD prune.**
  `select_independent` (FastIndep, deterministic greedy) is fed the per-chr cM-distance matrix, so
  it enforces a hard 0.1 cM minimum spacing and does no LD decorrelation. The full pipeline is now
  reconstructed as a traceable script — `agent/teonam_jlm_build.R` (FastIndep on the 51K pool +
  step-interpolated genotypes). `agent/teonam_jlm_verify_source.R` reproduces the prior
  47,750-pool hapmap **0-mismatch** (6,059 markers × 1,237 taxa), confirming the method.
  `analysis/marker-thinning.qmd` (written against the old 47,750 pool → 6,059) needs
  re-rendering to the 51K-pool set (6,049).
- **6,049 (`teonam_map_v5_jlm`) vs `chen2019_JLM` 4,578** = map length (1781 vs 1540 cM,
  ×1.156) × 0.1-cM-bin saturation (×1.14) = ×1.32. Not a thinning-algorithm artifact.
  See `analysis/marker-thinning.qmd`.
- **`teonam_map_v5_gwas_nr`** dedups markers sharing a cM (unique-cM grid, 47,750). It is
  **no longer in any active path** — JLM now thins directly from the 51K `_gwas` pool
  (cM-distance @0.1 subsumes the dedup: tied markers are 0 cM apart → one survivor per cluster).
  Retained only as the record of the duplicate-cM markers (~3,307 on the 51,004 base, centromeric
  low-recombination); see `agent/notes-redundant-markers.md`.

### What Chen 2019 has vs what we managed (reproduction status, 2026-07-05)

We are reproducing Chen 2019 on the **51,544-marker composite MAP panel** (obtained via EasiGP;
verified map-set-only). It's a strong basis and the JLM anchors reproduce on it. Chen's separate
**118,838-SNP GWAS set** is not yet in hand (sibling of the map set — same raw 955,690 ZeaGBSv2.7
SNPs, different filter chain); we keep working on the 50K and will fold in the 118K once obtained
(outreach: `agent/outreach_teonam_118k.md`).

| pipeline stage | `chen2019_*` (published) | `teonam_*` (this work) | status |
|---|---|---|---|
| raw GBS calls | 955,690 SNPs (ZeaGBSv2.7, AGPv2, HapMap) | not held (EasiGP ships post-filter 0/1/2 dosages) | to obtain |
| **map panel** | `chen2019_map` 51,544 — MAF<5% + 64-bp thin → FSFHap → composite | `teonam_map_v2_family_imputed` 51,482 / `teonam_map_v5_markers` 51,065 | ✓ have |
| genetic map | R/qtl `est.map`, 1540 cM, AGPv2 order | consensus Marey cM (Ed Coe composite), 1781 cM | different map |
| coordinates | AGPv2 → AGPv4 (CrossMap) | AGPv2 → AGPv5 (liftover) | v5 vs v4 |
| **JLM markers** | `chen2019_JLM` 4,578 (0.1 cM thin of map set) | `teonam_map_v5_jlm` 6,049 (FastIndep on cM-dist @0.1, 51K pool) | ✓ same op, denser/longer map |
| JLM engine | SAS `PROC GLMSELECT`, marker nested-in-family, permutation threshold | TASSEL5 StepwiseOLS, nested-in-family, default threshold | perm threshold TODO |
| JLM result | STAM 5 QTL, DTA 19 QTL | recovers STAM QTL (more enter w/o perm cutoff) | reproducing |
| **GWAS markers** | `chen2019_GWAS` 118,838 (MAF>0.01, unthinned → FSFHap) | `teonam_map_v5_gwas` 51,004 (full set; also the JLM thinning pool) | 118K to obtain |
| GWAS model | MLM: Q (5 PCs) + K (IBS) | OLS: Family + marker, 1-df F | MLM upgrade TODO |
| GWAS threshold | P < 1e-5 (LOD 5) | LOD 5 | ✓ |

Consequences to keep closing (not blockers): the GWAS Manhattan on `teonam_map_v5_gwas` is
sparser than `chen2019_GWAS` (some gaps, esp. low-recombination regions) and runs hotter (fixed
OLS vs Q+K MLM → peak inflation; λ_GC analysis pending). JLM QTL recovery is the primary
empirical-truth anchor and is unaffected by the GWAS-set gap.

## Still to freeze (plan B5 — blocks B2.3 / B2.4)

- [ ] **~400 paired BrB + skim cohort manifest** — IDs, donor species, per-source
      depth. Trace from `sample_metatada.csv` / `bzea-sample-donor-metadata` +
      the BRB and skim sample maps. Blocks `source-method-comparison.qmd`.
- [ ] Which MolBreeding target-seq samples for the calibration cross-check.
- [ ] B1 anthocyanin phenotype table + NIL panel. Blocks
      `b1-mapping-benchmark.qmd`.

## Migration from zealtiger (plan B3)

The sanity-paint sweep inputs currently read from the zealtiger working repo
(`results/sim_calibration/coverage_sweep_members.csv`, `data/rtiger_50K/`,
`results/sim_calibration/brbseq_ks_wideseq/counts/`). Only the **validated,
paper-bound** subsets migrate onto the mount / local `data/`; the `agent/`
scratch tree is not copied.
