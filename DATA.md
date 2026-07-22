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
  markers_v2.csv            51,544 markers; chromosome,name,start,end
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
- **Assembly = AGPv2, documented first-party by Qiuyue Chen (not by EasiGP).** The marker-name
  positions (`S[chr]_[pos]`) are B73 AGPv2, confirmed by three sources:

  1. **Chen's own genotype-release README** (`data/teonam/9250682/Genotype_files_ReadMe.docx`,
     shared directly by Qiuyue Chen from her Google Drive) — the authoritative, first-party
     statement:
     > "The genetic map was estimated separately for each subpopulation using est.map function
     > in R/qtl, and the physical position is represented in marker name with **B73 V2 position
     > (S[chr]_[pos])**."
  2. **Chen et al. 2019** (Genetics 213:1065):
     > "The genotypic data were uploaded with **AGPv2 position in the marker name**." (Data Availability)
     > "The **B73 reference genome v2** was used to determine marker order..." (Methods)

     AGPv4 appears in Chen 2019 **only** for *reporting* QTL positions (CrossMap-lifted from
     AGPv2); the genotype marker-name coordinates are AGPv2.
  3. **EasiGP is silent on assembly.** The EasiGP repo and paper (Tomura et al. 2025, *Plant
     Genome*, doi:10.1002/tpg2.70138) only cite "adjusted example data from the TeoNAM dataset
     (Chen et al., 2019)" and never mention a genome build — so the version documentation is
     entirely Chen's.

  **Marker-name identity checks (2026-07-07)** tying EasiGP back to Chen's documented panel:
  - Chen's `TeoNAM.composite.genetic.map.51544snp.txt` (51,544 names) **== EasiGP `marker_info.csv`
    (51,544 names)**: identical name sets, 0 unique to either side.
  - EasiGP `marker_info.csv` is byte-identical to our `markers_v2.csv` (md5 `6b29c715…`).
  - Our analyzed genotype matrix `TeoNAM_genotype_clean.csv` (51,482 columns) is a clean subset
    of Chen's names (0 unrecognized columns; 62 Chen markers absent = the known 51,544→51,482 drop).
  - Independent physical check: per-chromosome max marker position fills the AGPv2 length and
    falls ~6–7 Mb short of v4/v5 (e.g. chr1 max 301,252,447 < AGPv2 301,354,135 « v4 307,041,717).

  (EasiGP repo clone kept at `agent/EasiGP/`. Chen's Drive release, incl. the README, composite
  map, phenotypes, and — per the README — the 118,838-SNP GWAS HapMap + raw `.h5`, is under
  `data/teonam/9250682/`.)

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
| `teonam_map_v5_markers` | 51,065 | v5 | annotation (lift v2→v5 + consensus cM); `map_v5_coe2008.tsv` (`R/teonam_liftover.R`) |
| `teonam_map_v5_gwas` | 51,004 | v5 | **GWAS base** — per-marker, typed lines only; `stam_gwas_scan_family_imputed.csv` |
| `teonam_map_v5_gwas_nr` | 47,750 | v5 | non-redundant (one marker per unique cM). **Retained only as the record of the duplicate-cM markers (`agent/notes-redundant-markers.md`); NOT in any active path.** |
| `teonam_map_v5_jlm` | 6,049 | v5 | **JLM base** — FastIndep (**cM-distance** @0.1) directly on `_gwas` (51,004); genotypes step-interpolated; `data/teonam/tassel/geno.hmp.txt` (`scripts/teonam_jlm_build.R`) |
| `teonam_gwas118k_v2_hapmap` | 118,838 | v2 name / v4 pos | **authentic Chen 2019 GWAS panel** (the real Fig 4C set) — TASSEL HapMap, single-char IUPAC, 1,257 RILs; from Qiuyue Chen's Drive: `data/teonam/9250682/W22TILXX_Chr1-10.impute_filter_MR0.2_MAF0.05.hmp.txt` |
| `teonam_gwas118k_v5_markers` | 118,514 | v5 | v2→v5 liftover roster (chr from marker NAME, not the v4 `chrom` col); `markers_v5_gwas118k.tsv` (`scripts/teonam_lift_gwas118k.R`). 32,120/32,120 markers shared with `teonam_map_v5_markers` lift identically. |
| `teonam_gwas118k_v5_dosage` | 118,838 | v5 | additive dosage matrix [markers×1,257 lines], 2.72% NA; dosage = A2-allele count (TASSEL polarization, arbitrary vs W22/teo — fine for the flip-invariant STAM F-test); `teonam_gwas118k_dosage.rds` (`scripts/teonam_gwas118k_dosage.R`) |
| `teonam_gwas118k_v5_scan` | 118,514 | v5 | **STAM GWAS scan** on the 118K panel — `STAM ~ Family + marker`, 1-df F (same model as `stam_gwas_scan_family_imputed.csv`); `stam_gwas_scan_118k.csv` (`scripts/teonam_stam_gwas118k.R`) |
| `teonam_gwas118k_v5_dosage_polar` | 118,838 | v5 | **the coverage-sweep TRUTH** — dosage polarized to W22↔teo (0=W22, 2=teo; 0 flips, 99.6% concordant with 51K), from the distributed imputed-nucleotide genotypes. These are the *real* per-SNP genotypes (ancestry is choppy, ~889 flips/chr — real per-SNP imputation het). The sweep Poisson-samples reads from these; the **caller's ancestry inference + back-projection does the smoothing** (not a pre-smoothed input), which is why the truth stays the authentic per-SNP panel. `teonam_gwas118k_dosage_polar.rds` (`scripts/teonam_gwas118k_polarize.R`) |
| `teonam_gwas118k_v5_dosage_fsfhap` | 118,514 | v5 | the **ancestry-imputed panel** (Fig-4C panel C + the source of `stam_gwas_scan_118k_complete_baseline.csv`) — the 118K block-smoothed with TASSEL5 **FSFHap** as Chen (bc, Phet=0.03125, Fillgaps, per family), parental (A/C) output → clean mosaic (~9 breakpoints/chr) step-interpolated to the union. `teonam_gwas118k_dosage_fsfhap.rds` (`teonam_gwas118k_fsfhap.R` + `teonam_gwas118k_truth_assemble.R`). **NOT the sweep truth** (superseded — the sweep uses the authentic per-SNP `_dosage_polar`; the ancestry smoothing belongs in the caller, not the truth). |
| `teonam_gwas118k_v5_cm` | 118,514 | v5 | native-cM grid (31,411 est.map-placed + 87,103 Marey-spline); `markers_v5_gwas118k_cm.tsv` (`scripts/teonam_gwas118k_cm_grid.R`) |
| `teonam_gwas118k_v5_thin01` | 8,969 | v5 | **0.1 cM ancestry-inference grid** (the exact maximum distance-*d* independent set on the 1-D cM axis via O(n) greedy — Jena et al. 2018, `data/main.pdf`; no O(n²) matrix). The sweep HMMs infer segments on this grid, then back-project to the full 118K. `markers_v5_gwas118k_cm_thin01.tsv` (`scripts/teonam_gwas118k_thin01.R`) |
| `teonam_gwas118k_v5_<caller>_sweep` | 118,514 | v5 | **coverage-degradation sweeps** — Poisson reads from the authentic per-SNP truth (`_dosage_polar`), ancestry inferred on the **0.1 cM thin grid** (`_thin01`), segments **back-projected onto the full 118K** for the STAM `~Family+marker` scan, at λ∈{0.1,0.2,0.5,1,5,10,20,∞}. **λ=∞ is each caller at perfect coverage** (deterministic decisive counts) — its own ceiling, in-pipeline (no shared baseline; per-caller). Four callers: `teonam_{rtiger,nnil,lbimpute,control}_sweep_118k.R` → `results/sim/teonam/stam_gwas_<caller>_118k_sweep.csv` (+ `_sweep_manhattan.png` via `teonam_sweep_manhattan_png.R`). One replicate/cell; peak heights not monotone in coverage (caller smoothing × real per-SNP noise). |

Notes:
- **`teonam_gwas118k_*` is a NEW panel, a SIBLING of `map` — not derived from it.** Both were
  filtered from the same raw 955,690 ZeaGBSv2.7 SNPs via different chains (map: MAF<5% + 64-bp thin
  + FSFHap; GWAS: MAF>1%, MR<0.2, no thin, separate imputation). Name overlap with the 51K map panel
  is only 32,171 / 118,838 (86,667 GWAS-only, 19,373 map-only). This is the file whose absence
  made the earlier STAM Manhattan a non-faithful approximation — now resolved.
- **GWAS reproduction (2026-07-07):** 118K scan gives max −log10P 19.8 (vs 14.3 for the 51K approx)
  and 8,392 genome-wide-significant hits concentrated on chr1 (4,901; tb1 region, v5 ~270–282 Mb),
  chr3 (3,428), chr7 (58) — matching Chen 2019's STAM QTL architecture (STAM1.x/tb1, STAM3.x, STAM7.1).
  Decoding validated: all 32,158 shared markers correlate positively with the 51K coding (mean |r| 0.996).
- **Also in Chen's Drive release** (`9250682/`, per `Genotype_files_ReadMe.docx`), not yet all pulled:
  raw `W22TILXX_NoDepth.h5` (GBS in H5), `TeoNAM.filter0.1cM.snpInfor.txt` (4,578-SNP JLM subset),
  `TeoNAM_phenotype_07022018.txt` (4,578 SNPs + 22 traits for SAS JLM).
- **`teonam_map_v5_jlm` thinning is a cM-distance thin, not an r² LD prune.**
  `select_independent` (FastIndep, deterministic greedy) is fed the per-chr cM-distance matrix, so
  it enforces a hard 0.1 cM minimum spacing and does no LD decorrelation. The full pipeline is now
  reconstructed as a traceable script — `scripts/teonam_jlm_build.R` (FastIndep on the 51K pool +
  step-interpolated genotypes). `scripts/teonam_jlm_verify_source.R` reproduces the prior
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

- [ ] **~400 paired BrB + skim cohort manifest** — IDs, donor taxon, per-source
      depth. Trace from `sample_metatada.csv` / `bzea-sample-donor-metadata` +
      the BRB and skim sample maps. Blocks `source-method-comparison.qmd`.
- [ ] Which MolBreeding target-seq samples for the calibration cross-check.

The B1 anthocyanin positive control (formerly the deferred `b1-mapping-benchmark.qmd`,
now removed) was realized instead by the stem-pigment GWAS
`analysis/zeal-qtl-recovery-stpi-mlm-snp50k-fdr.qmd` — no dedicated B1 phenotype
table / NIL panel is required.

## §1 missing-data summaries (`data/missing_data/`, gitignored)

Inputs for the §1 foundations notebooks `analysis/missing-data-floor-model.qmd`
(part 2) and `analysis/missing-data-model-comparison.qmd` (part 3). Each is a tiny
per-sample table of coverage (`lambda`) and observed missingness (`missing_obs`):

- `wgs_per_sample.tsv`      (1437 samples) — from `picard CollectWgsMetrics`
- `wideseq_per_sample.tsv`  (1434 samples) — 27 M teosinte-variant binned counts
- `snp50k_per_sample.tsv`   (1439 samples) — 50 K panel per-site allelic counts

**Source tables (native HPC paths, `/rsstu` — the mount `/Volumes/rsstu/...` is the
same data but unstable):**

- `/rsstu/users/r/rrellan/BZea/bzeaseq/WGSmetrics_summary.tsv` (296 KB)
- `/rsstu/users/r/rrellan/BZea/bzeaseq/ancestry/all_samples_bin_genotypes.tsv` (309 MB)
- `/rsstu/users/r/rrellan/BZea/bzeaseq/50K/allelic_counts50K.tsv` (2.4 GB)

**How they were built (NCSU Hazel, LSF).** The heavy reduction (esp. the 2.4 GB 50K
read) runs on a compute node; only tiny summaries are pulled back. Scheduler is LSF.

1. Code (tracked): `scripts/make_missing_data_summaries.R` (reduces each source to a
   per-sample summary with `data.table::fread`, reading only the needed columns) and
   its submit wrapper `scripts/make_missing_data_summaries.lsf`.
2. R environment: the full `bzeaseq.yml` conda spec is no longer dependency-solvable,
   so a **minimal** env (`r-base` + `r-data.table`) is built on the **login node**
   (compute nodes have no internet) in fast scratch (ephemeral, ~1 month) by
   `scripts/build_missing_data_env.sh`:
   ```
   bash scripts/build_missing_data_env.sh   # -> $RDT_ENV (default /share/<grp>/<user>/conda/env/rdt)
   ```
3. Submit the reduction (native `/rsstu` reads, priority queue `sara`):
   ```
   cd /share/maize/frodrig4/zealhmm_jobs/missing_data
   # (make_missing_data_summaries.{R,lsf} staged here)
   bsub < make_missing_data_summaries.lsf     # -> out/{wgs,wideseq,snp50k}_per_sample.tsv
   ```
   The `.lsf` activates the `rdt` env, sets `BZEASEQ_DIR=/rsstu/users/r/rrellan/BZea/bzeaseq`
   and `OUT_DIR=$PWD/out`, and runs the R script (~20 s wall, incl. the 2.4 GB read).
4. Copy the three summaries back into `data/missing_data/` here. They are gitignored
   (local only); the notebooks read them and the executed render is shipped via the
   committed `docs/` HTML (Pages does not re-execute).

To reproduce anywhere with the mount instead of the HPC, run the same R script with
`BZEASEQ_DIR=/Volumes/rsstu/users/r/rrellan/BZea/bzeaseq`.

## ZEAL Inv4m tagging-SNP genotype comparison (`data/zeal/`, gitignored)

Inputs for `analysis/zeal-inv4m-rtiger-genotype.qmd` (genotype at PZE04175660223, chr4:181,637,780):

- `rtiger_200K_calls_introfinder.rds` — **200K (previous)** RTIGER segment calls, one table per NIL
  (`V1` chr, `V2` start, `V3` end, `V4` state ∈ {`B73`, `Introgression`}). Staged from
  `results_list_new_name.rds` in Nirwan Tandukar's app repo
  `github.com/nirwan1265/BZea_Introgression_Finder` (the calls behind the inversion paper's
  `inv4m_introgression` label). Names carry a `.B` bulk suffix → strip to match pedigree.
  **Provenance (Nirwan DM 2025-08-15; a lineage SEPARATE from this repo's BzeaSeq/wideseq work,
  and not reproducible here):** Nirwan ran GATK HaplotypeCaller → joint genotyping on the ZEAL BAMs
  (`/rsstu/.../DOE_CAREER/BZea/joint_genotype/all_samps/9_final_samples/more_filtered/BZea.vcf.gz`,
  raw, ~98% missing) → MAF ≥ 0.05 + quality filter + rename → **Beagle imputation by Asher (P.
  Balint-Kurti's lab)** → `BZea_MAF_0.05_qfiltered_newnames_imputed.vcf.gz` (196,120 SNPs) → RTIGER.
  It is **excluded from the release** (legacy, Beagle-imputed, no pipeline on hand; regenerate with
  recalibrated RTIGER on the SNP50K side if ever needed).
- `rtiger_50K_calls.csv` — **50K (current)** RTIGER 3-state segments (see `zeal_rtiger_mosaic.R`).
- `CLY25_ZEAL.csv`, `Bzea_metadata.csv` — spatially-corrected NC2025 phenotypes + donor metadata,
  staged from the inv4m Drive repo (`scripts/inversion_paper`).

Both call sets are RTIGER; they differ only in the reference variant set. See memory
`bzea-genotype-call-sets`.

## Migration from zealtiger (plan B3)

The sanity-paint sweep inputs currently read from the zealtiger working repo
(`results/sim_calibration/coverage_sweep_members.csv`, `data/rtiger_50K/`,
`results/sim_calibration/brbseq_ks_wideseq/counts/`). Only the **validated,
paper-bound** subsets migrate onto the mount / local `data/`; the `agent/`
scratch tree is not copied.

## Shareable release bundle (`release/`, gitignored → CyVerse)

`scripts/zeal_export_release.R` assembles the citable ZEAL genotype release into
`release/bzea_genotypes/` (bulk, gitignored; destined for the CyVerse Data Store with a
DOI — Phase E of `SHAREABLE_DATASET_PLAN.md`). Fully reproducible from the tracked script
+ `scripts/release_README.md` template; ~156 MB, 30 files.

All five objects share **one common callable panel** (intersection of every object's callable
lines ≈ 1,395; the 8 near-empty libraries RTIGER's per-chromosome QC can't call are dropped) and
ship in three encodings: PLINK `.bed/.bim/.fam` (canonical binary 012), tidy `_012.tsv.gz`, and
the native `.rds`. Four are **ancestry mosaics** (`rtiger/nnil/binhmm/lbimpute`); the fifth,
**`hwe_post_gt`**, is the **genotype** — the real `bcftools mpileup | call -mv` (HWE-prior MAP)
calls from `bzea_50K_cohort.vcf.gz`, extracted for the panel (values = the bcftools calls, not a
reconstruction). **No single-sample GL genotypes are shared.** Shared `markers/snp50k_markers.tsv`
+ `lines/snp50k_lines.tsv`; `MANIFEST.tsv` carries sha256 per file.
The legacy **200K** RTIGER introgression set is **not** in the release (regenerate with
recalibrated RTIGER if ever needed).

The cohort VCF is staged from `/Volumes/rsstu/.../bzeaseq/50K/results/joint/bzea_50K_cohort.vcf.gz`
into `data/zeal/` (gitignored) and consumed by `scripts/zeal_hwe_post_gt.R`. The README states the

**Teosinte reference panel on the SNP50K frame (for the ZEAL count-from-parents sim).** From the
same `50K/results/joint/`: `bzea_50K_ref_panel.vcf.gz` (**217 teosinte individuals**, Grzybowski
et al. 2023 resequencing panel, projected to the SNP50K sites) is the reference/donor pool the sim
resamples WITH replacement as taxon-matched stand-ins; `bzea_50K_cohort_ref.vcf.gz` (1,656 samples =
cohort + reference) and `bzea_50K_cohort_ref_metadata.csv` classify every sample (`is_reference`,
`maizegdb_prefix` = Zd/Zh/Zl/Zn/Zv/Zx/Zm, `donor_accession`, `donor_taxa`, `NIL_pedigree`). Staged
locally (mount unstable) into `data/zeal/reference_panel/` (gitignored) with a derived
`reference_ids_by_taxon.csv` (217 rows: Zv 88, Zx 79, Zn 14, Zl 9, Zh 9, Zd 5, Zm 1, NA 12). NOTE
(snp50k_genotype_identifiability.qmd §2.6): these are **taxon-level** panel common variants
(taxon, never "species" — the taxa span ranks; TERMINOLOGY.md), not the individual donor
accessions' own genotypes — hence the stand-in resampling. Skim coverage/missing
regime input: `data/missing_data/snp50k_per_sample.tsv`.

The README states the
wall (the genotype↔ancestry-mosaic distinction — `nilhmm/design/TERMINOLOGY.md`; ZEAL object names
in `TERMINOLOGY.md`): a `_mosaic` PLINK/012 file encodes **ancestry** dosage on the SNP's
ref/alt alleles for tooling compatibility — it is not a genotype and does not report the true
allele at invariant sites; only `hwe_post` does.

## RTIGER C++/Rcpp-vs-Julia benchmark inputs (`data/rtiger_shared3_input/`, `data/bench_ref/`)

Backs the optimization supplement (`nilhmm-paper/supplement_optimization.tex` §1) and
its companion `analysis/rtiger-optimization-benchmark.qmd`. Gitignored; staged from
`~/repos/rtiger-fork-assets/` (safety copy of the `faustovrz/RTIGER@optimize-julia-core`
fork workspace). The **paper reports the operating rigidity r=250**; the r=2 pair
(`scripts/bench_rtiger_cpp_vs_julia.R`, `data/bench_ref/`) is the earlier companion.

- `data/rtiger_shared3_input/sample{BN,Z,AU}.txt` — the **109,703-locus shared-3
  panel**: real *Arabidopsis* Col-0 x Ler allele counts (BioProject PRJNA720439 /
  SRA SRP313911), the loci covered in all three of the BN/Z/AU samples, one common
  grid (built by the fork's `27_build_shared_panel.R`). RTIGER count format:
  `chr pos refAllele refCount altAllele altCount`. The C++ core's benchmark input;
  thinned by odd index to 6,857 / 13,713 / 27,426 / 54,852 / 109,703 markers.
- `data/bench_ref/` — the **upstream-original Julia baseline**, NOT re-run (the
  110k convergence fit is ~3.4 h upstream). Fixed reference, from the fork:
  - `orig_conv/log_orig_conv_L{0..4}.log` — per-EM-iteration `delta`/`elapsed`.
  - `orig_conv/conv_orig_L{0..4}.txt` — converged params + Viterbi paths (the
    equivalence oracle the C++ paths are compared against).
  - `30_panel_sweep_orig.txt`, `33_panel_convergence.txt` — original per-iteration
    throughput and convergence iters/runtime per size.
  Provenance/method: `~/repos/rtiger-fork-assets/.../docs/optimization.md`.
- **r=250 (operating-point) benchmark** (the paper's §2). Reproduce the whole dataset
  + both figures end-to-end with the tracked driver `scripts/rtiger_r250_run_all.sh`
  (orig-Julia conv + materialize thinned panels + orig/C++ peak RSS + C++ bench +
  plot, in order; ~1.5 h, the full-panel orig fit being the pole). Component scripts:
  `scripts/bench_rtiger_r250.R`
  (C++ core + equivalence/throughput, all five sizes to convergence),
  `scripts/rtiger_julia_conv_worker.jl` (original Julia at r=250, all five sizes ->
  `results/bench/orig_conv_r250/`; the full 109,703 panel is ~43 min). Peak RSS:
  `scripts/materialize_thinned_panel.R` writes `data/rtiger_shared3_input/thin_L{0..4}/`
  so each measured process loads ONLY its marker set (clean per-level memory, no
  full-panel floor), then `scripts/{bench_rtiger_cpp_memory_markers.R (RTIGER_RIG=250),
  rtiger_julia_mem_worker.jl}` read those. `scripts/rtiger_r250_plot.R` (figures).
  Outputs: `results/bench/rtiger_{equiv,marker_scaling,conv_trace}_r250.csv`,
  `rtiger_{memory_markers,julia_memory_markers}_r250.csv`; figures
  `nilhmm-paper/figures/rtiger_{equivalence,marker_scaling}_r250.png`. Everything
  measured (no projection); equivalence 0 mismatches through the full panel; the
  original needs chains >= 2r = 500 markers (holds at every size).

## nNIL equivalence inputs (`data/zhong2025/`, `data/nnil_equiv/`)

Backs the optimization supplement Sec. 2 (nNIL reproduction) and `scripts/nnil_equiv/`.
Gitignored. Zhong et al. 2025, The Plant Journal (doi:10.1111/tpj.70228), nested
near-isogenic lines; Holland's HMM caller.

- `data/zhong2025/` -- the Wiley supplement (Files S1-S22; DOI above), PLUS the two
  >100 MB files not in the Wiley zip, sourced from the authors' Box folder:
  - `File_S01.nNIL_raw_SNPs_bgi_id_miss20.txt` -- raw GBS SNP calls, HapMap-numeric
    (901 samples x 93,540 markers, B73 v4), `{0=minor-homoz, 0.5=het, 1=major-homoz,
    NA}`. The caller INPUT.
  - `File_S18.nNIL_gbs_HMM_introgressionCalls_full_set.csv` -- Holland's PUBLISHED
    per-marker calls `{0,1,2}` (888 lines x 64,025 markers). The reproduction ORACLE.
  - Files S3/S4 (grid search), S9 (bgi->line translation), S11 (Holland's
    `call_intros`, hmmlearn), S16 (optimal params).
- `data/nnil_equiv/` -- the shared faithful input built by `scripts/nnil_equiv/01_build_input.py`
  (File S1 put through Holland's own `File_S10` filter, reproducing his filtered
  genotype table exactly: 888 lines x 64,025 markers, 0 cell differences; recoded
  `{0,1,2,3}`) and its compact PLINK `.bed` (`make_bed.py`, round-trip-asserted so
  `0=REF`). Consumed by `02_holland_calls.py` (hmmlearn) and `03_nilhmm_calls.R`
  (BEDMatrix stream); `04_compare.R` for the full-panel equivalence.
- Sweeps use a SPLIT-first design (like RTIGER): `materialize_thinned_bed.py` writes
  an odd-index-thinned `thin_L{0..6}/geno.bed` + `markers.csv` per density level, and
  every sweep worker reads ONLY its level's split (no in-script thinning). `06_sweep.R`
  (time + peak RSS -> `nnil_marker_scaling.png`; auto-materializes the split if absent)
  and `07_equiv_sweep.R` (`07_holland_level.py` per-size Holland calls vs nilHMM
  position-by-position -> `results/bench/nnil_equiv_sweep.csv`, supplement Table 5).
- Python: existing `~/anaconda3/envs/nilhmm` env (hmmlearn, bed-reader). R: BEDMatrix, genio.

## nNIL calibration foil (`data/nnil_foil/`, `scripts/nnil_foil/`)

Backs `analysis/nnil-calibration-foil.qmd` (the recombination rate `r` is a benign
hyperparameter; the emission `nir` carries the calibration signal). All markers on
maize **v5** (`S<chr>_<pos_v5>`); the v4-keyed `data/nnil_equiv/geno.bed` is joined
via the crosswalk and relabelled. Reproduction chain (`scripts/nnil_foil/`):
- `01_liftover_cm.R` -- lift the 64,025 nNIL markers AGPv4->v5 (vendored chain
  `data/ref/chain_files/B73_RefGen_v4_to_...NAM-5.0.chain`, `~/bio/bin/liftOver` via
  `lift_unique()`), keep 1:1 same-chr (63,904), interpolate cM from the TeoNAM native
  v5 map (`data/teonam/markers_v5_gwas118k_cm.tsv`) -> `markers_v5.tsv` (crosswalk).
- `02_chip_truth.py` -- reproduces Holland's chip calls by porting File S14 + calling
  File S11 (chip xlsx `File_S02`, V3->V4 beds `File_S12/S13` from `agent/nNIL/`);
  24 both-platform NILs x 11,310 shared v5 markers -> `chip_truth_projected.csv`.
- `03_chip_calibrate.R` -- sweep nnil `rrate` vs the chip calls (Holland's emission) ->
  `chip_rrate_sweep.csv`, `chip_calib.json` (mismatch flat over [1e-6, 8.2e-3]).
- `04_sim_calibrate.R` -- BC5S2 simcross on the 63,904 v5 markers + native cM,
  n=1500, degraded FORWARD through Holland's `emission_gt` at the measured biological
  `nir=0.594`; calibrate `rrate` on 300 vs dense truth -> `sim_rrate_sweep.csv`,
  `sim_calib.json`, `sim_truth_segments.csv` (single-locus BC5S2 check PASS).
- `05_foil_figure.R` -> `nilhmm-paper/figures/nnil_calibration_foil.png`.
- `06_fragsize_figure.R` -> `fragsize_ecdf.csv` + `agent/nnil_foil_fragsize.png`.
- `07_holland_sensitivity.R` -- variance decomposition of Holland's File S04 grid
  (nir 97.7%, r 0.03%) -> `holland_param_sensitivity.csv`.
- `08_nir_sweep.R` -- nir sweep at the map `r` (3 metrics) -> `nir_sweep.csv`,
  `nir_fragsize_ecdf.csv` (chip-supervised nir optimum ~0.9 vs founder f0 ~0.59).
- `09_sim_nir_sweep.R` -- sim-only nir sweep: data generated at `nir=0.594`, caller
  scored vs the simcross truth -> `sim_nir_sweep.csv`; recovers `nir* = 0.594` (the
  excess to 0.9 on the real chip-supervised fit is the GBS data-quality effect).
- Biological nir estimate: `agent/nnil_foil_estimate_nir.py` (chip NAM-founder donors).
- Env: `~/anaconda3/envs/nilhmm` (+ `openpyxl` for the chip xlsx). R: nilHMM, BEDMatrix,
  data.table, GenomicRanges/IRanges/rtracklayer, ggplot2, patchwork.
