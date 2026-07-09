# Shareable BZea genotype dataset & genotype browser — goals and plan

**Status:** proposal (2026-07-09) · **Owner:** GEMMA lab (Fausto R.) · **Priority:** this is the
current top priority for wrapping up `zealhmm`.

This document states the goals for turning the ZEAL/BZea ancestry work into a **shareable,
citable genotype dataset** plus a **Shiny genotype browser**, and lays out an implementation
plan. It also records the terminology / scripts / provenance standardization that should land
alongside it, and the still-open tasks in the repo (secondary to the dataset).

---

## 1. Goals

### G1 — A shareable, citable BZea genotype dataset
Package the BZea BC2S3 NIL ancestry calls into a documented, versioned release with a DOI
landing (CyVerse Data Store). Three components:

1. **50K dataset, in its caller variations** — the SNP50K ancestry mosaics produced by each
   nilHMM/RTIGER caller, as **binary 0/1/2 ancestry-state matrices** (0 = B73 hom,
   1 = het, 2 = teosinte hom).
2. **250K dataset** — the previous (inv4m-paper) RTIGER introgression call set
   (2-state B73 / introgression), as staged from Nirwan's introgression finder.
3. **Expanded wideseq-union dataset** — genotype/ancestry at the **union of (Poisson-QC'd
   wideseq 27 M positions) ∪ (SNP50K positions)**, derived from the low-coverage read counts.

### G2 — A Shiny genotype browser that serves genomic fragments
A **standalone repository** hosting a Shiny app that serves **up to 1 Mb fragments** of
genotype on demand, exported as **VCF or HapMap**. The 50K data ships **inside that repo** as
the binary 0/1/2 ancestry mosaic, so the app is self-contained and deployable to shinyapps.io.

### G3 — Standardize terminology, scripts, data provenance, and the app
Fix the naming overloads (esp. "mosaic"), settle a canonical vocabulary (taxa codes, caller
names, marker sets, call sets), and align the scripts, `DATA.md`, and the app to it.

---

## 2. What exists today (inventory)

All data lives under `data/` (gitignored; provenance in `DATA.md`). Caller outputs share the
schema `list(markers, state[marker × line, 0/1/2], lines)` — this **is** the "binary 012
ancestry mosaic" G1 calls for.

| Asset | Path (`data/zeal/`) | Form | Notes |
|-------|---------------------|------|-------|
| SNP50K dosage base | `zeal_snp50k_dosage.rds` | `n_ref,n_alt,dosage,cov` matrices | 51,991 HQ teo-vs-B73 sites, B73 v5, polarized (0=B73,2=teo) |
| SNP50K markers | `markers_snp50k_v5.tsv`, `markers_snp50k_cm.tsv` | marker,chr,pos(,cM) | v5 bp + native cM |
| **50K mosaics (variations)** | `zeal_{rtiger,nnil,binhmm,lbimpute}_mosaic.rds` | 0/1/2 marker×line | one per caller — the G1.1 payload |
| per-SNP genotype (`_gt`) | `zeal_ml_gt.rds` | 0/1/2 hardcall | `call_gt(prior="flat")` via `zeal_gt.R`; retired `persnp` |
| **250K calls** | `rtiger_250K_calls_introfinder.rds` | per-line segment tbls (2-state) | 1,077 lines; from Nirwan's finder |
| 50K RTIGER segments | `rtiger_50K_calls.csv` | 3-state segments | Fausto's calls |
| samplesheet | `samplesheet_3way.csv` | line metadata | `skim_id/brbseq_id/in_snp50k/gwas_nil/...` |
| kinship | `zeal_K_vanraden_*.rds` | K per caller | GWAS |

**Builders (tracked, `scripts/`):** `zeal_snp50k_dosage.R` (dosage base) →
`zeal_caller_mosaic.R` (nnil/binhmm/lbimpute), `zeal_rtiger_mosaic.R`, `zeal_binhmm_mosaic.R`,
`zeal_lbimpute_mosaic.R` (per-caller 012 mosaics).

**Not yet built:** the wideseq-union expansion (G1.3), any VCF/HapMap exporter, any Shiny app.

---

## 3. Terminology to standardize (G3)

The single worst overload:

- **"mosaic"** is used two incompatible ways: (a) the generic per-caller ancestry-state matrix
  file `zeal_<caller>_mosaic.rds`, and (b) *specifically the RTIGER caller* in the GWAS layer
  (`GENO=mosaic` → reads `zeal_rtiger_mosaic.rds`; Panel C in composites). This makes
  "the mosaic MLM" ambiguous.

**Proposed canonical vocabulary** (to apply across scripts, `DATA.md`, notebooks, app):

| Concept | Canonical term | Retire / disambiguate |
|---------|----------------|-----------------------|
| The 0/1/2 ancestry-state matrix (any caller) | **ancestry mosaic** (generic) | keep as file suffix `_mosaic.rds` but never as a caller name |
| The RTIGER-derived mosaic used as GWAS Panel C | **`rtiger`** | replace `GENO=mosaic` alias with `GENO=rtiger`; drop the "mosaic == rtiger" special-case |
| Per-site genotype from counts (no HMM) | **`<method>_gt`** via `call_gt` | `persnp` retired → `ml_gt` (`zeal_gt.R`) |
| Caller set | `rtiger`, `nnil`, `binhmm`, `lbimpute` | one spelling everywhere |
| 50K marker set | **SNP50K** | `50K` / `snp50k` / `SNP 50K` → pick one (`SNP50K`) |
| 250K call set | **250K (inv4m / introfinder)** | always tag "previous"; it is 2-state |
| Dense teosinte panel | **wideseq (~27.6 M)** | `27M` / `wideseq_ref` → "wideseq" |
| Germplasm units | **taxa** (Zv/Zx/Zl/Zd/Zh); accessions `Zx.####` | never "species"; see `bzea-taxa-naming` |
| States | 0 = B73(REF), 1 = HET, 2 = teosinte(ALT) | fixed order + palette (`R/plotting.R`) |

Deliverable: a short **GLOSSARY** section in `DATA.md` (or `TERMINOLOGY.md`) that the app and
scripts both cite, plus a mechanical rename pass (`GENO=mosaic` → `rtiger`).

---

## 4. Implementation plan

### Phase A — Standardize (do first; unblocks clean naming downstream)
- A1. Write the GLOSSARY (§3) into `DATA.md` / `TERMINOLOGY.md`.
- A2. Rename the `GENO=mosaic` alias to `GENO=rtiger` in `zeal_mlm_taxon.R`, `zeal_composite.R`
  and any callers; keep `mosaic` only as the generic file suffix. Re-render affected notebooks.
- A3. Normalize spellings (`SNP50K`, `wideseq`, caller names) across `scripts/`, `DATA.md`,
  `analysis/`.

### Phase B — Package the 50K + 250K datasets (G1.1, G1.2) — ✅ DONE (2026-07-09)
Built by `scripts/zeal_export_release.R` → `release/bzea_genotypes/` (gitignored, ~156 MB, 30
files, fully reproducible; awaits CyVerse upload = Phase E). Per 50K object (`rtiger/nnil/binhmm/
lbimpute`_mosaic + `ml_gt`): PLINK `.bed/.bim/.fam` + tidy `_012.tsv.gz` + `.rds`; shared
`markers/` (marker,chr,pos,ref,alt,cM) + `lines/` tables; 250K as 2-state segment TSV + rds;
`README` (from tracked `scripts/release_README.md`) + `MANIFEST.tsv` (sha256). 49,002 SNP50K sites
× ~1,400 lines, B73 v5. README carries the mosaic≠genotype wall (PLINK 0/1/2 for a `_mosaic` is
ancestry dosage, not true alleles).

- B1. Define a **release schema**: per caller, a plain-text 0/1/2 matrix + a marker table
  (`chr,pos,ref,alt,marker,cM`) + a line table (`skim_id,pedigree,taxon,donor_accession`).
  Store canonically as **PLINK binary (`.bed/.bim/.fam`)** for the union set (space-efficient,
  standard) **plus** the RDS mosaics for R users.
- B2. `scripts/zeal_export_release.R` — read `zeal_*_mosaic.rds` → write, per caller:
  `bzea_snp50k_<caller>.{bed,bim,fam}` and a tidy `bzea_snp50k_<caller>_012.tsv.gz`.
- B3. Stage the **250K** set into the same release layout (2-state, clearly labeled "previous").
- B4. `README` + `MANIFEST.tsv` (sha256 per file) for the release bundle.

### Phase C — Wideseq-union expansion with Poisson QC (G1.3)
Definition of the QC (per user): **per-position KS goodness-of-fit to a Poisson**, then
**discard positions whose normalized deviation z > 3**.
- C1. `scripts/zeal_wideseq_counts.R` — assemble per-sample allele-count tables at the wideseq
  (~27.6 M) positions from the mount (GATK `CollectAllelicCounts` outputs), analogous to
  `zeal_snp50k_dosage.R`.
- C2. `scripts/zeal_wideseq_poisson_qc.R` — for each position, fit the across-sample total-depth
  distribution to a Poisson (λ = per-position mean depth), run a **KS test**, compute a
  **normalized z of the KS statistic** across positions, and **drop z > 3** (repeat/paralog/CNV
  and mapping-artifact sites). Emit a QC table (`position, lambda, ks_D, z, keep`).
- C3. **Union** the surviving wideseq positions with the SNP50K positions; build the polarized
  teosinte dosage + a chosen caller's 0/1/2 mosaic on the union grid.
- C4. Export the union set to PLINK + VCF/HapMap slices (feeds the app and the archive).
- *Open param:* exact "normalized z" definition (z of KS-D across positions vs a per-position
  p-value → z). Pin it in the script header before running at scale.

### Phase D — Shiny genotype browser (G2, new standalone repo)
Proposed repo: **`bzea-genotype-browser`** (separate git repo; deploy to shinyapps.io).
- D1. Commit the **SNP50K 0/1/2 ancestry mosaic** (per caller) into the app repo as the served
  asset — compact binary (PLINK `.bed` or an indexed RDS/`fst`), sized to fit shinyapps.io.
- D2. UI: pick chromosome + start (bp), window ≤ **1 Mb** (hard cap), caller/call-set selector
  (rtiger/nnil/binhmm/lbimpute; 50K vs 250K), sample subset.
- D3. Export: assemble the requested fragment on the fly → **VCF** and **HapMap** download.
  Reuse the exporter from B2/C4 as a shared function.
- D4. Provenance panel + link to the CyVerse release and this repo.

### Phase E — Archive & cite (G1)
- E1. Push the full release bundle (50K variations + 250K + wideseq-union PLINK/VCF + manifest +
  README) to **CyVerse Data Store**; request a DOI landing.
- E2. Cross-link: `zealhmm` `DATA.md` → CyVerse DOI → app repo. Add a "Data availability" blurb.

---

## 5. Sequencing & priority

1. **Phase A** (fast, unblocks naming) →
2. **Phase B** (50K + 250K packaging — highest value, all inputs already on disk) →
3. **Phase D** (app; can start against B's 50K export while C runs) →
4. **Phase C** (wideseq-union; the heaviest compute; needs the mount) →
5. **Phase E** (archive + DOI once C lands).

The 50K/250K release (B) + app (D) can ship **before** the wideseq expansion (C); C then lands
as a versioned update to the same CyVerse collection.

---

## 6. Open decisions / risks
- **Poisson-z definition** (C2) — pin exact normalization before scaling.
- **App payload size** — shinyapps.io free tier is ~1 GB image; the SNP50K 012 mosaic across
  callers must fit (PLINK `.bed` of ~52 K sites × ~1,400 lines is small; fine). The wideseq
  union is NOT shipped in-app — served from a fragment index or kept archive-only.
- **250K ↔ 50K reconciliation** in the app (2-state vs 3-state) — surface both, don't merge.
- **Which caller is "the" release mosaic** vs shipping all four — default: ship all four, mark
  `rtiger` as reference.

---

## 7. Still-open repo tasks (secondary to the dataset)

From a repo-wide sweep of `agent/*.md` plans, `README`/`DATA.md` status language, `.qmd`
callouts, and source-tree TODO grep. **None of these block the dataset/app work above.**

**Genuine blockers (missing inputs, not code):**
- **plan-B5 ~400 paired BrB+skim cohort manifest** — gates the two paper-core notebooks
  `simulation-calibration.qmd` and `source-method-comparison.qmd` (both scaffolded/rendered but
  the §3.3 concordance is not real without it). Also needs the MolBreeding target-seq sample list.
- **B1 anthocyanin phenotype + NIL panel** — gates `b1-mapping-benchmark.qmd` (marked DEFERRED).
- **Authentic TeoNAM raw GBS / some Drive artifacts** — CyVerse `panzea` is access-gated; raw
  ZeaGBSv2.7 + a few Drive files not all pulled. Blocked on author access (outreach open).

**Actionable analysis work (not blocked):**
- TeoNAM 118K sweep: **nNIL** calib script written-but-unrun, **LB-Impute** calib not built →
  then regenerate the 4-caller MLM composite (`agent/HANDOVER.md`).
- Low-coverage **skim sim recalibration** — port fitted π/k/λ̄/error, re-run `02_calibrate.R`,
  re-render `simulation-calibration.qmd` (`agent/HANDOVER-low-coverage-sim.md`).
- TeoNAM JLM **permutation threshold** + **MLM Q+K upgrade** + a **λ_GC inflation** notebook
  (`DATA.md` status table; deferred in the OLS notebooks).
- `marker-thinning.qmd` re-render to the current 51K JLM pool (minor).
- ZEAL: fetch **Sanchez/Holland accession passport** (lat/long/elev) and merge; DTA Phase 5
  (coverage sweep + noisiness) deferred; lollipop candidate-label refinement (minor).

**Cleanups:** stale inline TODO in `teonam_gwas118k_dosage.R`; prune retired-map artifacts;
confirm zealtiger→zealhmm migration (B3) + `callers-and-methods.qmd` bundled fixture; spot-check
nilHMM package Part A (Python retirement, CITATION.cff).

**Plan files that are effectively DONE (candidates to archive):** `teonam-map-handover.md`,
`teonam-v5-genetic-map-plan.md`, `handover-51k-interpolated-mlm.md`,
`teonam-control-sweep-plan.md`, `teonam-rtiger-degradation-plan.md`, `hazel-cyverse-download.md`,
`notes-redundant-markers.md`. **Still-live:** `HANDOVER.md`, `zeal_dta_repro_plan.md`,
`zealhmm_paper_plan.md`, `HANDOVER-low-coverage-sim.md`, `teonam-qtl-recovery-plan.md`,
`outreach_teonam_118k.md`.