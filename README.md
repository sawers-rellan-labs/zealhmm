<p align="center">
  <img src="zeal.png" alt="ZEAL — Zea Exotic Allele Library" width="200">
</p>

<h1 align="center">ZEAL — Zea Exotic Allele Library</h1>

# zealhmm — simulation-calibrated HMM ancestry calling for NILs

Clean, reproducible analysis repo built on the [`nilHMM`](https://github.com/sawers-rellan-labs/nilhmm)
ancestry-caller package. It calibrates and validates caller parameters against Broman
`simcross` simulations, reproduces the TeoNAM QTL analysis (Chen et al. 2019) under
simulated low coverage, and applies the calibrated callers to the BZea/ZEAL
teosinte-introgression population.

📖 **Docs site:** <https://sawers-rellan-labs.github.io/zealhmm/>
📦 **Package:** <https://sawers-rellan-labs.github.io/nilhmm/>

## Papers this repo serves

The notebooks are organized on the docs site by **dataset** (TeoNAM, ZEAL/BZea) and,
within each, by the **paper** the work supports:

1. **nilHMM R package paper** (Oxford Bioinformatics) — **priority; written here.** A
   methods paper: the read-based genotype model, the shared HMM engine, simulation
   calibration, the TeoNAM reproduction, and a cross-modality application. See below.
2. **ZEAL population paper** — *second; the ZEAL trait analyses here are its
   exploratory work.* A reflection of the Chen 2019 TeoNAM analysis applied to the
   B73 × teosinte BC2S3 population (DTA / DTS / PH, and the binary stem traits
   StPi / StPu — modeled as empirical logits, no spatial correction).
3. **Inv4m inversion paper** — *already written, in review.*  This repo only
   supplies it updated genotypes (`analysis/zeal-inv4m-rtiger-genotype.qmd`); it is not
   authored here. https://github.com/sawers-rellan-labs/inv4m

The **nilHMM manuscript lives in its own repo** (RILAB LaTeX → Overleaf); this repo
produces the figures, tables, and rendered notes it consumes.

## The nilHMM paper — structure

| Part | Content | Notebook(s) | Status |
|------|---------|-------------|--------|
| §1 Motivation + foundations | Low coverage; single-read genotype calling; missing-data model (Poisson + floor) | *(to write — math currently in `R/simulate.R`)* | ⚠️ not started |
| §2 Shared HMM engine | One duration-aware 3-state engine; each caller (nNIL, RTIGER, binhmm, ATLAS) a different setup | `callers-and-methods.qmd`; engine fidelity `rtiger_vs_nilhmm_reproduction.qmd`; 2-source × caller demo `nilhmm_sanity_check_paint.qmd` | ready |
| §3 Calibration **+ TeoNAM reproduction** | simcross → params (two-stage search); **the comparison of GWAS profiles under simulated coverage is the results center** | `simulation-calibration.qmd`; `teonam-qtl-recovery{,-118k}.qmd` ⭐, `teonam-qtl-recovery-mlm-118k.qmd`, `teonam-qtl-recovery-dta-mlm-118k.qmd`, + Fig-4C reproduction and map/thinning notebooks | core done; nNIL/LB-Impute calib open |
| §4 Application | Cross-modality comparison on the narrow cohort of lines with **both** skim + BRB-seq | qualitative: `nilhmm_sanity_check_paint.qmd` (available now); quantitative: `source-method-comparison.qmd` | ⏸ blocked on the ~400 paired-cohort manifest (plan B5) |

The chromosome-painting figure (`nilhmm_sanity_check_paint`) ships two renderings: a
pared-down **main figure** (one representative chromosome × the 6 source×caller tracks
across a handful of NILs) and the **full genome-wide painting** as a supplement.

## Layout

```
zealhmm/
  README.md              # this file
  DESCRIPTION            # deps + Remotes: nilhmm, simcross (reproducibility)
  _quarto.yml            # renders analysis/*.qmd -> docs/  (Pages)
  DATA.md                # data provenance (data/ itself is gitignored)
  R/                     # shared analysis helpers + the sim generator (tracked code)
    staging.R            #   mount detection, file lookup, common call schema
    metrics.R            #   concordance + truth-based accuracy metrics
    plotting.R           #   REF/HET/ALT palette + chromosome-painting plot
    simulate.R           #   simcross NIL generator + per-source degradation
  analysis/              # paper notes -> docs/*.html (grouped by dataset -> paper)
  docs/                  # rendered site (Pages, .nojekyll) — HTML committed

  # gitignored (bulk inputs / generated outputs — see DATA.md):
  #   data/     inputs (live on the mount, referenced from DATA.md)
  #   sim/      generated simcross simulations (generator is R/simulate.R)
  #   results/  derived tables/figures
```

The docs landing page (`index.qmd`) is the authoritative index of notebooks, grouped
dataset → paper.

## Reproduce

```r
# 1. install the pinned package + deps
remotes::install_deps(dependencies = TRUE)   # honors DESCRIPTION Remotes:
# 2. (once) activate the styler pre-commit hook
#    git config core.hooksPath .githooks
# 3. render the analysis notes -> docs/
#    quarto render
```

`nilHMM` is a **dependency** (`library(nilHMM)`), not a vendored `load_all`.
Mount-dependent chunks are guarded by `R/staging.R::pick()` so the site still
renders off-mount. Seeds for `simcross` are pinned in `R/simulate.R`.

> Note: render notebooks **one at a time** (`quarto render analysis/<file>.qmd`).
> Passing several `.qmd` files to a single `quarto render` invocation can cross their
> metadata and fail at the pandoc step.

## Status

- **§2 (engine)** and the **TeoNAM reproduction** (§3 results center) are wired to the
  package and rendered.
- **§3 calibration** core is done; per-coverage **nNIL** and **LB-Impute** calibration
  remain open (`agent/HANDOVER.md`).
- **§1 (foundations)** must still be written.
- **§4 (paired cohort)** is scaffolded and **blocked** on the ~400 paired-cohort
  manifest (plan B5), tracked in [`DATA.md`](DATA.md).

## Relationship to other repos

- **`nilhmm`** — the R + Rcpp caller package this repo consumes.
- **`nilhmm-paper`** — the nilHMM manuscript (RILAB template, Bioinformatics
  structure). Its own standalone git repo; nested locally at `nilhmm-paper/` but
  gitignored here. This repo produces the figures it consumes.
- **`zealtiger`** — the exploratory lab notebook; stays as-is. Only validated,
  paper-bound notes are ported here; the `agent/` scratch is not.