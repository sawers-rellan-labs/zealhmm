# zealhmm — simulation-calibrated HMM ancestry calling for NILs

Clean, reproducible analysis repo backing the **nilHMM methods paper** (Oxford
Bioinformatics). It installs the [`nilHMM`](https://github.com/sawers-rellan-labs/nilhmm)
package, runs the ancestry-caller comparison across sequencing modalities on BZea
NILs, calibrates and validates caller parameters against Broman `simcross`
simulations, and renders the paper's figures and analysis notes.

📖 **Docs site:** <https://sawers-rellan-labs.github.io/zealhmm/>
📦 **Package:** <https://sawers-rellan-labs.github.io/nilhmm/>

## Central thesis

Ancestry-caller parameters (rigidity `r`, error `err`, emission means) can be
**calibrated from cheap `simcross` simulations** rather than from scarce
high-coverage truth samples — and the calibrated parameters transfer to real
data across the Skim (low-cov WGS) and BrB (3′ RNA-seq) modalities.

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
  analysis/              # paper notes -> docs/*.html
    callers-and-methods.qmd       # §2.1-2.2 engine + 4 callers  (ready)
    simulation-calibration.qmd    # §2.4-2.5, §3.1-3.2 CORE       (needs B5 sim design)
    source-method-comparison.qmd  # §3.3 400-cohort concordance   (needs B5 manifest)
    b1-mapping-benchmark.qmd      # B1 mapping                    (deferred)
    nilhmm_sanity_check_paint.qmd # supplementary QC painting     (ported)
  docs/                  # rendered site (Pages, .nojekyll) — HTML committed

  # gitignored (bulk inputs / generated outputs — see DATA.md):
  #   data/     inputs (live on the mount, referenced from DATA.md)
  #   sim/      generated simcross simulations (generator is R/simulate.R)
  #   results/  derived tables/figures
```

The **manuscript lives in its own repo** (RILAB LaTeX → Overleaf); this repo
produces the figures, tables, and rendered notes. Section → note → figure map:

| Manuscript section | Note | Emits |
|--------------------|------|-------|
| §2.1–2.2 engine + callers | `callers-and-methods.qmd` | Fig F1 |
| §2.4–2.5, §3.1–3.2 calibration + benchmark | `simulation-calibration.qmd` | Fig F2, F3 |
| §3.3 cross-modality concordance | `source-method-comparison.qmd` | Fig F4 |
| future application | `b1-mapping-benchmark.qmd` | (deferred) |
| supplementary QC | `nilhmm_sanity_check_paint.qmd` | Supp. painting |

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

## Status

The methods note (01) and the supplementary sanity-check paint are wired to the
package. The
core calibration (02), the 400-cohort comparison (03), and B1 mapping (04) are
scaffolded and **blocked on open inputs** tracked in [`DATA.md`](DATA.md) and
plan B5.

## Relationship to other repos

- **`nilhmm`** — the R + Rcpp caller package this repo consumes.
- **`nilhmm-paper`** — the manuscript (RILAB bioRxiv template, Bioinformatics
  structure). Its own standalone git repo; nested locally at `nilhmm-paper/` but
  gitignored here. This repo produces the figures it consumes.
- **`zealtiger`** — the exploratory lab notebook; stays as-is. Only validated,
  paper-bound notes are ported here (plan B3); the `agent/` scratch is not.
