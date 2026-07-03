# `sim/` — simcross simulation + truth generation

Generates simulated BZea NILs with known ancestry truth and degrades them to
each sequencing source's depth/error regime. These sims are the paper's
**calibration set** (fit `r`, `err`, emission means on truth) and **validation
set** (held-out accuracy/boundary/ROC metrics) — see
`analysis/02-simulation-calibration.qmd`.

- `simulate_nils.R` — generator. Uses `kbroman/simcross` + the map/fitting
  primitives in the `nilHMM` package (`load_map`, `expected_fragment_dist`,
  `calibrate_r`, `cm_to_mb`, `fit_design_gamma`). **Scaffold** until the BZea
  simcross design is fixed (plan B5).

Seeds are pinned (`DEFAULT_SEED`) so runs are reproducible (B4). Outputs land in
`results/sim/` (gitignored).

## Open before this runs (plan B5)

- generations: BC2S2 bulked skim vs BC2S3
- n per condition
- cM map + interference model to match BZea
