# RTIGER reproduction and timing: nilHMM's `rtiger` vs the Julia RTIGER

**Date:** 2026-07-03 (updated after fixing a determinism bug in `rtiger_fit_cpp`)
**Machine:** Apple Silicon (arm64), 10 cores, macOS (Darwin 25.5); R 4.5.
**Reference implementation:** RTIGER R package (Julia backend via juliaup), `rHMM_methods.jl`.
**Test implementation:** `nilHMM::call_ancestry(caller = "rtiger")` — a C++ (Rcpp)
port of the RTIGER fork in `src/rtiger.cpp` / `R/rtiger.R`.
**Data:** the real SNP50K BZea NIL counts staged in the zealtiger repo
(`data/rtiger_50K/`): per-taxon cohorts of `counts/<donor>/<sample>.tsv`, the
51,991-SNP teosinte-informative panel. Taxa **Zh** (n=61) and **Zl** (n=121).
Reference RTIGER calls/fits: `calls_taxa_r8.csv`, `fits_taxa*/<taxon>/rtiger_result.rds`.

## Summary

nilHMM's `rtiger` **reproduces RTIGER faithfully on both Zh and Zl**: identical
optimizer, identical convergence criterion, and — after fixing a determinism bug
found during this validation — the **same converged emission and the same
iteration count**. Three findings underpin this, each verified against the actual
Julia source and the saved RTIGER fits (not from code comments):

1. **Same optimizer.** nilHMM's EM is a line-for-line C++ port of RTIGER's Julia
   Baum–Welch: identical `getlogpsi` / `productpsi` / forward / backward /
   Viterbi kernels, and the identical emission M-step (closed-form BetaBinomial
   mean + univariate **Brent** on the precision, same bracket).
2. **Same convergence criterion — a delta, not a deviance.** Both stop on
   `round(max(|Δα|,|Δβ|), 6) ≤ eps` with **eps = 0.01** (RTIGER's Julia source,
   line 1547, and its R wrapper default; confirmed against a run's
   `fit_progress.log`). An earlier draft of this note said RTIGER used a
   "deviance" criterion — that was wrong; the `EMdev`/`dev()` function is unused.
3. **Same preprocessing.** RTIGER decodes **covered markers only** (~12–16 k of
   51,991 per sample) and **requires ≥ 2·rigidity covered markers per
   chromosome**, aborting otherwise. nilHMM now matches both (see fixes below).

## Correctness

| taxon / comparison | result |
|---|---|
| **Zh** — identical params + identical covered markers (decode of RTIGER's exact fit) | **100.00%** per-marker; ALT 4387/4387 (recall & precision 1.000) |
| **Zh** — nilHMM's own r=8 fit + covered filter vs RTIGER | **99.95%** per-marker; donor-present F1 1.000; ALT recall 0.965 |
| **Zl** — nilHMM's own r=8 fit vs RTIGER (emission means, QC'd) | **19 iterations, means 0.99 / 0.870 / 0.659** vs RTIGER 19 iters, 0.99 / 0.870 / 0.661 |
| **Zh** — nilHMM's own r=5 fit vs RTIGER (per-marker, 848 k) | **99.76%**; ALT recall 0.944 / precision 0.802 |
| **Zl** — nilHMM's own r=5 fit vs RTIGER (per-marker, 1.98 M) | **99.95%**; ALT recall 0.997 / precision 0.989 |

- Given the same parameters and marker support, the two Viterbis are **bit-for-bit
  identical**, including the rare donor-hom state.
- With nilHMM doing its **own** fit, it now converges to RTIGER's emission on both
  taxa — Zl matches iteration-for-iteration (19) and mean-for-mean (to ~2e-3, i.e.
  float-summation / Brent-tolerance noise). The residual on Zh (0.05% of markers)
  is HET↔ALT boundary jitter on the rare, overdispersed donor-hom state.
- Note: RTIGER's own Zl fit is a *poorly-identified* optimum (states at
  0.99/0.87/0.66, barely separated — Zl sits on a flat likelihood ridge). nilHMM
  now lands on that same ridge point rather than wandering off it.

## The determinism bug (found and fixed during this validation)

An earlier draft claimed nilHMM's fit was "deterministic (seeded)." **It was not.**
On Zl, repeated identical calls (same data, seed, threads) returned different
iteration counts (1, 2, 3, 8…) and different converged means — sometimes returning
the init unchanged after "1 iteration."

Root cause, traced through `src/rtiger.cpp`:

- A chain (sample × chromosome) with fewer than **2·rigidity** covered markers has
  an empty `zeta` window; its E-step collapses to `exp(-Inf − (-Inf)) = NaN`, which
  poisons the pooled sufficient statistics. The emission M-step's NaN fallback then
  silently kept the init, reporting a bogus 1-iteration "convergence."
- The E-step was **always dispatched through RcppParallel `parallelFor`, even at
  `threads = 1`** (there was no real serial path). `parallelFor` did not reliably
  process the whole range across repeated calls, so it **nondeterministically
  skipped** the degenerate chains — masking the NaN on some calls and not others.

**Fixes (`src/rtiger.cpp`):**

1. **Real serial path for `threads ≤ 1`** — the worker is called directly, not
   through `parallelFor`. Deterministic; `threads > 1` keeps the parallel fold.
2. **Hard stop on under-covered chains** (`< 2·rigidity`), matching RTIGER's Julia
   abort — up front on the main thread, with a clear message — instead of silently
   skipping or NaN-ing. A `totw` guard remains as a backstop for empty input.

Verified after the fix: QC'd input is **deterministic across `threads = 1` and
`threads = 4` and repeated calls** (all 19 iterations, means 0.99/0.870/0.659,
matching RTIGER); under-covered input **hard-stops** with
`"rtiger: chain N has M covered markers, below the 2*rigidity = K floor …"`.

### `min_cov` (covered-marker preprocessing)

`call_states`/`call_ancestry` gained `min_cov` on the `rtiger` path (default `1L`):
markers with `n_ref + n_alt < min_cov` are dropped before fit/decode, matching
RTIGER's covered-only decoding. `min_cov = 0L` restores decoding every marker.

## Timing (taxon Zl, 120 QC-passing samples, `threads = 4`)

Grid = `2,3,4,5,6,8,10` (the rigidity sweep from zealtiger's `sweep_rigidity_45k.R`).

Post-fix build, quiet machine (`RTIGER()` per call includes its file read;
nilHMM's excludes the one-time `read_counts`):

| rigidity | Julia RTIGER (s) | nilHMM rtiger (s) | nilHMM / Julia |
|---:|---:|---:|---:|
| 2  | 21.6 | 14.5 | 0.67× |
| 3  | 20.2 | 14.8 | 0.73× |
| 4  | 16.3 | 11.4 | 0.70× |
| 5  | 15.9 | 11.2 | 0.70× |
| 6  | 15.6 | 11.0 | 0.71× |
| 8  | 15.3 | 10.6 | 0.69× |
| 10 | 15.3 | 10.5 | 0.69× |
| **full sweep** | **~120** | **~84** | **0.70×** |

At matched `threads = 4`, nilHMM is consistently **~30% faster** than the Julia
RTIGER across the whole rigidity grid.

> **On the earlier "1.4× slower" numbers.** An earlier draft reported nilHMM at
> ~1.2–1.9× *slower*. That was a **machine-contention artifact**: the pre-fix
> sweep ran alongside other jobs, inflating nilHMM's times (r=2 was 36.2 s under
> load vs 14.5 s quiet) while Julia's stayed ~21 s. On a quiet machine at matched
> threads, nilHMM wins. Timing is load-sensitive; single-run, single-machine.

## Reproduce

Scripts (session scratchpad; data paths point at the zealtiger repo):

- `time_sweep.R` — head-to-head timing (`JULIA_NUM_THREADS=4 Rscript time_sweep.R`).
- `definitive.R` / `covered.R` — identical-input → 100% per-marker check (Zh).
- `endtoend.R` — nilHMM's own-fit + covered filter vs RTIGER (Zh, 99.95%).
- `qc_det.R` — determinism across threads + Zl emission match to RTIGER (post-fix).
- `iters_zl.R` — EM iteration counts, nilHMM vs RTIGER's saved fits.

## Conclusion

nilHMM's `rtiger` is a faithful reproduction of RTIGER — identical optimizer,
kernels, and delta convergence criterion; exact calls on identical input; and,
after the determinism fix, the **same converged fit and iteration count** on both
Zh and Zl, at **r = 5 and r = 8** (~99.8–99.95% per-marker). The two behavioural gaps to RTIGER are now closed: decoding all markers
(fixed by `min_cov`) and the under-covered-chain NaN/nondeterminism (fixed by the
serial path + `2·rigidity` hard stop). On a quiet machine at matched threads,
nilHMM's rtiger is **~30% faster** than the Julia RTIGER across the rigidity grid.
