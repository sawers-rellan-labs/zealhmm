# Parameter calibration against simulated truth: a two-stage search.
#   1. LOG SWEEP (this file): score a coarse log-spaced grid to bracket the optimum.
#   2. golden-ratio refine (next): shrink the bracket on the continuous parameter.
#
# Both stages score with R/metrics.R and get their segmentations from
# nilHMM::caller_sweep() -- ONE shared fit, decodes fanned across the grid -- so a
# whole log sweep costs about one fit (see caller_sweep's fit-once semantics).

#' Log-spaced parameter grid
#'
#' @param lo,hi Grid endpoints (`lo < hi`, both > 0).
#' @param n Number of points.
#' @param integer If `TRUE`, round to unique integers >= 1 (for `rtiger` rigidity).
#' @return Numeric (or integer) vector, ascending.
#' @export
log_grid <- function(lo, hi, n = 12L, integer = FALSE) {
  if (lo <= 0 || hi <= 0 || lo >= hi) stop("log_grid(): need 0 < lo < hi")
  v <- 10^seq(log10(lo), log10(hi), length.out = n)
  if (integer) v <- sort(unique(pmax(1L, as.integer(round(v)))))
  v
}

#' Coarse log sweep: score a caller's parameter grid against simulated truth
#'
#' Runs [nilHMM::caller_sweep()] over `values` (one shared fit, parallel decodes)
#' and scores each value's segmentation against `truth` with the marker- and
#' fragment-level metrics. The returned curve brackets the optimum for the
#' golden-ratio refine.
#'
#' @param data Marker input for the (degraded-sim) cohort.
#' @param truth Common-schema truth segments (from [simulate_source()]).
#' @param grid Shared marker grid `data.table(chr, pos)` for the marker metrics.
#' @param caller `"nnil"` (sweeps `rrate`) or `"rtiger"` (sweeps `rigidity`).
#' @param values Parameter grid (e.g. from [log_grid()]).
#' @param threads Fan-out width forwarded to `caller_sweep`.
#' @param refit `"none"` (fit once at `ref`; the fast scan default) or `"cold"`.
#' @param ... Forwarded to `caller_sweep` (`design`/`f_1`,`f_2`, `err`, `ref`, ...).
#' @return `data.table(param, value, donor_frag_F1, donor_frag_FDR,
#'   donor_marker_recall, marker_macroF1, n_breakpoints, truth_bp, ks_fragsize)`,
#'   ascending in `value`.
#' @export
sweep_calibrate <- function(data, truth, grid, caller = c("nnil", "rtiger"),
                            values, threads = 1L, refit = "none", ...) {
  caller <- match.arg(caller)
  pcol <- if (caller == "rtiger") "rigidity" else "rrate"
  segs <- data.table::as.data.table(
    nilHMM::caller_sweep(data,
      caller = caller, values = values,
      refit = refit, threads = threads, ...
    )
  )
  truth <- data.table::as.data.table(truth)
  tsz <- donor_block_sizes(truth)
  tbp <- breakpoint_count(truth)
  data.table::rbindlist(lapply(sort(values), function(v) {
    called <- segs[get(pcol) == v]
    mf <- marker_f1(called, truth, grid)
    ff <- donor_fragment_f1(called, truth)
    dr <- mf$per_class[class == "donor(>0)"]
    data.table::data.table(
      param = pcol, value = v,
      donor_frag_F1 = ff$f1, donor_frag_FDR = ff$fdr,
      donor_marker_recall = dr$recall, marker_macroF1 = mf$macro_f1,
      n_breakpoints = breakpoint_count(called), truth_bp = tbp,
      ks_fragsize = fragment_size_ks(donor_block_sizes(called), tsz)
    )
  }))
}

#' Pick the grid value optimizing an objective from a sweep curve
#'
#' @param scores A [sweep_calibrate()] table.
#' @param objective Column to optimize; maximized unless it is a "lower is better"
#'   metric (`ks_fragsize`, `donor_frag_FDR`), which is minimized.
#' @return List: `value` (the optimizer), `objective`, and the winning `row`.
#' @export
best_value <- function(scores, objective = "donor_frag_F1") {
  if (!objective %in% names(scores)) stop("best_value(): no column '", objective, "'")
  minimize <- objective %in% c("ks_fragsize", "donor_frag_FDR")
  i <- if (minimize) which.min(scores[[objective]]) else which.max(scores[[objective]])
  list(value = scores$value[i], objective = objective, row = scores[i])
}
