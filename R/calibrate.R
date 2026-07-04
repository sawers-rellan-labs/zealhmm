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

#' Keep only the rigidity values a cohort can support
#'
#' RTIGER requires more than `2 * rigidity` covered markers on every chromosome of
#' every sample (below that the E-step degenerates and `caller_sweep` hard-stops).
#' Drop grid values that violate this against the tightest `(sample, chromosome)`
#' in `data`, so a rigidity sweep never trips the floor.
#'
#' @param data Marker input (`name, chr, n_ref, n_alt`).
#' @param values Candidate rigidity grid (e.g. `2^(1:9)`).
#' @return The feasible subset of `values` (ascending); warns about any dropped.
#' @export
feasible_rigidity <- function(data, values) {
  d <- data.table::as.data.table(data)[n_ref + n_alt > 0L]
  min_cov <- d[, .N, by = list(name, chr)][, min(N)]
  keep <- sort(values[2L * values < min_cov])
  if (length(keep) < length(values)) {
    message(sprintf(
      "feasible_rigidity: min covered markers/chromosome = %d -> keeping rigidity <= %d (dropped %s)",
      min_cov, if (length(keep)) max(keep) else 0L,
      paste(setdiff(values, keep), collapse = ", ")
    ))
  }
  if (!length(keep)) stop("feasible_rigidity: no rigidity satisfies 2*r < ", min_cov)
  keep
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
#' @return `data.table(param, value, donor_frag_dice, donor_frag_FDR,
#'   donor_marker_recall, marker_macro_dice, n_breakpoints, truth_bp, ks_fragsize)`,
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
    mf <- marker_dice(called, truth, grid)
    ff <- donor_fragment_dice(called, truth)
    dr <- mf$per_class[class == "donor(>0)"]
    data.table::data.table(
      param = pcol, value = v,
      donor_frag_dice = ff$dice, donor_frag_FDR = ff$fdr,
      donor_marker_recall = dr$recall, marker_macro_dice = mf$macro_dice,
      n_breakpoints = breakpoint_count(called), truth_bp = tbp,
      ks_fragsize = fragment_size_ks(donor_block_sizes(called), tsz)
    )
  }))
}

#' Bracket the optimum from a log-sweep curve (for golden_refine)
#'
#' Returns the two grid values straddling the [best_value()] optimum -- the
#' interval to hand to [golden_refine()]. At a grid edge, uses the edge and its
#' one neighbour.
#'
#' @param scores A [sweep_calibrate()] table.
#' @param objective Column optimized (see [best_value()]).
#' @return Numeric `c(lo, hi)`.
#' @export
bracket_from_sweep <- function(scores, objective = "donor_frag_dice") {
  s <- scores[order(scores$value)]
  i <- match(best_value(s, objective)$value, s$value)
  lo <- s$value[max(1L, i - 1L)]
  hi <- s$value[min(nrow(s), i + 1L)]
  if (lo == hi) { # optimum at an edge: widen to the adjacent interior point
    if (i == 1L) hi <- s$value[min(nrow(s), 2L)] else lo <- s$value[max(1L, nrow(s) - 1L)]
  }
  c(lo = lo, hi = hi)
}

#' Golden-section refine of a caller's parameter within a bracket
#'
#' Stage 2 of calibration: golden-section search on `log10(parameter)` inside
#' `[lo, hi]` (from [bracket_from_sweep()]), evaluating the objective at each probe
#' via one [nilHMM::caller_sweep()] fit + truth scoring. Reuses one interior point
#' per iteration (one new evaluation per step) and caches by value, so an
#' integer parameter (rtiger rigidity) never re-fits the same rigidity.
#'
#' @param data,truth,grid,caller,threads,... As in [sweep_calibrate()].
#' @param lo,hi Bracket endpoints (`0 < lo < hi`).
#' @param objective Column to optimize (Dice maximized; `ks_fragsize`/`donor_frag_FDR` minimized).
#' @param tol Stop when the bracket width on the log10 scale is below this
#'   (default 0.05 ~ a 1.12x ratio).
#' @param max_iter Iteration cap.
#' @param integer Round probes to integers (default: `TRUE` for rtiger).
#' @return List: `value` (refined optimum), `objective`, `score` (its scored row),
#'   `trace` (per-iteration bracket), `evals` (all probed value/objective pairs).
#' @export
golden_refine <- function(data, truth, grid, caller = c("nnil", "rtiger"),
                          lo, hi, objective = "donor_frag_dice", threads = 1L,
                          tol = 0.05, max_iter = 20L, integer = NULL, ...) {
  caller <- match.arg(caller)
  if (lo <= 0 || hi <= 0 || lo >= hi) stop("golden_refine(): need 0 < lo < hi")
  if (is.null(integer)) integer <- caller == "rtiger"
  minimize <- objective %in% c("ks_fragsize", "donor_frag_FDR")
  cache <- new.env(parent = emptyenv())
  score_at <- function(logx) {
    v <- 10^logx
    if (integer) v <- max(1L, as.integer(round(v)))
    key <- format(v, digits = 15)
    if (!is.null(cache[[key]])) {
      return(cache[[key]])
    }
    sc <- sweep_calibrate(data, truth, grid,
      caller = caller, values = v,
      threads = threads, refit = "none", ...
    )
    y <- sc[[objective]][1]
    res <- list(value = v, y = if (minimize) -y else y, obj = y, row = sc)
    assign(key, res, cache)
    res
  }
  gr <- (sqrt(5) - 1) / 2 # 0.618...
  L <- log10(lo)
  H <- log10(hi)
  c1 <- H - gr * (H - L)
  c2 <- L + gr * (H - L)
  p1 <- score_at(c1) # probe at interior point 1 (not the golden-section objective)
  p2 <- score_at(c2)
  trace <- vector("list", max_iter)
  for (it in seq_len(max_iter)) {
    if (p1$y >= p2$y) {
      H <- c2
      c2 <- c1
      p2 <- p1
      c1 <- H - gr * (H - L)
      p1 <- score_at(c1)
    } else {
      L <- c1
      c1 <- c2
      p1 <- p2
      c2 <- L + gr * (H - L)
      p2 <- score_at(c2)
    }
    trace[[it]] <- data.table::data.table(iter = it, lo = 10^L, hi = 10^H)
    if ((H - L) < tol) break
    if (integer && (10^H - 10^L) < 1) break # bracket narrower than one integer
  }
  evals <- data.table::rbindlist(lapply(ls(cache), function(k) {
    e <- cache[[k]]
    data.table::data.table(value = e$value, objective = e$obj)
  }))[order(value)]
  best <- if (minimize) evals[which.min(objective)] else evals[which.max(objective)]
  win <- cache[[format(best$value, digits = 15)]]
  list(
    value = best$value, objective = objective, score = win$row,
    trace = data.table::rbindlist(trace), evals = evals
  )
}

#' Pick the grid value optimizing an objective from a sweep curve
#'
#' @param scores A [sweep_calibrate()] table.
#' @param objective Column to optimize; maximized unless it is a "lower is better"
#'   metric (`ks_fragsize`, `donor_frag_FDR`), which is minimized.
#' @return List: `value` (the optimizer), `objective`, and the winning `row`.
#' @export
best_value <- function(scores, objective = "donor_frag_dice") {
  if (!objective %in% names(scores)) stop("best_value(): no column '", objective, "'")
  minimize <- objective %in% c("ks_fragsize", "donor_frag_FDR")
  i <- if (minimize) which.min(scores[[objective]]) else which.max(scores[[objective]])
  list(value = scores$value[i], objective = objective, row = scores[i])
}
