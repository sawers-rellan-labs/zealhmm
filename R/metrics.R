# Concordance + accuracy metrics on the common call schema (see staging.R).
#
# Two families:
#   * truth-based (vs simcross ground truth) — per-bin accuracy, donor-footprint
#     ROC pieces, breakpoint-count error. Used by the simulation benchmark (B2.2).
#   * pairwise cross-track (no truth) — segment agreement, donor Jaccard,
#     breakpoint concordance, control cleanliness. Used by the 400-cohort
#     source x method comparison (B2.3), where the paired design means
#     disagreement is method, not sample.
#
# All functions take common-schema tables (name, chr, start_bp, end_bp, state)
# and evaluate on a shared genomic grid so segment tables with different
# breakpoints are directly comparable.

#' Rasterize a segment table onto a fixed bp grid, per chromosome
#'
#' @param seg Common-schema segments for ONE sample/track.
#' @param grid A `data.table(chr, pos)` of evaluation points (e.g. marker or
#'   bin midpoints), shared across the tracks being compared.
#' @return `grid` with an added integer `state` column (NA where uncovered).
#' @export
rasterize_states <- function(seg, grid) {
  seg <- data.table::as.data.table(seg)
  grid <- data.table::copy(data.table::as.data.table(grid))
  grid[, state := NA_integer_]
  for (cc in unique(grid$chr)) {
    s <- seg[chr == cc][order(start_bp)]
    if (!nrow(s)) next
    gi <- which(grid$chr == cc)
    # last segment whose start is <= pos, if pos also <= its end
    idx <- findInterval(grid$pos[gi], s$start_bp)
    ok <- idx >= 1L & grid$pos[gi] <= s$end_bp[pmax(idx, 1L)]
    grid$state[gi[ok]] <- s$state[idx[ok]]
  }
  grid[]
}

#' Per-bin state accuracy of a call set against truth (vs simcross)
#'
#' @param called,truth Common-schema segments for the same sample.
#' @param grid Shared evaluation grid (see [rasterize_states()]).
#' @return List: `accuracy`, `confusion` (3x3 REF/HET/ALT), `n` evaluated bins.
#' @export
state_accuracy <- function(called, truth, grid) {
  c1 <- rasterize_states(called, grid)$state
  t1 <- rasterize_states(truth, grid)$state
  ok <- !is.na(c1) & !is.na(t1)
  cm <- table(
    factor(t1[ok], 0:2, c("REF", "HET", "ALT")),
    factor(c1[ok], 0:2, c("REF", "HET", "ALT"))
  )
  list(
    accuracy = if (sum(cm)) sum(diag(cm)) / sum(cm) else NA_real_,
    confusion = cm, n = sum(ok)
  )
}

#' Donor footprint of a call set: bp with state in {HET, ALT} per chromosome
#'
#' @param seg Common-schema segments for one sample.
#' @return A `data.table(chr, donor_bp)`.
#' @export
donor_footprint <- function(seg) {
  seg <- data.table::as.data.table(seg)[state %in% c(1L, 2L)]
  if (!nrow(seg)) {
    return(data.table::data.table(chr = integer(), donor_bp = numeric()))
  }
  seg[, .(donor_bp = sum(end_bp - start_bp)), by = chr]
}

#' Donor-footprint Jaccard between two call sets on a shared grid
#'
#' Jaccard over grid points called donor (HET or ALT) by each track — the
#' cross-source/cross-caller footprint-overlap metric (B2.3).
#'
#' @param a,b Common-schema segments for the same sample, different tracks.
#' @param grid Shared evaluation grid.
#' @return Jaccard index in \[0, 1\] (NA if neither calls any donor bin).
#' @export
donor_jaccard <- function(a, b, grid) {
  da <- rasterize_states(a, grid)$state %in% c(1L, 2L)
  db <- rasterize_states(b, grid)$state %in% c(1L, 2L)
  inter <- sum(da & db)
  union <- sum(da | db)
  if (!union) NA_real_ else inter / union
}

#' Per-bin state agreement between two call sets on a shared grid
#'
#' @param a,b Common-schema segments for the same sample, different tracks.
#' @param grid Shared evaluation grid.
#' @return Fraction of jointly-covered bins with identical state.
#' @export
pairwise_agreement <- function(a, b, grid) {
  sa <- rasterize_states(a, grid)$state
  sb <- rasterize_states(b, grid)$state
  ok <- !is.na(sa) & !is.na(sb)
  if (!any(ok)) NA_real_ else mean(sa[ok] == sb[ok])
}

#' Breakpoint count (number of state changes) of a call set
#'
#' @param seg Common-schema segments for one sample.
#' @return Integer count of within-chromosome state transitions.
#' @export
breakpoint_count <- function(seg) {
  seg <- data.table::as.data.table(seg)[order(chr, start_bp)]
  seg[, sum(head(state, -1) != tail(state, -1)), by = chr][, sum(V1)]
}

#' Control cleanliness: donor-called fraction on a known-REF control (e.g. B73)
#'
#' Lower is better — a clean control should paint ~all REF. Complements the
#' truth-based metrics for tracks where the "truth" is that there is no donor.
#'
#' @param seg Common-schema segments for a control sample.
#' @param grid Shared evaluation grid.
#' @return Fraction of covered bins called HET or ALT.
#' @export
control_donor_rate <- function(seg, grid) {
  s <- rasterize_states(seg, grid)$state
  ok <- !is.na(s)
  if (!any(ok)) NA_real_ else mean(s[ok] %in% c(1L, 2L))
}
