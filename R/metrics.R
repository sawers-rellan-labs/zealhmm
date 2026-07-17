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

# ============================================================================
# Truth-based calibration/validation (vs simcross truth) — the paper's F2/F3.
# All operate on the common schema (name, chr, start_bp, end_bp, state).
# ============================================================================

#' Per-sample genome fractions and donor dosage (length-weighted)
#'
#' @param seg Common-schema segments (one or many samples).
#' @return `data.table(name, REF, HET, ALT, dosage)`; dosage = ALT + HET/2.
#' @export
genotype_fractions <- function(seg) {
  seg <- data.table::as.data.table(seg)
  seg[, len := end_bp - start_bp]
  out <- seg[, .(
    REF = sum(len[state == 0L]) / sum(len),
    HET = sum(len[state == 1L]) / sum(len),
    ALT = sum(len[state == 2L]) / sum(len)
  ), by = name]
  out[, dosage := ALT + HET / 2][]
}

#' Merge adjacent donor segments into maximal introgression blocks (per name, chr)
#' @keywords internal
.donor_blocks <- function(seg, states = c(1L, 2L)) {
  seg <- data.table::as.data.table(seg)[order(name, chr, start_bp)]
  seg[, don := as.integer(state %in% states)]
  seg[,
    {
      r <- rle(don)
      e <- cumsum(r$lengths)
      s <- c(1L, utils::head(e, -1L) + 1L)
      keep <- r$values == 1L
      list(start_bp = start_bp[s[keep]], end_bp = end_bp[e[keep]])
    },
    by = .(name, chr)
  ]
}

#' Donor introgression block sizes in Mb (across all samples)
#' @param seg Common-schema segments; `states` = donor states (default HET+ALT).
#' @return Numeric vector of block sizes (Mb).
#' @export
donor_block_sizes <- function(seg, states = c(1L, 2L)) {
  b <- .donor_blocks(seg, states)
  if (!nrow(b)) {
    return(numeric(0))
  }
  (b$end_bp - b$start_bp) / 1e6
}

#' KS distance between two donor-block-size distributions (fragment-size fit)
#'
#' `D = sup_x |F_a(x) - F_b(x)|`; lower = better fit to the simulated truth.
#' @param a,b Numeric size vectors (e.g. from [donor_block_sizes()]).
#' @return The KS D statistic (NA if either is empty).
#' @export
fragment_size_ks <- function(a, b) {
  if (!length(a) || !length(b)) {
    return(NA_real_)
  }
  suppressWarnings(as.numeric(stats::ks.test(a, b)$statistic))
}

#' Per-state marker precision/recall/Dice vs truth (pooled over samples)
#'
#' The marker-level score. Includes a binary `donor(>0)` row (introgression
#' present) — its **recall is the "marker true-positive rate" Holland maximized**;
#' its **Dice** is what we argue to optimize instead (penalizes the false donor
#' calls that drive over-fragmentation).
#'
#' @param called,truth Common-schema segments for the same samples.
#' @param grid Shared evaluation grid (`data.table(chr, pos)`).
#' @return List: `per_class` (REF/HET/ALT + donor(>0): precision, recall, dice,
#'   n_truth), `macro_dice` (mean over the 3 states), `accuracy`, `n`.
#' @export
marker_dice <- function(called, truth, grid) {
  called <- data.table::as.data.table(called)
  truth <- data.table::as.data.table(truth)
  nms <- intersect(unique(called$name), unique(truth$name))
  cc <- integer(0)
  tt <- integer(0)
  for (nm in nms) {
    c1 <- rasterize_states(called[name == nm], grid)$state
    t1 <- rasterize_states(truth[name == nm], grid)$state
    ok <- !is.na(c1) & !is.na(t1)
    cc <- c(cc, c1[ok])
    tt <- c(tt, t1[ok])
  }
  prf <- function(pred, tru) {
    tp <- sum(pred & tru)
    fp <- sum(pred & !tru)
    fn <- sum(!pred & tru)
    prec <- if (tp + fp) tp / (tp + fp) else NA_real_
    rec <- if (tp + fn) tp / (tp + fn) else NA_real_
    dice <- if (!is.na(prec) && !is.na(rec) && (prec + rec) > 0) 2 * prec * rec / (prec + rec) else NA_real_
    c(precision = prec, recall = rec, dice = dice)
  }
  lab <- c("REF", "HET", "ALT")
  per <- data.table::rbindlist(lapply(0:2, function(s) {
    v <- prf(cc == s, tt == s)
    data.table::data.table(
      class = lab[s + 1L], precision = v[["precision"]],
      recall = v[["recall"]], dice = v[["dice"]], n_truth = sum(tt == s)
    )
  }))
  vd <- prf(cc > 0L, tt > 0L)
  per <- rbind(per, data.table::data.table(
    class = "donor(>0)",
    precision = vd[["precision"]], recall = vd[["recall"]], dice = vd[["dice"]],
    n_truth = sum(tt > 0L)
  ))
  list(
    per_class = per, macro_dice = mean(per$dice[1:3], na.rm = TRUE),
    accuracy = if (length(cc)) mean(cc == tt) else NA_real_, n = length(cc)
  )
}

#' Donor-fragment precision/recall/Dice by reciprocal overlap (segment level)
#'
#' A called donor block matches a truth block when they **reciprocally** overlap
#' by at least `min_overlap` (overlap >= min_overlap of *each* block's length).
#' Recall = matched truth blocks / truth blocks; precision = matched called /
#' called. This is the block-level score; contrast with marker Dice to expose
#' over-fragmentation (many spurious blocks -> low precision, high FDR).
#'
#' @param called,truth Common-schema segments for the same samples.
#' @param states Donor states (default HET+ALT = "introgression present").
#' @param min_overlap Reciprocal-overlap threshold (default 0.5).
#' @return List: `precision, recall, dice, fdr, n_truth, n_called`.
#' @export
donor_fragment_dice <- function(called, truth, states = c(1L, 2L), min_overlap = 0.5) {
  cb <- .donor_blocks(called, states)
  tb <- .donor_blocks(truth, states)
  cb[, gk := paste(name, chr, sep = "\r")]
  tb[, gk := paste(name, chr, sep = "\r")]
  cb[, hit := FALSE]
  tb[, hit := FALSE]
  for (g in intersect(cb$gk, tb$gk)) {
    ci <- which(cb$gk == g)
    ti <- which(tb$gk == g)
    for (a in ci) {
      for (b in ti) {
        ov <- max(0, min(cb$end_bp[a], tb$end_bp[b]) - max(cb$start_bp[a], tb$start_bp[b]))
        # Guard the reciprocal-overlap denominators: a single-marker donor block
        # has end_bp == start_bp (0 bp span) -> max(span, 1) avoids the 0/0 = NaN
        # that would break the comparison. A point block has ov = 0 either way, so
        # it never reaches the 50% threshold (spurious 1-marker blocks stay unmatched).
        if (ov / max(cb$end_bp[a] - cb$start_bp[a], 1) >= min_overlap &&
          ov / max(tb$end_bp[b] - tb$start_bp[b], 1) >= min_overlap) {
          cb$hit[a] <- TRUE
          tb$hit[b] <- TRUE
        }
      }
    }
  }
  n_truth <- nrow(tb)
  n_called <- nrow(cb)
  recall <- if (n_truth) sum(tb$hit) / n_truth else NA_real_
  precision <- if (n_called) sum(cb$hit) / n_called else NA_real_
  dice <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0) {
    2 * precision * recall / (precision + recall)
  } else {
    NA_real_
  }
  list(
    precision = precision, recall = recall, dice = dice, fdr = 1 - precision,
    n_truth = n_truth, n_called = n_called
  )
}

#' Calibration sweep: score a caller across a duration-parameter grid
#'
#' Runs `nilHMM::call_ancestry(data, caller, <param> = value)` for each value,
#' scoring vs `truth`. Returns the curve used for the F2 calibration panel and
#' for picking the Dice optimum vs the marker-true-positive optimum.
#'
#' @param data Marker input for the (degraded-sim) cohort.
#' @param truth Common-schema truth segments (BC2S2 simcross).
#' @param grid Shared marker/bin grid for the marker metrics.
#' @param caller "nnil" / "rtiger" / etc.
#' @param param The swept argument name ("rrate" or "rigidity").
#' @param values Values to sweep.
#' @param ... Forwarded to [nilHMM::call_ancestry()] (e.g. design, err).
#' @return `data.table(param, value, marker_macro_dice, donor_marker_recall,
#'   donor_marker_dice, donor_frag_dice, donor_frag_FDR, n_breakpoints, ks_fragsize)`.
#' @export
calibrate_sweep <- function(data, truth, grid, caller, param, values, ...) {
  truth_sizes <- donor_block_sizes(truth)
  data.table::rbindlist(lapply(values, function(v) {
    args <- c(list(data = data, caller = caller), stats::setNames(list(v), param), list(...))
    called <- data.table::as.data.table(do.call(nilHMM::call_ancestry, args))
    mf <- marker_dice(called, truth, grid)
    ff <- donor_fragment_dice(called, truth)
    donor_row <- mf$per_class[class == "donor(>0)"]
    data.table::data.table(
      param = param, value = v,
      marker_macro_dice = mf$macro_dice,
      donor_marker_recall = donor_row$recall, # the Holland "true-positive" objective
      donor_marker_dice = donor_row$dice,
      donor_frag_dice = ff$dice, donor_frag_FDR = ff$fdr,
      n_breakpoints = breakpoint_count(called), # over-fragmentation proxy
      ks_fragsize = fragment_size_ks(donor_block_sizes(called), truth_sizes)
    )
  }))
}
