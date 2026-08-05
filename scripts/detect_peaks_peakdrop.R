#!/usr/bin/env Rscript
# Port of R/qtl2's find_peaks peakdrop algorithm (Broman) onto our R/qtl (v1)
# bcsft scanone LOD curves. find_peaks operates purely on a LOD-vs-position
# curve, so it is cross-type agnostic: this reproduces Broman's single-pass
# peakdrop peak detector (src/find_peaks.cpp: find_peaks_valleys) plus the
# LOD-drop support interval (lod_int_contained), giving a principled middle
# between R/qtl refine (over-segments) and one peak per chromosome (too coarse).
#
# peakdrop rule: on a chromosome, two local maxima above `threshold` are counted
# as separate peaks only if the LOD between them drops by at least `peakdrop`
# below both. peakdrop = Inf collapses to one peak per chromosome.
#
# Parameters are FUNCTION ARGUMENTS with defaults (R-package idiom), not env
# vars. `source()` this file to get the functions; run it with Rscript to
# execute run_peakdrop() over the standard trait set (guarded at the bottom).

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(qtl)
})

logf <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")

# ---- Broman find_peaks_valleys: single forward pass over one chromosome's LOD ----
# Returns 1-based peak indices and the valley indices bounding them
# (length n_peaks + 1: [start, valley_1_2, ..., end]).
find_peaks_valleys <- function(lod, threshold, peakdrop) {
  n <- length(lod)
  peaks <- integer(0)
  valleys <- 1L # start
  n_peaks <- 0L
  last_peak <- 0
  min_since <- 0
  min_since_loc <- 1L
  for (i in seq_len(n)) {
    if (lod[i] < min_since) {
      min_since <- lod[i]
      min_since_loc <- i
    }
    if (lod[i] > threshold) { # possible peak
      if (n_peaks == 0L) { # first one
        peaks <- c(peaks, i)
        last_peak <- lod[i]
        min_since <- lod[i]
        n_peaks <- 1L
      } else if (last_peak > min_since + peakdrop && lod[i] > min_since + peakdrop) { # new peak
        peaks <- c(peaks, i)
        valleys <- c(valleys, min_since_loc)
        min_since_loc <- i
        last_peak <- lod[i]
        min_since <- lod[i]
        n_peaks <- n_peaks + 1L
      } else if (lod[i] > last_peak) { # move the current peak up the same hill
        peaks[n_peaks] <- i
        last_peak <- lod[i]
        min_since <- lod[i]
        min_since_loc <- i
      }
    }
  }
  valleys <- c(valleys, n) # end
  list(peaks = peaks, valleys = valleys)
}

# ---- Broman lod_int_contained: LOD-drop interval, search bounded by adjacent valleys ----
lod_int_contained <- function(lod, peak, drop, left_bound, right_bound) {
  target <- lod[peak] - drop
  lo <- peak
  j <- peak - 1L
  while (j >= left_bound && lod[j] >= target) {
    lo <- j
    j <- j - 1L
  }
  hi <- peak
  j <- peak + 1L
  while (j <= right_bound && lod[j] >= target) {
    hi <- j
    j <- j + 1L
  }
  c(lo, hi)
}

marker_bp <- function(m) as.numeric(sub("^S[0-9]+_", "", m)) # v5 bp from "S<chr>_<bp>"

# v5 bp for any grid position (incl. pseudomarkers) by interpolating cM -> bp
# from the real (S-named) markers on that chromosome. Returns a function of cM.
bp_interp <- function(sub) {
  real <- sub[grepl("^S[0-9]+_", marker)]
  bp <- marker_bp(real$marker)
  ok <- is.finite(real$pos) & is.finite(bp)
  if (sum(ok) < 2) {
    return(function(cm) rep(NA_real_, length(cm)))
  }
  approxfun(real$pos[ok], bp[ok], rule = 2, ties = mean)
}

#' Peakdrop peaks for one scanone LOD table.
#'
#' @param so scanone data.table with columns marker, chr, pos (cM), lod.
#' @param threshold genome-wide LOD threshold (a peak must exceed it).
#' @param peakdrop LOD valley depth required to split two peaks on a chromosome.
#'   Inf collapses to one peak per chromosome (Broman default).
#' @param drop LOD-drop support interval width (must be <= peakdrop).
#' @return data.table, one row per peak, columns matching the existing peak
#'   tables plus pos_mb (v5 Mb of the peak).
peakdrop_peaks <- function(so, threshold, peakdrop = 5, drop = 1.5) {
  stopifnot(drop <= peakdrop)
  rows <- list()
  for (cc in sort(unique(so$chr))) {
    sub <- so[chr == cc][order(pos)]
    bpf <- bp_interp(sub)
    pv <- find_peaks_valleys(sub$lod, threshold, peakdrop)
    for (k in seq_along(pv$peaks)) {
      pk <- pv$peaks[k]
      iv <- lod_int_contained(sub$lod, pk, drop, pv$valleys[k], pv$valleys[k + 1L])
      bpl <- bpf(sub$pos[iv[1]])
      bpr <- bpf(sub$pos[iv[2]])
      rows[[length(rows) + 1L]] <- data.table(
        name = sprintf("%d@%.1f", cc, sub$pos[pk]),
        chr = cc, pos = sub$pos[pk], lod = sub$lod[pk],
        thresh = threshold, peakdrop = peakdrop, drop = drop,
        ci.low = sub$pos[iv[1]], ci.high = sub$pos[iv[2]],
        marker = sub$marker[pk], pos_mb = bpf(sub$pos[pk]) / 1e6,
        ci_left = min(bpl, bpr), ci_right = max(bpl, bpr),
        width_mb = abs(bpr - bpl) / 1e6
      )
    }
  }
  if (length(rows)) {
    rbindlist(rows)
  } else {
    data.table( # typed 0-row table so fwrite still writes a header
      name = character(), chr = integer(), pos = numeric(), lod = numeric(),
      thresh = numeric(), peakdrop = numeric(), drop = numeric(),
      ci.low = numeric(), ci.high = numeric(), marker = character(),
      pos_mb = numeric(), ci_left = numeric(), ci_right = numeric(), width_mb = numeric()
    )
  }
}

#' Peakdrop peaks for one trait, reading its scanone and perm files.
#'
#' @param trait trait stem (lowercase, e.g. "dta").
#' @param tag scan model tag (default "taxon", the taxon-covariate joint scan).
#' @param peakdrop,drop see [peakdrop_peaks()].
#' @param alpha genome-wide significance for the permutation threshold.
#' @param rq directory holding the R/qtl scan outputs.
#' @return peak data.table, or NULL if inputs are missing.
peakdrop_trait <- function(trait, tag = "taxon", peakdrop = 5, drop = 1.5,
                           alpha = 0.05, rq = here("results/sim/zeal/rqtl")) {
  sf <- file.path(rq, sprintf("zeal_%s_scanone_%s.csv", trait, tag))
  pf <- file.path(rq, sprintf("zeal_%s_perms_%s.rds", trait, tag))
  if (!file.exists(sf) || !file.exists(pf)) {
    logf("SKIP %s: missing scanone or perms", trait)
    return(NULL)
  }
  thr <- as.numeric(summary(readRDS(pf), alpha = alpha))[1] # genome-wide 1000-perm threshold
  peakdrop_peaks(fread(sf), thr, peakdrop = peakdrop, drop = drop)
}

#' Driver: write per-trait peakdrop peak tables and a peakdrop-sweep count summary.
#'
#' @param traits trait stems; NULL auto-detects every scanone_<tag> file (minus "rqtl").
#' @param peakdrop peakdrop used for the WRITTEN per-trait tables.
#' @param peakdrops peakdrop values swept in the count-comparison summary.
#' @param write if TRUE, write CSVs; always returns the summary data.table.
run_peakdrop <- function(traits = NULL, tag = "taxon", peakdrop = 5,
                         peakdrops = c(1.5, 3, 5, Inf), drop = 1.5, alpha = 0.05,
                         rq = here("results/sim/zeal/rqtl"), write = TRUE) {
  if (is.null(traits)) {
    traits <- setdiff(
      sub(
        sprintf("^zeal_(.*)_scanone_%s\\.csv$", tag), "\\1",
        list.files(rq, pattern = sprintf("_scanone_%s\\.csv$", tag))
      ),
      "rqtl"
    )
  }
  logf(
    "primary peakdrop=%.2f  sweep={%s}  drop=%.2f  alpha=%.2f  tag=%s",
    peakdrop, paste(peakdrops, collapse = ","), drop, alpha, tag
  )
  summ_rows <- list()
  for (tr in traits) {
    sf <- file.path(rq, sprintf("zeal_%s_scanone_%s.csv", tr, tag))
    pf <- file.path(rq, sprintf("zeal_%s_perms_%s.rds", tr, tag))
    if (!file.exists(sf) || !file.exists(pf)) {
      logf("SKIP %s: missing scanone or perms", tr)
      next
    }
    so <- fread(sf)
    thr <- as.numeric(summary(readRDS(pf), alpha = alpha))[1]
    if (write) {
      fwrite(
        peakdrop_peaks(so, thr, peakdrop = peakdrop, drop = drop),
        file.path(rq, sprintf("zeal_%s_peaks_peakdrop_%s.csv", tr, tag))
      )
    }
    counts <- vapply(peakdrops, function(pd) nrow(peakdrop_peaks(so, thr, pd, drop)), integer(1))
    srow <- data.table(trait = tr, threshold = round(thr, 2))
    for (j in seq_along(peakdrops)) srow[[sprintf("pd_%s", peakdrops[j])]] <- counts[j]
    summ_rows[[length(summ_rows) + 1L]] <- srow
    logf("%-10s thr=%.2f  peaks@pd={%s}", tr, thr, paste(counts, collapse = ","))
  }
  summ <- rbindlist(summ_rows, fill = TRUE)
  if (write) fwrite(summ, file.path(rq, sprintf("zeal_peaks_peakdrop_%s_summary.csv", tag)))
  summ
}

# Run when invoked as a script (Rscript), stay quiet when source()'d for the functions.
if (sys.nframe() == 0L) {
  summ <- run_peakdrop()
  cat("\n==== peak counts by peakdrop (pd_Inf = one peak per chromosome) ====\n")
  print(summ)
}
