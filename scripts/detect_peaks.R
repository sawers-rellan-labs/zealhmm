# =============================================================================
# detect_peaks.R — multi-peak QTL detection, ported from airmine (F. Rodríguez,
# scripts/detect_peaks.R + map_single_marker.R).
#
# R/qtl's summary.scanone reports at most ONE peak per chromosome. This module
# adds two things on top of a scanone LOD profile:
#
#   get_peak_table()  — the per-chromosome significant peaks WITH lodint 1.5-LOD
#                       confidence intervals, plus the peak marker and the CI
#                       bounds expressed in v5 bp + width (Mb).
#   refine_peaks()    — the "better" search: per chromosome, Akima-interpolate the
#                       LOD curve (intensity x), run pracma::findpeaks to get up to
#                       `npeaks` sub-peaks above threshold, take each sub-peak's
#                       1.5-LOD-drop interval, then resolve OVERLAPPING sub-peaks
#                       into a maximal set of NON-overlapping QTL via an interval
#                       graph (igraph independent vertex sets), preferring the set
#                       that contains the global max-LOD peak / has the greatest
#                       summed LOD. This recovers multiple linked QTL on one chr.
#
# ZEAL adaptation vs the airmine original:
#   * single trait — the multi-trait "trait : pseudomarker" rowname parsing and
#     dplyr/tidyr row-binding are dropped; tabByChr for one trait is a length-1
#     list holding one data.frame spanning all significant chromosomes.
#   * bp bounds come from the `mk` map (marker, chr, pos = v5 bp, cm). ZEAL markers
#     are renamed to S<chr_v5>_<pos_v5> upstream (zeal_rqtl_dta.R), so the name and
#     the bp agree; we still take bp from `mk$pos` as the authoritative source.
#   * the interval-graph sub-peak search is factored into find_subpeaks() so the
#     whole-genome (scanone) and per-taxon (flat LOD table) paths share it.
# =============================================================================

suppressMessages({
  library(qtl)
  library(pracma)
  library(igraph)
})

# ---- map a cM position on a chromosome to the nearest marker in `mk` ----------
# `mk` columns: marker, chr, pos (v5 bp), cm. Returns a one-row list(marker, bp).
# Uses vector/position indexing only — `mk` may be a data.table, so `mk[mk$chr ==
# chr, ]` would resolve `chr` to the COLUMN (NSE) and never filter by chromosome.
.nearest_marker <- function(mk, chr, cm) {
  w <- which(mk$chr == chr)
  i <- w[which.min(abs(mk$cm[w] - cm))]
  list(marker = mk$marker[i], bp = mk$pos[i], cm = mk$cm[i])
}

# ---- physical (bp) span of a cM interval on a chromosome ----------------------
# Honest under cM-degeneracy: in the pericentromere many markers share a cM, so a
# single nearest-marker endpoint is arbitrary (can exclude a gene by <1 Mb). Return
# the min/max bp of ALL markers whose cM falls in the interval. Returns c(lo_bp, hi_bp).
.ci_bp <- function(mk, chr, lo_cm, hi_cm) {
  w <- which(mk$chr == chr & mk$cm >= min(lo_cm, hi_cm) & mk$cm <= max(lo_cm, hi_cm))
  if (!length(w)) {
    b <- .nearest_marker(mk, chr, lo_cm)$bp
    return(c(b, b))
  }
  c(min(mk$pos[w]), max(mk$pos[w]))
}

# ---- core: multiple sub-peaks in ONE chromosome's LOD profile -----------------
# pos, lod : cM positions and LOD along a single chromosome (any spacing).
# thresh   : minimum peak height (LOD threshold).
# min_trough : two maxima are reported as DISTINCT QTL only if the valley between
#   them drops at least this many LOD below the lower of the two (topographic
#   prominence). Shallower separations are one broad peak (a shoulder/plateau wiggle,
#   e.g. the chr9 pericentromere) and are merged into the taller. Stops the
#   interval-graph step from over-splitting a broad peak.
# smooth_k : optional running-median window (odd, in markers) on the LOD before
#   interpolation. DISABLED by default (smooth_k = 1): thin-miscall LOD dives are
#   corrected upstream in the genotypes (zeal_rqtl_dta.R realistic calc.genoprob
#   error.prob), not smoothed post-scan here. Kept as an optional fallback.
# ci_on : where the reported peak + 1.5-LOD support interval are measured.
#   "markers" (default) -- on the REAL (pos, lod) markers within the peak's basin:
#     peak = the real argmax, support = real markers within `drop` of it. This matches
#     R/qtl lodint and avoids the interpolation grid excluding a near-max marker by one
#     grid step (the ZmCCT10-at-52.4-vs-CI-low-52.5 case).
#   "grid" -- on the Akima-interpolated grid (the earlier behaviour), kept as an option.
# Peak DETECTION and the trough-merge always use the interpolated grid (findpeaks needs
# a regular grid); only the reported peak/CI switch with ci_on.
# The interval is the OUTERMOST-crossing 1.5-LOD support (Lander-Botstein): it spans
# shallow internal dips, and is bounded by the troughs to neighbouring kept peaks (its
# basin) so it never swallows another peak.
# Returns a data.frame(lod, pos_cm, ci_low_cm, ci_high_cm), one row per retained
# sub-peak, or NULL if the profile is too short to interpolate.
find_subpeaks <- function(pos, lod, thresh,
                          intensity = 10, minpeakdist_cm = 5,
                          drop = 1.5, npeaks = 5, smooth_k = 1, min_trough = 2,
                          ci_on = c("markers", "grid")) {
  ci_on <- match.arg(ci_on)
  ok <- is.finite(pos) & is.finite(lod)
  pos <- pos[ok]
  lod <- lod[ok]
  if (length(unique(pos)) < 3) {
    return(NULL)
  }
  o <- order(pos)
  pos <- pos[o]
  lod <- lod[o]
  if (smooth_k >= 3 && length(lod) > smooth_k) {
    lod <- as.numeric(runmed(lod, smooth_k, endrule = "median"))
  }

  grid <- seq(from = min(pos), to = max(pos), by = 1 / intensity)
  akima_lod <- as.numeric(akimaInterp(pos, lod, grid))
  L <- length(akima_lod)
  grid_cm <- function(idx) min(pos) + idx / intensity # airmine index->cM convention

  subpeaks <- pracma::findpeaks(akima_lod,
    minpeakdistance = minpeakdist_cm * intensity,
    minpeakheight = thresh, npeaks = npeaks
  )
  if (is.null(subpeaks) || nrow(subpeaks) == 0) {
    idx <- which.max(akima_lod) # fallback: single peak at the interpolated max
    subpeaks <- matrix(c(akima_lod[idx], idx, idx, idx), nrow = 1)
  }

  # narrow (nearest-crossing) 1.5-LOD interval -- used ONLY to let the interval graph
  # pick a non-overlapping candidate set, NOT as the reported CI.
  in_drop <- t(apply(subpeaks, 1, function(x) {
    d <- round(which((x[1] - akima_lod) > drop) - x[2], 0)
    lo <- if (any(d < 0)) x[2] + max(d[d < 0]) else 1L
    hi <- if (any(d > 0)) x[2] + min(d[d > 0]) else L
    c(lo, hi)
  }))
  if (is.null(dim(in_drop))) in_drop <- matrix(in_drop, ncol = 2, byrow = TRUE)

  # resolve overlapping sub-peaks into a maximal non-overlapping candidate set
  if (nrow(subpeaks) > 1) {
    edge_list <- as.matrix(IRanges::findOverlaps(IRanges::IRanges(in_drop[, 1], in_drop[, 2])))
    g <- graph_from_edgelist(edge_list, directed = FALSE)
    ivs_list <- ivs(g, min = 2)
    if (length(ivs_list) == 0) ivs_list <- ivs(g, min = 1)
    vsize <- sapply(ivs_list, length)
    vs <- ivs_list[vsize == max(vsize)] # largest independent sets
    is_max <- sapply(vs, function(v) max(subpeaks[v, 1]) == max(subpeaks[, 1]))
    keep_sets <- vs[is_max] # ... that contain the top peak
    if (length(keep_sets) > 1) { # tie-break: greatest LOD sum
      sum_lod <- sapply(keep_sets, function(v) sum(subpeaks[v, 1]))
      sel <- as.integer(keep_sets[[which.max(sum_lod)]])
    } else {
      sel <- as.integer(keep_sets[[1]])
    }
  } else {
    sel <- 1L
  }

  # candidate peaks (height, grid index), sorted along the chromosome
  pk <- subpeaks[sel, c(1, 2), drop = FALSE]
  pk <- pk[order(pk[, 2]), , drop = FALSE]

  # merge peaks separated by a trough shallower than min_trough (prominence): one
  # broad peak, not two. Keep the taller of the merged pair.
  repeat {
    if (nrow(pk) < 2) break
    merged <- FALSE
    for (k in seq_len(nrow(pk) - 1L)) {
      valley <- min(akima_lod[pk[k, 2]:pk[k + 1L, 2]])
      if (min(pk[k, 1], pk[k + 1L, 1]) - valley < min_trough) {
        pk <- pk[-(if (pk[k, 1] >= pk[k + 1L, 1]) k + 1L else k), , drop = FALSE]
        merged <- TRUE
        break
      }
    }
    if (!merged) break
  }

  # basins: grid-index boundaries = the troughs between adjacent kept peaks, so a
  # peak's support interval spans internal dips but never crosses into another peak.
  npk <- nrow(pk)
  bnd <- integer(npk + 1L)
  bnd[1] <- 1L
  bnd[npk + 1L] <- L
  if (npk > 1L) {
    for (k in seq_len(npk - 1L)) {
      seg <- pk[k, 2]:pk[k + 1L, 2]
      bnd[k + 1L] <- seg[which.min(akima_lod[seg])]
    }
  }

  out <- t(vapply(seq_len(npk), function(k) {
    if (ci_on == "grid") {
      # interval measured on the interpolated grid (optional)
      basin <- bnd[k]:bnd[k + 1L]
      above <- basin[akima_lod[basin] >= pk[k, 1] - drop]
      if (!length(above)) above <- pk[k, 2]
      c(pk[k, 1], grid_cm(pk[k, 2]), grid_cm(min(above)), grid_cm(max(above)))
    } else {
      # interval measured on the REAL markers within the peak's basin (default):
      # peak = real argmax; support = real markers within `drop` of that max, then
      # EXPANDED to the flanking markers (matches R/qtl lodint expandtomarkers=TRUE),
      # so the interval brackets the 1.5-LOD drop rather than stopping at the last
      # supra-threshold marker (which would exclude a gene lying in the gap to the next).
      ib <- which(pos >= grid_cm(bnd[k]) & pos <= grid_cm(bnd[k + 1L]))
      if (!length(ib)) ib <- which.min(abs(pos - grid_cm(pk[k, 2])))
      pmax_i <- ib[which.max(lod[ib])]
      above <- ib[lod[ib] >= lod[pmax_i] - drop]
      lo_i <- min(above)
      hi_i <- max(above)
      if (lo_i > 1L) lo_i <- lo_i - 1L
      if (hi_i < length(pos)) hi_i <- hi_i + 1L
      c(lod[pmax_i], pos[pmax_i], pos[lo_i], pos[hi_i])
    }
  }, numeric(4)))

  data.frame(lod = out[, 1], pos_cm = out[, 2], ci_low_cm = out[, 3], ci_high_cm = out[, 4])
}

# ---- per-chromosome significant peaks + lodint CIs (one-peak-per-chr) ----------
# single_scan : scanone object (LOD in column 3); mk : marker map.
# Returns a data.frame or NULL if nothing is significant.
get_peak_table <- function(single_scan, perms, mk, alpha = 0.05) {
  s <- summary(single_scan,
    perms = perms, alpha = alpha, format = "tabByChr",
    ci.function = "lodint", expandtomarkers = TRUE, pvalues = TRUE
  )
  df <- if (length(s) >= 1 && !is.null(s[[1]])) as.data.frame(s[[1]]) else NULL
  if (is.null(df) || nrow(df) == 0) {
    return(NULL)
  }
  df$name <- rownames(df)
  df$chr <- as.integer(as.character(df$chr))
  df$thresh <- as.numeric(summary(perms, alpha = alpha))[1]
  # peak marker (nearest) + CI bounds as the physical span of the cM interval
  pk <- Map(.nearest_marker, list(mk), df$chr, df$pos)
  cibp <- t(mapply(function(ch, lo, hi) .ci_bp(mk, ch, lo, hi), df$chr, df$ci.low, df$ci.high))
  df$marker <- vapply(pk, `[[`, "", "marker")
  df$ci_left <- cibp[, 1]
  df$ci_right <- cibp[, 2]
  df$width_mb <- (df$ci_right - df$ci_left) / 1e6
  rownames(df) <- NULL
  df[, c(
    "name", "chr", "pos", "lod", "pval", "thresh", "ci.low", "ci.high",
    "marker", "ci_left", "ci_right", "width_mb"
  )]
}

# ---- multi-peak refinement over the significant chromosomes -------------------
# peaks : output of get_peak_table (gives the significant chromosomes + thresh).
# single_scan : scanone object (LOD in column 3).
# Returns a data.frame with one row per refined QTL (may be >1 per chromosome).
refine_peaks <- function(peaks, single_scan, mk) {
  if (is.null(peaks) || nrow(peaks) == 0) {
    return(NULL)
  }
  lodcol <- names(single_scan)[3] # scanone convention: chr, pos, <lod>
  out <- lapply(unique(peaks$chr), function(ch) {
    thr <- peaks$thresh[match(ch, peaks$chr)]
    d <- single_scan[as.character(single_scan$chr) == as.character(ch), ]
    sp <- find_subpeaks(d$pos, d[[lodcol]], thresh = thr)
    if (is.null(sp) || nrow(sp) == 0) {
      return(NULL)
    }
    pkm <- Map(.nearest_marker, list(mk), ch, sp$pos_cm)
    cibp <- t(mapply(function(lo, hi) .ci_bp(mk, ch, lo, hi), sp$ci_low_cm, sp$ci_high_cm))
    data.frame(
      chr = ch,
      pos = sp$pos_cm,
      lod = sp$lod,
      thresh = thr,
      ci.low = sp$ci_low_cm,
      ci.high = sp$ci_high_cm,
      marker = vapply(pkm, `[[`, "", "marker"),
      ci_left = cibp[, 1],
      ci_right = cibp[, 2]
    )
  })
  out <- do.call(rbind, out)
  if (is.null(out)) {
    return(NULL)
  }
  out <- out[out$lod > out$thresh, ]
  out$width_mb <- (out$ci_right - out$ci_left) / 1e6
  out$name <- sprintf("%d@%.1f", out$chr, out$pos)
  rownames(out) <- NULL
  out[, c(
    "name", "chr", "pos", "lod", "thresh", "ci.low", "ci.high",
    "marker", "ci_left", "ci_right", "width_mb"
  )]
}
