# map_tools.R — shared helpers for the TeoNAM genetic-map build
# (sourced by teonam_qtl_permap.R and teonam_qtl_map.R).

# --- Canonical TeoNAM map ------------------------------------------------------
# Default genetic map for downstream analyses = the NATIVE TeoNAM v5 est.map
# (replaces the Ed Coe consensus `map_v5_coe2008.tsv`, which was the old default).
# NOTE: in this file the native cM is column `cm`; the Ed Coe consensus cM is
# retained alongside it as `cm_coe2008` for comparison. (The standalone consensus
# map `map_v5_coe2008.tsv` keeps its own column `cm` — the filename disambiguates.)
DEFAULT_TEONAM_MAP <- "data/teonam/teonam_v5_native.tsv"

#' Flag "quirky" markers as isolated small connected-components on an ordered map.
#'
#' Generalizes the singleton both-adjacent-gap test to isolated CLUSTERS. On a
#' per-chromosome ordered cM map, connected components under a distance threshold
#' are a 1-D segmentation (no graph needed): split the sorted positions wherever a
#' gap exceeds `thr`; each run is a component. A component is an ISLAND (quirky) iff
#'   (1) it is isolated from BOTH neighbours by > `thr` (chromosome ends count as
#'       isolated), AND
#'   (2) it is small: <= `max_n` markers AND more compact than the gaps isolating
#'       it (internal span < min(left_gap, right_gap)).
#' The compactness test (2b) is scale-free — it flags a cluster only when the group
#' is tighter than the empty map around it, so a genuinely spread-out isolated
#' segment is left alone. The singleton both-adjacent rule is the n==1 special case,
#' so this is a strict superset: it flags everything the old rule did, plus the
#' multi-marker islands the old rule missed (e.g. the chr7 v2->v5 displaced block).
#'
#' @param pos  named numeric cM vector for ONE chromosome (names = markers; any order).
#' @param thr  gap-outlier threshold in cM (e.g. the data-driven `gap_out_thr`).
#' @param max_n safety cap on island size (default 20); components larger than this
#'   are never flagged even if compact+isolated.
#' @return character vector of marker names in flagged island components.
#'
#' NOTE on threshold choice: `thr` is the ISOLATION gap. For the fine, data-driven
#' singleton threshold (e.g. `gap_out_thr` ~0.2 cM) a size-1 cap reproduces the old
#' both-adjacent singleton rule; but at that fine threshold a size cap > 1 over-flags
#' normal dense-map clumping (tiny groups separated by ~0.2 cM look "isolated"). To
#' catch true displaced CLUSTERS you need a COARSE `thr` (cM-scale break). Use
#' [find_quirky()] which combines both regimes.
find_quirky_islands <- function(pos, thr, max_n = 20L) {
  n <- length(pos)
  if (n == 0L) {
    return(character(0))
  }
  o <- order(pos)
  pos <- pos[o]
  nm <- names(pos)
  seg <- cumsum(c(TRUE, diff(pos) > thr)) # component id, O(n) single pass
  quirky <- character(0)
  for (s in unique(seg)) {
    idx <- which(seg == s)
    lo <- idx[1L]
    hi <- idx[length(idx)]
    left_gap <- if (lo == 1L) Inf else pos[lo] - pos[lo - 1L] # chr end = isolated
    right_gap <- if (hi == n) Inf else pos[hi + 1L] - pos[hi]
    span <- pos[hi] - pos[lo]
    iso <- left_gap > thr && right_gap > thr
    small <- length(idx) <= max_n && span < min(left_gap, right_gap)
    if (iso && small) quirky <- c(quirky, nm[idx])
  }
  quirky
}

#' Combined quirky-marker finder: fine singletons + coarse isolated clusters.
#'
#' Two regimes, unioned:
#'   (1) SINGLETONS at the fine, data-driven `fine_thr` (e.g. `gap_out_thr`), size
#'       cap 1 — reproduces the original both-adjacent singleton rule.
#'   (2) ISLANDS (clusters) at a COARSE `island_thr` (a genetically implausible
#'       break; the map's 99.99th-percentile gap is ~1 cM, so a few cM is a real
#'       break), size cap `island_max_n` — catches displaced blocks the singleton
#'       rule misses (e.g. the chr7 v2->v5 relocated block, isolated by ~5 cM gaps).
#'
#' @param pos named cM vector for ONE chromosome (names = markers).
#' @param fine_thr fine isolation threshold for singletons (cM).
#' @param island_thr coarse isolation threshold for clusters (cM, default 2).
#' @param island_max_n max markers in a flagged cluster (default 20).
#' @return character vector of quirky marker names (singletons + island clusters).
find_quirky <- function(pos, fine_thr, island_thr = 2, island_max_n = 20L) {
  union(
    find_quirky_islands(pos, thr = fine_thr, max_n = 1L),
    find_quirky_islands(pos, thr = island_thr, max_n = island_max_n)
  )
}
