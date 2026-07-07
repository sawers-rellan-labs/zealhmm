#!/usr/bin/env Rscript
# =============================================================================
# Per-subpopulation TeoNAM genetic maps (R/qtl bcsft(1,4)) -- airmine pipeline.
# Markers ordered on B73 v5 (chr_v5, pos_v5). Validation vs Chen Table 1
# (per-subpop avg 1461 cM, range 1348-1596; ~13,733 SNPs/subpop).
#
# Cross = bcsft(BC.gen=1, F.gen=4) -- R/qtl's BC-then-self model = Shannon 2012,
# the "modified R/qtl" Chen cites. Reproducible bad-marker removal (replacing
# Chen's "visual inspection"), ported from airmine/scripts/make_map_qtl_cross.R:
#   (1) SEGREGATION DISTORTION: per-marker chi-square of the observed (AA,Aa,aa)
#       counts vs the expected BC1S4 freqs (0.734, 0.031, 0.234, from the
#       transition matrices). Flagged as OUTLIERS by the empirical-CDF -> qnorm
#       z renormalization (z > 1.96), per chromosome. Drop, then est.map.
#   (2) QUIRKY / isolated markers: on the estimated map, a marker whose BOTH
#       adjacent cM gaps are DISTRIBUTIONAL OUTLIERS is dropped. The gap threshold
#       is NOT fixed -- it's set by the SAME ecdf->qnorm z>1.96 rule applied to the
#       observed gap distribution (outlier detection), then the both-adjacent rule.
#       Drop, then re-est.map.
# est.map: error.prob=0.001, Haldane, per chromosome.
#
# Run: Rscript scripts/teonam_qtl_permap.R
# =============================================================================
suppressMessages({
  library(qtl)
  library(data.table)
  library(parallel)
  library(logger)
  if (requireNamespace("nilHMM", quietly = TRUE)) {
    library(nilHMM)
  } else {
    devtools::load_all("/Users/fvrodriguez/repos/nilhmm", quiet = TRUE)
  }
})
t0 <- Sys.time()
# logger: timestamped INFO/WARN/ERROR, sprintf-style messages ("%s"/"%d"/"%.1f")
log_layout(layout_glue_generator(format = '[{format(time, "%H:%M:%S")}] {level}: {msg}'))
log_formatter(formatter_sprintf)
log_threshold(INFO)
source("scripts/map_tools.R") # find_quirky_islands (isolated-cluster quirky finder)
FAMILIES <- c("W22TIL01", "W22TIL03", "W22TIL11", "W22TIL14", "W22TIL25")
GENO_DIR <- "data/teonam"
INFO <- fread(file.path(GENO_DIR, "map_v5_coe2008.tsv"))
N_CLUSTER <- min(10L, parallel::detectCores())
ISLAND_MAX_N <- 20L # quirky finder: max markers in an isolated cluster to flag
ISLAND_GAP_CM <- 2 # quirky finder: coarse isolation gap (cM) for clusters (99.99%ile gap ~1 cM)
dir.create("results/sim/teonam", showWarnings = FALSE, recursive = TRUE)

# ---- expected BC1S4 genotype freqs (AA,Aa,aa) from transition matrices (airmine)
AA <- matrix(c(1, 1 / 2, 0, 0, 1 / 2, 1, 0, 0, 0), 3, byrow = TRUE) # backcross to AA
S <- matrix(c(1, 1 / 4, 0, 0, 1 / 2, 0, 0, 1 / 4, 1), 3, byrow = TRUE) # selfing
bc1 <- AA %*% c(0, 1, 0) # F1 (het) backcrossed to W22 -> BC1
bc1s4 <- as.vector(S %*% S %*% S %*% S %*% bc1) # -> (0.734, 0.031, 0.234); codes 1/2/3 = AA/AB/BB
log_info("expected BC1S4 (AA,Aa,aa) = %s", paste(round(bc1s4, 4), collapse = ", "))

# ---- empirical-CDF -> qnorm z renormalization; upper-tail outliers z>1.96 [airmine]
# maps an arbitrary right-skewed statistic onto standard-normal quantiles and flags
# the ~upper 2.5% as distributional outliers -- distribution-free, RELATIVE to the
# observed spread (so a uniform shift, e.g. genome-wide het excess, is not flagged).
renorm_z <- function(x) {
  z <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (sum(ok) < 5L) {
    return(z)
  }
  d <- ecdf(x[ok])
  u <- suppressWarnings(predict(smooth.spline(x[ok], d(x[ok])), x[ok])$y)
  u <- pmin(pmax(u, 1e-6), 1 - 1 / sum(ok))
  z[ok] <- qnorm(u)
  z
}
is_outlier <- function(x) {
  z <- renorm_z(x)
  !is.na(z) & z > 1.96
}

# ---- bcsft(1,4) cross from a lines x marker integer matrix (1/2/3), v5-ordered ----
build_bcsft <- function(ord, Gr, ids) {
  cross <- list(geno = list())
  for (ch in sort(unique(ord$chr_v5))) {
    idx <- which(ord$chr_v5 == ch)
    dat <- Gr[, idx, drop = FALSE]
    colnames(dat) <- ord$marker[idx]
    mp <- ord$pos_v5[idx] / 1e6
    names(mp) <- ord$marker[idx] # placeholder; only ORDER used
    cross$geno[[as.character(ch)]] <- structure(list(data = dat, map = mp), class = "A")
  }
  cross$pheno <- data.frame(id = ids, stringsAsFactors = FALSE)
  class(cross) <- c("f2", "cross")
  convert2bcsft(cross, BC.gen = 1, F.gen = 4, estimate.map = FALSE) # Shannon 2012 BC1S4
}
run_estmap <- function(cross, cores = N_CLUSTER, tag = "") {
  chrs <- names(cross$geno)
  res <- parallel::mclapply(chrs, function(ch) {
    tc <- Sys.time()
    m <- est.map(subset(cross, chr = ch),
      error.prob = 0.001, map.function = "haldane",
      maxit = 10000, tol = 1e-6, n.cluster = 1
    )[[1]]
    log_info(
      "%s  chr %s done - %d markers, %.1f cM (%.1f min)", tag, ch, length(m),
      max(m) - min(m), as.numeric(difftime(Sys.time(), tc, units = "mins"))
    )
    m
  }, mc.cores = cores, mc.preschedule = FALSE)
  names(res) <- chrs
  bad <- vapply(res, function(x) inherits(x, "try-error") || is.null(x), logical(1))
  if (any(bad)) {
    log_error("est.map failed on chr: %s", paste(chrs[bad], collapse = ","))
    stop("est.map failed")
  }
  res
}
chr_len <- function(m) sapply(m, function(x) max(x) - min(x))

# ---- one family: distortion drop -> est.map -> quirky drop -> re-est.map ----------
per_family_map <- function(fam) {
  g <- fread(file.path(GENO_DIR, paste0(fam, "_genotype.csv")))
  # Genotype files are line x phenotype-rep long format: each RIL appears once per
  # `factor` level (col 3), the two rows exact genotype copies. Dedup to one row per
  # unique line ID (col 1) -> Chen Table 1 RIL counts. Exact dups leave est.map cM
  # unchanged (rf MLE is scale-invariant) but halve compute and fix XO/line counts.
  g <- g[!duplicated(g[[1]])]
  mk <- names(g)[-(1:3)]
  mk <- mk[mk %in% INFO$marker]
  ord <- INFO[match(mk, marker), .(marker, chr_v5, pos_v5)][order(chr_v5, pos_v5)]
  G <- as.matrix(g[, ord$marker, with = FALSE])
  storage.mode(G) <- "integer" # lines x markers, 0/1/2
  Gr <- G + 1L
  Gr[!(G %in% 0:2)] <- NA_integer_ # -> 1/2/3
  ids <- paste0(fam, ":", g[[1]])
  n0 <- nrow(ord)
  log_info("%s: start - %d markers x %d lines", fam, n0, length(ids))

  # (1) segregation distortion: per-marker chi vs expected BC1S4, outlier per chr
  cnt <- vapply(1:3, function(k) colSums(Gr == k, na.rm = TRUE), numeric(n0)) # markers x (AA,AB,BB)
  e <- rowSums(cnt) %o% bc1s4
  chi <- rowSums((cnt - e)^2 / e)
  dist_out <- logical(n0)
  for (ch in unique(ord$chr_v5)) {
    i <- which(ord$chr_v5 == ch)
    dist_out[i] <- is_outlier(chi[i])
  }
  keep1 <- !dist_out
  log_info(
    "%s: distortion drop - %d flagged, %d kept; ROUND 1 (PRELIMINARY) est.map, chr-parallel...",
    fam, sum(dist_out), sum(keep1)
  )

  # (2) PRELIMINARY map (round 1) -- used only to locate quirky/isolated markers
  s1 <- system.time(m1 <- run_estmap(build_bcsft(ord[keep1], Gr[, keep1, drop = FALSE], ids),
    tag = sprintf("%s R1", fam)
  ))["elapsed"]
  log_info("%s: ROUND 1 (PRELIMINARY) done - %d markers, %.1f cM (%.1f min)", fam, sum(keep1), sum(chr_len(m1)), s1 / 60)

  # (3) quirky: data-driven gap-outlier threshold, isolated-cluster rule
  # find_quirky_islands (scripts/map_tools.R) generalizes the old singleton
  # both-adjacent test to isolated small CLUSTERS (connected components).
  gaps_all <- unlist(lapply(m1, function(v) diff(v[order(v)])), use.names = FALSE)
  gap_out_thr <- {
    o <- is_outlier(gaps_all)
    if (any(o)) min(gaps_all[o]) else Inf # smallest OUTLIER gap
  }
  quirky <- unlist(
    lapply(m1, find_quirky,
      fine_thr = gap_out_thr,
      island_thr = ISLAND_GAP_CM, island_max_n = ISLAND_MAX_N
    ),
    use.names = FALSE
  )
  keep2 <- keep1 & !(ord$marker %in% quirky)
  log_info(
    "%s: quirky drop - gap-outlier thr %.2f cM, %d flagged, %d kept; ROUND 2 (REFINED, bad markers excluded) est.map...",
    fam, gap_out_thr, length(quirky), sum(keep2)
  )

  # (4) REFINED map (round 2) -- distortion + quirky markers excluded
  s2 <- system.time(m2 <- run_estmap(build_bcsft(ord[keep2], Gr[, keep2, drop = FALSE], ids),
    tag = sprintf("%s R2", fam)
  ))["elapsed"]
  cm2 <- unlist(lapply(m2, function(v) v - min(v)))
  names(cm2) <- unlist(lapply(m2, names))
  L2 <- chr_len(m2)
  log_info(
    "%s: ROUND 2 (REFINED) done - markers %d -> %d (dist -%d, quirky -%d) | total %.1f cM (%.1f min)",
    fam, n0, sum(keep2), sum(dist_out), length(quirky), sum(L2), s2 / 60
  )

  # QC warnings (Chen Table 1: per-subpop 1348-1596 cM; flag heavy drop or off-range length)
  drop_rate <- (n0 - sum(keep2)) / n0
  if (drop_rate > 0.10) {
    log_warn(
      "%s: distortion+quirky drop rate %.1f%% exceeds 10%% (%d of %d markers dropped)",
      fam, 100 * drop_rate, n0 - sum(keep2), n0
    )
  }
  if (sum(L2) < 1348 || sum(L2) > 1596) {
    log_warn(
      "%s: refined length %.1f cM outside Chen per-subpop range 1348-1596 cM",
      fam, sum(L2)
    )
  }
  list(
    long = data.table(
      family = fam, marker = names(cm2),
      chr_v5 = INFO$chr_v5[match(names(cm2), INFO$marker)],
      pos_v5 = INFO$pos_v5[match(names(cm2), INFO$marker)], cm = unname(cm2)
    ),
    summ = data.table(
      family = fam, n_lines = length(ids), n_markers_in = n0,
      n_markers_map = sum(keep2), dist_dropped = sum(dist_out),
      quirky_dropped = length(quirky), gap_outlier_thr = round(gap_out_thr, 2),
      total_cM = round(sum(L2), 1)
    )
  )
}

log_info(
  "=== per-family bcsft(1,4) maps (airmine pipeline): %d families, 2 est.map rounds each ===",
  length(FAMILIES)
)
rows <- vector("list", length(FAMILIES))
for (i in seq_along(FAMILIES)) {
  rows[[i]] <- per_family_map(FAMILIES[i])
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  avg <- el / i
  log_info(
    ">>> %d/%d families done | elapsed %.1f min | avg %.1f min/family | ETA ~%.1f min remaining",
    i, length(FAMILIES), el, avg, avg * (length(FAMILIES) - i)
  )
}
per_family <- rbindlist(lapply(rows, `[[`, "long"))
fam_tab <- rbindlist(lapply(rows, `[[`, "summ"))
fwrite(per_family, "results/sim/teonam/teonam_v5_native_perfam.csv")
fwrite(fam_tab, "results/sim/teonam/teonam_v5_native_family_summary.csv")

cat("\n==================== PER-FAMILY SUMMARY ====================\n")
print(fam_tab)
cat(sprintf(
  "\nOurs: avg %.0f cM (range %.0f-%.0f), avg %.0f markers\n",
  mean(fam_tab$total_cM), min(fam_tab$total_cM), max(fam_tab$total_cM),
  mean(fam_tab$n_markers_map)
))
cat("Chen: avg 1461 cM (range 1348-1596), ~13,733 markers/subpop\n")
cat(sprintf("Elapsed: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
