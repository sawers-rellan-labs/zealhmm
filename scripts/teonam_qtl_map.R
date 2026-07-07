#!/usr/bin/env Rscript
# =============================================================================
# Native TeoNAM COMPOSITE genetic map on B73 v5 marker order (R/qtl) -- the
# DELIVERABLE: a native cM per marker (`cm`) to replace the borrowed Ed Coe
# consensus cM in data/teonam/map_v5_coe2008.tsv.
# Plan: agent/teonam-v5-genetic-map-plan.md  Handover: agent/teonam-map-handover.md
# -----------------------------------------------------------------------------
# FAITHFUL to Chen 2019 (Methods, "Genetic map construction and marker
# imputation", line 76). Depends on the per-family run (teonam_qtl_permap.R):
#
#   1. UNION = union of the per-family KEPT (post distortion+quirky filter)
#      markers, read from results/sim/teonam/teonam_v5_native_perfam.csv -- Chen's
#      "51,544 high-quality SNPs" analog, NOT the raw union.
#   2. FLANKING-MARKER IMPUTATION exactly per Chen: "If the flanking markers had
#      same genotypes, the missing genotype was imputed as the same with flanking
#      markers, or otherwise left as missing." Implemented by nilHMM's
#      interpolate_genotype(mode = "chen2019", coord = "bp") -- concordant flanks
#      fill, discordant flanks or chromosome ends -> NA, order-only (distance is
#      NOT used for the fill decision). This is the SAME code path we use for the
#      Tian 2011 JLM/GWAS densification (continuous/step modes) -- one function,
#      Chen's map rule is just its own mode. Result is partially-missing (correct;
#      est.map's HMM handles NA).
#   3. ONE bcsft(1,4) est.map per chromosome (error.prob=0.001, Haldane) on the
#      imputed union -> PRELIMINARY composite map.
#   4. JOINT QUIRKY pass (same data-driven gap-outlier rule as teonam_qtl_permap.R)
#      on the preliminary map; drop quirky markers; re-est.map -> REFINED composite
#      map = the final cm.
#   5. Write data/teonam/teonam_v5_native.tsv (+ results/sim/teonam/teonam_v5_native_qc.csv).
#
# NON-DESTRUCTIVE: does not overwrite map_v5_coe2008.tsv; touches no sweep/JLM/
# notebook/calibrate code. Does not commit.
#
# Run:  Rscript scripts/teonam_qtl_map.R   (after teonam_qtl_permap.R has finished)
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
if (!"chen2019" %in% eval(formals(interpolate_genotype)$mode)) {
  stop("nilHMM lacks the 'chen2019' interpolate_genotype mode -- reinstall from repos/nilhmm")
}
source("scripts/map_tools.R") # find_quirky_islands (isolated-cluster quirky finder)

t0 <- Sys.time()
log_layout(layout_glue_generator(format = '[{format(time, "%H:%M:%S")}] {level}: {msg}'))
log_formatter(formatter_sprintf)
log_threshold(INFO)

FAMILIES <- c("W22TIL01", "W22TIL03", "W22TIL11", "W22TIL14", "W22TIL25")
GENO_DIR <- "data/teonam"
INFO_PATH <- file.path(GENO_DIR, "map_v5_coe2008.tsv")
PERFAM_CSV <- "results/sim/teonam/teonam_v5_native_perfam.csv"
N_CLUSTER <- min(10L, parallel::detectCores())
ISLAND_MAX_N <- 20L # quirky finder: max markers in an isolated cluster to flag
ISLAND_GAP_CM <- 2 # quirky finder: coarse isolation gap (cM) for clusters (99.99%ile gap ~1 cM)
dir.create("results/sim/teonam", showWarnings = FALSE, recursive = TRUE)

# ---- expected BC1S4 genotype freqs (AA,Aa,aa) from transition matrices ------
AA <- matrix(c(1, 1 / 2, 0, 0, 1 / 2, 1, 0, 0, 0), 3, byrow = TRUE)
S <- matrix(c(1, 1 / 4, 0, 0, 1 / 2, 0, 0, 1 / 4, 1), 3, byrow = TRUE)
bc1s4 <- as.vector(S %*% S %*% S %*% S %*% (AA %*% c(0, 1, 0))) # (0.734,0.031,0.234)

# ---- empirical-CDF -> qnorm z renormalization; upper-tail outliers z>1.96 ----
# (identical rule to teonam_qtl_permap.R -- flags the ~upper 2.5% of a right-skewed
# statistic RELATIVE to its observed spread; distribution-free.)
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

# ---- marker annotation, ordered on v5 --------------------------------------
info <- fread(INFO_PATH)
setnames(info, "cm", "cm_coe2008") # consensus (Ed Coe 2008) cM; native est.map cM is `cm`
setorder(info, chr_v5, pos_v5)
info[, rank_v5 := seq_len(.N)]
CHRS <- sort(unique(info$chr_v5))
log_info("marker_info: %d markers, %d chromosomes (v5)", nrow(info), length(CHRS))

# QC watchlist: v2->v5 chromosome changers & local rank inversions (non-dropping)
info[, chr_change := as.integer(chr_v2 != chr_v5)]
info[, inversion := 0L]
for (ch in CHRS) {
  ix <- info[chr_v5 == ch, which = TRUE] # already ordered by pos_v5
  sub <- info[ix]
  co <- sub$chr_v2 == ch
  p2 <- sub$pos_v2
  inv <- rep(0L, nrow(sub))
  if (nrow(sub) > 1L) {
    for (i in 2:nrow(sub)) if (co[i] && co[i - 1] && p2[i] < p2[i - 1]) inv[i] <- 1L
  }
  info[ix, inversion := inv]
}
log_info(
  "QC watchlist: %d chr-changers, %d local inversions",
  sum(info$chr_change), sum(info$inversion)
)

# ---- UNION = per-family KEPT markers (post distortion+quirky filter) --------
if (!file.exists(PERFAM_CSV)) {
  stop("missing ", PERFAM_CSV, " -- run scripts/teonam_qtl_permap.R first")
}
perfam <- fread(PERFAM_CSV)
kept_by_fam <- lapply(FAMILIES, function(f) perfam[family == f, marker])
names(kept_by_fam) <- FAMILIES
n_fam_marker <- table(perfam$marker) # families keeping each marker
union_mk <- info[marker %in% names(n_fam_marker)]
setorder(union_mk, chr_v5, pos_v5) # target must be (chr, bp)-sorted
union_mk[, n_fam := as.integer(n_fam_marker[marker])]
log_info(
  "UNION (per-family kept): %d markers | shared(>=2 fam) %d | private(1 fam) %d",
  nrow(union_mk), sum(union_mk$n_fam >= 2L), sum(union_mk$n_fam == 1L)
)

# target grid = the v5-ordered union; pos_v5 is strictly increasing within chr
# (verified: 0 tied positions), so coord="bp" is well-posed for interpolate_genotype.
target <- data.frame(chr = union_mk$chr_v5, bp = union_mk$pos_v5)
rownames(target) <- union_mk$marker

# ---- Chen 2019 flanking imputation of one family onto the union ------------
# obs = family's KEPT markers (complete genotypes); target = full union.
# interpolate_genotype(mode="chen2019") fills a gap union marker from the family's
# nearest kept flanks iff they agree, else NA; chromosome ends -> NA. Union markers
# the family itself kept fall on an obs bp -> returned as the observed genotype.
impute_family <- function(fam) {
  g <- fread(file.path(GENO_DIR, paste0(fam, "_genotype.csv")))
  g <- g[!duplicated(g[[1]])] # line x phenotype-rep long format -> one row per RIL (exact dups)
  ids <- paste0(fam, ":", g[[1]])
  ok <- union_mk[marker %in% kept_by_fam[[fam]]] # kept markers (subset of union)
  setorder(ok, chr_v5, pos_v5)
  G <- t(as.matrix(g[, ok$marker, with = FALSE])) # markers x lines, dosage 0/1/2
  storage.mode(G) <- "double"
  colnames(G) <- ids
  if (anyNA(G)) stop(fam, ": kept-marker genotypes contain NA (interpolate needs complete)")
  obs <- data.frame(chr = ok$chr_v5, bp = ok$pos_v5)
  out <- interpolate_genotype(G,
    obs = obs, target = target,
    mode = "chen2019", coord = "bp"
  ) # union x lines, NA where discordant/end
  filled <- nrow(out) - nrow(ok)
  log_info(
    "%s: kept %d markers -> union %d rows x %d lines | %.1f%% NA in gap rows",
    fam, nrow(ok), nrow(out), ncol(out),
    100 * mean(is.na(out[!(union_mk$marker %in% ok$marker), , drop = FALSE]))
  )
  list(mat = out, ids = ids)
}

# ---- bcsft(1,4) cross from a lines x marker integer matrix (1/2/3, NA ok) ----
build_bcsft <- function(ord_dt, Gr, ids) {
  cross <- list(geno = list())
  for (ch in sort(unique(ord_dt$chr_v5))) {
    idx <- which(ord_dt$chr_v5 == ch)
    dat <- Gr[, idx, drop = FALSE]
    colnames(dat) <- ord_dt$marker[idx]
    mp <- ord_dt$pos_v5[idx] / 1e6
    names(mp) <- ord_dt$marker[idx] # placeholder; only ORDER used
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

# =============================================================================
# Build the imputed union, estimate the composite map, joint-quirky refine.
# =============================================================================
log_info("=== COMPOSITE map: Chen-2019 flanking imputation -> bcsft(1,4) est.map ===")
blocks <- lapply(FAMILIES, impute_family)
big <- do.call(cbind, lapply(blocks, `[[`, "mat")) # union x all lines (0/1/2, NA)
ids_comp <- unlist(lapply(blocks, `[[`, "ids"))
stopifnot(nrow(big) == nrow(union_mk), all(rownames(big) == union_mk$marker))
rm(blocks)
gc()
log_info(
  "imputed union: %d markers x %d lines | %.1f%% NA (discordant/chr-end -> HMM handles)",
  nrow(big), ncol(big), 100 * mean(is.na(big))
)

# pooled genotype-class counts (from the imputed union) for the seg-distortion QC flag
cnt <- vapply(0:2, function(k) rowSums(big == k, na.rm = TRUE), numeric(nrow(big))) # markers x (AA,AB,BB)

# recode 0/1/2 -> 1/2/3, ONE bcsft(1,4) cross, PRELIMINARY est.map per chromosome
Gr <- t(big)
rm(big)
gc()
storage.mode(Gr) <- "integer"
Gr <- Gr + 1L
ord_u <- union_mk[, .(marker, chr_v5, pos_v5)]
log_info("ROUND 1 (PRELIMINARY) bcsft(1,4) est.map per chromosome (the heavy step)...")
s1 <- system.time(m1 <- run_estmap(build_bcsft(ord_u, Gr, ids_comp), tag = "R1"))["elapsed"]
log_info(
  "ROUND 1 (PRELIMINARY) done - %d markers, %.1f cM (%.1f min)",
  nrow(union_mk), sum(chr_len(m1)), s1 / 60
)

# ---- JOINT QUIRKY pass: data-driven gap-outlier threshold, isolated-cluster rule
# find_quirky_islands (scripts/map_tools.R) generalizes the old singleton
# both-adjacent test to isolated small CLUSTERS (connected components), so an
# 8-marker v2->v5 displaced block cut off by outlier gaps on both ends is caught,
# not just single isolated markers.
gaps_all <- unlist(lapply(m1, function(v) diff(v[order(v)])), use.names = FALSE)
gap_out_thr <- {
  o <- is_outlier(gaps_all)
  if (any(o)) min(gaps_all[o]) else Inf
}
quirky <- unlist(
  lapply(m1, find_quirky,
    fine_thr = gap_out_thr,
    island_thr = ISLAND_GAP_CM, island_max_n = ISLAND_MAX_N
  ),
  use.names = FALSE
)
keep2 <- !(union_mk$marker %in% quirky)
log_info(
  "joint quirky drop - gap-outlier thr %.2f cM, %d flagged, %d kept; ROUND 2 (REFINED) est.map...",
  gap_out_thr, length(quirky), sum(keep2)
)

# ---- ROUND 2 (REFINED) est.map on the quirky-filtered union -----------------
s2 <- system.time(m2 <- run_estmap(build_bcsft(ord_u[keep2], Gr[, keep2, drop = FALSE], ids_comp), tag = "R2"))["elapsed"]
rm(Gr)
gc()
cm2 <- unlist(lapply(m2, function(v) v - min(v)))
names(cm2) <- unlist(lapply(m2, names))
Lc <- setNames(chr_len(m2)[as.character(CHRS)], as.character(CHRS))
union_mk[, cm := unname(cm2[union_mk$marker])]
log_info("ROUND 2 (REFINED) done - %d markers, %.1f cM (%.1f min)", sum(keep2), sum(Lc), s2 / 60)
if (sum(Lc) < 1400 || sum(Lc) > 1700) {
  log_warn("composite length %.1f cM far from Chen's 1540 cM", sum(Lc))
}

# ---- QC: pooled segregation distortion (vs BC1S4 expected) + large gaps -----
Nm <- rowSums(cnt)
Exp <- Nm %o% bc1s4
pv <- pchisq(rowSums((cnt - Exp)^2 / Exp), df = 2, lower.tail = FALSE)
bonf <- 0.05 / nrow(cnt)
union_mk[, seg_distort := as.integer(is.finite(pv) & pv < bonf)]
union_mk[, quirky_drop := as.integer(marker %in% quirky)]
union_mk[, gap_prev := c(NA_real_, diff(cm)), by = chr_v5]
gap_thr_report <- if (is.finite(gap_out_thr)) gap_out_thr else 30
union_mk[, big_gap := as.integer(!is.na(gap_prev) & gap_prev > gap_thr_report)]
log_info(
  "QC: %d seg-distorted (P<%.1e), %d gaps > %.2f cM, %d quirky-dropped",
  sum(union_mk$seg_distort), bonf, sum(union_mk$big_gap, na.rm = TRUE),
  gap_thr_report, length(quirky)
)

# ---- outputs ----------------------------------------------------------------
# Output schema: `cm` = native est.map cM; `cm_coe2008` = Ed Coe consensus cM.
qc <- merge(info[, .(marker, chr_v2, pos_v2, chr_v5, pos_v5, cm_coe2008, rank_v5, chr_change, inversion)],
  union_mk[, .(marker, n_fam, cm, gap_prev, big_gap, seg_distort, quirky_drop)],
  by = "marker", all.x = TRUE
)
setorder(qc, chr_v5, pos_v5)
info_out <- copy(info)[, cm := union_mk[match(info$marker, marker), cm]]
info_out <- info_out[, .(marker, chr_v2, pos_v2, chr_v5, pos_v5, cm, cm_coe2008)]
fwrite(qc, "results/sim/teonam/teonam_v5_native_qc.csv")
fwrite(info_out, "data/teonam/teonam_v5_native.tsv", sep = "\t")

sp <- cor(info_out$cm, info_out$cm_coe2008, method = "spearman", use = "complete.obs")
cat("\n==================== COMPOSITE SUMMARY ====================\n")
cat(sprintf("Native composite map: %.1f cM  (vs Chen 1540 / Ed Coe consensus ~1781)\n", sum(Lc)))
print(data.table(chr = names(Lc), cM = round(Lc, 2)))
cat(sprintf(
  "union markers: %d | with cm (post-quirky): %d | quirky-dropped: %d\n",
  nrow(union_mk), sum(!is.na(union_mk$cm)), length(quirky)
))
cat(sprintf(
  "shared (>=2 fam): %d | private (1 fam): %d\n",
  sum(union_mk$n_fam >= 2L), sum(union_mk$n_fam == 1L)
))
cat(sprintf("Spearman cor(cm [native], cm_coe2008 [consensus]): %.4f\n", sp))
cat(sprintf(
  "QC: chr-changers=%d, inversions=%d, seg-distorted=%d, big-gaps(>%.2f)=%d\n",
  sum(qc$chr_change), sum(qc$inversion), sum(qc$seg_distort, na.rm = TRUE),
  gap_thr_report, sum(qc$big_gap, na.rm = TRUE)
))
cat("Chen composite: 51,544 SNPs, 1540 cM\n")
cat("Wrote: data/teonam/teonam_v5_native.tsv, results/sim/teonam/teonam_v5_native_qc.csv\n")
cat(sprintf("Elapsed: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
