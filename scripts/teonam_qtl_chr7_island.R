#!/usr/bin/env Rscript
# Recalculate LG7 with the isolated island removed, using the GENERIC cluster
# detector (find_quirky_islands, scripts/map_tools.R) rather than hardcoded marker
# names. Reuses teonam_qtl_map.R's flanking-imputation -> bcsft(1,4) est.map logic,
# restricted to chr7. NON-DESTRUCTIVE: writes only a chr7 comparison CSV.
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
setwd("/Users/fvrodriguez/repos/zealhmm")
source("scripts/map_tools.R")
log_layout(layout_glue_generator(format = '[{format(time,"%H:%M:%S")}] {level}: {msg}'))
log_formatter(formatter_sprintf)
log_threshold(INFO)

FAMILIES <- c("W22TIL01", "W22TIL03", "W22TIL11", "W22TIL14", "W22TIL25")
GENO_DIR <- "data/teonam"
PERFAM <- "results/sim/teonam/teonam_v5_native_perfam.csv"
ISLAND_MAX_N <- 20L
ISLAND_GAP_CM <- 2
CH <- 7L

# empirical-CDF -> qnorm z; upper-tail outlier gaps (verbatim from the pipeline)
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

# ---- detect islands GENERICALLY on the existing composite map ---------------
mp <- fread("data/teonam/teonam_v5_native.tsv")[!is.na(cm)]
gaps_all <- unlist(lapply(split(mp, mp$chr_v5), function(d) {
  v <- sort(d$cm)
  diff(v)
}), use.names = FALSE)
gap_out_thr <- {
  o <- is_outlier(gaps_all)
  if (any(o)) min(gaps_all[o]) else Inf
}
per_chr <- split(mp, mp$chr_v5)
islands_all <- unlist(lapply(per_chr, function(d) {
  find_quirky(setNames(d$cm, d$marker),
    fine_thr = gap_out_thr,
    island_thr = ISLAND_GAP_CM, island_max_n = ISLAND_MAX_N
  )
}), use.names = FALSE)
d7 <- per_chr[["7"]]
isl7 <- find_quirky(setNames(d7$cm, d7$marker),
  fine_thr = gap_out_thr,
  island_thr = ISLAND_GAP_CM, island_max_n = ISLAND_MAX_N
)
log_info(
  "gap-outlier thr = %.2f cM | islands genome-wide: %d markers | chr7 island: %d markers",
  gap_out_thr, length(islands_all), length(isl7)
)
log_info("chr7 island markers: %s", paste(isl7, collapse = ", "))

# ---- rebuild chr7 union WITHOUT the detected island, re-est.map -------------
info <- fread(file.path(GENO_DIR, "map_v5_coe2008.tsv"))
setorder(info, chr_v5, pos_v5)
perfam <- fread(PERFAM)
kept_by_fam <- lapply(FAMILIES, function(f) perfam[family == f, marker])
names(kept_by_fam) <- FAMILIES
n_fam_marker <- table(perfam$marker)
union_mk <- info[marker %in% names(n_fam_marker) & chr_v5 == CH & !(marker %in% isl7)]
setorder(union_mk, chr_v5, pos_v5)
target <- data.frame(chr = union_mk$chr_v5, bp = union_mk$pos_v5)
rownames(target) <- union_mk$marker

impute_family <- function(fam) {
  g <- fread(file.path(GENO_DIR, paste0(fam, "_genotype.csv")))
  g <- g[!duplicated(g[[1]])]
  ids <- paste0(fam, ":", g[[1]])
  ok <- union_mk[marker %in% kept_by_fam[[fam]]]
  setorder(ok, chr_v5, pos_v5)
  G <- t(as.matrix(g[, ok$marker, with = FALSE]))
  storage.mode(G) <- "double"
  colnames(G) <- ids
  if (anyNA(G)) stop(fam, ": kept genotypes contain NA")
  list(mat = interpolate_genotype(G,
    obs = data.frame(chr = ok$chr_v5, bp = ok$pos_v5), target = target,
    mode = "chen2019", coord = "bp"
  ), ids = ids)
}
build_bcsft <- function(ord_dt, Gr, ids) {
  cross <- list(geno = list())
  for (ch in sort(unique(ord_dt$chr_v5))) {
    idx <- which(ord_dt$chr_v5 == ch)
    dat <- Gr[, idx, drop = FALSE]
    colnames(dat) <- ord_dt$marker[idx]
    mpv <- ord_dt$pos_v5[idx] / 1e6
    names(mpv) <- ord_dt$marker[idx]
    cross$geno[[as.character(ch)]] <- structure(list(data = dat, map = mpv), class = "A")
  }
  cross$pheno <- data.frame(id = ids, stringsAsFactors = FALSE)
  class(cross) <- c("f2", "cross")
  convert2bcsft(cross, BC.gen = 1, F.gen = 4, estimate.map = FALSE)
}
blocks <- lapply(FAMILIES, impute_family)
big <- do.call(cbind, lapply(blocks, `[[`, "mat"))
ids <- unlist(lapply(blocks, `[[`, "ids"))
Gr <- t(big) + 1L
storage.mode(Gr) <- "integer"
ord <- union_mk[, .(marker, chr_v5, pos_v5)]
log_info("est.map chr7 (island removed, %d markers x %d lines)...", nrow(big), ncol(big))
m7 <- est.map(subset(build_bcsft(ord, Gr, ids), chr = CH),
  error.prob = 0.001,
  map.function = "haldane", maxit = 10000, tol = 1e-6, n.cluster = 1
)[[1]]
m7 <- m7 - min(m7)
log_info("chr7 NEW length: %.2f cM (%d markers) | was 154.5 cM WITH island", max(m7), length(m7))

below <- union_mk[pos_v5 < 143.9e6][.N]
above <- union_mk[pos_v5 > 144.2e6][1]
log_info(
  "FLANK GAP now: %s (v5 %.3f) -> %s (v5 %.3f) = %.3f cM  (was 10.22 cM WITH island)",
  below$marker, below$pos_v5 / 1e6, above$marker, above$pos_v5 / 1e6,
  m7[above$marker] - m7[below$marker]
)
out <- merge(data.table(marker = names(m7), cm_noisland = as.numeric(m7)),
  union_mk[, .(marker, pos_v5)],
  by = "marker"
)[order(cm_noisland)]
fwrite(out, "results/sim/teonam/teonam_v5_native_dropped.csv")
cat("wrote results/sim/teonam/teonam_v5_native_dropped.csv\n")
