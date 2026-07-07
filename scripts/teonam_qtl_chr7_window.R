#!/usr/bin/env Rscript
# FAST local estimate of the chr7 flank gap: re-est.map only a cM window around the
# island (island removed via the generic detector), instead of the whole 154 cM
# chromosome. For a LOCAL gap this is ~identical to the full-chr est.map but runs in
# minutes. Reuses teonam_qtl_map.R's flanking-imputation -> bcsft(1,4) logic.
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
WIN <- c(80, 96) # cM window around the gap (~87.7)

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

mp <- fread("data/teonam/teonam_v5_native.tsv")[!is.na(cm)]
per <- split(mp, mp$chr_v5)
gaps_all <- unlist(lapply(per, function(d) {
  v <- sort(d$cm)
  diff(v)
}), use.names = FALSE)
gap_out_thr <- {
  o <- is_outlier(gaps_all)
  if (any(o)) min(gaps_all[o]) else Inf
}
isl7 <- find_quirky(setNames(per[["7"]]$cm, per[["7"]]$marker),
  fine_thr = gap_out_thr, island_thr = ISLAND_GAP_CM, island_max_n = ISLAND_MAX_N
)
win <- per[["7"]][cm >= WIN[1] & cm <= WIN[2]]
bp_lo <- min(win$pos_v5)
bp_hi <- max(win$pos_v5)
log_info(
  "window cM [%g,%g] = v5 %.2f-%.2f Mb | island removed: %d markers", WIN[1], WIN[2],
  bp_lo / 1e6, bp_hi / 1e6, length(isl7)
)

info <- fread(file.path(GENO_DIR, "map_v5_coe2008.tsv"))
setorder(info, chr_v5, pos_v5)
perfam <- fread(PERFAM)
kept_by_fam <- lapply(FAMILIES, function(f) perfam[family == f, marker])
names(kept_by_fam) <- FAMILIES
n_fam_marker <- table(perfam$marker)
union_mk <- info[chr_v5 == CH & pos_v5 >= bp_lo & pos_v5 <= bp_hi &
  marker %in% names(n_fam_marker) & !(marker %in% isl7)]
setorder(union_mk, chr_v5, pos_v5)
target <- data.frame(chr = union_mk$chr_v5, bp = union_mk$pos_v5)
rownames(target) <- union_mk$marker
log_info("window union: %d markers", nrow(union_mk))

impute_family <- function(fam) {
  g <- fread(file.path(GENO_DIR, paste0(fam, "_genotype.csv")))
  g <- g[!duplicated(g[[1]])]
  ids <- paste0(fam, ":", g[[1]])
  ok <- union_mk[marker %in% kept_by_fam[[fam]]]
  setorder(ok, chr_v5, pos_v5)
  G <- t(as.matrix(g[, ok$marker, with = FALSE]))
  storage.mode(G) <- "double"
  colnames(G) <- ids
  list(mat = interpolate_genotype(G,
    obs = data.frame(chr = ok$chr_v5, bp = ok$pos_v5), target = target,
    mode = "chen2019", coord = "bp"
  ), ids = ids)
}
build_bcsft <- function(ord_dt, Gr, ids) {
  cross <- list(geno = list())
  idx <- seq_len(nrow(ord_dt))
  dat <- Gr[, idx, drop = FALSE]
  colnames(dat) <- ord_dt$marker
  mpv <- ord_dt$pos_v5 / 1e6
  names(mpv) <- ord_dt$marker
  cross$geno[["7"]] <- structure(list(data = dat, map = mpv), class = "A")
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
log_info("est.map window (%d markers x %d lines)...", nrow(big), ncol(big))
mw <- est.map(subset(build_bcsft(ord, Gr, ids), chr = CH),
  error.prob = 0.001,
  map.function = "haldane", maxit = 10000, tol = 1e-6, n.cluster = 1
)[[1]]
mw <- mw - min(mw)
below <- union_mk[pos_v5 < 143.9e6][.N]
above <- union_mk[pos_v5 > 144.2e6][1]
log_info(
  "WINDOW len %.2f cM (%d markers) | FLANK GAP: %s (%.3f Mb) -> %s (%.3f Mb) = %.3f cM (was 10.22 WITH island)",
  max(mw), length(mw), below$marker, below$pos_v5 / 1e6, above$marker, above$pos_v5 / 1e6,
  mw[above$marker] - mw[below$marker]
)
out <- merge(data.table(marker = names(mw), cm_window = as.numeric(mw)), union_mk[, .(marker, pos_v5)], by = "marker")[order(cm_window)]
fwrite(out, "results/sim/teonam/teonam_v5_native_chr7_window.csv")
cat("wrote results/sim/teonam/teonam_v5_native_chr7_window.csv\n")
