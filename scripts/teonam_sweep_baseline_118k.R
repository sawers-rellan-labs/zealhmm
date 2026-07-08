#!/usr/bin/env Rscript
# n=infinity (perfect-caller) baseline for the 118K coverage sweeps.
# The lambda=Inf reference must be the GWAS on the COMPLETE truth matrix -- the
# dense polarized 118K truth genotypes interpolated onto the union grid (exactly
# what an infinite-coverage caller recovers), NOT the authentic per-marker scan
# (which carries per-marker NA/scatter and breaks the geom_line geometry against
# the smooth block-interpolated degraded curves). Mirrors the 51K design, where
# the sweeps' lambda=Inf baseline was the interpolated-truth scan.
#
# Caller-INDEPENDENT (it is just the truth), so computed once here and shared.
# Output: data/teonam/stam_gwas_scan_118k_complete_baseline.csv (SNP,CHR,BP,P)
# Side effect: patches the coverage==Inf rows of the four
#   results/sim/teonam/stam_gwas_<caller>_118k_sweep.csv files in place.
# Run: Rscript scripts/teonam_sweep_baseline_118k.R
suppressMessages({
  library(data.table)
  library(parallel)
  library(readxl)
  devtools::load_all("/Users/fvrodriguez/repos/nilhmm", quiet = TRUE)
})
ROOT <- "/Users/fvrodriguez/repos/zealhmm"
setwd(ROOT)
source(file.path(ROOT, "scripts/logging.R"))
OUTDIR <- file.path(ROOT, "results/sim/teonam")
THREADS <- max(1L, detectCores() - 2L)
t0 <- Sys.time()

# --- 118K cM grid + dense polarized truth (identical to the sweeps) -----------
mc <- fread(file.path(ROOT, "data/teonam/markers_v5_gwas118k_cm.tsv"))
setnames(mc, "pos_v5", "pos")
setorder(mc, chr, cm)
mt_all <- mc[, .SD[!duplicated(cm)], by = chr]
setorder(mt_all, chr, cm)
u <- copy(mc)
setorder(u, chr, cm)
target_df <- data.frame(chr = as.integer(u$chr), cm = as.numeric(u$cm))
union_markers <- u$marker
union_pos <- as.integer(u$pos)
union_chr <- as.integer(u$chr)

g118 <- readRDS(file.path(ROOT, "data/teonam/teonam_gwas118k_dosage_fsfhap.rds"))
dos <- g118$dos
FAMS <- c("TIL01", "TIL03", "TIL11", "TIL14", "TIL25")

load_family <- function(fam) {
  keys <- colnames(dos)[substr(colnames(dos), 1, 5) == fam]
  D <- dos[mt_all$marker, keys, drop = FALSE]
  storage.mode(D) <- "double"
  if (anyNA(D)) {
    rm <- apply(D, 1, function(z) {
      z <- z[!is.na(z)]
      if (!length(z)) 0 else as.numeric(names(which.max(table(z))))
    })
    for (j in seq_len(ncol(D))) {
      na <- is.na(D[, j])
      if (any(na)) D[na, j] <- rm[na]
    }
  }
  list(mt = mt_all, D = D, keys = keys)
}
log_info("assembling complete truth matrix ...")
Gtruth <- do.call(cbind, lapply(FAMS, function(fam) {
  fd <- load_family(fam)
  b <- interpolate_genotype(fd$D, data.frame(chr = fd$mt$chr, cm = fd$mt$cm), target_df, mode = "step")
  colnames(b) <- fd$keys
  b
}))
rownames(Gtruth) <- union_markers

# --- STAM ~ Family + marker scan (identical model to the sweeps) --------------
ph <- as.data.frame(read_excel(file.path(ROOT, "data/teonam/9250682/TeoNAM_1257RILs_22traits_phenotype_data.xlsx")))
names(ph)[1] <- "line"
stam_by <- setNames(ph$STAM, ph$line)
y <- suppressWarnings(as.numeric(stam_by[colnames(Gtruth)]))
fam <- factor(substr(colnames(Gtruth), 1, 5))
ok <- !is.na(y)
y <- y[ok]
fam <- droplevels(fam[ok])
Gm <- Gtruth[, ok, drop = FALSE]
Xr <- model.matrix(~fam)
n <- length(y)
RSS0 <- sum(lm.fit(Xr, y)$residuals^2)
P <- unlist(mclapply(seq_len(nrow(Gm)), function(i) {
  g <- Gm[i, ]
  if (sd(g) == 0) {
    return(NA_real_)
  }
  fit <- lm.fit(cbind(Xr, g), y)
  RSS1 <- sum(fit$residuals^2)
  df2 <- n - fit$rank
  pf(((RSS0 - RSS1) / 1) / (RSS1 / df2), 1, df2, lower.tail = FALSE)
}, mc.cores = THREADS))
baseline <- data.table(SNP = union_markers, CHR = union_chr, BP = union_pos, P = P)[order(CHR, BP)]

BL_PATH <- file.path(ROOT, "data/teonam/stam_gwas_scan_118k_complete_baseline.csv")
fwrite(baseline, BL_PATH)
tb1 <- 272330564L
w <- baseline[CHR == 1 & abs(BP - tb1) <= 5e5 & is.finite(P) & P > 0]
log_info(
  "complete-truth baseline: %d markers, tb1 peak -log10P = %.2f -> %s",
  nrow(baseline), max(-log10(w$P)), BL_PATH
)

# --- patch the coverage==Inf rows of each sweep CSV in place ------------------
bl_inf <- copy(baseline)[, coverage := Inf]
callers <- c("rtiger", "nnil", "lbimpute", "control")
for (ci in seq_along(callers)) {
  caller <- callers[ci]
  f <- file.path(OUTDIR, sprintf("stam_gwas_%s_118k_sweep.csv", caller))
  if (!file.exists(f)) {
    log_info("%s", paste0("  (skip, not found: ", basename(f), ")"))
    next
  }
  sw <- fread(f)
  sw <- sw[is.finite(coverage)] # drop old (scattery) Inf rows
  sw <- rbindlist(list(sw, bl_inf), use.names = TRUE)
  fwrite(sw, f)
  log_info("  patched %s (%d rows)", basename(f), nrow(sw))
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  log_info(
    ">>> %d/%d done | elapsed %.1f min | avg %.1f min | ETA ~%.1f min remaining",
    ci, length(callers), el, el / ci, (el / ci) * (length(callers) - ci)
  )
}
