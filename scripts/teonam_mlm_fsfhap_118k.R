#!/usr/bin/env Rscript
# =============================================================================
# STAM MLM (FAMILY + K, P3D, joint 2-df additive+dominance F) on the FSFHap
# ANCESTRY MOSAIC (data/teonam/teonam_gwas118k_dosage_fsfhap.rds) -- the
# ancestry-imputed complete-truth panel. This is the MLM analog of the OLS
# ancestry-imputed panel (teonam-qtl-recovery-118k.qmd panel C): same clean
# W22<->teosinte block mosaic, scanned with the OLS model's FIXED part (family factor,
# 5 levels) PLUS a K random effect -- i.e. exactly "OLS (Family + marker) + K", so the
# OLS<->MLM contrast is the K term alone. Family is used instead of 5 PCs because in
# TeoNAM the structure IS the 5 families (Q approx Family; verified interchangeable in
# scripts/teonam_mlm_family_118k.R: candidate peaks within 0.2, lambda_GC 1.09 vs 1.12).
#
# K is estimated from the FSFHap genotypes THEMSELVES (full data, no read simulation).
# R-EMMAX (VanRaden K, REML delta, whitening) -- validated vs TASSEL MLM (scripts/
# emmax_qk.R header). NA cells are per-marker mean-imputed (mosaic is near-complete).
#
# Output: data/teonam/stam_gwas_mlm_fsfhap_118k.csv  (SNP, CHR, BP[v5], P)
# Run: Rscript scripts/teonam_mlm_fsfhap_118k.R
# =============================================================================
suppressMessages({
  library(data.table)
  library(readxl)
})
ROOT <- "/Users/fvrodriguez/repos/zealhmm"
setwd(ROOT)
source(file.path(ROOT, "scripts/logging.R"))
source(file.path(ROOT, "scripts/emmax_qk.R")) # emmax_qk_scan + mem headroom

# --- FSFHap ancestry mosaic + v5 roster --------------------------------------
g <- readRDS("data/teonam/teonam_gwas118k_dosage_fsfhap.rds")
dos <- g$dos # markers x lines, 0/1/2 (0=W22, 2=teo), some NA
mc <- fread("data/teonam/markers_v5_gwas118k.tsv") # marker, chr_v5, pos_v5
chr_by <- setNames(mc$chr_v5, mc$marker)
pos_by <- setNames(mc$pos_v5, mc$marker)
mk <- intersect(rownames(dos), mc$marker) # lifted markers only
G <- dos[mk, , drop = FALSE]
CHR <- as.integer(chr_by[mk])
BP <- as.integer(pos_by[mk])
o <- order(CHR, BP)
G <- G[o, , drop = FALSE]
CHR <- CHR[o]
BP <- BP[o]
mk <- mk[o]

# --- phenotype (STAM) ---------------------------------------------------------
ph <- as.data.frame(read_excel("data/teonam/9250682/TeoNAM_1257RILs_22traits_phenotype_data.xlsx"))
names(ph)[1] <- "line"
TRAIT <- toupper(Sys.getenv("TRAIT", "STAM"))
TTAG <- tolower(TRAIT) # phenotype col; STAM default
if (!TRAIT %in% names(ph)) stop("TRAIT '", TRAIT, "' is not a phenotype column")
stam <- suppressWarnings(as.numeric(ph[[TRAIT]]))
names(stam) <- ph$line
lines <- intersect(colnames(G), names(stam)[is.finite(stam)])
G <- G[, lines, drop = FALSE]
storage.mode(G) <- "double"

# per-marker mean-impute NA (near-complete mosaic; keeps K/Q/scan defined)
if (anyNA(G)) {
  rmn <- rowMeans(G, na.rm = TRUE)
  rmn[!is.finite(rmn)] <- 0
  na <- which(is.na(G), arr.ind = TRUE)
  G[na] <- rmn[na[, 1]]
}
log_info("FSFHap panel: %d markers x %d lines (STAM)", nrow(G), ncol(G))

# --- build the full-data Q+K null (R-EMMAX), then scan ------------------------
vanraden <- function(Ms) { # Ms = lines x markers
  pr <- colMeans(Ms) / 2
  Z <- sweep(Ms, 2, 2 * pr, "-")
  tcrossprod(Z) / (2 * sum(pr * (1 - pr)))
}
M <- t(G) # lines x markers
y <- stam[lines]
n <- length(y)
vpc <- (colMeans(M^2) - colMeans(M)^2) > 1e-12
fam <- factor(substr(lines, 1, 5))
X <- model.matrix(~fam) # fixed = family factor (the OLS/JLM structure), not 5 PCs
p <- ncol(X)
K <- vanraden(M[, vpc, drop = FALSE])
eig <- eigen(K, symmetric = TRUE)
U <- eig$vectors
xi <- pmax(eig$values, 1e-8)
yt <- as.numeric(crossprod(U, y))
Xt <- crossprod(U, X)
nll <- function(ld) {
  w <- 1 / (xi + exp(ld))
  Xs <- Xt * sqrt(w)
  ys <- yt * sqrt(w)
  XtX <- crossprod(Xs)
  b <- solve(XtX, crossprod(Xs, ys))
  0.5 * ((n - p) * log(sum((ys - Xs %*% b)^2) / (n - p)) + sum(log(xi + exp(ld))) + as.numeric(determinant(XtX, log = TRUE)$modulus))
}
delta <- exp(optimize(nll, c(-10, 10))$minimum)
ws <- sqrt(1 / (xi + delta))
Xs <- Xt * ws
XtXinv <- solve(crossprod(Xs))
ry <- as.numeric((yt * ws) - Xs %*% (XtXinv %*% crossprod(Xs, yt * ws)))
null <- list(lines = lines, U = U, ws = ws, Xs = Xs, XtXinv = XtXinv, ry = ry, n = n, p = p, delta = delta)
log_info("  Family+K null: n=%d, p=%d (family levels + int), delta=%.3g", n, p, delta)

scan <- emmax_qk_scan(G, null, CHR, BP)[order(CHR, BP)]
lgc <- with(scan[is.finite(P) & P > 0], median(qchisq(P, 1, lower.tail = FALSE)) / qchisq(0.5, 1))
fwrite(scan, sprintf("data/teonam/%s_gwas_mlm_fsfhap_118k.csv", TTAG))
log_info(
  "wrote data/teonam/stam_gwas_mlm_fsfhap_118k.csv  (%d markers, lambda_GC=%.3f, max -log10P=%.2f)",
  nrow(scan), lgc, max(-log10(scan[is.finite(P) & P > 0, P]))
)
