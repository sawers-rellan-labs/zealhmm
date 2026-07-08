#!/usr/bin/env Rscript
# =============================================================================
# STAM MLM structure-specification test on the AUTHENTIC per-SNP 118K panel:
#   Q+K       : fixed = intercept + 5 PCs         (Chen Fig-4C structure)
#   Family+K  : fixed = family factor (5 levels)  (the OLS/JLM structure)
# both with the SAME VanRaden K (random) and the joint 2-df additive+dominance F,
# R-EMMAX (method-matched). Question: is the tb1/domestication-peak suppression driven
# by K (the ancestry confound -> stays under both) or by the 5 PCs (-> tb1 returns under
# Family)? Reports candidate-gene peaks + lambda_GC for both; writes the Family+K scan.
#
# Output: data/teonam/stam_gwas_mlm_family_118k.csv, results/.../stam_gwas_mlm_family_manhattan_118k.png
# Run: Rscript scripts/teonam_mlm_family_118k.R
# =============================================================================
suppressMessages({
  library(data.table)
  library(readxl)
})
ROOT <- "/Users/fvrodriguez/repos/zealhmm"
setwd(ROOT)
source(file.path(ROOT, "scripts/logging.R"))
source(file.path(ROOT, "scripts/emmax_qk.R"))
NPC <- 5L

g <- readRDS("data/teonam/teonam_gwas118k_dosage_polar.rds")
dos <- g$dos # AUTHENTIC per-SNP, markers x lines, 0/1/2, ~2.7% NA
mc <- fread("data/teonam/markers_v5_gwas118k.tsv")
chr_by <- setNames(mc$chr_v5, mc$marker)
pos_by <- setNames(mc$pos_v5, mc$marker)
mk <- intersect(rownames(dos), mc$marker)
G <- dos[mk, , drop = FALSE]
CHR <- as.integer(chr_by[mk])
BP <- as.integer(pos_by[mk])
o <- order(CHR, BP)
G <- G[o, , drop = FALSE]
CHR <- CHR[o]
BP <- BP[o]

ph <- as.data.frame(read_excel("data/teonam/9250682/TeoNAM_1257RILs_22traits_phenotype_data.xlsx"))
names(ph)[1] <- "line"
TRAIT <- toupper(Sys.getenv("TRAIT", "STAM"))
TTAG <- tolower(TRAIT) # phenotype col; STAM default, e.g. DTA
if (!TRAIT %in% names(ph)) stop("TRAIT '", TRAIT, "' is not a phenotype column")
stam <- suppressWarnings(as.numeric(ph[[TRAIT]]))
names(stam) <- ph$line
lines <- intersect(colnames(G), names(stam)[is.finite(stam)])
G <- G[, lines, drop = FALSE]
storage.mode(G) <- "double"
if (anyNA(G)) {
  rmn <- rowMeans(G, na.rm = TRUE)
  rmn[!is.finite(rmn)] <- 0
  na <- which(is.na(G), arr.ind = TRUE)
  G[na] <- rmn[na[, 1]]
}
fam <- factor(substr(lines, 1, 5))
log_info("authentic panel: %d markers x %d lines, %d families", nrow(G), ncol(G), nlevels(fam))

vanraden <- function(Ms) {
  pr <- colMeans(Ms) / 2
  Z <- sweep(Ms, 2, 2 * pr, "-")
  tcrossprod(Z) / (2 * sum(pr * (1 - pr)))
}
lambda_gc <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  round(median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1), 3)
}
M <- t(G)
y <- stam[lines]
n <- length(y)
vpc <- (colMeans(M^2) - colMeans(M)^2) > 1e-12
K <- vanraden(M[, vpc, drop = FALSE]) # SAME K for both models
eig <- eigen(K, symmetric = TRUE)
U <- eig$vectors
xi <- pmax(eig$values, 1e-8)
yt <- as.numeric(crossprod(U, y))

build_null <- function(X) {
  p <- ncol(X)
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
  list(lines = lines, U = U, ws = ws, Xs = Xs, XtXinv = XtXinv, ry = ry, n = n, p = p, delta = delta)
}
X_q <- cbind(1, prcomp(M[, vpc, drop = FALSE], center = TRUE, scale. = FALSE, rank. = NPC)$x[, seq_len(NPC)])
X_fam <- model.matrix(~fam)
scan_q <- emmax_qk_scan(G, build_null(X_q), CHR, BP)
scan_fam <- emmax_qk_scan(G, build_null(X_fam), CHR, BP)[order(CHR, BP)]

# candidate-gene peaks for both
ov <- fread(sprintf("results/sim/teonam/%s_candidate_overlap.csv", TTAG))
peak <- function(s, ch, st) {
  w <- s[CHR == ch & abs(BP - st) <= 5e5 & is.finite(P) & P > 0]
  if (!nrow(w)) NA_real_ else round(max(-log10(w$P)), 2)
}
cmp <- ov[, .(
  gene = symbol, chr,
  `Q+K` = mapply(function(c, s) peak(scan_q, c, s), chr, start),
  `Family+K` = mapply(function(c, s) peak(scan_fam, c, s), chr, start)
)]
setorder(cmp, -`Q+K`)
log_info("lambda_GC:  Q+K = %.3f   Family+K = %.3f", lambda_gc(scan_q$P), lambda_gc(scan_fam$P))
log_info("candidate peak -log10P (+/-500 kb):")
print(cmp)

fwrite(scan_fam, sprintf("data/teonam/%s_gwas_mlm_family_118k.csv", TTAG))
log_info("wrote data/teonam/%s_gwas_mlm_family_118k.csv (max -log10P = %.2f)", TTAG, max(-log10(scan_fam[is.finite(P) & P > 0, P])))
