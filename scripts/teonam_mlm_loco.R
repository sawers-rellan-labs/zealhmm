#!/usr/bin/env Rscript
# LOCO-MLM test: is the MLM's loss of tb1/tsh4/zmm16 real deflation or proximal
# contamination (kinship absorbing each QTL's own-chromosome signal)?
# Approach: EMMAX in R (VanRaden/Centered-IBS kinship + REML for delta + P3D GLS),
#   (1) FULL-K  -> must reproduce the TASSEL MLM (lambda~1.36, tb1 dead)  [validation]
#   (2) LOCO    -> kinship built from all markers EXCEPT the tested marker's chr.
# If peaks RETURN under LOCO, the MLM loss was proximal contamination (over-correction).
# Uses the interpolated complete matrix data/teonam/tassel/geno_gwas_interp.rds.
# Run: Rscript agent/teonam_mlm_loco.R
suppressMessages({
  library(data.table)
  library(readxl)
})
setwd("/Users/fvrodriguez/repos/zealhmm")

d <- readRDS("data/teonam/tassel/geno_gwas_interp.rds") # G (markers x lines), markers(marker,chr,cm)
G <- d$G
mk <- as.data.table(d$markers)
mc <- fread("data/teonam/map_v5_coe2008.tsv")
pos_by <- setNames(mc$pos_v5, mc$marker)
mk[, bp := as.integer(pos_by[marker])]
M <- t(G) # lines x markers, 0/1/2, complete

ph <- as.data.frame(read_excel("data/teonam/9250682/TeoNAM_1257RILs_22traits_phenotype_data.xlsx"))
names(ph)[1] <- "line"
stam <- suppressWarnings(as.numeric(ph$STAM))
names(stam) <- ph$line
keep <- intersect(rownames(M), names(stam)[is.finite(stam)])
M <- M[keep, ]
y <- stam[keep]
n <- length(y)
Mv <- M[, (colMeans(M^2) - colMeans(M)^2) > 0]
Q <- prcomp(Mv, center = TRUE, scale. = FALSE)$x[, 1:5]
X <- cbind(1, Q)
p <- ncol(X)
cat(sprintf("lines=%d markers=%d PCs=5\n", n, ncol(M)))

vanraden <- function(Ms) {
  pr <- colMeans(Ms) / 2
  Z <- sweep(Ms, 2, 2 * pr, "-")
  tcrossprod(Z) / (2 * sum(pr * (1 - pr)))
}
lam <- function(pv) {
  pv <- pv[is.finite(pv) & pv > 0 & pv <= 1]
  round(median(qchisq(pv, 1, lower.tail = FALSE)) / qchisq(.5, 1), 3)
}

# EMMAX: null delta via REML on (y,X,K); then vectorized P3D GLS on markers Gc (lines x mk)
emmax <- function(y, X, K, Gc) {
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
    rss <- sum((ys - Xs %*% b)^2)
    0.5 * ((n - p) * log(rss / (n - p)) + sum(log(xi + exp(ld))) + as.numeric(determinant(XtX, log = TRUE)$modulus))
  }
  delta <- exp(optimize(nll, c(-10, 10))$minimum)
  ws <- sqrt(1 / (xi + delta))
  ys <- yt * ws
  Xs <- Xt * ws
  XtXinv <- solve(crossprod(Xs))
  ry <- as.numeric(ys - Xs %*% (XtXinv %*% crossprod(Xs, ys))) # y residualized on X (whitened)
  ronx <- function(Z) Z - Xs %*% (XtXinv %*% crossprod(Xs, Z))
  A <- ronx(crossprod(U, Gc) * ws)
  D <- ronx(crossprod(U, (Gc == 1) * 1) * ws) # additive + dominance(het)
  aa <- colSums(A^2)
  dd <- colSums(D^2)
  ad <- colSums(A * D)
  ay <- colSums(A * ry)
  dy <- colSums(D * ry)
  rss0 <- sum(ry^2)
  det <- aa * dd - ad^2
  expl <- ifelse(dd > 1e-9 & det > 1e-9, (dd * ay^2 - 2 * ad * ay * dy + aa * dy^2) / det, ay^2 / aa) # TASSEL default: 2-df add+dom
  k2 <- ifelse(dd > 1e-9 & det > 1e-9, 2, 1)
  Fst <- (expl / k2) / ((rss0 - expl) / (n - p - k2))
  pf(Fst, k2, n - p - k2, lower.tail = FALSE) # joint F (matches TASSEL 'p')
}

# (1) FULL-K validation
Kfull <- vanraden(M)
p_full <- emmax(y, X, Kfull, M)
# (2) LOCO
p_loco <- rep(NA_real_, ncol(M))
names(p_loco) <- colnames(M)
for (c in sort(unique(mk$chr))) {
  Kc <- vanraden(M[, mk$chr != c])
  ci <- which(mk$chr == c)
  p_loco[ci] <- emmax(y, X, Kc, M[, ci, drop = FALSE])
  cat(sprintf("  LOCO chr%d done (%d markers)\n", c, length(ci)))
}

res <- data.table(marker = colnames(M), chr = mk$chr, bp = mk$bp, p_full = p_full, p_loco = p_loco)
fwrite(res, "results/sim/teonam/mlm_loco_scan.csv")
cat(sprintf("\nlambda_GC:  EMMAX full-K = %.3f (cf. TASSEL MLM 1.36)   LOCO = %.3f\n", lam(p_full), lam(p_loco)))

# candidate-gene peaks: OLS, TASSEL-MLM, EMMAX-full, LOCO
ov <- fread("results/sim/teonam/stam_candidate_overlap.csv")
ols <- fread("data/teonam/stam_gwas_scan_interpolated.csv")[is.finite(P) & P > 0]
tmlm <- fread("data/teonam/tassel/mlm_interp2.txt")[Marker != "None", .(chr = as.integer(Chr), bp = as.integer(Pos), p = as.numeric(p))][is.finite(p) & p > 0]
pk <- function(dt, pc, ch, ps) {
  w <- dt[get(names(dt)[grep("chr|CHR", names(dt))[1]]) == ch]
  NULL
}
peak <- function(chr, bp, p, ch, ps) {
  w <- which(chr == ch & abs(bp - ps) <= 5e5 & is.finite(p) & p > 0)
  if (length(w)) round(max(-log10(p[w])), 2) else NA
}
tab <- ov[, .(
  gene = symbol, chr,
  OLS = mapply(function(c, s) peak(ols$CHR, ols$BP, ols$P, c, s), chr, start),
  MLM = mapply(function(c, s) peak(tmlm$chr, tmlm$bp, tmlm$p, c, s), chr, start),
  EMMAXfull = mapply(function(c, s) peak(res$chr, res$bp, res$p_full, c, s), chr, start),
  LOCO = mapply(function(c, s) peak(res$chr, res$bp, res$p_loco, c, s), chr, start)
)]
cat("\ncandidate-gene peak -log10P:\n")
print(tab[order(-OLS)])
