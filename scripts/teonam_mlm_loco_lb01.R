#!/usr/bin/env Rscript
# MLM (Q+K) and LOCO on the 0.1x LB-Impute genotypes — the low-coverage analogue
# of agent/teonam_mlm_loco.R (which runs on the interpolated-truth baseline).
# Purpose: test whether the chr3-pericentromere "ghost" (a structure/LD confound
# that full-K MLM removes at baseline) survives Q+K / resurfaces under LOCO once
# low-coverage LB-Impute block imputation inflates it. Same EMMAX as the baseline
# script (VanRaden K + 5 PCs + REML delta + P3D GLS, 2-df add+dom to match TASSEL).
#
# Input : results/sim/teonam/cache/geno_lb0.1.rds — the 51,004 x 1,237 union matrix
#         at coverage lambda=0.1, calibrated recombdist_star (drp=TRUE). Regenerate
#         with `Rscript scripts/teonam_lbimpute_sweep.R --generate --save-matrix`.
# Output: results/sim/teonam/mlm_loco_scan_lb01.csv  (marker, chr, bp, p_full, p_loco)
# Run   : Rscript agent/teonam_mlm_loco_lb01.R
suppressMessages({
  library(data.table)
  library(readxl)
})
setwd("/Users/fvrodriguez/repos/zealhmm")

d <- readRDS("results/sim/teonam/cache/geno_lb0.1.rds") # G (markers x lines), markers, chr, bp
G <- d$G
mk <- data.table(marker = d$markers, chr = d$chr, bp = d$bp)
M <- t(G) # lines x markers, 0/1/2

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
cat(sprintf("0.1x LB-Impute: lines=%d markers=%d PCs=5\n", n, ncol(M)))

vanraden <- function(Ms) {
  pr <- colMeans(Ms) / 2
  Z <- sweep(Ms, 2, 2 * pr, "-")
  tcrossprod(Z) / (2 * sum(pr * (1 - pr)))
}
lam <- function(pv) {
  pv <- pv[is.finite(pv) & pv > 0 & pv <= 1]
  round(median(qchisq(pv, 1, lower.tail = FALSE)) / qchisq(.5, 1), 3)
}
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
  ry <- as.numeric(ys - Xs %*% (XtXinv %*% crossprod(Xs, ys)))
  ronx <- function(Z) Z - Xs %*% (XtXinv %*% crossprod(Xs, Z))
  A <- ronx(crossprod(U, Gc) * ws)
  D <- ronx(crossprod(U, (Gc == 1) * 1) * ws)
  aa <- colSums(A^2)
  dd <- colSums(D^2)
  ad <- colSums(A * D)
  ay <- colSums(A * ry)
  dy <- colSums(D * ry)
  rss0 <- sum(ry^2)
  det <- aa * dd - ad^2
  expl <- ifelse(dd > 1e-9 & det > 1e-9, (dd * ay^2 - 2 * ad * ay * dy + aa * dy^2) / det, ay^2 / aa)
  k2 <- ifelse(dd > 1e-9 & det > 1e-9, 2, 1)
  Fst <- (expl / k2) / ((rss0 - expl) / (n - p - k2))
  pf(Fst, k2, n - p - k2, lower.tail = FALSE)
}

# (1) FULL-K
Kfull <- vanraden(M)
p_full <- emmax(y, X, Kfull, M)
# (2) LOCO — kinship from all markers except the tested marker's chromosome
p_loco <- rep(NA_real_, ncol(M))
names(p_loco) <- colnames(M)
for (c in sort(unique(mk$chr))) {
  Kc <- vanraden(M[, mk$chr != c])
  ci <- which(mk$chr == c)
  p_loco[ci] <- emmax(y, X, Kc, M[, ci, drop = FALSE])
  cat(sprintf("  LOCO chr%d done (%d markers)\n", c, length(ci)))
}

res <- data.table(marker = colnames(M), chr = mk$chr, bp = mk$bp, p_full = p_full, p_loco = p_loco)
fwrite(res, "results/sim/teonam/mlm_loco_scan_lb01.csv")
cat(sprintf("\nlambda_GC (0.1x):  full-K = %.3f   LOCO = %.3f\n", lam(p_full), lam(p_loco)))
