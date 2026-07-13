#!/usr/bin/env Rscript
# =============================================================================
# GWAS genome-wide FWER permutation threshold for the ZEAL SNP50K scans.
# Threshold = the alpha-quantile of the genome-wide max(-log10 P) under permutation
# (alpha = 0.05 and 0.10 stored). Sourced by zeal_ols_taxon.R / zeal_mlm_taxon.R;
# each writes its own (trait, model, geno) row into data/zeal/gwas_perm_thresholds.csv
# (the same table the JLM uses).
#
# OLS  (y ~ taxon + marker, 1-df F): residualise y and every marker against the taxon
#   fixed part, then permute the residualised phenotype and recompute all-marker F.
# MLM  (EMMAX Q/Family + K, 2-df add+dom F): the null already gives the whitened,
#   residualised phenotype `ry` and the machinery to residualise markers. The whitened
#   markers A (additive) and D (dominance) are permutation-invariant, so we permute `ry`
#   and recompute A'ry / D'ry per marker (rss0 = sum(ry^2) is permutation-invariant).
#   This is the standard EMMAX rotation permutation — the VC/eigendecomposition is
#   estimated ONCE (in build_null), not per permutation.
#
# Seeded (SEED, default 1) so the stored threshold is reproducible.
# =============================================================================
suppressMessages(library(data.table))

# --- max(-log10 P) null for the OLS (taxon + marker) 1-df scan ---------------
# G = markers x lines (mean-imputed), fam = factor over lines, y = phenotype over lines.
perm_max_ols <- function(G, fam, y, nperm = 1000L, seed = 1L) {
  ok <- is.finite(y)
  G <- G[, ok, drop = FALSE]
  y <- y[ok]
  fam <- droplevels(fam[ok])
  X <- if (nlevels(fam) > 1) model.matrix(~fam) else matrix(1, length(y), 1)
  n <- length(y)
  p <- ncol(X)
  XtXinv <- solve(crossprod(X))
  resid <- function(z) z - X %*% (XtXinv %*% crossprod(X, z)) # project taxon out
  yr <- as.numeric(resid(y))
  Gr <- t(resid(t(G))) # residualise each marker (row) against taxon: markers x lines
  ggr <- rowSums(Gr^2)
  rss0 <- sum(yr^2)
  df2 <- n - p - 1
  set.seed(seed)
  YRP <- vapply(seq_len(nperm), function(k) yr[sample.int(n)], numeric(n)) # n x nperm
  GY <- Gr %*% YRP # markers x nperm
  expl <- (GY^2) / ggr # 1-df explained SS
  Fst <- (expl / 1) / ((rss0 - expl) / df2)
  Pv <- pf(Fst, 1, df2, lower.tail = FALSE)
  apply(-log10(Pv), 2, max, na.rm = TRUE)
}

# --- max(-log10 P) null for the EMMAX MLM (add+dom 2-df) scan -----------------
# G = markers x lines (0/1/2), null = build_null(...) from zeal_mlm_taxon.R, CHR aligned
# to rows of G. Mirrors emmax_qk_scan() with `ry` permuted.
perm_max_mlm <- function(G, null, CHR, nperm = 1000L, seed = 1L) {
  U <- null$U
  ws <- null$ws
  Xs <- null$Xs
  XtXinv <- null$XtXinv
  ry <- null$ry
  n <- null$n
  p <- null$p
  M <- t(G[, null$lines, drop = FALSE])
  storage.mode(M) <- "double"
  rss0 <- sum(ry^2)
  ron <- function(Z) Z - Xs %*% (XtXinv %*% crossprod(Xs, Z))
  set.seed(seed)
  RYP <- vapply(seq_len(nperm), function(k) ry[sample.int(n)], numeric(n)) # n x nperm
  maxlp <- rep(0, nperm)
  for (ch in unique(CHR)) {
    cols <- which(CHR == ch)
    Mc <- M[, cols, drop = FALSE]
    A <- ron(crossprod(U, Mc) * ws) # n x m_chr, additive whitened+residualised
    D <- ron(crossprod(U, (Mc == 1) * 1) * ws) # dominance (het indicator)
    aa <- colSums(A^2)
    dd <- colSums(D^2)
    ad <- colSums(A * D)
    det <- aa * dd - ad^2
    use2 <- dd > 1e-9 & det > 1e-9 # 2-df where dominance is estimable, else 1-df additive
    AY <- crossprod(A, RYP) # m_chr x nperm
    DY <- crossprod(D, RYP)
    expl <- (dd * AY^2 - 2 * ad * AY * DY + aa * DY^2) / det # 2-df (recycles the per-marker vectors down columns)
    if (any(!use2)) expl[!use2, ] <- (AY[!use2, , drop = FALSE]^2) / aa[!use2]
    k2 <- ifelse(use2, 2L, 1L)
    df2 <- n - p - k2
    Fst <- (expl / k2) / ((rss0 - expl) / df2)
    Pv <- pf(Fst, k2, df2, lower.tail = FALSE) # k2/df2 recycle down columns
    cmax <- apply(-log10(Pv), 2, max, na.rm = TRUE)
    maxlp <- pmax(maxlp, cmax, na.rm = TRUE)
  }
  maxlp
}

# --- store one (trait, model, geno) pair of rows (alpha 0.05 + 0.10) ---------
upsert_gwas_threshold <- function(csv, trait, model, geno, famcol, maxnull, nperm, seed = 1L) {
  q <- quantile(maxnull[is.finite(maxnull)], c(0.95, 0.90))
  new <- data.table(
    trait = trait, model = model, geno = geno, famcol = famcol,
    alpha = c(0.05, 0.10), thr_neglog10p = as.numeric(q),
    enter_p = NA_real_, nperm = nperm, note = sprintf("seed=%d", seed)
  )
  cur <- if (file.exists(csv)) fread(csv) else data.table()
  if (nrow(cur)) {
    tr <- trait
    mdl <- model
    gn <- geno
    # locals renamed off the column names so the i-filter resolves them from scope
    cur <- cur[!(trait == tr & model == mdl & geno == gn)]
  }
  fwrite(rbind(cur, new, fill = TRUE), csv)
  q
}
