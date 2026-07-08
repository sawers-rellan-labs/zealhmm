# R-EMMAX MLM (Q+K) per-marker scan against a PRECOMPUTED null. Chen 2019's Fig-4C
# model (Q = 5 PCs + K kinship, P3D, joint 2-df additive+dominance F; verified
# lambda_GC 1.36 vs TASSEL 1.363, agent/teonam_emmax_adddom.R).
#
# The null (Q, K, eigendecomposition, REML delta, whitening) is built ONCE PER COVERAGE
# from the GL-dosage genotypes of the downsampled reads -- BEFORE ancestry calling and
# from the same reads the ancestry caller sees -- so the structure/kinship correction
# is the honest low-coverage estimate, not the recovered ancestry blocks (which absorb
# the QTL). See scripts/teonam_mlm_null_bycov_118k.R. This function only runs the
# per-marker test on the coverage-recovered (ancestry) genotypes.
#
# NOTE on the "Q" / "q+K" name: emmax_qk_scan is null-agnostic -- it runs whatever
# fixed effect X the null carries. Chen's ORIGINAL Fig-4C null uses Q = 5 PCs (above),
# but the 118K COVERAGE-SWEEP null (teonam_mlm_null_bycov_118k.R) instead uses the
# 5-FAMILY FACTOR as X (PCs degrade with depth; the family structure is known a priori
# and coverage-independent). So wherever the 118K sweeps say "MLM (Q+K)", the "Q" IS
# the family factor -- it is exactly MLM (Family + K). The q+K name is kept for
# continuity across the codebase; read it as Family+K in the 118K coverage sweeps.
#
# emmax_qk_scan(G, null, CHR, BP): G = markers x lines recovered genotypes (0/1/2);
#   null = readRDS(mlm_null_118k_l<lambda>.rds); CHR/BP aligned to rows of G.
suppressWarnings(try(if (mem.maxVSize() < 22000) mem.maxVSize(22000), silent = TRUE))

emmax_qk_scan <- function(G, null, CHR, BP) {
  U <- null$U
  ws <- null$ws
  Xs <- null$Xs
  XtXinv <- null$XtXinv
  ry <- null$ry
  n <- null$n
  p <- null$p
  if (!all(null$lines %in% colnames(G))) { # MLM needs the full panel the null was built on (e.g. skipped in 1-family smoke)
    return(data.table::data.table(SNP = rownames(G), CHR = CHR, BP = BP, P = NA_real_))
  }
  M <- t(G[, null$lines, drop = FALSE]) # lines (null order) x markers
  storage.mode(M) <- "double"
  rss0 <- sum(ry^2)
  ron <- function(Z) Z - Xs %*% (XtXinv %*% crossprod(Xs, Z))
  # per-marker test is independent given the null; CHROMOSOME chunks keep the
  # n x m_chr A/D matrices ~150 MB instead of forming full n x 118K (>16 GB OOM).
  Pv <- rep(NA_real_, ncol(M))
  for (ch in unique(CHR)) {
    cols <- which(CHR == ch)
    Mc <- M[, cols, drop = FALSE]
    A <- ron(crossprod(U, Mc) * ws) # additive, residualized (n x m_chr)
    D <- ron(crossprod(U, (Mc == 1) * 1) * ws) # dominance = het indicator
    aa <- colSums(A^2)
    dd <- colSums(D^2)
    ad <- colSums(A * D)
    ay <- colSums(A * ry)
    dy <- colSums(D * ry)
    det <- aa * dd - ad^2
    expl <- ifelse(dd > 1e-9 & det > 1e-9, (dd * ay^2 - 2 * ad * ay * dy + aa * dy^2) / det, ay^2 / aa)
    k2 <- ifelse(dd > 1e-9 & det > 1e-9, 2, 1)
    Fst <- (expl / k2) / ((rss0 - expl) / (n - p - k2))
    Pv[cols] <- pf(Fst, k2, n - p - k2, lower.tail = FALSE)
  }
  data.table::data.table(SNP = rownames(G), CHR = CHR, BP = BP, P = Pv)
}
