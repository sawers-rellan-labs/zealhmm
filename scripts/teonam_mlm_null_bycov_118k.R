#!/usr/bin/env Rscript
# Per-coverage FIXED MLM (Family + K) null models for the 118K sweeps -- the honest
# low-coverage homologue of the structure/kinship correction.
#
# Fixed structure = the 5-family FACTOR (the OLS/JLM structure), which is known a priori
# and COVERAGE-INDEPENDENT -- a cleaner homologue than 5 PCs, which degrade with depth.
# The RANDOM K must still be estimated from the downsampled reads: per-SNP GL DOSAGE
# (posterior expected 0/1/2 from the counts), computed BEFORE ancestry calling from the
# SAME reads the callers see (identical seeds); GL dosage (not hard calls) so K does not
# fall apart at 0.1x. K is caller-independent (depends only on the reads), so one null
# per coverage is shared by all ancestry sweeps. Only the TESTED genotypes differ.
#
# Output: data/teonam/mlm_null_118k_l<lambda>.rds  (one per coverage)
# Run: Rscript scripts/teonam_mlm_null_bycov_118k.R
suppressMessages({
  library(data.table)
  library(readxl)
})
ROOT <- "/Users/fvrodriguez/repos/zealhmm"
setwd(ROOT)
source(file.path(ROOT, "scripts/logging.R"))
source(file.path(ROOT, "R/simulate.R")) # .draw_counts()
t0 <- Sys.time()
suppressWarnings(try(if (mem.maxVSize() < 22000) mem.maxVSize(22000), silent = TRUE))

LAMBDAS <- c(0.1, 0.2, 0.5, 1, 5, 10, 20, Inf) # MUST match the sweeps (seed alignment)
READ_PARS <- list(pi_floor = 0, k_decay = 1, error = 0.01)
NPC <- 5L
FAMS <- c("TIL01", "TIL03", "TIL11", "TIL14", "TIL25")

mt_thin <- fread("data/teonam/markers_v5_gwas118k_cm_thin01.tsv") # 0.1 cM inference grid
g <- readRDS("data/teonam/teonam_gwas118k_dosage_polar.rds")
dos <- g$dos
ph <- as.data.frame(read_excel("data/teonam/9250682/TeoNAM_1257RILs_22traits_phenotype_data.xlsx"))
names(ph)[1] <- "line"
stam <- suppressWarnings(as.numeric(ph$STAM))
names(stam) <- ph$line

# per-family truth on the thin grid (+ per-marker teosinte AF for the GL prior)
load_family <- function(fam) {
  keys <- colnames(dos)[substr(colnames(dos), 1, 5) == fam]
  D <- dos[mt_thin$marker, keys, drop = FALSE]
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
  list(D = D, keys = keys, af = rowMeans(D) / 2)
}
fam_data <- lapply(FAMS, load_family)
names(fam_data) <- FAMS

# GL posterior expected dosage E[g | counts], HWE(af) prior; 0-depth -> prior mean 2af
gl_dosage <- function(n_alt, depth, af, err) {
  L0 <- dbinom(n_alt, depth, err)
  L1 <- dbinom(n_alt, depth, 0.5)
  L2 <- dbinom(n_alt, depth, 1 - err)
  p0 <- L0 * (1 - af)^2
  p1 <- L1 * 2 * af * (1 - af)
  p2 <- L2 * af^2
  (p1 + 2 * p2) / (p0 + p1 + p2)
}

vanraden <- function(Ms) {
  pr <- colMeans(Ms) / 2
  Z <- sweep(Ms, 2, 2 * pr, "-")
  tcrossprod(Z) / (2 * sum(pr * (1 - pr)))
}

for (li in seq_along(LAMBDAS)) {
  lambda <- LAMBDAS[li]
  # assemble the GL-dosage genotype matrix (thin markers x all lines) at this coverage
  Dose <- do.call(cbind, lapply(FAMS, function(fam) {
    fi <- match(fam, FAMS)
    fd <- fam_data[[fam]]
    M <- nrow(fd$D)
    N <- length(fd$keys)
    set.seed(1000L + 100L * fi + li) # SAME seed as the ancestry sweeps
    if (is.infinite(lambda)) {
      p_alt <- c(0, 0.5, 1)[as.vector(fd$D) + 1L]
      p_eff <- p_alt * (1 - READ_PARS$error) + (1 - p_alt) * READ_PARS$error
      n_alt <- as.integer(round(100L * p_eff))
      depth <- rep(100L, length(n_alt))
    } else {
      ac <- .draw_counts(as.vector(fd$D), lambda,
        pi_floor = READ_PARS$pi_floor, k_decay = READ_PARS$k_decay, error = READ_PARS$error
      )
      n_alt <- as.integer(ac$alt)
      depth <- as.integer(ac$ref + ac$alt)
    }
    dose <- gl_dosage(n_alt, depth, fd$af, READ_PARS$error) # af recycled per marker
    matrix(dose, M, N, dimnames = list(mt_thin$marker, fd$keys))
  }))

  M <- t(Dose) # lines x thin-markers
  keep <- intersect(rownames(M), names(stam)[is.finite(stam)])
  M <- M[keep, , drop = FALSE]
  y <- stam[keep]
  n <- length(y)
  vpc <- (colMeans(M^2) - colMeans(M)^2) > 1e-12
  X <- model.matrix(~ factor(substr(keep, 1, 5))) # fixed = family factor (coverage-independent), not 5 PCs
  p <- ncol(X)
  K <- vanraden(M)
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
  ys <- yt * ws
  Xs <- Xt * ws
  XtXinv <- solve(crossprod(Xs))
  ry <- as.numeric(ys - Xs %*% (XtXinv %*% crossprod(Xs, ys)))
  saveRDS(
    list(lines = keep, U = U, ws = ws, Xs = Xs, XtXinv = XtXinv, ry = ry, n = n, p = p, delta = delta, lambda = lambda),
    sprintf("data/teonam/mlm_null_118k_l%s.rds", lambda)
  )
  log_info("  lambda=%-4g null: n=%d, thin markers=%d, delta=%.3g", lambda, n, ncol(M), delta)
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  log_info(">>> %d/%d done | elapsed %.1f min | avg %.1f min | ETA ~%.1f min remaining", li, length(LAMBDAS), el, el / li, (el / li) * (length(LAMBDAS) - li))
}
log_info("wrote per-coverage MLM nulls: data/teonam/mlm_null_118k_l*.rds")
