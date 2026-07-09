#!/usr/bin/env Rscript
# =============================================================================
# ZEAL/BZea Phase 4 — MLM (Taxon + K) GWAS on the RTIGER ancestry mosaic.
# The ZEAL analog of teonam_mlm_family_118k.R: fixed = the taxon factor (Zx/Zv/Zd/Zl/Zh,
# the 5-family analog), random = VanRaden K, joint 2-df additive+dominance F, R-EMMAX
# (emmax_qk.R). Genotype = the RTIGER mosaic state (0/1/2), since the raw per-SNP dosage
# is 70% missing at 0.4x. 102 invariant-REF markers dropped (zero variance). No thinning.
#
# TRAIT via env (default DTA): phenotype = data/zeal/pheno_<trait>_blue.csv (<TRAIT>_mean BLUE).
# Output: data/zeal/<trait>_gwas_mlm_taxon_snp50k.csv (SNP, CHR, BP, P), data/zeal/zeal_K_vanraden.rds
# =============================================================================
suppressMessages({
  library(here)
  library(data.table)
})
source(here("scripts/logging.R"))
source(here("scripts/emmax_qk.R"))

TRAIT <- toupper(Sys.getenv("TRAIT", "DTA"))
TTAG <- tolower(TRAIT)
GENO <- Sys.getenv("GENO", "mosaic") # mosaic = RTIGER ancestry (Panel C) | persnp = per-SNP dosage (Panel B)

# --- genotype (0/1/2), drop invariant markers ---------------------------------
inv <- tryCatch(readRDS(here("data/zeal/snp50k_invariant_markers.rds")), error = function(e) character(0))
if (GENO %in% c("mosaic", "rtiger", "nnil", "lbimpute", "binhmm")) {
  caller_file <- if (GENO == "mosaic") "rtiger" else GENO # "mosaic" == the RTIGER reference mosaic
  M0 <- readRDS(here(sprintf("data/zeal/zeal_%s_mosaic.rds", caller_file)))
  state <- M0$state
  mk <- M0$markers
} else { # per-SNP authentic: teosinte dosage 2*alt/cov rounded to 0/1/2 (70% missing -> mean-imputed below)
  D <- readRDS(here("data/zeal/zeal_snp50k_dosage.rds"))
  ss0 <- fread(here("data/zeal/samplesheet_3way.csv"))[gwas_nil == TRUE & !is.na(skim_id)]
  sk <- intersect(ss0$skim_id, colnames(D$n_alt))
  dose <- 2 * D$n_alt[, sk] / pmax(D$cov[, sk], 1)
  dose[D$cov[, sk] == 0] <- NA_real_
  state <- round(dose)
  colnames(state) <- ss0$pedigree[match(sk, ss0$skim_id)]
  mk <- copy(D$markers)[, .(marker, chr = as.integer(chr), pos = as.integer(pos))]
}
keep <- !(mk$marker %in% inv)
state <- state[keep, , drop = FALSE]
mk <- mk[keep]
if (any(duplicated(colnames(state)))) state <- state[, !duplicated(colnames(state)), drop = FALSE]

# --- phenotype BLUE + taxon (family) ------------------------------------------
ph <- fread(here(sprintf("data/zeal/pheno_%s_blue.csv", TTAG)))
mcol <- paste0(TRAIT, "_mean")
y_all <- setNames(ph[[mcol]], ph$Genotype)
ss <- fread(here("data/zeal/samplesheet_3way.csv"))[gwas_nil == TRUE]
FAMCOL <- Sys.getenv("FAMILY_COL", "taxon") # taxon (5) or donor_accession (82) fixed factor
taxon_by <- setNames(ss[[FAMCOL]], ss$pedigree)

lines <- intersect(colnames(state), names(y_all)[is.finite(y_all)])
lines <- lines[!is.na(taxon_by[lines])]
G <- state[, lines, drop = FALSE]
storage.mode(G) <- "double"
CHR <- mk$chr
BP <- mk$pos
o <- order(CHR, BP)
G <- G[o, , drop = FALSE]
CHR <- CHR[o]
BP <- BP[o]
if (anyNA(G)) {
  rmn <- rowMeans(G, na.rm = TRUE)
  rmn[!is.finite(rmn)] <- 0
  na <- which(is.na(G), arr.ind = TRUE)
  G[na] <- rmn[na[, 1]]
}
fam <- factor(taxon_by[lines])
y <- y_all[lines]
log_info(
  "MLM panel: %d markers x %d lines | %d taxa (%s)", nrow(G), ncol(G), nlevels(fam),
  paste(levels(fam), table(fam), sep = ":", collapse = " ")
)

# --- VanRaden K + R-EMMAX null (ported from teonam_mlm_family_118k.R) ----------
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
n <- length(y)
vpc <- (colMeans(M^2) - colMeans(M)^2) > 1e-12
K <- vanraden(M[, vpc, drop = FALSE])
saveRDS(list(K = K, lines = lines), here(sprintf("data/zeal/zeal_K_vanraden_%s.rds", GENO)))
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
scan_fam <- emmax_qk_scan(G, build_null(model.matrix(~fam)), CHR, BP)[order(CHR, BP)]

# --- QC + write ---------------------------------------------------------------
log_info(
  "lambda_GC (Taxon+K) = %.3f | max -log10P = %.2f", lambda_gc(scan_fam$P),
  max(-log10(scan_fam[is.finite(P) & P > 0, P]))
)
gg <- fread(here("data/teonam/dta_candidate_genes.tsv")) # same maize v5 flowering candidates
gg[, chr := as.integer(chr)]
peak <- function(ch, st) {
  w <- scan_fam[CHR == ch & abs(BP - st) <= 5e5 & is.finite(P) & P > 0]
  if (!nrow(w)) NA_real_ else round(max(-log10(w$P)), 2)
}
cand <- gg[, .(gene = symbol, chr, `MLM peak -log10P` = mapply(peak, chr, start))][order(-`MLM peak -log10P`)]
print(cand)
fwrite(scan_fam, here(sprintf("data/zeal/%s_gwas_mlm_%s_%s_snp50k.csv", TTAG, GENO, FAMCOL)))
log_info("wrote data/zeal/%s_gwas_mlm_%s_%s_snp50k.csv", TTAG, GENO, FAMCOL)
