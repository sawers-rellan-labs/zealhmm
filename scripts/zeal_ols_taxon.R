#!/usr/bin/env Rscript
# ZEAL Phase 4 — OLS (Taxon + marker) GWAS on the RTIGER mosaic (no K).
# Analog of teonam_stam_gwas118k.R: per-marker y ~ taxon + marker(additive), 1-df F.
# The OLS<->MLM contrast (this vs zeal_mlm_taxon.R, same taxon fixed part) isolates the
# K term = the ancestry/flowering confound. TRAIT via env (default DTA).
# Output: data/zeal/<trait>_gwas_ols_taxon_snp50k.csv (SNP, CHR, BP, P, n)
suppressMessages({
  library(here)
  library(data.table)
  library(parallel)
})
source(here("scripts/logging.R"))
source(here("scripts/zeal_gwas_perm.R"))
TRAIT <- toupper(Sys.getenv("TRAIT", "DTA"))
TTAG <- tolower(TRAIT)

M0 <- readRDS(here("data/zeal/zeal_rtiger_mosaic.rds"))
state <- M0$state
mk <- M0$markers
inv <- tryCatch(readRDS(here("data/zeal/snp50k_invariant_markers.rds")), error = function(e) character(0))
keep <- !(mk$marker %in% inv)
state <- state[keep, , drop = FALSE]
mk <- mk[keep]
if (any(duplicated(colnames(state)))) state <- state[, !duplicated(colnames(state)), drop = FALSE]

PHENO <- Sys.getenv("PHENO", "blue") # blue = SpATS BLUE (default) | direct = raw genotype mean (StPi)
ph <- fread(here(sprintf("data/zeal/pheno_%s_%s.csv", TTAG, PHENO)))
# phenotype columns keep the trait's native case (e.g. StPi_mean); TRAIT is upper-cased -> match case-insensitively
mcol <- names(ph)[tolower(names(ph)) == tolower(paste0(TRAIT, "_mean"))][1]
stopifnot(!is.na(mcol))
y_all <- setNames(ph[[mcol]], ph$Genotype)
ss <- fread(here("data/zeal/samplesheet_3way.csv"))[gwas_nil == TRUE]
FAMCOL <- Sys.getenv("FAMILY_COL", "taxon")
taxon_by <- setNames(ss[[FAMCOL]], ss$pedigree)
lines <- intersect(colnames(state), names(y_all)[is.finite(y_all)])
lines <- lines[!is.na(taxon_by[lines])]
G <- state[, lines, drop = FALSE]
storage.mode(G) <- "double"
CHR <- mk$chr
BP <- mk$pos
o <- order(CHR, BP)
G <- G[o, ]
CHR <- CHR[o]
BP <- BP[o]
y <- y_all[lines]
fam <- factor(taxon_by[lines])
log_info("OLS panel: %d markers x %d lines, %d taxa", nrow(G), ncol(G), nlevels(fam))

scan1 <- function(i) {
  g <- G[i, ]
  ok <- is.finite(g) & is.finite(y)
  n <- sum(ok)
  if (n < 20 || sd(g[ok]) == 0) {
    return(c(NA_real_, n))
  }
  yy <- y[ok]
  ff <- droplevels(fam[ok])
  Xr <- if (nlevels(ff) > 1) model.matrix(~ff) else matrix(1, n, 1)
  RSS0 <- sum(lm.fit(Xr, yy)$residuals^2)
  fit <- lm.fit(cbind(Xr, g[ok]), yy)
  RSS1 <- sum(fit$residuals^2)
  df2 <- n - fit$rank
  if (df2 <= 0 || RSS1 <= 0) {
    return(c(NA_real_, n))
  }
  c(pf(((RSS0 - RSS1) / 1) / (RSS1 / df2), 1, df2, lower.tail = FALSE), n)
}
out <- mclapply(seq_len(nrow(G)), scan1, mc.cores = max(1L, detectCores() - 2L))
P <- vapply(out, `[`, numeric(1), 1)
N <- vapply(out, `[`, numeric(1), 2)
scan <- data.table(SNP = rownames(G), CHR = CHR, BP = BP, P = P, n = as.integer(N))[order(CHR, BP)]
lambda_gc <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  round(median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1), 3)
}
fwrite(scan, here(sprintf("data/zeal/%s_gwas_ols_%s_snp50k.csv", TTAG, FAMCOL)))
log_info(
  "OLS lambda_GC = %.3f | max -log10P = %.2f | wrote %s_gwas_ols_%s_snp50k.csv",
  lambda_gc(scan$P), max(-log10(scan[is.finite(P) & P > 0, P])), TTAG, FAMCOL
)
cand_file <- here(sprintf("data/teonam/%s_candidate_genes.tsv", TTAG))
if (!file.exists(cand_file)) cand_file <- here("data/teonam/dta_candidate_genes.tsv")
gg <- fread(cand_file)
gg[, chr := as.integer(chr)]
peak <- function(ch, st) {
  w <- scan[CHR == ch & abs(BP - st) <= 5e5 & is.finite(P) & P > 0]
  if (!nrow(w)) NA_real_ else round(max(-log10(w$P)), 2)
}
print(gg[, .(gene = symbol, chr, `OLS peak` = mapply(peak, chr, start))][order(-`OLS peak`)])

# --- genome-wide FWER permutation threshold (Phase 1) -------------------------
NPERM <- as.integer(Sys.getenv("NPERM", "1000"))
if (NPERM > 0) {
  Gimp <- G
  if (anyNA(Gimp)) {
    rmn <- rowMeans(Gimp, na.rm = TRUE)
    rmn[!is.finite(rmn)] <- 0
    na <- which(is.na(Gimp), arr.ind = TRUE)
    Gimp[na] <- rmn[na[, 1]]
  }
  t0 <- Sys.time()
  mn <- perm_max_ols(Gimp, fam, y, nperm = NPERM, seed = 1L)
  thr <- upsert_gwas_threshold(
    here("data/zeal/gwas_perm_thresholds.csv"), TRAIT, "ols", "rtiger_mosaic", FAMCOL, mn, NPERM
  )
  log_info(
    "OLS FWER perm threshold (%d perm, %.1fs): 5%%=%.2f 10%%=%.2f",
    NPERM, as.numeric(Sys.time() - t0, units = "secs"), thr[1], thr[2]
  )
}
