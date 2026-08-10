#!/usr/bin/env Rscript
# =============================================================================
# ZEAL flowering photoperiod scan — marker x macroenvironment (G x ME), R/qtl
# Approach B (interactive covariate).
# -----------------------------------------------------------------------------
# ME (macroenvironment) = day length: long-day temperate (Clayton BLUE, from the
# zeal_dta_cross.rds phenotypes) vs short-day subtropical (CUCBA BLUE). Same NILs
# in both MEs -> a doubled cross with an environment covariate:
#   out.full = scanone(BLUE, addcovar=[taxon,env], intcovar=env)   # marker x env
#   out.add  = scanone(BLUE, addcovar=[taxon,env])
#   gxe LOD  = out.full - out.add                                   # interaction
# Threshold = genome-wide 1000-perm, permuted at the LINE level (breaks
# genotype<->phenotype while keeping the within-line env pairing) to respect the
# two-records-per-line non-independence of the doubled cross.
# Faithful to scripts/zeal_rqtl_scan.R: bcsft, calc.genoprob step=1
# error.prob=0.01 haldane, method=hk, taxon dummy covariate.
#
# CAVEATS (provisional): short-day side = CUCBA only (La Soledad flowering is a
# derived copy, excluded); one unreplicated environment; CUCBA<->pedigree join via
# the Genealogia crosswalk (pending the definitive three-way sheet).
#
# In:  results/sim/zeal/rqtl/zeal_dta_cross.rds              (Clayton cross: genos + DTA/DTS phenos)
#      results/sim/zeal/rqtl/zeal_cucba_flowering_blues.csv  (from scripts/zeal_cucba_pheno.R)
#      data/zeal/samplesheet_3way.csv                        (taxon covariate)
# Out: results/sim/zeal/rqtl/mexg_<trait>_scan.csv, mexg_<trait>_result.rds
# Env: TRAIT (DTA|DTS|ASI; default DTA), NPERM (default 1000), MEXG_OUT
# Run: TRAIT=DTA NPERM=1000 Rscript scripts/zeal_mexg_scan.R
#      TRAIT=DTS NPERM=1000 Rscript scripts/zeal_mexg_scan.R
# =============================================================================
suppressPackageStartupMessages({
  library(qtl)
  library(data.table)
  library(here)
})
set.seed(1234567890)
RQ <- here("results/sim/zeal/rqtl")
OUT <- Sys.getenv("MEXG_OUT", RQ)
NPERM <- as.integer(Sys.getenv("NPERM", "1000"))
TRAIT <- Sys.getenv("TRAIT", "DTA")

# ---- Clayton cross (light rds: genotypes + ALL Clayton trait phenos) --------
cr <- readRDS(file.path(RQ, "zeal_dta_cross.rds"))
ids <- as.character(cr$pheno$id)
dta_cly <- as.numeric(cr$pheno[[TRAIT]]) # long-day (Clayton) BLUE

# ---- CUCBA short-day BLUE, joined by line_id (= cross id) -------------------
cu <- fread(file.path(RQ, "zeal_cucba_flowering_blues.csv"))
cu <- cu[!is.na(line_id) & line_id != "" & is.finite(get(TRAIT))]
dta_cu <- setNames(cu[[TRAIT]], cu$line_id)[ids] # short-day (CUCBA) BLUE

# ---- taxon covariate (same source as the additive pipeline) ----------------
ss <- fread(here("data/zeal/samplesheet_3way.csv"))[, .(pedigree, taxon)]
tax <- ss[match(ids, pedigree), taxon]

keep <- which(is.finite(dta_cly) & is.finite(dta_cu) & !is.na(tax) & tax != "")
cat(sprintf("ME x G lines (both MEs + taxon): %d\n", length(keep)))
crs <- cr[, keep]
n <- nind(crs)
idk <- as.character(crs$pheno$id)
dcly <- dta_cly[keep]
dcu <- dta_cu[keep]
txk <- tax[keep]

# ---- doubled cross: 2 records/line (env 0=CLY long, 1=CUCBA short) ---------
# double the genotype DATA, then compute genoprob ONCE on the doubled cross.
cr2 <- crs
for (c in seq_along(cr2$geno)) {
  cr2$geno[[c]]$data <- crs$geno[[c]]$data[c(1:n, 1:n), , drop = FALSE]
  cr2$geno[[c]]$prob <- NULL
}
cr2$pheno <- data.frame(
  id = rep(idk, 2), BLUE = c(dcly, dcu),
  env = c(rep(0L, n), rep(1L, n)), taxon = rep(txk, 2),
  stringsAsFactors = FALSE
)
cr2 <- calc.genoprob(cr2, step = 1, error.prob = 0.01, map.function = "haldane")

envv <- cr2$pheno$env
scan1 <- function(cx, addc, intc = NULL) {
  scanone(cx, pheno.col = "BLUE", addcovar = addc, intcovar = intc, method = "hk")
}

addcov <- cbind(model.matrix(~ factor(cr2$pheno$taxon))[, -1, drop = FALSE], env = envv)
out.full <- scan1(cr2, addcov, intc = envv)
out.add <- scan1(cr2, addcov)
gxe <- out.add
gxe$lod <- out.full$lod - out.add$lod # marker x ME interaction

fwrite(
  data.table(
    marker = rownames(gxe), chr = gxe$chr, pos = gxe$pos,
    lod_full = out.full$lod, lod_add = out.add$lod, lod_gxe = gxe$lod
  ),
  file.path(OUT, sprintf("mexg_%s_scan.csv", tolower(TRAIT)))
)

gtab <- as.data.table(gxe, keep.rownames = "marker")
gtab[, is_real := grepl("^S[0-9]+_", marker)]
gtab[, bp := ifelse(is_real, as.integer(sub("^S[0-9]+_", "", marker)), NA_integer_)]
cat(sprintf("\nmax interaction LOD = %.2f\n", max(gxe$lod, na.rm = TRUE)))
top <- gtab[is_real == TRUE][, .SD[which.max(lod)], by = chr][order(-lod)]
top[, Mb := round(bp / 1e6, 2)]
cat("top interaction peaks (per chr, max at a real marker):\n")
print(head(top[, .(chr, marker, cM = round(pos, 2), Mb, lod = round(lod, 2))], 10))

# ---- LINE-LEVEL permutation for the interaction threshold ------------------
cat(sprintf("\npermuting interaction (line-level, NPERM=%d) ...\n", NPERM))
permmax <- numeric(NPERM)
for (b in seq_len(NPERM)) {
  p <- sample(n)
  cr2$pheno$BLUE <- c(dcly[p], dcu[p])
  addcov_b <- cbind(model.matrix(~ factor(c(txk[p], txk[p])))[, -1, drop = FALSE], env = envv)
  f <- scan1(cr2, addcov_b, intc = envv)$lod
  a <- scan1(cr2, addcov_b)$lod
  permmax[b] <- max(f - a, na.rm = TRUE)
}
cr2$pheno$BLUE <- c(dcly, dcu)
thr <- quantile(permmax, c(0.90, 0.95), na.rm = TRUE)
cat(sprintf("interaction LOD threshold: 10%%=%.2f  5%%=%.2f  (NPERM=%d)\n", thr[1], thr[2], NPERM))
saveRDS(
  list(gxe = gxe, full = out.full, add = out.add, permmax = permmax, thr = thr, n = n),
  file.path(OUT, sprintf("mexg_%s_result.rds", tolower(TRAIT)))
)
sig <- top[lod >= thr[2]]
cat("\nchromosomes with interaction LOD >= 5% threshold:\n")
print(sig[, .(chr, Mb, lod = round(lod, 2))])
cat("DONE.\n")
