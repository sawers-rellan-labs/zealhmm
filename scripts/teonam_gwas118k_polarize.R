#!/usr/bin/env Rscript
# Polarize the 118K GWAS dosage to W22<->teosinte ancestry: 0 = W22 hom,
# 1 = het, 2 = teosinte hom -- required to use the 118K panel as SIMULATION TRUTH
# in the coverage sweep (read simulation is polarization-sensitive; the STAM GWAS
# F-test was not, so teonam_gwas118k_dosage.rds left polarization arbitrary).
#
# Principle: TeoNAM is a BC1S4 backcross to the common W22 recurrent parent, so
# teosinte contributes ~25% genome-wide and the TEOSINTE allele is the MINOR
# allele at essentially every marker. We recode each marker so the minor allele
# counts as 2 (teosinte). Validated below against the definitive 51K W22/teo
# coding (Chen's AA=maize/CC=teo R/qtl calls) on the 32K shared markers.
#
# Output: data/teonam/teonam_gwas118k_dosage_polar.rds (list: dos, markers, lines)
# Run: Rscript scripts/teonam_gwas118k_polarize.R
suppressMessages({
  library(data.table)
})

g <- readRDS("data/teonam/teonam_gwas118k_dosage.rds")
dos <- g$dos # markers x lines, dosage = TASSEL A2-allele count (arbitrary sign)

# per-marker frequency of the "2" allele; flip so the MINOR allele = teosinte = 2
p2 <- rowMeans(dos, na.rm = TRUE) / 2
flip <- p2 > 0.5
cat(sprintf(
  "markers flipped (A2 was major -> W22): %d / %d (%.1f%%)\n",
  sum(flip), length(flip), 100 * mean(flip)
))
dos[flip, ] <- 2L - dos[flip, ]
storage.mode(dos) <- "integer"

# --- validate against the 51K definitive W22/teo coding on shared markers ------
m51 <- fread("data/teonam/TeoNAM_genotype_clean.csv")
key51 <- paste0(sub("^W22", "", m51[[2]]), sub("^.*Line_", "", m51[[1]]))
sl <- intersect(colnames(dos), key51)
smk <- intersect(rownames(dos), names(m51))
Bi <- match(sl, key51)
set.seed(1)
samp <- sample(smk, 4000)
A <- t(dos[samp, sl, drop = FALSE]) # lines x markers (polarized 118K)
B <- as.matrix(m51[Bi, ..samp]) # lines x markers (51K W22/teo)
conc <- r <- numeric(length(samp))
for (j in seq_along(samp)) {
  a <- A[, j]
  b <- B[, j]
  ok <- !is.na(a) & !is.na(b)
  conc[j] <- if (sum(ok) >= 10) mean(a[ok] == b[ok]) else NA_real_
  r[j] <- if (sum(ok) >= 10 && sd(a[ok]) > 0 && sd(b[ok]) > 0) cor(a[ok], b[ok]) else NA_real_
}
r <- r[!is.na(r)]
conc <- conc[!is.na(conc)]
cat(sprintf(
  "vs 51K W22/teo (n=%d shared markers sampled): mean call-concordance %.3f | cor +1 %.1f%% / -1 %.1f%% (want ~all +1)\n",
  length(r), mean(conc), 100 * mean(r > 0.99), 100 * mean(r < -0.99)
))

saveRDS(
  list(dos = dos, markers = g$markers, lines = g$lines),
  "data/teonam/teonam_gwas118k_dosage_polar.rds"
)
cat("wrote data/teonam/teonam_gwas118k_dosage_polar.rds\n")
