#!/usr/bin/env Rscript
# Convert the Chen 2019 118,838-SNP GWAS HapMap (single-char IUPAC, TASSEL) to an
# additive dosage matrix for the STAM GWAS scan.
#
# Encoding: per marker the `alleles` column gives A1/A2; dosage = number of A2
# copies in the diploid IUPAC call:
#   hom A1 -> 0, het -> 1, hom A2 -> 2, N -> NA.
# Polarization (which allele is A2) is TASSEL's, i.e. ARBITRARY w.r.t.
# W22/teosinte. This is fine for the STAM GLM F-test (1 df), which is invariant
# to allele flip. True W22<->teosinte polarization is deferred to the
# coverage/ancestry sweep step (see scripts/teonam_gwas118k_polarize.R, TODO).
#
# Output: data/teonam/teonam_gwas118k_dosage.rds  (list: dos = integer matrix
#   [markers x lines], markers, lines)  -- lines keyed TIL<pop><id> e.g. TIL01A001.
# Run: Rscript scripts/teonam_gwas118k_dosage.R
suppressMessages({
  library(data.table)
  library(parallel)
})

HMP <- "data/teonam/9250682/W22TILXX_Chr1-10.impute_filter_MR0.2_MAF0.05.hmp.txt"
OUT <- "data/teonam/teonam_gwas118k_dosage.rds"

cat("reading HapMap ...\n")
dt <- fread(HMP, colClasses = "character")
markers <- dt[["rs#"]]
alleles <- dt[["alleles"]]
a1 <- substr(alleles, 1, 1)
a2 <- substr(alleles, 3, 3)
lines <- names(dt)[-(1:11)] # 1257 sample IDs, already TIL01A001 form
cat("markers:", length(markers), " lines:", length(lines), "\n")

# IUPAC single-char -> the two diploid bases
base1 <- c(A = "A", C = "C", G = "G", T = "T", R = "A", Y = "C", S = "C", K = "G", M = "A", W = "A")
base2 <- c(A = "A", C = "C", G = "G", T = "T", R = "G", Y = "T", S = "G", K = "T", M = "C", W = "T")

geno <- as.matrix(dt[, -(1:11)]) # character [markers x lines]
rm(dt)
gc()

# dosage per column = (base1(call)==a2) + (base2(call)==a2); N -> NA
cat("recoding to dosage ...\n")
cols <- mclapply(seq_len(ncol(geno)), function(j) {
  call <- geno[, j]
  b1 <- base1[call]
  b2 <- base2[call]
  d <- (b1 == a2) + (b2 == a2) # 0/1/2
  d[call == "N" | is.na(b1)] <- NA_integer_
  as.integer(d)
}, mc.cores = max(1L, detectCores() - 2L))

dos <- do.call(cbind, cols)
dimnames(dos) <- list(markers, lines)
rm(geno, cols)
gc()

miss <- mean(is.na(dos))
cat(sprintf("dosage matrix: %d x %d | missing: %.3f%%\n", nrow(dos), ncol(dos), 100 * miss))
saveRDS(list(dos = dos, markers = markers, lines = lines), OUT)
cat("wrote", OUT, "\n")
