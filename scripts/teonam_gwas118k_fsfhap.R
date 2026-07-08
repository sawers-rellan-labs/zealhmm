#!/usr/bin/env Rscript
# Block-smooth the 118K GWAS truth with TASSEL5 FSFHap, exactly as Chen 2019
# (backcross, Phet=0.03125, Fillgaps=TRUE, per family). The distributed 118K hmp
# is the imputed NUCLEOTIDE genotypes (alleles T/C, G/C, ...), which carry ~7%
# scattered single-marker het (imputation residual) -> a choppy ancestry profile
# (~800 flips/chr vs the ~3 breakpoints of a real BC1S4 RIL). Feeding the ancestry
# coding to FSFHap and taking its PARENTAL (A/C) output collapses that scatter into
# clean haplotype blocks -- the biologically-sane ancestry mosaic used as sweep truth.
# Mirrors agent/teonam_fsfhap.R (the 51K completion run), on the 118K panel.
#
# Run one family (pilot):  Rscript scripts/teonam_gwas118k_fsfhap.R TIL01
# Output: data/teonam/tassel/fsfhap118k/imputed_<fam><chr>.hmp.txt (A/C parental)
suppressMessages(library(data.table))
ROOT <- "/Users/fvrodriguez/repos/zealhmm"
setwd(ROOT)
TASSEL <- "/Applications/TASSEL 5/run_pipeline.pl"
TDIR <- file.path(ROOT, "data/teonam/tassel/fsfhap118k")
dir.create(TDIR, recursive = TRUE, showWarnings = FALSE)
fam <- commandArgs(TRUE)[1]
if (is.na(fam)) fam <- "TIL01"

# --- family ancestry skeleton from the polarized 118K dosage (0=W22,1=het,2=teo) ---
g <- readRDS("data/teonam/teonam_gwas118k_dosage_polar.rds")
dos <- g$dos # markers x lines, polarized ancestry dosage
mc <- fread("data/teonam/markers_v5_gwas118k_cm.tsv") # marker, chr, pos_v5, cm
chr_by <- setNames(mc$chr, mc$marker)
pos_by <- setNames(mc$pos_v5, mc$marker)
mk <- intersect(rownames(dos), mc$marker) # 118,514 lifted markers
mk <- mk[order(as.integer(chr_by[mk]), as.integer(pos_by[mk]))] # v5 order
keys <- colnames(dos)[substr(colnames(dos), 1, 5) == fam]
cat(sprintf("family %s: %d lines; markers %d\n", fam, length(keys), length(mk)))

M <- dos[mk, keys, drop = FALSE] # markers x famlines (ancestry 0/1/2, some NA)
cn <- matrix("N", nrow(M), ncol(M))
cn[M == 0] <- "A"
cn[M == 1] <- "M"
cn[M == 2] <- "C"
cat(sprintf("  skeleton: %.1f%% observed (het %.1f%%)\n", 100 * mean(cn != "N"), 100 * mean(cn == "M")))
hd <- data.table(
  `rs#` = mk, alleles = "A/C", chrom = as.integer(chr_by[mk]),
  pos = as.integer(pos_by[mk]), strand = "+", `assembly#` = NA, center = NA,
  protLSID = NA, assayLSID = NA, panelLSID = NA, QCcode = NA
)
hmp <- cbind(hd, as.data.table(cn))
setnames(hmp, c(names(hd), keys))
HMPF <- file.path(TDIR, sprintf("skeleton_%s.hmp.txt", fam))
fwrite(hmp, HMPF, sep = "\t", quote = FALSE, na = "N")

# --- pedigree: BC1S4 backcross to W22 (0.75 W22 / 0.25 teo; F~0.9375) ---------
PED <- file.path(TDIR, sprintf("pedigree_%s.txt", fam))
ped <- data.table(
  family = fam, taxon = keys, parent1 = "W22",
  parent2 = paste0(fam, "_teo"), p1 = 0.75, p2 = 0.25, F = 0.9375
)
fwrite(ped, PED, sep = "\t", quote = FALSE, col.names = TRUE)

# --- FSFHap (Chen's params) ---------------------------------------------------
OUT <- file.path(TDIR, sprintf("imputed_%s", fam))
LOG <- file.path(TDIR, sprintf("fsfhap_%s.log", fam))
cmd <- sprintf(
  '"%s" -Xmx12g -h "%s" -FSFHapImputationPlugin -pedigrees "%s" -bc true -phet 0.03125 -fillgaps true -maxMissing 1.0 -minMaf 0.0 -minR 0.0 -logfile "%s" -endPlugin -export "%s" -exportType Hapmap',
  TASSEL, HMPF, PED, LOG, OUT
)
cat("running FSFHap ...\n")
t <- system.time(rc <- system(cmd))
cat(sprintf("FSFHap %s: rc=%d, %.1f s (%.2f min)\n", fam, rc, t["elapsed"], t["elapsed"] / 60))
cat("outputs:", paste(list.files(TDIR, pattern = sprintf("imputed_%s", fam)), collapse = ", "), "\n")
