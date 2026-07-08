#!/usr/bin/env Rscript
# Lift the Chen 2019 118,838-SNP TeoNAM GWAS panel AGPv2 -> v5, mirroring the
# 51K map-panel liftover (map-neutral: roster + positions only, no cM).
#
# Source: data/teonam/9250682/W22TILXX_Chr1-10.impute_filter_MR0.2_MAF0.05.hmp.txt
#   (Qiuyue Chen's Drive release; HapMap format, 1257 RILs). The marker NAME is
#   the AGPv2 identifier S<chr>_<pos>; the file's own `pos` column is AGPv4
#   (CrossMap-lifted, per Chen 2019 Methods) -- we IGNORE the v4 pos col and lift
#   the v2 name-embedded positions v2->v4->v5 with the same chains as the 51K set,
#   so the GWAS lives in the repo's v5 world.
# Output: data/teonam/markers_v5_gwas118k.tsv (marker, chr_v2, pos_v2, chr_v5, pos_v5)
# Run: Rscript scripts/teonam_lift_gwas118k.R
suppressMessages({
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
})
source("R/teonam_liftover.R")

HMP <- "data/teonam/9250682/W22TILXX_Chr1-10.impute_filter_MR0.2_MAF0.05.hmp.txt"
V2_INFO <- "data/teonam/markers_v2_gwas118k.csv" # intermediate marker_info (v2)
OUT <- "data/teonam/markers_v5_gwas118k.tsv"

# --- read only the marker names (rs#) from the 300 MB HapMap -------------------
# NOTE: the HapMap `chrom` column is v4-aligned; the AGPv2 identifier lives in the
# NAME (S<chr>_<pos>, per Chen's README). We drive the v2->v5 lift entirely from
# the name -- chromosome AND position -- exactly as the 51K map lift does. This
# correctly drops v2-unplaced `S0_` scaffolds (no v2 chr0 in the chains).
info <- fread(HMP, select = "rs#", colClasses = "character")
setnames(info, "rs#", "name")
cat("HapMap markers read:", nrow(info), "\n")

# chromosome parsed from the NAME, not the (v4) chrom column
name_chr <- sub("_.*$", "", sub("^S", "", info$name))
n_s0 <- sum(name_chr == "0")
cat("v2-unplaced (S0_) markers, will be dropped:", n_s0, "\n")

# marker_info shape liftover_teonam expects: chromosome, name, start, end
mi <- data.frame(
  chromosome = name_chr, name = info$name,
  start = 0L, end = 0L, stringsAsFactors = FALSE
)
fwrite(mi, V2_INFO)

# --- v2 -> v4 -> v5, keep unique 1:1, same-chromosome, chr 1-10 ---------------
res <- liftover_teonam(V2_INFO)
n_in <- attr(res, "n_in")

fwrite(res, OUT, sep = "\t")
cat(sprintf(
  "lifted %d / %d markers to v5 (%.1f%%); dropped %d (multimap / chr-change / off 1-10)\n",
  nrow(res), n_in, 100 * nrow(res) / n_in, n_in - nrow(res)
))
cat("wrote", OUT, "\n")
