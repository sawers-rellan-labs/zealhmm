#!/usr/bin/env Rscript
# Assign native-map cM to all 118,514 lifted 118K-GWAS markers, mirroring the cM
# logic in scripts/teonam_rtiger_sweep.R exactly: native TeoNAM est.map cM for the
# markers it placed, and for the rest a cM interpolated ON THE NATIVE MAP via its
# monotone bp->cM Marey spline (nilHMM::bp_to_cm, Hyman, clamped), fit per chr to
# the placed markers. This is the union cM grid the 118K coverage sweep runs on.
#
# Output: data/teonam/markers_v5_gwas118k_cm.tsv (marker, chr, pos_v5, cm)
# Run: Rscript scripts/teonam_gwas118k_cm_grid.R
suppressMessages({
  library(nilHMM)
  library(data.table)
})
ROOT <- "/Users/fvrodriguez/repos/zealhmm"
setwd(ROOT)
source("scripts/map_tools.R") # DEFAULT_TEONAM_MAP

mc <- fread("data/teonam/markers_v5_gwas118k.tsv") # 118,514 markers: marker, chr_v2, pos_v2, chr_v5, pos_v5
setnames(mc, "chr_v5", "chr")
nat_cm <- fread(DEFAULT_TEONAM_MAP) # native est.map cM (placed 51K markers)

mc[, cm := nat_cm$cm[match(marker, nat_cm$marker)]] # native cM; NA where not placed
n_placed <- sum(!is.na(mc$cm))
# place est.map-unplaced markers on the NATIVE cM scale via a per-chr Marey spline
# fit on the placed markers (bp_to_cm splits by chr internally; a chr with
# <2 placed markers has no spline and its unplaced markers stay NA).
placed <- mc[!is.na(cm), .(chr, bp = pos_v5, cm)]
fit_chr <- placed[, .N, by = chr][N >= 2L, chr]
to_cm <- bp_to_cm(placed[chr %in% fit_chr])
mc[is.na(cm) & chr %in% fit_chr, cm := to_cm(chr, pos_v5)]

out <- mc[, .(marker, chr, pos_v5, cm)][order(chr, pos_v5)]
fwrite(out, "data/teonam/markers_v5_gwas118k_cm.tsv", sep = "\t")
cat(sprintf(
  "118K cM grid: %d markers | %d native-placed + %d Marey-spline | cM range %.1f-%.1f\n",
  nrow(out), n_placed, nrow(out) - n_placed, min(out$cm), max(out$cm)
))
cat("per-chr max cM:\n")
print(out[, .(max_cm = round(max(cm), 1), n = .N), by = chr][order(chr)])
cat("wrote data/teonam/markers_v5_gwas118k_cm.tsv\n")
