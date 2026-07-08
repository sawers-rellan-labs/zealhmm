#!/usr/bin/env Rscript
# Assign native-map cM to all 118,514 lifted 118K-GWAS markers, mirroring the cM
# logic in scripts/teonam_rtiger_sweep.R exactly: native TeoNAM est.map cM for the
# markers it placed, and for the rest a cM interpolated ON THE NATIVE MAP via its
# monotone bp->cM Marey spline (.bp_to_cm_fun, Hyman, clamped), fit per chr to the
# placed markers. This is the union cM grid the 118K coverage sweep runs on.
#
# Output: data/teonam/markers_v5_gwas118k_cm.tsv (marker, chr, pos_v5, cm)
# Run: Rscript scripts/teonam_gwas118k_cm_grid.R
suppressMessages({
  library(data.table)
})
ROOT <- "/Users/fvrodriguez/repos/zealhmm"
setwd(ROOT)
source("R/simulate.R") # .bp_to_cm_fun
source("scripts/map_tools.R") # DEFAULT_TEONAM_MAP

mc <- fread("data/teonam/markers_v5_gwas118k.tsv") # 118,514 markers: marker, chr_v2, pos_v2, chr_v5, pos_v5
setnames(mc, "chr_v5", "chr")
nat_cm <- fread(DEFAULT_TEONAM_MAP) # native est.map cM (placed 51K markers)

mc[, cm := nat_cm$cm[match(marker, nat_cm$marker)]] # native cM; NA where not placed
n_placed <- sum(!is.na(mc$cm))
mc[, cm := {
  ok <- !is.na(cm)
  if (any(!ok) && sum(ok) >= 2L) {
    f <- .bp_to_cm_fun(data.table(bp = pos_v5[ok], cm = cm[ok])) # native Marey spline
    cm[!ok] <- f(pos_v5[!ok])
  }
  cm
}, by = chr]

out <- mc[, .(marker, chr, pos_v5, cm)][order(chr, pos_v5)]
fwrite(out, "data/teonam/markers_v5_gwas118k_cm.tsv", sep = "\t")
cat(sprintf(
  "118K cM grid: %d markers | %d native-placed + %d Marey-spline | cM range %.1f-%.1f\n",
  nrow(out), n_placed, nrow(out) - n_placed, min(out$cm), max(out$cm)
))
cat("per-chr max cM:\n")
print(out[, .(max_cm = round(max(cm), 1), n = .N), by = chr][order(chr)])
cat("wrote data/teonam/markers_v5_gwas118k_cm.tsv\n")
