#!/usr/bin/env Rscript
# ZEAL/BZea Phase 3 — cM grid for the FULL SNP50K panel, on the TeoNAM native map
# (analog of teonam_gwas118k_cm_grid.R). The GWAS is NEVER thinned — it runs on the
# full marker set. (Thinning is only an ancestry-inference speed-up for some coverage-
# sweep callers, back-projected to the full set; that belongs to the deferred Phase 5,
# not here.) This grid just attaches native-map cM to every SNP50K marker for the
# map-aware ancestry callers (LB-Impute cM transitions, etc.).
#
# The SNP50K sites are teosinte-informative positions (not TeoNAM marker ids), so every
# marker's cM is interpolated from bp via the native map's monotone bp->cM Marey spline
# (.bp_to_cm_fun, Hyman, clamped), fit per chromosome to the native placed markers.
# Genetic-distance reference = the TeoNAM native v5 map (user-specified), already local.
#
# Input : data/zeal/markers_snp50k_v5.tsv (marker, chr, pos)
#         data/teonam/teonam_v5_native.tsv (native est.map: chr_v5, pos_v5, cm)
# Output: data/zeal/markers_snp50k_cm.tsv  (marker, chr, pos, cm) -- FULL grid, no thinning

suppressMessages({
  library(here)
  library(data.table)
})
source(here("scripts/logging.R"))
source(here("R/simulate.R")) # .bp_to_cm_fun

mk <- fread(here("data/zeal/markers_snp50k_v5.tsv")) # marker, chr, pos
nat <- fread(here("data/teonam/teonam_v5_native.tsv"))[, .(chr = chr_v5, bp = pos_v5, cm)]
log_info(
  "SNP50K markers=%d | native map placed=%d (chr %s)", nrow(mk), nrow(nat),
  paste(range(nat$chr), collapse = "-")
)

# per-chromosome bp -> cM via the native Marey spline
mk[, cm := {
  nm <- nat[chr == .BY$chr]
  if (nrow(nm) >= 2L) .bp_to_cm_fun(nm[, .(bp, cm)])(pos) else NA_real_
}, by = chr]
n_na <- mk[is.na(cm), .N]
if (n_na) log_warn("%d markers have no cM (chr not in native map?)", n_na)

grid <- mk[!is.na(cm), .(marker, chr, pos, cm)][order(chr, pos)]
fwrite(grid, here("data/zeal/markers_snp50k_cm.tsv"), sep = "\t")
log_info("cM grid (FULL, no thinning): %d markers | cM range %.1f-%.1f", nrow(grid), min(grid$cm), max(grid$cm))
print(grid[, .(n = .N, max_cm = round(max(cm), 1)), by = chr][order(chr)])
