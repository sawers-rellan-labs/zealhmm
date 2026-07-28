#!/usr/bin/env Rscript
# MolBreeding informative-site cM grid (SNP_wsfilt: 9,157 donor-informative,
# teosinte-vs-B73 sites) for the ZEAL MolBreeding calibration (port of
# molb_sim_calibration.qmd). Native TeoNAM cM via the native Marey spline on
# pos_v5 -- the analog of scripts/zeal_snp50k_cm_grid.R (native map, NOT consensus).
#
# Marker IDs are v5-concatenations S{chr_v5}_{pos_v5} so the id always matches the
# v5 coordinate it is placed on (avoids the version-mismatch trap where a marker
# named in one assembly carries positions from another). The original v3
# MolBreeding id is retained in `marker_v3` to keep tabs on the source.
#
# Input : data/zeal/molbreeding/wideseq_keep_v5.tsv (chr, pos_v5; 9,157 sites)
#         data/zeal/molbreeding/sites_v5_SNP.tsv     (original v3 id <-> v5 map)
#         data/teonam/teonam_v5_native.tsv           (native est.map: chr_v5, pos_v5, cm)
# Output: data/zeal/molbreeding/markers_molbreeding_cm.tsv
#         (marker [S{chr}_{pos_v5}], chr, pos, cm, marker_v3)

suppressMessages({
  library(nilHMM)
  library(here)
  library(data.table)
})
source(here("scripts/logging.R"))

sites <- fread(here("data/zeal/molbreeding/wideseq_keep_v5.tsv"),
  header = FALSE, col.names = c("contig", "pos")
)
sites[, chr := as.integer(sub("^chr", "", contig))]

# original v3 id (keep tabs on the source), matched on v5 position
smap <- fread(here("data/zeal/molbreeding/sites_v5_SNP.tsv")) # marker(v3), chr_v5, pos_v5, ...
smap2 <- smap[, .(chr = chr_v5, pos = pos_v5, marker_v3 = marker)]
sites <- merge(sites, smap2, by = c("chr", "pos"), all.x = TRUE, sort = FALSE)

# v5-concatenation primary id (id matches its v5 coordinate)
sites[, marker := paste0("S", chr, "_", pos)]

# native TeoNAM cM via the native Marey spline (position-based, per chr)
nat <- fread(here("data/teonam/teonam_v5_native.tsv"))[, .(chr = chr_v5, bp = pos_v5, cm)]
fit_chr <- nat[, .N, by = chr][N >= 2L, chr]
to_cm <- bp_to_cm(nat[chr %in% fit_chr])
sites[, cm := NA_real_]
sites[chr %in% fit_chr, cm := to_cm(chr, pos)]
n_na <- sites[is.na(cm), .N]
if (n_na) log_warn("%d MolBreeding sites have no cM (chr not in native map?)", n_na)

out <- sites[!is.na(cm), .(marker, chr, pos, cm, marker_v3)][order(chr, pos)]
fwrite(out, here("data/zeal/molbreeding/markers_molbreeding_cm.tsv"), sep = "\t")
log_info(
  "MolBreeding informative cM grid: %d sites | %d with original v3 id | cM range %.1f-%.1f",
  nrow(out), out[!is.na(marker_v3), .N], min(out$cm), max(out$cm)
)
print(out[, .(n = .N, max_cm = round(max(cm), 1)), by = chr][order(chr)])
