#!/usr/bin/env Rscript
# Calibration metrics for a REFERENCE nnil config that is off the sweep grid:
# Holland's published nir=0.9 emission at r = 1.67e-3 (the molbreeding-panel map r),
# ML skim caller vs the MolBreeding truth. Same metric computation as
# scripts/molbreeding_skim_calibration.R; single config, so a targeted eval.
# Out: data/zeal/molbreeding_skim_holland_config.csv (one row, sweep-CSV schema).
suppressMessages({
  library(devtools)
  library(data.table)
})
devtools::load_all(path.expand("~/repos/nilhmm"), quiet = TRUE)
root <- here::here()
for (f in list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)) source(f)
N_CORES <- min(parallel::detectCores() - 2L, 8L)
DESIGN <- "BC2S3"
MB <- file.path(root, "data/zeal/molbreeding")
SK <- file.path(root, "data/zeal/skim")
KEEP <- c("name", "chr", "start_bp", "end_bp", "state")
CFG <- list(nir = 0.9, germ = 1e-3, gert = 1e-4, p = 0.9, r = 1.67e-3) # Holland's nir=0.9 at map-panel r
pair <- fread(file.path(root, "data/zeal/correspondence/calibration_pairing.csv"))[in_calibration == TRUE]

# truth: nnil on molbreeding hardcalls (fixed truth config), relabel to skim id
hc <- fread(file.path(MB, "molbreeding_hardcalls_wsfilt.tsv"))[name %in% pair$truth_sample]
hc[is.na(g), g := 3L]
grid_eval <- unique(hc[, .(chr, pos)])[order(chr, pos)]
truth <- as.data.table(caller_grid(hc[, .(name, chr, pos, g)],
  caller = "nnil",
  emission_grid = data.frame(nir = 0.7, germ = 1e-3, gert = 1e-4, p = 0.9, mr = 0),
  rrate = 1.67e-3, design = DESIGN, threads = N_CORES
))[, ..KEEP]
truth[, name := pair$test_sample[match(name, pair$truth_sample)]]
tr_sizes <- donor_block_sizes(truth)
truth_bp <- breakpoint_count(truth)

# skim ML (flat) hard call
skim <- rbindlist(lapply(pair$test_sample, function(s) {
  cf <- fread(file.path(SK, "counts", paste0(s, ".tsv")),
    header = FALSE,
    col.names = c("contig", "pos", "rb", "rc", "ab", "ac")
  )
  data.table(
    name = s, chr = as.integer(sub("chr", "", cf$contig)), pos = cf$pos,
    n_ref = as.integer(cf$rc), n_alt = as.integer(cf$ac)
  )
}))
mr_skim <- skim[, mean(n_ref + n_alt == 0)]
skim[, g := {
  gg <- call_gt(n_ref, n_alt, prior = "flat", error = 0.01)
  gg[is.na(gg)] <- 3L
  as.integer(gg)
}]
setorder(skim, name, chr, pos)
gdat <- skim[, .(name, chr, pos, g)]

seg <- as.data.table(caller_grid(gdat,
  caller = "nnil",
  emission_grid = data.frame(nir = CFG$nir, germ = CFG$germ, gert = CFG$gert, p = CFG$p, mr = mr_skim),
  rrate = CFG$r, design = DESIGN, threads = N_CORES
))[, ..KEEP]
md <- marker_dice(seg, truth, grid_eval)
pc <- md$per_class
ff <- donor_fragment_dice(seg, truth)
row <- data.table(
  prior = "ml", nir = CFG$nir, germ = CFG$germ, gert = CFG$gert, p = CFG$p, r = CFG$r,
  mismatch = 1 - md$accuracy,
  donor_marker_recall = pc[class == "donor(>0)"]$recall,
  het_recall = pc[class == "HET"]$recall, alt_recall = pc[class == "ALT"]$recall,
  donor_frag_dice = ff$dice, donor_frag_FDR = ff$fdr,
  ks_fragsize = fragment_size_ks(donor_block_sizes(seg), tr_sizes),
  breakpoint_ratio = breakpoint_count(seg) / truth_bp
)
fwrite(row, file.path(root, "data/zeal/molbreeding_skim_holland_config.csv"))
print(row)
