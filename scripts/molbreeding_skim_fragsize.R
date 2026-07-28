#!/usr/bin/env Rscript
# Donor-fragment-size distributions for the ECDF panels of the ZEAL MolBreeding->skim
# calibration notebook (foil Section 4 analog). Reuses the calibration setup -- truth =
# nnil on the MolBreeding hard calls; test = skim ML (flat/argmax-GL) hard call -> nnil --
# and saves pooled donor-block sizes for reference configs: an r ladder at the ML
# operating-point emission, and a nir ladder at the map-derived r, plus the truth.
# Targeted (a handful of configs), NOT the full grid.
# Out: data/zeal/molbreeding_skim_fragsize.csv (series, size_mb)
suppressMessages({
  library(devtools)
  library(data.table)
})
devtools::load_all(path.expand("~/repos/nilhmm"), quiet = TRUE)
root <- here::here()
for (f in list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)) source(f)
source(file.path(root, "scripts/logging.R"))
N_CORES <- min(parallel::detectCores() - 2L, 8L)
DESIGN <- "BC2S3"
MB <- file.path(root, "data/zeal/molbreeding")
SK <- file.path(root, "data/zeal/skim")
R_MOLB <- 1.67e-3 # molbreeding-panel map-derived r (truth caller)
R_MAP_SKIM <- 3.2e-4 # SNP50K map-derived r (skim; the nir ladder holds r here)
# ML operating-point emission (best-by-mismatch skim config)
ML <- list(nir = 0.7, germ = 1e-4, gert = 1e-2, p = 0.1)
KEEP <- c("name", "chr", "start_bp", "end_bp", "state")

pair <- fread(file.path(root, "data/zeal/correspondence/calibration_pairing.csv"))[in_calibration == TRUE]

# --- truth: nnil on MolBreeding hard calls (fixed config), relabel to skim id ---
hc <- fread(file.path(MB, "molbreeding_hardcalls_wsfilt.tsv"))[name %in% pair$truth_sample]
hc[is.na(g), g := 3L]
truth_seg <- as.data.table(caller_grid(hc[, .(name, chr, pos, g)],
  caller = "nnil",
  emission_grid = data.frame(nir = 0.7, germ = 1e-3, gert = 1e-4, p = 0.9, mr = 0),
  rrate = R_MOLB, design = DESIGN, threads = N_CORES
))
truth_seg[, name := pair$test_sample[match(name, pair$truth_sample)]]
truth_sizes <- donor_block_sizes(truth_seg[, ..KEEP])

# --- skim ML (flat = argmax-GL = maximum-likelihood) hard call ---
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

# pooled donor-block sizes for one skim nnil config
skim_sizes <- function(nir, germ, gert, p, r) {
  seg <- as.data.table(caller_grid(gdat,
    caller = "nnil",
    emission_grid = data.frame(nir = nir, germ = germ, gert = gert, p = p, mr = mr_skim),
    rrate = r, design = DESIGN, threads = N_CORES
  ))
  donor_block_sizes(seg[, ..KEEP])
}

out <- list(data.table(series = "truth", panel = "both", size_mb = truth_sizes))
# r ladder at the ML emission: low / map / over-fragmenting
for (r in c(1e-6, R_MAP_SKIM, 1e-2)) {
  lab <- sprintf("r = %.0e%s", r, if (isTRUE(all.equal(r, R_MAP_SKIM))) " (map)" else "")
  out[[paste0("r_", r)]] <- data.table(
    series = lab, panel = "r",
    size_mb = skim_sizes(ML$nir, ML$germ, ML$gert, ML$p, r)
  )
}
# nir ladder at the map-derived r
for (nir in c(0.1, 0.3, 0.5, 0.7, 0.9)) {
  out[[paste0("nir_", nir)]] <- data.table(
    series = sprintf("nir = %.1f", nir), panel = "nir",
    size_mb = skim_sizes(nir, ML$germ, ML$gert, ML$p, R_MAP_SKIM)
  )
}
res <- rbindlist(out)[is.finite(size_mb) & size_mb > 0]
fwrite(res, file.path(root, "data/zeal/molbreeding_skim_fragsize.csv"))
log_info("wrote molbreeding_skim_fragsize.csv: %d donor blocks across %d series", nrow(res), uniqueN(res$series))
print(res[, .(n_blocks = .N, median_mb = round(median(size_mb), 3)), by = series])
