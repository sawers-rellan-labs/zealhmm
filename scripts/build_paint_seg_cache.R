#!/usr/bin/env Rscript
# Build paint_seg_cache/{nnil,bbnil,rtiger}.rds for the coverage-sweep chromosome painting
# (scripts/fig_coverage_sweep_chr_paint.R) by fitting each SKIM-COUNT caller as a SINGLE POOL on the
# 330-NIL dnarna cohort (+ the B73 check, so the control lane decodes) at its operating point, then
# selecting the 11 coverage-sweep samples. This is the reproducible replacement for the old basis,
# which fit each caller on just the 11 (coverage-spanning, unrepresentative) sweep samples.
#
# Counts come from the unified GATK store (data/zeal/snp50k_counts). Operating points from
# reference/zeal_skim_calibration_operating_points.csv (rtiger rigidity; bbnil rrate). nnil is not in
# that CSV; its map-derived rrate + nir=0.70 are the values used throughout (fig_zeal_fragment_size_cm_
# validation.R). Samples below the r per-chromosome floor (2*rigidity covered markers) are dropped.
#
# NOT produced here (different modality/input; handled by their own caches):
#   binhmm  -> 1 Mb bins (data/skimsweep/skim/bins), not SNP50K counts
#   googa/atlas -> BrB RNA-seq pangene (googa_nir_seg_cache/, atlas_nir_seg_cache/)
#
#   Rscript scripts/build_paint_seg_cache.R    # sandbox DISABLED (mount not needed; store is local)
suppressMessages({
  library(data.table)
  library(here)
})
devtools::load_all(path.expand("~/repos/nilhmm"), quiet = TRUE)
source(here::here("scripts/logging.R"))

KEEP <- c("name", "chr", "start_bp", "end_bp", "state")
SEG <- here::here("results/sim/zeal_nil/paint_seg_cache")
dir.create(SEG, showWarnings = FALSE, recursive = TRUE)
OPS <- fread(here::here("reference/zeal_skim_calibration_operating_points.csv"))
R_RTIGER <- as.integer(OPS[caller == "rtiger", value])
RRATE <- 6.36e-4 # map-derived skim rrate (nnil/bbnil); matches fig_zeal_fragment_size_cm_validation.R
threads <- min(parallel::detectCores() - 2L, 8L)

roster <- fread(here::here("data/zeal/snp50k_count_roster.tsv"))
pool_ids <- roster[in_dnarna == TRUE | in_dnasweep == TRUE, sample] # dnarna NILs + B73 check
sweep_ids <- roster[in_dnasweep == TRUE, sample]
long <- rbindlist(lapply(pool_ids, function(s) {
  cf <- fread(file.path(here::here("data/zeal/snp50k_counts"), paste0(s, ".tsv")),
    header = FALSE, col.names = c("contig", "pos", "rb", "rc", "ab", "ac")
  )
  data.table(
    name = s, chr = as.integer(sub("chr", "", cf$contig)), pos = as.integer(cf$pos),
    n_ref = as.integer(cf$rc), n_alt = as.integer(cf$ac)
  )
}))
# r=rigidity per-chromosome floor (rtiger needs 2*rigidity covered markers on every chr)
cov <- long[n_ref + n_alt > 0, .N, by = .(name, chr)][, .(minc = min(N), nchr = .N), by = name]
keep <- cov[minc >= 2L * R_RTIGER & nchr == 10L, name]
long <- long[name %in% keep]
log_info(
  "[paint-cache] pool %d -> fittable %d (r=%d floor); sweep in fit %d/%d",
  length(pool_ids), length(keep), R_RTIGER, length(intersect(sweep_ids, keep)), length(sweep_ids)
)

pn <- setNames(roster$skim_prefix, roster$sample)
save_sweep <- function(seg, tag) {
  s <- as.data.table(seg)[, ..KEEP][name %in% sweep_ids]
  s[, name := pn[name]]
  saveRDS(as.data.frame(s), file.path(SEG, paste0(tag, ".rds")))
  log_info("[paint-cache] %s: %d sweep samples, %d segments", tag, uniqueN(s$name), nrow(s))
}

# ---- fit each skim-count caller as a single pool, then select the sweep samples ----
log_info("[paint-cache] nnil (design-prior hard call -> HMM, per-sample; pool = select)")
g <- call_gt(long$n_ref, long$n_alt, prior = breeding_prior("BC2S3"), error = 0.01, return = "call")
g[is.na(g)] <- 3L
save_sweep(call_ancestry(data.frame(name = long$name, chr = long$chr, pos = long$pos, g = as.integer(g)),
  caller = "nnil", rrate = RRATE, nir = 0.70, germ = 1e-4, gert = 1e-2, p = 0.1,
  design = "BC2S3", threads = threads
), "nnil")

log_info("[paint-cache] bbnil (fit_means emission, pooled)")
save_sweep(call_ancestry(as.data.frame(long),
  caller = "bbnil", rrate = RRATE, fit_means = TRUE, conc = 20, err = 0.01,
  design = "BC2S3", parallel = TRUE, threads = threads
), "bbnil")

log_info("[paint-cache] rtiger (native BetaBinomial EM, pooled, r=%d)", R_RTIGER)
save_sweep(call_ancestry(as.data.frame(long),
  caller = "rtiger", rigidity = R_RTIGER, design = "BC2S3", threads = threads
), "rtiger")

log_info("[paint-cache] done. binhmm/googa/atlas caches unchanged (different modality).")
