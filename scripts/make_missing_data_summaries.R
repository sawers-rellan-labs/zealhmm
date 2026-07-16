#!/usr/bin/env Rscript
# Build small per-sample missing-data summaries for the §1 foundations notebooks
# (analysis/missing-data-floor-model.qmd, analysis/missing-data-model-comparison.qmd).
#
# Reduces three heavy ZEAL tables to tiny per-sample tables (~1,400 rows) of
# coverage (lambda) and observed missingness (missing_obs). Designed to run on the
# HPC over NATIVE paths (submitted via scripts/make_missing_data_summaries.lsf),
# but is path-agnostic: set BZEASEQ_DIR / OUT_DIR to run anywhere (e.g. the mount).
#
#   BZEASEQ_DIR  bzeaseq root (default: /rsstu/users/r/rrellan/BZea/bzeaseq)
#   OUT_DIR      where to write the summaries (default: current directory)
#
# Outputs (SAMPLE + lambda + missing_obs; wideseq/snp50k also carry the raw sums):
#   $OUT_DIR/wgs_per_sample.tsv
#   $OUT_DIR/wideseq_per_sample.tsv
#   $OUT_DIR/snp50k_per_sample.tsv

suppressMessages(library(data.table))

say <- function(...) message(format(Sys.time(), "%H:%M:%S"), "  ", sprintf(...))

bzeaseq_dir <- Sys.getenv("BZEASEQ_DIR", "/rsstu/users/r/rrellan/BZea/bzeaseq")
out_dir <- Sys.getenv("OUT_DIR", ".")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

wgs_file <- file.path(bzeaseq_dir, "WGSmetrics_summary.tsv")
wideseq_file <- file.path(bzeaseq_dir, "ancestry", "all_samples_bin_genotypes.tsv")
snp50k_file <- file.path(bzeaseq_dir, "50K", "allelic_counts50K.tsv")
for (f in c(wgs_file, wideseq_file, snp50k_file)) {
  if (!file.exists(f)) stop("Input not found: ", f)
}
say("BZEASEQ_DIR = %s", bzeaseq_dir)
say("OUT_DIR     = %s", out_dir)

# --- Dataset 1: WGS metrics (small) ------------------------------------------
say("WGS metrics: %s", wgs_file)
wgs <- fread(wgs_file, skip = "SAMPLE", select = c("SAMPLE", "MEAN_COVERAGE", "PCT_1X"))
wgs_per_sample <- wgs[
  , .(SAMPLE, lambda = MEAN_COVERAGE, missing_obs = 1 - PCT_1X)
][is.finite(lambda) & is.finite(missing_obs)]
fwrite(wgs_per_sample, file.path(out_dir, "wgs_per_sample.tsv"), sep = "\t")
say("  -> %d samples", nrow(wgs_per_sample))

# --- Dataset 2: Wideseq 1 Mb binned counts (~309 MB) -------------------------
say("Wideseq bins: %s", wideseq_file)
ws <- fread(
  wideseq_file,
  skip = "SAMPLE",
  select = c("SAMPLE", "VARIANT_COUNT", "INFORMATIVE_VARIANT_COUNT", "DEPTH_SUM")
)
wideseq_per_sample <- ws[, .(
  VARIANT_COUNT             = sum(VARIANT_COUNT, na.rm = TRUE),
  INFORMATIVE_VARIANT_COUNT = sum(INFORMATIVE_VARIANT_COUNT, na.rm = TRUE),
  DEPTH_SUM                 = sum(DEPTH_SUM, na.rm = TRUE)
), by = SAMPLE][
  , `:=`(
    lambda = DEPTH_SUM / VARIANT_COUNT,
    missing_obs = 1 - INFORMATIVE_VARIANT_COUNT / VARIANT_COUNT
  )
]
fwrite(wideseq_per_sample, file.path(out_dir, "wideseq_per_sample.tsv"), sep = "\t")
say("  -> %d samples", nrow(wideseq_per_sample))

# --- Dataset 3: 50K per-site allelic counts (~2.4 GB) ------------------------
say("SNP50K allelic counts: %s (large — reading 3 columns)", snp50k_file)
sk <- fread(snp50k_file, skip = "SAMPLE", select = c("SAMPLE", "REF_COUNT", "ALT_COUNT"))
snp50k_per_sample <- sk[, .(
  VARIANT_COUNT             = .N,
  INFORMATIVE_VARIANT_COUNT = sum((REF_COUNT + ALT_COUNT) > 0, na.rm = TRUE),
  DEPTH_SUM                 = sum(REF_COUNT + ALT_COUNT, na.rm = TRUE)
), by = SAMPLE][
  , `:=`(
    lambda = DEPTH_SUM / VARIANT_COUNT,
    missing_obs = 1 - INFORMATIVE_VARIANT_COUNT / VARIANT_COUNT
  )
]
fwrite(snp50k_per_sample, file.path(out_dir, "snp50k_per_sample.tsv"), sep = "\t")
say("  -> %d samples", nrow(snp50k_per_sample))

say("Done. Summaries written to %s", normalizePath(out_dir))
