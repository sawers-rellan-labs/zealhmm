#!/usr/bin/env Rscript
# Stage the MolBreeding calibration cohort for panel A of the two-panel fragment-size cM figure
# ("calibration against molb truth").
#
# WHY THIS EXISTS: `data/zeal/skim/counts/` already holds skim counts for these samples, but they
# came via the zealtiger repo and are on the WRONG FRAME: 51,991 sites (vs the mount's 49,002) with
# a SYNTHETIC PLACEHOLDER ALT column (ALT = "A" unless REF = "A", then "C"). Panel A must sit on the
# same frame as panel B, so counts are re-staged from the mount. See [[snp50k-count-source-mount-frame]].
#
# Cohort = the 14 in_calibration pairs of calibration_pairing.csv. `truth_sample` is the MolBreeding
# well and `test_sample` the skim sample; they differ for the one MolBreeding well swap (truth
# PN4_SID330's DNA is actually PN4_SID322), so skim inputs are keyed on TEST_SAMPLE.
#
# Files are named by PEDIGREE (the NIL ID), and the in-file bins `SAMPLE` column is rewritten to it,
# for the same namespace reason as the paired-cohort staging: PN#_SID# is platform-local.
#
# NOTE: no BrB staging here. Only 4 of the 14 molb NILs have BrB, and all 4 are already inside the
# staged 332 paired cohort, so panel A reuses their googa/atlas segments from that decode.
#
#   Rscript scripts/stage_zeal_molb_cohort.R          # sandbox DISABLED (rsstu automount)
#   ZEAL_STAGE_FORCE=1 Rscript scripts/stage_zeal_molb_cohort.R

suppressMessages({
  library(data.table)
  library(here)
})
source(here::here("scripts/logging.R"))

BZ <- "/Volumes/rsstu/users/r/rrellan/BZea/bzeaseq"
AC <- file.path(BZ, "50K/results/allelic_counts")
CORR <- here::here("data/zeal/correspondence")
OUT <- here::here("data/zeal/molb_cohort")
FORCE <- Sys.getenv("ZEAL_STAGE_FORCE") != ""
canon <- function(x) sub("\\.B$", "", trimws(x))

for (d in c("skim/counts_50k", "skim/bins")) dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)

# ---- cohort: the 14 calibration pairs, keyed on pedigree --------------------
pr <- fread(file.path(CORR, "calibration_pairing.csv"))
log_info("[molb] calibration_pairing.csv: %d rows, %d in_calibration", nrow(pr), sum(pr$in_calibration))
coh <- pr[in_calibration == TRUE, .(
  pedigree = canon(true_pedigree), truth_sample = trimws(truth_sample), skim_prefix = trimws(test_sample), note
)]
swapped <- coh[truth_sample != skim_prefix]
if (nrow(swapped)) {
  log_info("[molb] %d pair(s) where the MolBreeding well label != the skim sample (well swap):", nrow(swapped))
  print(swapped[, .(pedigree, truth_sample, skim_prefix)])
}
if (uniqueN(coh$pedigree) != nrow(coh)) stop("duplicate pedigree in the molb calibration cohort")
log_info("[molb] cohort: %d NILs", nrow(coh))

# ---- sources: MOUNT ONLY ----------------------------------------------------
# skim counts: the BAM read tally (GATK CollectAllelicCounts) via the shared extractor + cache,
# already filtered to the 49,002-site biallelic panel. NOT 50K/results/allelic_counts/, which is
# `AD` from the joint-called cohort.vcf.gz. One counter everywhere; see DATA.md.
GATK_CACHE <- here::here("data/zeal/gatk_counts_50k")
dir.create(GATK_CACHE, recursive = TRUE, showWarnings = FALSE)
lst <- tempfile()
writeLines(sort(unique(coh$skim_prefix)), lst)
rc <- system2("bash", c(here::here("scripts/extract_gatk_counts_50k.sh"), lst, GATK_CACHE))
if (rc != 0) stop("extract_gatk_counts_50k.sh failed (rc=", rc, ")")
coh[, `:=`(
  src_counts = file.path(GATK_CACHE, sprintf("%s.tsv", skim_prefix)),
  src_bins = file.path(BZ, "ancestry", sprintf("%s_bin_genotypes.tsv", skim_prefix))
)]
coh[, `:=`(has_counts = file.exists(src_counts), has_bins = file.exists(src_bins))]
log_info("[molb] mount availability: counts %d/%d | bins %d/%d", sum(coh$has_counts), nrow(coh), sum(coh$has_bins), nrow(coh))
if (coh[!(has_counts & has_bins), .N]) {
  log_warn("[molb] dropping %d NILs missing a mount input:", coh[!(has_counts & has_bins), .N])
  print(coh[!(has_counts & has_bins), .(pedigree, skim_prefix, has_counts, has_bins)])
}
coh <- coh[has_counts & has_bins]
setorder(coh, pedigree)

# ---- stage ------------------------------------------------------------------
coh[, `:=`(
  dst_counts = file.path(OUT, "skim/counts_50k", sprintf("%s.tsv", pedigree)),
  dst_bins = file.path(OUT, "skim/bins", sprintf("%s.tsv", pedigree))
)]
n_written <- 0L
for (i in seq_len(nrow(coh))) {
  r <- coh[i]
  if (FORCE || !file.exists(r$dst_counts)) {
    file.copy(r$src_counts, r$dst_counts, overwrite = TRUE)
    n_written <- n_written + 1L
  }
  if (FORCE || !file.exists(r$dst_bins)) {
    b <- fread(r$src_bins)
    b[, SAMPLE := r$pedigree]
    fwrite(b, r$dst_bins, sep = "\t")
    n_written <- n_written + 1L
  }
  log_info("[molb] %d/%d %s (skim %s)", i, nrow(coh), r$pedigree, r$skim_prefix)
}

# ---- verify the staged frame: 49,002 biallelic panel sites, BAM-tally counts ----
# The ALT base column is FILLER on this source (GATK writes N where the panel defines no ALT and
# the extractor substitutes A, or C when REF is A), so a skewed ALT composition is EXPECTED here
# and is NOT a defect. Assert the site count, and report the composition only as a provenance
# fingerprint distinguishing this source from the bcftools-AD files, which carry a real %ALT.
chk <- rbindlist(lapply(coh$dst_counts, function(f) {
  d <- fread(f, header = FALSE)
  data.table(file = basename(f), n_sites = nrow(d), alt_A = sum(d$V5 == "A"), alt_G = sum(d$V5 == "G"))
}))
log_info("[molb] staged frame: site counts = %s (expect 49002, the biallelic panel)", paste(unique(chk$n_sites), collapse = ", "))
log_info(
  "[molb] ALT-column fingerprint: median A = %d, median G = %d (G near 0 = the expected filler column of the BAM-tally source)",
  median(chk$alt_A), median(chk$alt_G)
)
if (any(chk$n_sites != 49002L)) log_error("[molb] staged counts are NOT 49,002 biallelic panel sites")

man <- coh[, .(pedigree, truth_sample, skim_prefix, counts_50k = basename(dst_counts), bins = basename(dst_bins))]
fwrite(man, file.path(OUT, "molb_cohort_manifest.csv"))
log_info("[molb] wrote %s (%d NILs, %d files written)", file.path(OUT, "molb_cohort_manifest.csv"), nrow(man), n_written)
