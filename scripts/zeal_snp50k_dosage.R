#!/usr/bin/env Rscript
# ZEAL/BZea Phase 3 — SNP50K teosinte-ancestry dosage base (analog of TeoNAM
# teonam_gwas118k_dosage.R + _polarize.R). The SNP50K hard GT is unusable (70%
# missing, 0% ALT-hom at ~0.4x), so — like the bzeaseq callers — we build from the
# per-sample allele-count tables. Sites are the 51,991 HQ teosinte-vs-B73 panel;
# ref = B73 allele, alt = teosinte-informative allele -> teo dosage = 2*alt/(ref+alt),
# already polarized (0 = B73, 2 = teosinte). No liftover (already B73 v5).
#
# Input : samplesheet (in_snp50k skim_ids) + Hazel allelic_counts/<skim_id>_allele_counts.tsv
# Output: data/zeal/zeal_snp50k_dosage.rds  list(markers, n_ref, n_alt, dosage, cov, samples)
#         data/zeal/markers_snp50k_v5.tsv   (marker, chr, pos)
# QC    : B73 check samples must show ~0 teosinte dosage (validates polarization).

suppressMessages({
  library(here)
  library(data.table)
})
source(here("scripts/logging.R"))

AC_DIR <- Sys.getenv(
  "SNP50K_AC_DIR",
  "/Volumes/rsstu/users/r/rrellan/BZea/bzeaseq/50K/results/allelic_counts"
)
OUT_RDS <- here("data/zeal/zeal_snp50k_dosage.rds")
OUT_MARKERS <- here("data/zeal/markers_snp50k_v5.tsv")

ss <- fread(here("data/zeal/samplesheet_3way.csv"))
samp <- ss[in_snp50k == TRUE & !is.na(skim_id), unique(skim_id)]
log_info("panel samples (in_snp50k): %d", length(samp))
files <- file.path(AC_DIR, paste0(samp, "_allele_counts.tsv"))
ok <- file.exists(files)
if (any(!ok)) log_warn("%d count files missing (e.g. %s)", sum(!ok), basename(files[!ok][1]))
samp <- samp[ok]
files <- files[ok]

# roster from the first file (all files are bcftools-query'd over the same cohort
# positions in the same order); assert row-count consistency per file below
cn <- c("chr", "pos", "ref", "n_ref", "alt", "n_alt")
first <- fread(files[1], header = FALSE, col.names = cn)
markers <- first[, .(chr = as.integer(sub("^chr", "", chr)), pos = as.integer(pos), ref, alt)]
markers <- markers[chr %in% 1:10]
keep <- first[, chr %in% paste0("chr", 1:10)]
M <- sum(keep)
log_info("panel sites (chr1-10): %d", M)
markers[, marker := sprintf("S%d_%d", chr, pos)]

n_ref <- matrix(0L, nrow = M, ncol = length(samp), dimnames = list(markers$marker, samp))
n_alt <- matrix(0L, nrow = M, ncol = length(samp), dimnames = list(markers$marker, samp))
t0 <- Sys.time()
for (j in seq_along(files)) {
  d <- fread(files[j], header = FALSE, select = c(4, 6), col.names = c("nr", "na"))
  if (nrow(d) != nrow(first)) {
    log_warn("row mismatch in %s (%d vs %d) - skipped", basename(files[j]), nrow(d), nrow(first))
    next
  }
  n_ref[, j] <- d$nr[keep]
  n_alt[, j] <- d$na[keep]
  if (j %% 200 == 0) {
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    log_info(">>> %d/%d | elapsed %.1f min | ETA ~%.1f min", j, length(files), el, el / j * (length(files) - j))
  }
}

cov <- n_ref + n_alt
dosage <- ifelse(cov > 0, 2 * n_alt / cov, NA_real_) # teosinte-allele dosage (0=B73, 2=teo)
log_info(
  "built dosage %d x %d | mean coverage %.2f | cells covered %.1f%%",
  nrow(dosage), ncol(dosage), mean(cov), 100 * mean(cov > 0)
)

# QC: B73 checks must be ~0 teosinte dosage
b73 <- ss[in_snp50k == TRUE & is_B73 == TRUE, skim_id]
b73 <- intersect(b73, colnames(dosage))
if (length(b73)) {
  log_info(
    "QC polarization: B73 checks (n=%d) mean teo dosage = %.3f (expect ~0)",
    length(b73), mean(dosage[, b73], na.rm = TRUE)
  )
}
nil <- ss[gwas_nil == TRUE, skim_id]
nil <- intersect(nil, colnames(dosage))
log_info(
  "QC: NIL panel (n=%d) mean teo dosage = %.3f (expect low, BC2S3 ~0.06-0.12)",
  length(nil), mean(dosage[, nil], na.rm = TRUE)
)

fwrite(markers[, .(marker, chr, pos)], OUT_MARKERS, sep = "\t")
saveRDS(
  list(markers = markers, n_ref = n_ref, n_alt = n_alt, dosage = dosage, cov = cov, samples = samp),
  OUT_RDS
)
log_info("wrote %s and %s", OUT_MARKERS, OUT_RDS)
