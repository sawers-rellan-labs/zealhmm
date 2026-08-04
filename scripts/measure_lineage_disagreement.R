#!/usr/bin/env Rscript
# Measure the REAL count disagreement between the two SNP50K count lineages, keyed on
# (sample, chr, pos). This replaces the previous agent's un-derived "~21%" assertion.
#
#   GATK lineage  : data/skimsweep/skim/counts_50k/<S>.tsv  (CollectAllelicCounts, filler ALT)
#   bcftools lineage: $BZ/50K/results/allelic_counts/<S>_allele_counts.tsv
#                     (get_allelic_counts.sh: bcftools query AD{0}/AD{1} off cohort.vcf.gz)
#
# Both are 49,002 sites, PN-named, same 6-col format (chr pos ref n_ref alt n_alt).
# We compare the numeric counts (n_ref, n_alt); the ALT-base column is filler on the GATK side
# so it is NOT part of the disagreement test.
#
#   Rscript scripts/measure_lineage_disagreement.R   # run with sandbox DISABLED (rsstu automount)

suppressMessages(library(data.table))
here <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
BZ <- "/Volumes/rsstu/users/r/rrellan/BZea/bzeaseq/50K/results"
GATK_DIR <- file.path(here, "data/skimsweep/skim/counts_50k")
CC <- c("contig", "pos", "ref", "n_ref", "alt", "n_alt")

rd <- function(p) {
  d <- fread(p, header = FALSE, col.names = CC)
  d[, .(
    chr = sub("chr", "", contig), pos = as.integer(pos),
    n_ref = as.integer(n_ref), n_alt = as.integer(n_alt)
  )]
}

samples <- sub("\\.tsv$", "", list.files(GATK_DIR, pattern = "\\.tsv$"))
res <- rbindlist(lapply(samples, function(s) {
  g <- rd(file.path(GATK_DIR, sprintf("%s.tsv", s)))
  b <- rd(file.path(BZ, "allelic_counts", sprintf("%s_allele_counts.tsv", s)))
  j <- merge(g, b, by = c("chr", "pos"), suffixes = c("_g", "_b"))
  j[, `:=`(dep_g = n_ref_g + n_alt_g, dep_b = n_ref_b + n_alt_b)]
  cov_either <- j[dep_g > 0 | dep_b > 0]
  cov_both <- j[dep_g > 0 & dep_b > 0]
  data.table(
    sample = s,
    n_sites = nrow(j),
    cov_either = nrow(cov_either),
    cov_both = nrow(cov_both),
    # disagreement among sites covered by EITHER counter
    dep_diff_either = cov_either[dep_g != dep_b, .N],
    ref_diff_either = cov_either[n_ref_g != n_ref_b, .N],
    alt_diff_either = cov_either[n_alt_g != n_alt_b, .N],
    any_diff_either = cov_either[n_ref_g != n_ref_b | n_alt_g != n_alt_b, .N],
    # disagreement among sites covered by BOTH counters
    any_diff_both = cov_both[n_ref_g != n_ref_b | n_alt_g != n_alt_b, .N]
  )
}))

res[, pct_any_either := round(100 * any_diff_either / cov_either, 1)]
res[, pct_any_both := round(100 * any_diff_both / cov_both, 1)]
print(res)

cat("\n==== POOLED across", length(samples), "shared samples ====\n")
tot_either <- sum(res$cov_either)
tot_both <- sum(res$cov_both)
cat(sprintf("sites covered by EITHER counter: %d\n", tot_either))
cat(sprintf("sites covered by BOTH counters : %d\n", tot_both))
cat(sprintf(
  "(n_ref or n_alt) disagree, of covered-EITHER: %d = %.1f%%\n",
  sum(res$any_diff_either), 100 * sum(res$any_diff_either) / tot_either
))
cat(sprintf(
  "(n_ref or n_alt) disagree, of covered-BOTH  : %d = %.1f%%\n",
  sum(res$any_diff_both), 100 * sum(res$any_diff_both) / tot_both
))
cat(sprintf(
  "total-depth disagree, of covered-EITHER     : %d = %.1f%%\n",
  sum(res$dep_diff_either), 100 * sum(res$dep_diff_either) / tot_either
))
