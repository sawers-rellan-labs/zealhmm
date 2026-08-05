#!/usr/bin/env Rscript
# Marker-coverage QC across the staged 332-NIL paired cohort.
#
# MOTIVATION (Fausto): RTIGER's internal QC refusing a line is a REAL finding, not an inconvenience.
# It rejects any (sample, chromosome) chain with < 2*rigidity covered markers, and it is the ONLY
# caller that refuses: nnil, bbnil, googa, atlas and binhmm all emit fragments for such a line, from
# essentially no evidence. So the question is not "how do we get rtiger to run" but "which lines are
# too thin to call AT ALL, for every caller".
#
# Covered marker = a SNP50K site with n_ref + n_alt > 0 (any read). Reports the per-line, per-chromosome
# distribution, flags lines below a ladder of thresholds, and writes a per-line QC table so the figure
# script can exclude thin lines consistently across ALL SIX callers instead of per caller.
#
#   Rscript scripts/zeal_paired_cohort_coverage_qc.R

suppressMessages({
  library(data.table)
  library(here)
})
source(here::here("scripts/logging.R"))

COH <- here::here("data/zeal/paired_cohort")
OUT <- here::here("results/sim/zeal_nil/paired_cohort")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
RIGIDITY <- 2L # rtiger operating point; its QC floor is 2*rigidity covered markers per chain

man <- fread(file.path(COH, "cohort_manifest.csv"))
setorder(man, pedigree)
log_info("[qc] scanning marker coverage for %d NILs", nrow(man))

t0 <- Sys.time()
per_chr <- rbindlist(lapply(seq_len(nrow(man)), function(i) {
  p <- man$pedigree[i]
  d <- fread(file.path(COH, "skim/counts_50k", sprintf("%s.tsv", p)),
    header = FALSE, col.names = c("contig", "pos", "rb", "rc", "ab", "ac")
  )
  d[, chr := as.integer(sub("chr", "", contig))]
  r <- d[, .(
    n_sites = .N,
    n_covered = sum(rc + ac > 0L),
    n_donor_reads = sum(ac),
    depth = sum(rc + ac)
  ), by = chr]
  r[, pedigree := p]
  if (i %% 50 == 0 || i == nrow(man)) {
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    log_info(">>> %d/%d NILs | elapsed %.1f min | ETA ~%.1f min", i, nrow(man), el, (el / i) * (nrow(man) - i))
  }
  r[]
}))

# ---- per-line summary -------------------------------------------------------
per_line <- per_chr[, .(
  n_chr = .N,
  covered_total = sum(n_covered),
  covered_frac = sum(n_covered) / sum(n_sites),
  min_chr_covered = min(n_covered),
  chr_below_rtiger = sum(n_covered < 2L * RIGIDITY),
  chr_below_10 = sum(n_covered < 10L),
  chr_below_100 = sum(n_covered < 100L),
  mean_depth = sum(depth) / sum(n_sites)
), by = pedigree]
setorder(per_line, covered_total)

log_info(
  "[qc] cohort-wide covered-marker fraction: median %.3f | min %.4f | max %.3f",
  median(per_line$covered_frac), min(per_line$covered_frac), max(per_line$covered_frac)
)
log_info("[qc] quantiles of covered markers per line:")
print(round(quantile(per_line$covered_total, c(0, .01, .05, .25, .5, .75, 1))))

log_info("[qc] lines with ANY chromosome below the rtiger floor (%d covered markers):", 2L * RIGIDITY)
print(per_line[chr_below_rtiger > 0, .(pedigree, covered_total, covered_frac = round(covered_frac, 4), min_chr_covered, chr_below_rtiger)])

log_info("[qc] 15 thinnest lines overall:")
print(head(per_line[, .(pedigree, covered_total, covered_frac = round(covered_frac, 4), min_chr_covered, chr_below_10, chr_below_100, mean_depth = round(mean_depth, 4))], 15))

log_info(
  "[qc] lines with any chromosome below 100 covered markers: %d | below 10: %d | below rtiger floor: %d",
  per_line[chr_below_100 > 0, .N], per_line[chr_below_10 > 0, .N], per_line[chr_below_rtiger > 0, .N]
)

fwrite(per_chr, file.path(OUT, "coverage_qc_per_chr.csv"))
fwrite(per_line, file.path(OUT, "coverage_qc_per_line.csv"))
log_info("[qc] wrote coverage_qc_per_chr.csv and coverage_qc_per_line.csv to %s", OUT)
