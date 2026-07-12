#!/usr/bin/env Rscript
# ZEAL/BZea — SPAD leaf greenness at the two measurement dates (20DAS, 36DAS), USED AS-IS.
# Source: agent/introfinder/Bzea_merged.csv (the relabeled phenotype file the zealbrowser reads).
# These two columns are taken as provided — the only aggregation is a per-genotype MEAN over a
# line's rows (needed for one value per genotype in GWAS); NO spatial correction is applied here
# (the merged file has no plot/row/range; the values may already have been adjusted upstream by HH).
# This is deliberately separate from the SpATS single-SPAD analysis (zeal_spats_blues.R), which
# spatially corrects the fieldbook SPAD column.
# Writes pheno_<ttag>_direct.csv (Genotype, <TRAIT>_mean) + tassel/pheno_<ttag>_all.txt.
# Traits: SPAD20DAS, SPAD36DAS. GWAS: TRAIT=SPAD20DAS/SPAD36DAS PHENO=direct.
suppressMessages({
  library(here)
  library(data.table)
})
source(here("scripts/logging.R"))

d <- fread(here("agent/introfinder/Bzea_merged.csv"), na.strings = c("#N/A", "NA", ""))
setnames(d, "Lines", "Genotype")
d <- d[!is.na(Genotype) & Genotype != ""]
ss <- fread(here("data/zeal/samplesheet_3way.csv"))

spec <- c(SPAD20DAS = "SPAD_leaf_greenness_20DAS", SPAD36DAS = "SPAD_leaf_greenness_36DAS")
for (tr in names(spec)) {
  col <- spec[[tr]]
  v <- suppressWarnings(as.numeric(d[[col]]))
  m <- data.table(Genotype = d$Genotype, v = v)[is.finite(v), .(mn = mean(v)), by = Genotype]
  setnames(m, "mn", paste0(tr, "_mean"))
  fwrite(m, here(sprintf("data/zeal/pheno_%s_direct.csv", tolower(tr))))
  tass <- merge(ss[gwas_nil == TRUE, .(pedigree, taxon)],
    m[, .(pedigree = Genotype, y = get(paste0(tr, "_mean")))],
    by = "pedigree"
  )[is.finite(y)]
  ph_out <- here(sprintf("data/zeal/tassel/pheno_%s_all.txt", tolower(tr)))
  writeLines(c("<Phenotype>", "taxa\tdata\tfactor", sprintf("Taxa\t%s\tFamily", tr)), ph_out)
  fwrite(tass[, .(pedigree, round(y, 4), taxon)], ph_out, sep = "\t", append = TRUE, col.names = FALSE)
  log_info(
    "%s (as-is, per-line mean): %d lines (range %.1f-%.1f) | %d gwas_nil in TASSEL",
    tr, nrow(m), min(m[[2]]), max(m[[2]]), nrow(tass)
  )
}
