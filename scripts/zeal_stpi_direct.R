#!/usr/bin/env Rscript
# ZEAL/BZea — StPi (stem anthocyanin) phenotype used DIRECTLY, no spatial correction.
# StPi is a 0-1 pigment score; a SpATS genotype-fixed P-spline is the wrong model for it
# (the SpATS BLUEs ran out of range, e.g. 13 genotypes > 1). Here the phenotype is the raw
# per-genotype mean of the plot scores: mean within each field, then mean across fields
# (same field combination as the BLUE pipeline, only SpATS -> raw mean).
# Reuses the exact field manifests from zeal_spats_blues.R (plot_id -> Genotype mapping).
# Outputs: data/zeal/pheno_stpi_direct.csv (Genotype, StPi_cly23, StPi_cly25, StPi_mean)
#          data/zeal/tassel/pheno_stpi_all.txt  (TASSEL: Taxa | StPi | Family=taxon, direct)
suppressMessages({
  library(here)
  library(data.table)
  library(readxl)
})
source(here("scripts/logging.R"))
EXCEL_1970 <- 25569L
canon_ped <- function(x) sub("\\.B$", "", x)

manifest_cly25 <- function() {
  fm <- fread(here("data/zeal/cly25_b5_fieldmap.csv"))
  ph <- as.data.table(read_excel(here("data/zeal/CLY25-Fieldbook.xlsx"), sheet = "B5_BZea_eval"))
  setnames(ph, 1, "plot_id")
  ph[, plot_id := suppressWarnings(as.integer(plot_id))]
  ph[, `:=`(StPi = as.numeric(StPi), Genotype = canon_ped(Description))]
  merge(fm, ph[, .(plot_id, Genotype, StPi)], by = "plot_id")
}

manifest_cly23 <- function() {
  fm <- fread(here("data/zeal/cly23_d4_fieldmap.csv"))
  ph <- as.data.table(read_excel(here("data/zeal/CLY23_D4_FieldBook.xlsx"), sheet = "UPDATED_CLY23_D4_FieldBook"))
  setnames(ph, "CLY23_D4", "plot_id")
  ph[, plot_id := suppressWarnings(as.integer(plot_id))]
  gc <- as.data.table(read_excel(here("data/zeal/CLY23_D4_FieldBook.xlsx"), sheet = "GENOTYPE-CONVERSION"))
  o2n <- unique(rbind(
    gc[!is.na(oldold_genotype), .(old = oldold_genotype, new = new_genotype)],
    gc[!is.na(old_genotype), .(old = old_genotype, new = new_genotype)]
  ))[, .SD[1], by = old]
  ph <- merge(ph, o2n, by.x = "old_genotype", by.y = "old", all.x = TRUE)
  ph[, Genotype := fcase(
    Species == "B73" | old_genotype == "B73", "B73",
    Species == "Check", "Purple",
    !is.na(new), canon_ped(new),
    default = old_genotype
  )]
  ph[, StPi := as.numeric(StPi)]
  merge(fm, ph[, .(plot_id, Genotype, StPi)], by = "plot_id")
}

man23 <- manifest_cly23()
man25 <- manifest_cly25()

# raw per-genotype mean within each field
field_mean <- function(man) man[is.finite(StPi) & !is.na(Genotype), .(v = mean(StPi)), by = Genotype]
c23 <- field_mean(man23)
c25 <- field_mean(man25)
setnames(c23, "v", "StPi_cly23")
setnames(c25, "v", "StPi_cly25")
m <- merge(c23, c25, by = "Genotype", all = TRUE)
m[, StPi_mean := rowMeans(.SD, na.rm = TRUE), .SDcols = c("StPi_cly23", "StPi_cly25")]
fwrite(m, here("data/zeal/pheno_stpi_direct.csv"))
log_info(
  "direct StPi: %d genotypes | mean %.3f sd %.3f range [%.3f, %.3f]",
  nrow(m), mean(m$StPi_mean, na.rm = TRUE), sd(m$StPi_mean, na.rm = TRUE),
  min(m$StPi_mean, na.rm = TRUE), max(m$StPi_mean, na.rm = TRUE)
)

# --- empirical-logit phenotype -----------------------------------------------
# StPi is a binary 0/1 plot score; the per-genotype phenotype is the proportion of
# pigmented plots, k/n. Model it on the logit scale via the Haldane-Anscombe empirical
# logit  log((k+0.5)/(n-k+0.5))  (finite at k=0 and k=n), pooling plots across fields.
plots <- rbind(man23, man25, fill = TRUE)[is.finite(StPi) & !is.na(Genotype)]
el <- plots[, .(k = sum(StPi), n = .N), by = Genotype]
el[, prop := k / n][, StPi_mean := log((k + 0.5) / (n - k + 0.5))]
fwrite(el[, .(Genotype, k, n, prop, StPi_mean)], here("data/zeal/pheno_stpi_elogit.csv"))
log_info(
  "elogit StPi: %d genotypes | %d ever-pigmented (k>0) | elogit range [%.2f, %.2f]",
  nrow(el), sum(el$k > 0), min(el$StPi_mean), max(el$StPi_mean)
)

# TASSEL phenotype (gwas_nil lines, Family=taxon). PHENO env picks direct or elogit;
# elogit is the modeled StPi phenotype, so it is the default written to the TASSEL file.
PHENO <- Sys.getenv("PHENO", "elogit")
src <- if (PHENO == "elogit") el[, .(Genotype, StPi_mean)] else m[, .(Genotype, StPi_mean)]
ss <- fread(here("data/zeal/samplesheet_3way.csv"))
tass <- merge(ss[gwas_nil == TRUE, .(pedigree, taxon)], src[, .(pedigree = Genotype, y = StPi_mean)], by = "pedigree")[is.finite(y)]
ph_out <- here("data/zeal/tassel/pheno_stpi_all.txt")
writeLines(c("<Phenotype>", "taxa\tdata\tfactor", "Taxa\tStPi\tFamily"), ph_out)
fwrite(tass[, .(pedigree, round(y, 4), taxon)], ph_out, sep = "\t", append = TRUE, col.names = FALSE)
log_info("wrote %s (%s phenotype, %d gwas_nil lines)", ph_out, PHENO, nrow(tass))
