#!/usr/bin/env Rscript
# ZEAL — StPu (stem pubescence / macrohairs) phenotype used DIRECTLY, no spatial correction.
# StPu is a binary 0/1 macrohair-presence plot score (CLY23 56/4389, CLY25 28/2898 positive); a
# SpATS genotype-fixed P-spline is the wrong model for it — the current SpATS pipeline returns
# "no fit" and the stale BLUEs ran out of [0,1] range (min -0.034). Here, exactly as for StPi
# (zeal_stpi_direct.R), the phenotype is the empirical logit of the per-genotype pubescent-plot
# proportion. Reuses the exact field manifests from zeal_spats_blues.R (plot_id -> Genotype).
# Outputs: data/zeal/pheno_stpu_direct.csv (Genotype, StPu_cly23, StPu_cly25, StPu_mean)
#          data/zeal/pheno_stpu_elogit.csv (Genotype, k, n, prop, StPu_mean)
#          data/zeal/tassel/pheno_stpu_all.txt  (TASSEL: Taxa | StPu | Family=taxon)
suppressMessages({
  library(here)
  library(data.table)
  library(readxl)
})
source(here("scripts/logging.R"))
canon_ped <- function(x) sub("\\.B$", "", x)

manifest_cly25 <- function() {
  fm <- fread(here("data/zeal/cly25_b5_fieldmap.csv"))
  ph <- as.data.table(read_excel(here("data/zeal/CLY25-Fieldbook.xlsx"), sheet = "B5_BZea_eval"))
  setnames(ph, 1, "plot_id")
  ph[, plot_id := suppressWarnings(as.integer(plot_id))]
  ph[, `:=`(StPu = as.numeric(StPu), Genotype = canon_ped(Description))]
  merge(fm, ph[, .(plot_id, Genotype, StPu)], by = "plot_id")
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
  ph[, StPu := as.numeric(StPu)]
  merge(fm, ph[, .(plot_id, Genotype, StPu)], by = "plot_id")
}

man23 <- manifest_cly23()
man25 <- manifest_cly25()

# raw per-genotype mean within each field
field_mean <- function(man) man[is.finite(StPu) & !is.na(Genotype), .(v = mean(StPu)), by = Genotype]
c23 <- field_mean(man23)
c25 <- field_mean(man25)
setnames(c23, "v", "StPu_cly23")
setnames(c25, "v", "StPu_cly25")
m <- merge(c23, c25, by = "Genotype", all = TRUE)
m[, StPu_mean := rowMeans(.SD, na.rm = TRUE), .SDcols = c("StPu_cly23", "StPu_cly25")]
fwrite(m, here("data/zeal/pheno_stpu_direct.csv"))
log_info(
  "direct StPu: %d genotypes | mean %.3f sd %.3f range [%.3f, %.3f]",
  nrow(m), mean(m$StPu_mean, na.rm = TRUE), sd(m$StPu_mean, na.rm = TRUE),
  min(m$StPu_mean, na.rm = TRUE), max(m$StPu_mean, na.rm = TRUE)
)

# --- empirical-logit phenotype -----------------------------------------------
# StPu is a binary 0/1 plot score; the per-genotype phenotype is the proportion of
# pubescent plots, k/n. Model it on the logit scale via the Haldane-Anscombe empirical
# logit  log((k+0.5)/(n-k+0.5))  (finite at k=0 and k=n), pooling plots across fields.
plots <- rbind(man23, man25, fill = TRUE)[is.finite(StPu) & !is.na(Genotype)]
el <- plots[, .(k = sum(StPu), n = .N), by = Genotype]
el[, prop := k / n][, StPu_mean := log((k + 0.5) / (n - k + 0.5))]
fwrite(el[, .(Genotype, k, n, prop, StPu_mean)], here("data/zeal/pheno_stpu_elogit.csv"))
log_info(
  "elogit StPu: %d genotypes | %d ever-pubescent (k>0) | elogit range [%.2f, %.2f]",
  nrow(el), sum(el$k > 0), min(el$StPu_mean), max(el$StPu_mean)
)

# TASSEL phenotype (gwas_nil lines, Family=taxon). PHENO env picks direct or elogit;
# elogit is the modeled StPu phenotype, so it is the default written to the TASSEL file.
PHENO <- Sys.getenv("PHENO", "elogit")
src <- if (PHENO == "elogit") el[, .(Genotype, StPu_mean)] else m[, .(Genotype, StPu_mean)]
ss <- fread(here("data/zeal/samplesheet_3way.csv"))
tass <- merge(ss[gwas_nil == TRUE, .(pedigree, taxon)], src[, .(pedigree = Genotype, y = StPu_mean)], by = "pedigree")[is.finite(y)]
ph_out <- here("data/zeal/tassel/pheno_stpu_all.txt")
writeLines(c("<Phenotype>", "taxa\tdata\tfactor", "Taxa\tStPu\tFamily"), ph_out)
fwrite(tass[, .(pedigree, round(y, 4), taxon)], ph_out, sep = "\t", append = TRUE, col.names = FALSE)
log_info("wrote %s (%s phenotype, %d gwas_nil lines)", ph_out, PHENO, nrow(tass))
