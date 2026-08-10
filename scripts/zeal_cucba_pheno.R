#!/usr/bin/env Rscript
# =============================================================================
# ZEAL CUCBA (short-day) flowering BLUEs -> R/qtl phenotype files
# -----------------------------------------------------------------------------
# Guadalajara CUCBA (subtropical, short-day) is the only real short-day flowering
# environment in the 2024 GDL evaluation: the La Soledad flowering columns are a
# derived copy of CUCBA (FM = CUCBA-3 exactly), so they are excluded here. CUCBA
# is an unreplicated augmented design (each NIL in one incomplete block; only the
# 3 checks replicated), so BLUEs come from a mixed model with block random:
#     trait ~ genotype (fixed) + (1 | Bloque)      # FM->DTA, FF->DTS, ASI=FF-FM
# BLUEs are keyed by Origen (BZeaPV23-...) and joined to the BC2S3 pedigree via
# the Genealogia crosswalk (PROVISIONAL, pending the definitive three-way sheet).
#
# In:  data/zeal/gdl/Libro 2024T CUCBA.xlsx        (raw, staged; see DATA.md)
#      data/zeal/gdl/gdl_lineid_attempt.csv        (Origen -> line_id, from the
#                                                    agent GDL harmonization)
# Out: results/sim/zeal/rqtl/zeal_cucba_flowering_blues.csv  (origen,DTA,DTS,ASI,line_id)
#      data/zeal/pheno_cucbadta_blue.csv  (Genotype, CUCBADTA_mean)
#      data/zeal/pheno_cucbadts_blue.csv  (Genotype, CUCBADTS_mean)
# Next: TRAITS="CUCBADTA,CUCBADTS" NPERM=1000 Rscript scripts/zeal_rqtl_scan.R
# Run:  Rscript scripts/zeal_cucba_pheno.R
# =============================================================================
suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
  library(here)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})
GDL <- here("data/zeal/gdl")
RQ <- here("results/sim/zeal/rqtl")
dir.create(RQ, recursive = TRUE, showWarnings = FALSE)
num <- function(x) suppressWarnings(as.numeric(as.character(x)))

d <- as.data.table(read_excel(file.path(GDL, "Libro 2024T CUCBA.xlsx"), "CUCBA",
  .name_repair = "minimal", col_types = "text"
))
d[, ASI := num(FF) - num(FM)]

blue_one <- function(trait, col) {
  dd <- data.table(genotype = as.factor(d$Origen), block = as.factor(d$Bloque), y = num(d[[col]]))[!is.na(y)]
  m <- suppressMessages(lmer(y ~ genotype + (1 | block),
    data = dd,
    control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
  ))
  em <- as.data.table(emmeans(m, ~genotype))
  cat(sprintf("  %-4s obs=%d geno=%d\n", trait, nrow(dd), uniqueN(dd$genotype)))
  data.table(origen = as.character(em$genotype), trait = trait, BLUE = em$emmean)
}
cat("CUCBA BLUEs: trait ~ genotype (fixed) + (1|Bloque)\n")
blues <- rbindlist(list(blue_one("DTA", "FM"), blue_one("DTS", "FF"), blue_one("ASI", "ASI")))
wide <- dcast(blues, origen ~ trait, value.var = "BLUE")

# attach line_id via the Genealogia crosswalk (provisional)
xw <- fread(file.path(GDL, "gdl_lineid_attempt.csv"))
xw <- unique(xw[, .(origen, line_id = line_id_from_genealogia)])
wide <- merge(wide, xw, by = "origen", all.x = TRUE)
fwrite(wide, file.path(RQ, "zeal_cucba_flowering_blues.csv"))
cat(
  "wrote zeal_cucba_flowering_blues.csv | line_id attached:",
  sum(!is.na(wide$line_id) & wide$line_id != ""), "/", nrow(wide), "\n"
)

# R/qtl pheno files, keyed by pedigree (= line_id = cross id), deduped by mean
ph <- wide[!is.na(line_id) & line_id != ""]
agg <- ph[, .(
  CUCBADTA_mean = mean(DTA, na.rm = TRUE),
  CUCBADTS_mean = mean(DTS, na.rm = TRUE)
), by = .(Genotype = line_id)]
fwrite(agg[, .(Genotype, CUCBADTA_mean)], here("data/zeal/pheno_cucbadta_blue.csv"))
fwrite(agg[, .(Genotype, CUCBADTS_mean)], here("data/zeal/pheno_cucbadts_blue.csv"))
cat("wrote pheno_cucbadta_blue.csv & pheno_cucbadts_blue.csv | lines:", nrow(agg), "\n")
