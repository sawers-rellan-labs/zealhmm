#!/usr/bin/env Rscript
# =============================================================================
# ZEAL Guadalajara (GDL, subtropical) two-field BLUEs -> R/qtl phenotype files
# -----------------------------------------------------------------------------
# The 2024 GDL evaluation was grown at two fields, CUCBA and La Soledad, each a
# single-replicate augmented layout (39 incomplete blocks; only the 3 checks
# repeated per block). Treating the two FIELDS as two reps of ONE Guadalajara
# (subtropical) environment gives each test line 2 plots, so a two-field mixed
# model with genotype fixed yields a BLUE per line over the environment:
#     trait ~ genotype (fixed) + (1|field) + (1|field:Bloque)
# (checks anchor the field/block terms; the 2-level field variance is weakly
# identified and a singular fit there is expected and harmless for the BLUE).
#
# TRAIT SCOPE / EXCEPTIONS (verified from the data, not assumed):
#  - FLOWERING (FM/FF) is NOT done here: La Soledad FM/FF = CUCBA - 3 days exactly
#    (a constant 3-day offset from CUCBA, so no independent second rep; cause not yet
#    established, to be confirmed with the Guadalajara team). DTA/DTS stay CUCBA-only and
#    are owned by scripts/zeal_cucba_pheno.R; this script only *reuses* that file's
#    DTA/DTS columns to complete the wide table (no recomputation, no divergence).
#  - DISEASE (incidence/severity, 3 dates) is La Soledad ONLY (single field), so it
#    is fit block-random single-field: trait ~ genotype (fixed) + (1|Bloque).
#    NOTE: INC* and SEV* are ORDINAL SCORES (INC 0-4, SEV 0-5 in the raw data), not
#    0/1 proportions, so they are modeled AS-MEASURED (Gaussian) -- NOT empirical
#    logit (which would require proportions/counts). Dates 2-3 incidence are largely
#    at ceiling (median = max = 5); interpret those BLUEs with care.
#  - All other traits (PH, EH, EN, lodging, ear damage, husk cover, weights,
#    moisture) are independent between the two fields, so two-field BLUEs are valid.
#
# BLUEs are keyed by Origen and joined to the BC2S3 pedigree (= line_id = cross id)
# via the Genealogia crosswalk (PROVISIONAL, pending the definitive three-way sheet).
#
# In:  data/zeal/gdl/Libro 2024T CUCBA.xlsx        (sheet "CUCBA")
#      data/zeal/gdl/Libro 2024T La Soledad.xlsx   (sheet "La Soledad", + disease)
#      data/zeal/gdl/gdl_lineid_attempt.csv        (Origen -> line_id, provisional)
#      results/sim/zeal/rqtl/zeal_jal_flowering_blues.csv (optional; for DTA_JAL/
#                                                    DTS_JAL columns in the wide table)
# Out: results/sim/zeal/rqtl/zeal_jal_blues.csv    (origen, all JAL traits, line_id)
#      data/zeal/pheno_<trait>_jal_blue.csv        (Genotype, <TRAIT>_JAL_mean) per trait
# Next: TRAITS="PH_JAL,EH_JAL,EN_JAL,..." NPERM=1000 Rscript scripts/zeal_rqtl_scan.R
# Run:  Rscript scripts/zeal_gdl_pheno.R
# =============================================================================
suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
  library(here)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})
source(here::here("scripts/logging.R"))
emm_options(msg.interaction = FALSE)
GDL <- here("data/zeal/gdl")
RQ <- here("results/sim/zeal/rqtl")
dir.create(RQ, recursive = TRUE, showWarnings = FALSE)
# global fit counter + timer for a running ETA across all mixed-model fits
NFIT <- 10L + 6L # length(TWOFIELD) + length(DISEASE), for the ETA denominator
T0 <- Sys.time()
FI <- 0L
tick <- function(code, col, obs, geno, extra = "") {
  FI <<- FI + 1L
  el <- as.numeric(difftime(Sys.time(), T0, units = "mins"))
  log_info(
    ">>> %2d/%2d %-9s <- %-12s obs=%4d geno=%4d %s| elapsed %.1f min | avg %.2f | ETA ~%.1f min",
    FI, NFIT, code, col, obs, geno, extra, el, el / FI, (el / FI) * (NFIT - FI)
  )
}
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
rd <- function(f, s) {
  as.data.table(read_excel(file.path(GDL, f), s,
    .name_repair = "minimal", col_types = "text"
  ))
}

cu <- rd("Libro 2024T CUCBA.xlsx", "CUCBA")
ls <- rd("Libro 2024T La Soledad.xlsx", "La Soledad")
cu[, field := "CUCBA"]
ls[, field := "LaSoledad"]

# ---- two-field independent traits: <TRAIT>_JAL <- raw column ----------------
# _JAL = Guadalajara (subtropical macro-environment); these BLUEs pool over the two
# JAL sites (CUCBA + La Soledad). Base tokens match the Clayton _CLY traits so the
# two environments line up (PH_JAL vs PH_CLY, EH_JAL vs EH_CLY, EN_JAL vs EN_CLY).
TWOFIELD <- list(
  PH_JAL    = "AP", # altura de planta  = plant height (m)      <-> CLY PH
  EH_JAL    = "AM", # altura de mazorca = ear height (m)        <-> CLY EH
  EN_JAL    = "NM", # numero de mazorcas = ear number           <-> CLY EN
  AT_JAL    = "AT", # acame de tallo    = stalk lodging
  AR_JAL    = "AR", # acame de raiz     = root lodging
  MD_JAL    = "MD", # mazorcas danadas  = ear damage
  CM_JAL    = "CM", # cobertura mazorca = husk cover
  FW_JAL    = "PC", # peso de campo     = field/ear weight (fieldbook col PC)
  GW_JAL    = "PG", # peso de grano     = grain weight (fieldbook col PG)
  MOIST_JAL = "%H" # % humedad         = grain moisture
)
common <- c("Parcela", "Bloque", "Origen", "field")
stack <- rbind(
  cu[, c(common, unlist(TWOFIELD)), with = FALSE],
  ls[, c(common, unlist(TWOFIELD)), with = FALSE]
)

blue_twofield <- function(code, col) {
  dd <- data.table(
    genotype = as.factor(stack$Origen),
    field    = as.factor(stack$field),
    fblock   = as.factor(paste(stack$field, stack$Bloque, sep = ":")),
    y        = num(stack[[col]])
  )[!is.na(y)]
  m <- suppressWarnings(suppressMessages(lmer(
    y ~ genotype + (1 | field) + (1 | fblock),
    data = dd,
    control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
  )))
  em <- as.data.table(emmeans(m, ~genotype))
  tick(code, col, nrow(dd), uniqueN(dd$genotype))
  data.table(origen = as.character(em$genotype), trait = code, BLUE = em$emmean)
}
log_info("Two-field GDL BLUEs: trait ~ genotype (fixed) + (1|field) + (1|field:Bloque)")
blues <- rbindlist(lapply(names(TWOFIELD), function(k) blue_twofield(k, TWOFIELD[[k]])))

# ---- disease (La Soledad only): trait ~ genotype (fixed) + (1|Bloque) -------
DISEASE <- list(
  INC1_JAL = "INC 18-09-24", SEV1_JAL = "SEV1",
  INC2_JAL = "INC 26-09-24", SEV2_JAL = "SEV2",
  INC3_JAL = "INC 03-10-24", SEV3_JAL = "SEV3"
)
blue_singlefield <- function(code, col) {
  dd <- data.table(
    genotype = as.factor(ls$Origen), block = as.factor(ls$Bloque),
    y = num(ls[[col]])
  )[!is.na(y)]
  m <- suppressWarnings(suppressMessages(lmer(
    y ~ genotype + (1 | block),
    data = dd,
    control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
  )))
  em <- as.data.table(emmeans(m, ~genotype))
  tick(code, col, nrow(dd), uniqueN(dd$genotype), extra = "(LaSoledad only) ")
  data.table(origen = as.character(em$genotype), trait = code, BLUE = em$emmean)
}
log_info("Disease BLUEs (La Soledad only, as-measured ordinal): trait ~ genotype + (1|Bloque)")
dblues <- rbindlist(lapply(names(DISEASE), function(k) blue_singlefield(k, DISEASE[[k]])))
# across-date mean incidence / severity (simple mean of the per-date BLUEs)
disw <- dcast(dblues, origen ~ trait, value.var = "BLUE")
disw[, INC_JAL := rowMeans(.SD, na.rm = TRUE), .SDcols = c("INC1_JAL", "INC2_JAL", "INC3_JAL")]
disw[, SEV_JAL := rowMeans(.SD, na.rm = TRUE), .SDcols = c("SEV1_JAL", "SEV2_JAL", "SEV3_JAL")]

# ---- assemble wide table (origen x trait), + CUCBA flowering + line_id ------
wide <- dcast(blues, origen ~ trait, value.var = "BLUE")
wide <- merge(wide, disw, by = "origen", all = TRUE)

flw <- file.path(RQ, "zeal_jal_flowering_blues.csv")
if (file.exists(flw)) {
  fl <- fread(flw)[, .(origen, DTA_JAL, DTS_JAL)]
  wide <- merge(wide, fl, by = "origen", all.x = TRUE)
  cat("reused JAL (CUCBA) flowering DTA_JAL/DTS_JAL for the wide table\n")
} else {
  cat("NOTE: zeal_jal_flowering_blues.csv not found; run zeal_cucba_pheno.R for DTA_JAL/DTS_JAL\n")
}

xw <- fread(file.path(GDL, "gdl_lineid_attempt.csv"))
xw <- unique(xw[, .(origen, line_id = line_id_from_genealogia)])
wide <- merge(wide, xw, by = "origen", all.x = TRUE)
fwrite(wide, file.path(RQ, "zeal_jal_blues.csv"))
cat(
  "wrote zeal_jal_blues.csv | rows:", nrow(wide),
  "| line_id attached:", sum(!is.na(wide$line_id) & wide$line_id != ""), "\n"
)

# ---- per-trait R/qtl pheno files, keyed by pedigree (line_id), deduped ------
ph <- wide[!is.na(line_id) & line_id != ""]
ALL <- c(names(TWOFIELD), names(DISEASE), "INC_JAL", "SEV_JAL")
for (code in ALL) {
  agg <- ph[!is.na(get(code)), .(m = mean(get(code), na.rm = TRUE)), by = .(Genotype = line_id)]
  setnames(agg, "m", paste0(code, "_mean"))
  fwrite(agg, here(sprintf("data/zeal/pheno_%s_blue.csv", tolower(code))))
}
cat(
  "wrote", length(ALL), "per-trait pheno files (data/zeal/pheno_*_jal_blue.csv) |",
  "lines with line_id:", uniqueN(ph$line_id), "\n"
)
cat("traits:", paste(ALL, collapse = ","), "\n")
