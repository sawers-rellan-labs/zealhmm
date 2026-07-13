#!/usr/bin/env Rscript
# ZEAL/BZea Phase 2b — spatially-corrected phenotype BLUEs (BZeaPheno recipe).
#
# Per field, per contiguous map grid: SpATS 2-D P-spline surface PSANOVA(Range,Row) with
# genotype FIXED (-> BLUEs) and B73 as check (statgenHTP engine="SpATS", what="fixed").
# Genotype BLUE per field = mean over that field's grids; overall = mean across fields
# (per-field correction, union mean; GxE deferred). Traits: DTA, DTS, PH, StPi, StPu.
#
# Inputs : data/zeal/{cly25_b5,cly23_d4}_fieldmap.csv (plot_id, block, range, col)
#          data/zeal/CLY25-Fieldbook.xlsx :: B5_BZea_eval          (new-form pedigree)
#          data/zeal/CLY23_D4_FieldBook.xlsx :: UPDATED_... + GENOTYPE-CONVERSION (old->new)
#          data/zeal/samplesheet_3way.csv    (pedigree -> taxon, gwas_nil)
# Outputs: data/zeal/pheno_<trait>_blue.csv  (genotype, per-field + mean BLUE)
#          data/zeal/pheno_blues_all.csv     (wide, all traits)
#          data/zeal/tassel/pheno_dta_all.txt (TASSEL: Taxa | DTA | Family=taxon, gwas_nil lines)

suppressMessages({
  library(here)
  library(data.table)
  library(readxl)
  library(SpATS)
})
source(here("scripts/logging.R"))

IQR_MULT <- as.numeric(Sys.getenv("IQR_MULT", "3")) # per field x taxa upper Tukey fence on raw plots (pre-BLUE) — the ONLY outlier step

TRAITS <- c("DTA", "DTS", "PH", "EH", "SPAD", "EN", "Prolif", "LAE", "NBR", "StPi", "StPu")
CONT_TRAITS <- c("DTA", "DTS", "PH", "EH", "SPAD", "LAE") # genuinely continuous; count traits (EN/Prolif/NBR) and binary (StPi/StPu) excluded
EXCEL_1970 <- 25569L # Excel(1900) serial for 1970-01-01
plant_serial <- function(d) as.integer(as.Date(d)) + EXCEL_1970
canon_ped <- function(x) sub("\\.B$", "", x)

# ---- per-field standardized manifests: Plot_id, Genotype, Rep, block, range, col, traits ----
manifest_cly25 <- function() {
  fm <- fread(here("data/zeal/cly25_b5_fieldmap.csv"))
  ph <- as.data.table(read_excel(here("data/zeal/CLY25-Fieldbook.xlsx"), sheet = "B5_BZea_eval"))
  setnames(ph, 1, "plot_id")
  ph[, plot_id := suppressWarnings(as.integer(plot_id))]
  p0 <- plant_serial("2025-04-03")
  ph[, `:=`(
    DTA = as.numeric(DOA) - p0, DTS = as.numeric(DOS) - p0,
    PH = as.numeric(PH), EH = as.numeric(EH), StPi = as.numeric(StPi), StPu = as.numeric(StPu),
    SPAD = as.numeric(SPAD), EN = as.numeric(EN), Prolif = as.numeric(Prolif),
    LAE = as.numeric(LAE), NBR = as.numeric(NBR),
    Genotype = canon_ped(Description),
    Rep = fifelse(grepl("Rep2", `Who/What`), 2L, 1L)
  )]
  merge(fm, ph[, .(plot_id, Genotype, Rep, DTA, DTS, PH, EH, StPi, StPu, SPAD, EN, Prolif, LAE, NBR)], by = "plot_id")
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
  ph[, `:=`(
    DTA = as.numeric(DTA), DTS = as.numeric(DTS), PH = as.numeric(PH), EH = as.numeric(EH),
    StPi = as.numeric(StPi), StPu = as.numeric(StPu), Rep = as.integer(Rep),
    EN = as.numeric(EN), Prolif = as.numeric(Prolif), NBR = as.numeric(NBR)
  )]
  merge(fm, ph[, .(plot_id, Genotype, Rep, DTA, DTS, PH, EH, StPi, StPu, EN, Prolif, NBR)], by = "plot_id")
}

# ---- fit one grid, one trait -> genotype BLUEs (SpATS, genotype FIXED) -------
# SpATS fit (single pass, genotype fixed). Outlier removal is done ONCE upstream as the
# per field x taxa IQR/Tukey fence on the RAW plots (see below) — there is no second,
# residual-based removal inside the fit.
fit_grid_trait <- function(dat, trait, tag = "") {
  d0 <- dat[is.finite(get(trait)) & !is.na(Genotype)]
  if (uniqueN(d0$Genotype) < 5 || nrow(d0) < 20) {
    return(NULL)
  }
  d <- data.frame(
    Genotype = factor(as.character(d0$Genotype)), Rep = factor(d0$Rep),
    Range = as.numeric(d0$range), Row = as.numeric(d0$col), y = as.numeric(d0[[trait]])
  )
  d$Rf <- factor(d$Range)
  d$Cf <- factor(d$Row)
  fixed <- if (nlevels(d$Rep) > 1) ~Rep else NULL
  fit_one <- function(dat) {
    tryCatch(SpATS(
      response = "y", genotype = "Genotype", genotype.as.random = FALSE,
      spatial = ~ PSANOVA(Row, Range), fixed = fixed, random = ~ Rf + Cf,
      data = dat, control = list(monitoring = 0, maxit = 50)
    ), error = function(e) NULL)
  }

  m <- fit_one(d)
  if (is.null(m)) {
    return(NULL)
  }
  pr <- tryCatch(predict(m, which = "Genotype"), error = function(e) NULL)
  if (is.null(pr)) {
    return(NULL)
  }
  pr <- as.data.table(pr)
  vc <- grep("predicted", names(pr), ignore.case = TRUE, value = TRUE)[1]
  pr[!is.na(get(vc)), .(Genotype = as.character(Genotype), blue = get(vc))]
}

# ---- run all fields x grids x traits ----------------------------------------
fields <- list(cly23 = manifest_cly23(), cly25 = manifest_cly25())
# Per field x taxa upper Tukey fence on RAW plots (pre-BLUE): set continuous-trait
# plot values > Q3 + IQR_MULT*IQR (within field x taxon) to NA before SpATS. This is
# the per-trait cleaning that feeds the BLUEs; the multivariate outlier DETECTION
# (a flag list, no removal) is separate in agent/mv_outlier_detect.R.
taxon_map <- unique(fread(here("data/zeal/samplesheet_3way.csv"))[, .(Genotype = pedigree, taxon)], by = "Genotype")
per_field <- list()
for (fld in names(fields)) {
  man <- fields[[fld]]
  man[taxon_map, taxon := i.taxon, on = "Genotype"]
  for (tr in intersect(CONT_TRAITS, names(man))) {
    n0 <- man[is.finite(get(tr)), .N]
    man[, (tr) := {
      v <- .SD[[1]]
      ok <- is.finite(v)
      if (sum(ok) >= 4 && IQR(v[ok]) > 0) {
        fen <- quantile(v[ok], 0.75) + IQR_MULT * IQR(v[ok])
        v[ok & v > fen] <- NA
      }
      v
    }, by = taxon, .SDcols = tr]
    nf <- n0 - man[is.finite(get(tr)), .N]
    if (nf) log_info("  %s %-4s: per-taxa %gxIQR upper fence -> %d plots set NA", fld, tr, IQR_MULT, nf)
  }
  log_info(
    "=== %s: %d plots, %d grids, %d genotypes ===", fld, nrow(man),
    uniqueN(man$block), uniqueN(man$Genotype)
  )
  ftrait <- list()
  for (tr in TRAITS) {
    if (!tr %in% names(man)) next # trait not scored in this field (e.g. SPAD/LAE are CLY25-only)
    grids <- lapply(unique(man$block), function(b) fit_grid_trait(man[block == b], tr, tag = fld))
    grids <- rbindlist(Filter(Negate(is.null), grids))
    if (!nrow(grids)) {
      log_warn("  %-4s: no fit", tr)
      next
    }
    bl <- grids[, .(blue = mean(blue)), by = Genotype] # mean over grids within field
    setnames(bl, "blue", tr)
    ftrait[[tr]] <- bl
    log_info(
      "  %-4s: %d genotype BLUEs (mean %.2f, sd %.2f)", tr, nrow(bl),
      mean(bl[[tr]], na.rm = TRUE), sd(bl[[tr]], na.rm = TRUE)
    )
  }
  per_field[[fld]] <- Reduce(function(a, b) merge(a, b, by = "Genotype", all = TRUE), ftrait)
}

# ---- combine fields: per-trait per-field cols + union mean ------------------
alltraits <- data.table(Genotype = character())
for (tr in TRAITS) {
  cols <- lapply(names(per_field), function(fld) {
    x <- per_field[[fld]]
    if (is.null(x) || !tr %in% names(x)) {
      return(NULL)
    }
    setnames(x[, .(Genotype, get(tr))], c("Genotype", paste0(tr, "_", fld)))
  })
  cols <- Filter(Negate(is.null), cols)
  if (!length(cols)) next
  m <- Reduce(function(a, b) merge(a, b, by = "Genotype", all = TRUE), cols)
  fcols <- setdiff(names(m), "Genotype")
  m[, (paste0(tr, "_mean")) := rowMeans(.SD, na.rm = TRUE), .SDcols = fcols]
  fwrite(m, here(sprintf("data/zeal/pheno_%s_blue.csv", tolower(tr))))
  alltraits <- merge(alltraits, m[, .(Genotype, get(paste0(tr, "_mean")))], by = "Genotype", all = TRUE)
  setnames(alltraits, "V2", paste0(tr, "_mean"))
}
fwrite(alltraits, here("data/zeal/pheno_blues_all.csv"))
log_info("wrote per-trait BLUEs + pheno_blues_all.csv (%d genotypes)", nrow(alltraits))

# ---- TASSEL files (Taxa | <trait> | Family=taxon), gwas_nil lines only --------
# One pheno_<trait>_all.txt per trait — JLM (zeal_jlm_run.R) reads pheno_<TTAG>_all.txt.
ss <- fread(here("data/zeal/samplesheet_3way.csv"))
dir.create(here("data/zeal/tassel"), showWarnings = FALSE)
for (tr in TRAITS) {
  bl <- fread(here(sprintf("data/zeal/pheno_%s_blue.csv", tolower(tr))))[
    , .(pedigree = Genotype, y = get(paste0(tr, "_mean")))
  ]
  tass <- merge(ss[gwas_nil == TRUE, .(pedigree, taxon)], bl, by = "pedigree")[is.finite(y)]
  log_info(
    "TASSEL %s: %d gwas_nil lines with a BLUE (of %d panel lines)",
    tr, nrow(tass), ss[gwas_nil == TRUE, .N]
  )
  ph_out <- here(sprintf("data/zeal/tassel/pheno_%s_all.txt", tolower(tr)))
  writeLines(c("<Phenotype>", "taxa\tdata\tfactor", sprintf("Taxa\t%s\tFamily", tr)), ph_out)
  fwrite(tass[, .(pedigree, round(y, 4), taxon)], ph_out, sep = "\t", append = TRUE, col.names = FALSE)
  log_info("wrote %s", ph_out)
}
