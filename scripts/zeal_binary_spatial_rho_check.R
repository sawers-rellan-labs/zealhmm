#!/usr/bin/env Rscript
# "Does spatial correction matter?" check for the BINARY traits (StPi, StPu) that SpATS
# can't fit — the justification for using the empirical logit with NO spatial correction
# (cited in the StPi/StPu notebook phenotype callouts).
# Per field, fit a plot-level binomial GAM WITH vs
# WITHOUT a 2-D spatial smooth (genotype as a random effect in both), and correlate the
# per-genotype logit BLUPs. High rho => spatial correction does not move the genotype
# ranking => the empirical-logit phenotype (no spatial correction) is defensible.
# Also reports rho(spatial BLUP vs the empirical logit currently used in the GWAS).
suppressMessages({
  library(here)
  library(data.table)
  library(readxl)
  library(mgcv)
})
source(here("scripts/logging.R"))
canon_ped <- function(x) sub("\\.B$", "", x)

man25 <- function() {
  fm <- fread(here("data/zeal/cly25_b5_fieldmap.csv"))
  ph <- as.data.table(read_excel(here("data/zeal/CLY25-Fieldbook.xlsx"), sheet = "B5_BZea_eval"))
  setnames(ph, 1, "plot_id")
  ph[, plot_id := suppressWarnings(as.integer(plot_id))]
  ph[, `:=`(StPi = as.numeric(StPi), StPu = as.numeric(StPu), Genotype = canon_ped(Description))]
  merge(fm, ph[, .(plot_id, Genotype, StPi, StPu)], by = "plot_id")
}
man23 <- function() {
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
  ph[, `:=`(StPi = as.numeric(StPi), StPu = as.numeric(StPu))]
  merge(fm, ph[, .(plot_id, Genotype, StPi, StPu)], by = "plot_id")
}

fields <- list(cly23 = man23(), cly25 = man25())
NTH <- max(1L, parallel::detectCores() - 2L)

rho_check <- function(dat, trait) {
  d <- dat[is.finite(get(trait)) & !is.na(Genotype)]
  d <- data.table(
    y = as.integer(d[[trait]]), Genotype = factor(d$Genotype),
    Range = as.numeric(d$range), Row = as.numeric(d$col)
  )
  if (sum(d$y) < 5) {
    return(NULL)
  }
  m0 <- tryCatch(bam(y ~ s(Genotype, bs = "re"), family = binomial, data = d, discrete = TRUE, nthreads = NTH), error = function(e) NULL)
  m1 <- tryCatch(bam(y ~ s(Genotype, bs = "re") + te(Range, Row), family = binomial, data = d, discrete = TRUE, nthreads = NTH), error = function(e) NULL)
  if (is.null(m0) || is.null(m1)) {
    return(NULL)
  }
  glev <- levels(d$Genotype)
  re <- function(m) coef(m)[paste0("s(Genotype).", seq_along(glev))]
  b0 <- re(m0)
  b1 <- re(m1)
  el <- d[, .(k = sum(y), n = .N), by = Genotype][, e := log((k + 0.5) / (n - k + 0.5))]
  ev <- el$e[match(glev, el$Genotype)]
  list(
    trait = trait, n = nrow(d), npos = sum(d$y), ngeno = length(glev),
    rho_spatial_vs_nospatial = cor(b0, b1),
    rho_spatial_vs_elogit = cor(b1, ev),
    spatial_edf = sum(m1$edf[grep("Range|Row", names(m1$edf))])
  )
}

for (tr in c("StPi", "StPu")) {
  for (fld in names(fields)) {
    r <- rho_check(fields[[fld]], tr)
    if (is.null(r)) {
      log_warn("%s %s: skipped (too few positives)", fld, tr)
      next
    }
    log_info(
      "%s %-4s | plots=%d pos=%d geno=%d | rho(spatial vs no-spatial)=%.4f | rho(spatial vs elogit)=%.4f | spatial edf=%.1f",
      fld, tr, r$n, r$npos, r$ngeno, r$rho_spatial_vs_nospatial, r$rho_spatial_vs_elogit, r$spatial_edf
    )
  }
}
