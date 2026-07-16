#!/usr/bin/env Rscript
# ZEAL Phase 4 — JLM via TASSEL5 StepwiseOLS (analog of teonam_jlm_run.R).
# Family (=taxon) fitted first; marker effects NESTED within Family enter/leave the
# model on marginal-F P. Model-independent QTL analysis (Panel A). TRAIT via env
# (default DTA).
#
# Significance = PERMUTATION empirical alpha (Chen 2019 NAM JLM), via TASSEL's native
# -nperm: the phenotype is permuted NPERM times, the null min-P per permutation is
# collected (exported "Permuted_Pvalues"), and the enter limit = the alpha quantile of
# that null (FWER). We keep -enter/-exit permissive (1e-4) only as the forward-search
# start; TASSEL then applies the permuted enter limit. The empirical threshold is stored
# to data/zeal/gwas_perm_thresholds.csv so the notebooks (plot_lollipop) draw + filter at
# the per-trait permutation cutoff instead of a fixed global LOD.
#
# Inputs : data/zeal/tassel/geno_jlm.hmp.txt (zeal_jlm_build.R), data/zeal/tassel/pheno_<trait>_all.txt
# Output : results/sim/zeal/tassel/<trait>_jlm_native{1..5}.txt  (native1 = ANOVA_stepwise QTL;
#          native5 = Permuted_Pvalues null) + a row in data/zeal/gwas_perm_thresholds.csv
# Run    : TRAIT=DTA NPERM=1000 Rscript scripts/zeal_jlm_run.R
suppressMessages({
  library(here)
  library(data.table)
})
source(here("scripts/logging.R"))
TRAIT <- toupper(Sys.getenv("TRAIT", "DTA"))
TTAG <- tolower(TRAIT)
NPERM <- as.integer(Sys.getenv("NPERM", "1000"))
TDIR <- here("data/zeal/tassel") # JLM inputs (built HapMap + phenotype)
RDIR <- here("results/sim/zeal/tassel") # JLM outputs (results)
dir.create(RDIR, showWarnings = FALSE, recursive = TRUE)
TASSEL <- "/Applications/TASSEL 5/run_pipeline.pl"
HMP <- file.path(TDIR, "geno_jlm.hmp.txt")
PHENO <- file.path(TDIR, sprintf("pheno_%s_all.txt", TTAG))
OUTBASE <- file.path(RDIR, sprintf("%s_jlm_native", TTAG))
CONSOLE <- file.path(RDIR, sprintf("%s_jlm_console.log", TTAG))
THRCSV <- here("data/zeal/gwas_perm_thresholds.csv") # shared perm-threshold table
ENTER <- "0.0001" # forward-search start; permutation sets the real enter limit
EXIT <- "0.0001"
MAXMK <- "100"
stopifnot(file.exists(HMP), file.exists(PHENO), file.exists(TASSEL))

cmd <- sprintf(
  paste(
    '"%s" -Xmx10g',
    "-fork1 -h %s -fork2 -r %s -combine3 -input1 -input2 -intersect",
    "-StepwiseOLSModelFitterPlugin -modelType pvalue -enter %s -exit %s -nperm %d",
    "-nestMarkers true -nestFactor Family -maxMarkers %s -endPlugin",
    "-export %s -runfork1 -runfork2 > %s 2>&1"
  ),
  TASSEL, shQuote(HMP), shQuote(PHENO), ENTER, EXIT, NPERM, MAXMK,
  shQuote(OUTBASE), shQuote(CONSOLE)
)
log_info("TRAIT=%s NPERM=%d | %s", TRAIT, NPERM, cmd)
t0 <- Sys.time()
rc <- system(cmd)
secs <- as.numeric(Sys.time() - t0, units = "secs")
log_info(
  "StepwiseOLS rc=%d, %.1f s | outputs: %s", rc, secs,
  paste(list.files(RDIR, pattern = sprintf("^%s_jlm_native[0-9]", TTAG)), collapse = ", ")
)

# --- empirical enter limit -------------------------------------------------
# alpha = 0.05: TASSEL prints "the Enter limit is <p>" (its default alpha). alpha = 0.10:
# from the exported Permuted_Pvalues null with TASSEL's index convention
# (idx = floor(alpha*N), 0-based -> the (idx+1)-th smallest permuted min-P).
# DEGENERACY GUARD: near-binary/elogit traits (StPu) hit perfect-separation under
# permutation, so the null min-P collapses to ~1e-28 (LOD >> 10) — a numerical artifact,
# not a real cutoff. Above DEGEN_LOD we discard the permutation, re-run -nperm 0 at the
# fixed enter (restoring a sensible native1), and store a FLAGGED fixed fallback.
DEGEN_LOD <- as.numeric(Sys.getenv("DEGEN_LOD", "10"))
upsert_threshold <- function(al, enter_p, note = "") {
  # `al`/`TRAIT` deliberately avoid the column names (trait, alpha) so data.table
  # resolves them from scope inside the i-filter (the `..` prefix only works in j).
  row <- data.table(
    trait = TRAIT, model = "jlm", geno = NA_character_, famcol = "Family",
    alpha = al, thr_neglog10p = -log10(enter_p), enter_p = enter_p, nperm = NPERM, note = note
  )
  cur <- if (file.exists(THRCSV)) fread(THRCSV) else data.table()
  if (nrow(cur)) cur <- cur[!(trait == TRAIT & model == "jlm" & alpha == al)]
  fwrite(rbind(cur, row, fill = TRUE), THRCSV)
}
# Permuted_Pvalues export index varies (4 or 5 files depending on the CI-scan tables);
# find it by its "P-value" header rather than assuming a fixed suffix.
perm_null_file <- function() {
  fs <- list.files(RDIR, pattern = sprintf("^%s_jlm_native[0-9]+\\.txt$", TTAG), full.names = TRUE)
  for (f in fs) {
    h <- readLines(f, n = 1, warn = FALSE)
    if (length(h) && grepl("^P-value", h)) {
      return(f)
    }
  }
  NA_character_
}

con <- readLines(CONSOLE, warn = FALSE)
lim_line <- grep("the Enter limit is", con, value = TRUE)
if (NPERM > 0 && length(lim_line)) {
  enter05 <- as.numeric(sub(".*the Enter limit is\\s*", "", lim_line[1]))
  lod05 <- -log10(enter05)
  if (is.finite(lod05) && lod05 > DEGEN_LOD) {
    log_info(
      "DEGENERATE perm null for %s (enter=%.2e, LOD %.1f > %g) — binary/separable trait; re-running -nperm 0 at fixed enter=%s, storing flagged fixed fallback",
      TRAIT, enter05, lod05, DEGEN_LOD, ENTER
    )
    system(sub("-nperm [0-9]+", "-nperm 0", cmd)) # regenerate native1..4 at the fixed enter
    upsert_threshold(0.05, as.numeric(ENTER), note = sprintf("degenerate_perm(LOD%.1f)_fixed_fallback", lod05))
  } else {
    upsert_threshold(0.05, enter05)
    log_info("perm enter limit (alpha=0.05) = %.3e -> LOD %.2f", enter05, lod05)
    pf <- perm_null_file()
    if (!is.na(pf)) {
      nullp <- sort(as.numeric(fread(pf)[[1]]))
      enter10 <- nullp[floor(0.10 * NPERM) + 1]
      upsert_threshold(0.10, enter10)
      log_info("perm enter limit (alpha=0.10) = %.3e -> LOD %.2f", enter10, -log10(enter10))
    } else {
      log_info("no Permuted_Pvalues export found; alpha=0.10 not stored")
    }
  }
} else {
  log_info("NPERM=0 or no permuted enter limit found; no threshold stored")
}
