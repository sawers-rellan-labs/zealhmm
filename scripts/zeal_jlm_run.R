#!/usr/bin/env Rscript
# ZEAL/BZea Phase 4 — JLM via TASSEL5 StepwiseOLS (analog of teonam_jlm_run.R).
# Family (=taxon) fitted first; marker effects NESTED within Family enter/leave by
# marginal-F P at Chen's cutoff (LOD 5 == P 1e-5). Model-independent QTL analysis
# (Panel A). TRAIT via env (default DTA).
# Inputs : data/zeal/tassel/geno_jlm.hmp.txt (zeal_jlm_build.R), data/zeal/tassel/pheno_<trait>_all.txt
# Output : results/sim/zeal/tassel/<trait>_jlm_native{1..}.txt  (JLM QTL = analysis result)
suppressMessages(library(here))
source(here("scripts/logging.R"))
TRAIT <- toupper(Sys.getenv("TRAIT", "DTA"))
TTAG <- tolower(TRAIT)
TDIR <- here("data/zeal/tassel") # JLM inputs (built HapMap + phenotype)
RDIR <- here("results/sim/zeal/tassel") # JLM outputs (results)
dir.create(RDIR, showWarnings = FALSE, recursive = TRUE)
TASSEL <- "/Applications/TASSEL 5/run_pipeline.pl"
HMP <- file.path(TDIR, "geno_jlm.hmp.txt")
PHENO <- file.path(TDIR, sprintf("pheno_%s_all.txt", TTAG))
OUTBASE <- file.path(RDIR, sprintf("%s_jlm_native", TTAG))
ENTER <- "0.0001"
EXIT <- "0.0001"
MAXMK <- "100" # LOD 4 == P 1e-4 (Chen/JLM convention)
stopifnot(file.exists(HMP), file.exists(PHENO), file.exists(TASSEL))

cmd <- sprintf(
  paste(
    '"%s" -Xmx10g',
    "-fork1 -h %s -fork2 -r %s -combine3 -input1 -input2 -intersect",
    "-StepwiseOLSModelFitterPlugin -modelType pvalue -enter %s -exit %s -nperm 0",
    "-nestMarkers true -nestFactor Family -maxMarkers %s -endPlugin",
    "-export %s -runfork1 -runfork2"
  ),
  TASSEL, shQuote(HMP), shQuote(PHENO), ENTER, EXIT, MAXMK, shQuote(OUTBASE)
)
log_info("%s", cmd)
t0 <- Sys.time()
rc <- system(cmd)
log_info(
  "StepwiseOLS rc=%d, %.1f s | outputs: %s", rc, as.numeric(Sys.time() - t0, units = "secs"),
  paste(list.files(RDIR, pattern = sprintf("^%s_jlm_native", TTAG)), collapse = ", ")
)
