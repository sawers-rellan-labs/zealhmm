#!/usr/bin/env Rscript
# =============================================================================
# TeoNAM JLM (joint linkage mapping) — TASSEL5 StepwiseOLSModelFitterPlugin.
# The committed, reproducible JLM runner. Previously the JLM stepwise fit was run
# ad-hoc (GUI/inline) and never saved, so it could not be reproduced when the JLM
# marker pool changed — this script closes that traceability hole.
#
# FAITHFUL to Chen 2019 (Methods; after Buckler 2009 / Tian 2011): the Family main
# effect is fitted first, then marker effects NESTED WITHIN family enter/leave the
# model by the marginal-F P-value; the plugin's backward step is the paper's
# drop-one refit. TASSEL5's StepwiseOLSModelFitter implements this NAM stepwise.
#
# Significance cutoff — PERMUTATION empirical alpha (Chen 2019 NAM JLM), via TASSEL's
#   native -nperm: the phenotype is permuted NPERM times, the null min-P per permutation
#   is collected (exported "Permuted_Pvalues"), and the enter limit = the alpha quantile
#   of that null (genome-wide FWER). We keep -enter/-exit permissive (1e-5) only as the
#   forward-search start; TASSEL applies the permuted enter limit. The empirical threshold
#   is stored to data/teonam/jlm_perm_thresholds.csv so the notebooks' Fig-4A lollipop
#   draws + filters at it instead of the fixed LOD 5. (Was: fixed LOD 5 == P 1e-5,
#   -nperm 0 — the "unify search & display at LOD 5" owner decision, now superseded to
#   match the ZEAL JLM + GWAS + R/qtl permutation approach.)
#
# Inputs : data/teonam/tassel/geno.hmp.txt        (native JLM pool, built by
#                                                  scripts/teonam_jlm_build.R)
#          data/teonam/tassel/pheno_<trait>_all.txt  (Taxa | <TRAIT> data | Family factor)
# Output : data/teonam/tassel/<trait>_jlm_native{1..5}.txt (native5 = Permuted_Pvalues) +
#          a row in data/teonam/jlm_perm_thresholds.csv
# Run    : TRAIT=STAM NPERM=1000 Rscript scripts/teonam_jlm_run.R
# =============================================================================
suppressMessages(library(data.table))
setwd("/Users/fvrodriguez/repos/zealhmm")

TASSEL <- "/Applications/TASSEL 5/run_pipeline.pl"
TDIR <- "data/teonam/tassel"
HMP <- file.path(TDIR, "geno.hmp.txt") # native 9,063-marker JLM pool
TRAIT <- toupper(Sys.getenv("TRAIT", "STAM"))
TTAG <- tolower(TRAIT) # phenotype; STAM default, e.g. DTA
NPERM <- as.integer(Sys.getenv("NPERM", "1000"))
PHENO <- file.path(TDIR, sprintf("pheno_%s_all.txt", TTAG)) # <TRAIT> (data) + Family (factor)
OUTBASE <- file.path(TDIR, sprintf("%s_jlm_native", TTAG))
CONSOLE <- file.path(TDIR, sprintf("%s_jlm_console.log", TTAG))
THRCSV <- "data/teonam/jlm_perm_thresholds.csv"
DEGEN_LOD <- as.numeric(Sys.getenv("DEGEN_LOD", "10"))
ENTER <- "0.00001" # forward-search start; permutation sets the real enter limit
EXIT <- "0.00001"
MAXMK <- "100" # safety cap only; the P-value cutoff is the real gate (STAM ~ 5 QTL)
stopifnot(file.exists(HMP), file.exists(PHENO), file.exists(TASSEL))

# StepwiseOLS on genotype+phenotype (intersected taxa); Family fitted first,
# markers nested within Family, enter/exit by marginal-F P; -nperm sets the empirical alpha.
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
cat(cmd, "\n\n")
t0 <- Sys.time()
rc <- system(cmd)
cat(sprintf(
  "\nStepwiseOLS rc=%d, %.1f s\noutputs: %s\n", rc,
  as.numeric(Sys.time() - t0, units = "secs"),
  paste(list.files(TDIR, pattern = sprintf("^%s_jlm_native[0-9]", TTAG)), collapse = ", ")
))

# --- empirical enter limit -> data/teonam/jlm_perm_thresholds.csv --------------
# Mirrors scripts/zeal_jlm_run.R. Degeneracy guard: if the perm null collapses
# (LOD > DEGEN_LOD, a binary/separable artifact) re-run -nperm 0 at the fixed enter
# and store a flagged fixed fallback.
upsert_threshold <- function(al, enter_p, note = "") {
  row <- data.table(
    trait = TRAIT, model = "jlm", geno = NA_character_, famcol = "Family",
    alpha = al, thr_neglog10p = -log10(enter_p), enter_p = enter_p, nperm = NPERM, note = note
  )
  cur <- if (file.exists(THRCSV)) fread(THRCSV) else data.table()
  if (nrow(cur)) cur <- cur[!(trait == TRAIT & model == "jlm" & alpha == al)]
  fwrite(rbind(cur, row, fill = TRUE), THRCSV)
}
perm_null_file <- function() {
  fs <- list.files(TDIR, pattern = sprintf("^%s_jlm_native[0-9]+\\.txt$", TTAG), full.names = TRUE)
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
    cat(sprintf("DEGENERATE perm null for %s (enter=%.2e, LOD %.1f > %g) — re-running -nperm 0 at fixed enter, storing flagged fallback\n", TRAIT, enter05, lod05, DEGEN_LOD))
    system(sub("-nperm [0-9]+", "-nperm 0", cmd))
    upsert_threshold(0.05, as.numeric(ENTER), note = sprintf("degenerate_perm(LOD%.1f)_fixed_fallback", lod05))
  } else {
    upsert_threshold(0.05, enter05)
    cat(sprintf("perm enter limit (alpha=0.05) = %.3e -> LOD %.2f\n", enter05, lod05))
    pf <- perm_null_file()
    if (!is.na(pf)) {
      nullp <- sort(as.numeric(fread(pf)[[1]]))
      enter10 <- nullp[floor(0.10 * NPERM) + 1]
      upsert_threshold(0.10, enter10)
      cat(sprintf("perm enter limit (alpha=0.10) = %.3e -> LOD %.2f\n", enter10, -log10(enter10)))
    }
  }
}
