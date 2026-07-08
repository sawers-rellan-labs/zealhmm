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
# Significance cutoff — from the paper, NO permutation here (-nperm 0):
#   Chen 2019 state the LOD/P convention "P < 0.00001 (LOD = 5)", i.e.
#   LOD = -log10(P). We use LOD 5 == P = 1e-5 as the stepwise enter/exit
#   marginal-F cutoff, matching the -log10P >= 5 filter applied to the Fig 4A JLM
#   QTL display (owner decision: unify the JLM search and display at LOD 5). Note
#   this is stricter than Chen's own JLM permutation cutoff (~LOD 4.0 at P<0.05).
#
# Inputs : data/teonam/tassel/geno.hmp.txt        (native JLM pool, built by
#                                                  scripts/teonam_jlm_build.R)
#          data/teonam/tassel/pheno_stam_all.txt  (Taxa | STAM data | Family factor)
# Output : data/teonam/tassel/stam_jlm_native*.txt (leaves the old stam_jlm_*.txt
#          consensus-pool outputs intact for comparison)
# Run    : Rscript scripts/teonam_jlm_run.R
# =============================================================================
suppressMessages(library(data.table))
setwd("/Users/fvrodriguez/repos/zealhmm")

TASSEL <- "/Applications/TASSEL 5/run_pipeline.pl"
TDIR <- "data/teonam/tassel"
HMP <- file.path(TDIR, "geno.hmp.txt") # native 9,063-marker JLM pool
TRAIT <- toupper(Sys.getenv("TRAIT", "STAM"))
TTAG <- tolower(TRAIT) # phenotype; STAM default, e.g. DTA
PHENO <- file.path(TDIR, sprintf("pheno_%s_all.txt", TTAG)) # <TRAIT> (data) + Family (factor)
OUTBASE <- file.path(TDIR, sprintf("%s_jlm_native", TTAG))
ENTER <- "0.00001" # LOD 5 == P 1e-5 (Chen 2019 convention: LOD = -log10 P)
EXIT <- "0.00001"
MAXMK <- "100" # safety cap only; the P-value cutoff is the real gate (STAM ~ 5 QTL)
stopifnot(file.exists(HMP), file.exists(PHENO), file.exists(TASSEL))

# StepwiseOLS on genotype+phenotype (intersected taxa); Family fitted first,
# markers nested within Family, enter/exit by marginal-F P at the paper cutoff.
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
cat(cmd, "\n\n")
t0 <- Sys.time()
rc <- system(cmd)
cat(sprintf(
  "\nStepwiseOLS rc=%d, %.1f s\noutputs: %s\n", rc,
  as.numeric(Sys.time() - t0, units = "secs"),
  paste(list.files(TDIR, pattern = sprintf("^%s_jlm_native", TTAG)), collapse = ", ")
))
