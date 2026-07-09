#!/usr/bin/env Rscript
# ZEAL/BZea — per-site GENOTYPE (_gt) object from call_gt on the real SNP50K counts.
# The GENOTYPE LAYER (input/evidence), parallel to the ancestry mosaics (_mosaic) but
# NOT an HMM: each (marker, line) is decided independently from its own (n_ref, n_alt)
# — no linkage. Replaces the retired ad-hoc "persnp" (round(2*alt/cov)) with the
# principled call_gt caller. See TERMINOLOGY.md.
#
# METHOD env (default ml):
#   ml       -> call_gt(prior = "flat")            maximum-likelihood (flat prior; het-blind; ~ old persnp)
#   gphwe    -> call_gt(prior = "hwe", af = NULL)  MAP, HWE (self-estimated AF; het-excess)
#   gpdesign -> call_gt(prior = design_prior("BC2S3"))  MAP, Mendelian design prior
#
# Output: data/zeal/zeal_<method>_gt.rds  list(markers, state[marker x line, 0/1/2], lines)
suppressMessages({
  library(here)
  library(data.table)
})
devtools::load_all("/Users/fvrodriguez/repos/nilhmm", quiet = TRUE)
source(here("scripts/logging.R"))
METHOD <- Sys.getenv("METHOD", "ml")
ERR <- 0.01

D <- readRDS(here("data/zeal/zeal_snp50k_dosage.rds"))
mk <- copy(D$markers)[, .(marker, chr = as.integer(chr), pos = as.integer(pos))]
ss <- fread(here("data/zeal/samplesheet_3way.csv"))[gwas_nil == TRUE & !is.na(skim_id)]
panel <- ss[skim_id %in% colnames(D$n_ref), unique(skim_id)]
ped <- ss[match(panel, skim_id), pedigree]
nref <- D$n_ref[mk$marker, panel]
nalt <- D$n_alt[mk$marker, panel]

prior <- switch(METHOD,
  ml = "flat",
  gphwe = "hwe",
  gpdesign = design_prior("BC2S3"),
  stop("unknown METHOD '", METHOD, "': expected ml | gphwe | gpdesign")
)
log_info(
  "method=%s (prior=%s) | %d markers x %d lines",
  METHOD, if (is.character(prior)) prior else "BC2S3-vector", nrow(mk), length(panel)
)

# per-site genotype call: 0/1/2 hardcall (argmax; MAP under a non-flat prior), NA at 0 depth
state <- call_gt(nref, nalt, prior = prior, error = ERR, return = "call")
rownames(state) <- mk$marker
colnames(state) <- ped

comp <- prop.table(table(factor(state, levels = 0:2)))
log_info(
  "%s_gt: %d x %d | NA=%.1f%% | B73=%.1f%% het=%.1f%% teo=%.1f%% | presence=%.3f (BC2S3 ~0.11)",
  METHOD, nrow(state), ncol(state), 100 * mean(is.na(state)),
  100 * comp["0"], 100 * comp["1"], 100 * comp["2"], sum(comp[c("1", "2")])
)
saveRDS(
  list(markers = mk, state = state, lines = data.table(skim_id = panel, pedigree = ped)),
  here(sprintf("data/zeal/zeal_%s_gt.rds", METHOD))
)
log_info("wrote data/zeal/zeal_%s_gt.rds", METHOD)
