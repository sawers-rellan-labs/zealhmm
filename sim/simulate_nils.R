# sim/simulate_nils.R — simcross NIL generator + per-source degradation.
#
# Emits, per (design, source): a truth mosaic and degraded allele counts that
# look like a Skim or BRB sample, for calibration + validation (analysis note
# 02-simulation-calibration.qmd). Seeds are pinned for reproducibility (B4).
#
# STATUS: scaffold. The concrete BZea simcross design (generations, n, cM map,
# interference) is open (plan B5) and must be fixed with the user before the
# numbers are final. This file records the intended interface; the map + fitting
# primitives already live in the nilHMM package (load_map, expected_fragment_dist,
# calibrate_r, cm_to_mb, fit_design_gamma).

suppressMessages({
  library(nilHMM)
  library(data.table)
})
stopifnot(requireNamespace("simcross", quietly = TRUE))

# Reproducibility: fixed seed (B4). Vary only via the `seed` argument.
DEFAULT_SEED <- 1L

#' Simulate NILs for one design and degrade to one sequencing source
#'
#' @param design "BC2S2" (bulked skim) or "BC2S3".
#' @param source "skim", "brb", or "target".
#' @param n Number of NILs to simulate.
#' @param seed RNG seed (pinned).
#' @return list(truth = <common-schema segments>, counts = <per-marker counts>,
#'   grid = <chr,pos evaluation grid>). NOTE: body is a stub pending B5.
simulate_source <- function(design = c("BC2S3", "BC2S2"),
                            source = c("skim", "brb", "target"),
                            n = 200L, seed = DEFAULT_SEED) {
  design <- match.arg(design)
  source <- match.arg(source)
  set.seed(seed)

  # 1. map + pedigree ---------------------------------------------------------
  # map <- nilHMM::load_map()                       # consensus cM<->bp Marey map
  # ped <- <build BC2S2/BC2S3 pedigree>             # simcross::sim_from_pedigree
  #   check_pedigree(ped, ignore_sex = TRUE)        # selfing => mom == dad

  # 2. simulate mosaics -> per-locus truth genotypes at the panel markers -----
  #   donor = teosinte = ALT, recurrent = B73 = REF.

  # 3. degrade to `source` regime --------------------------------------------
  #   coverage-driven missingness: missing(lambda) = pi + (1-pi) e^{-k lambda}
  #   + per-read error; source-specific depth (skim ~0.4x, brb expression-driven).

  stop("simulate_source(): stub — fix the BZea simcross design (plan B5) first.")
}

if (sys.nframe() == 0L) {
  # CLI entry point (once implemented):
  #   Rscript sim/simulate_nils.R BC2S3 skim 200
  args <- commandArgs(trailingOnly = TRUE)
  message("sim/simulate_nils.R is a scaffold; see plan B5 before running.")
}
