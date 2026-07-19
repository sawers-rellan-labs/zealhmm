#!/usr/bin/env Rscript
# Compare the three per-marker call sets on the identical shared input:
#   nilHMM nnil   vs  Holland File_S11   -> the EQUIVALENCE (target: 0 mismatches)
#   each          vs  published File S18 -> reproduction cross-check
# Reports overall + per-line concordance and donor-state (state>0) agreement.
#
#   Rscript scripts/nnil_equiv/04_compare.R

ROOT <- "/Users/fvrodriguez/repos/zealhmm"
OUT <- file.path(ROOT, "data/nnil_equiv")
rd <- function(f) as.matrix(read.csv(file.path(OUT, f), row.names = 1, check.names = FALSE))
log_info <- function(...) cat(sprintf("[04_compare] %s\n", sprintf(...)))

nh <- readRDS(file.path(OUT, "nilhmm_calls.rds")) # compact integer matrix
ho <- rd("holland_calls.csv")
s18 <- rd("s18_aligned.csv")
g <- rd("geno_recoded.csv") # {0,1,2,3}; 3 = missing
# align (same lines/markers, same order)
stopifnot(identical(dimnames(nh), dimnames(ho)))
s18 <- s18[rownames(nh), colnames(nh), drop = FALSE]
g <- g[rownames(nh), colnames(nh), drop = FALSE]

concord <- function(a, b) {
  ok <- a == b
  list(
    pos = mean(ok), n = length(ok), mism = sum(!ok),
    donor = mean((a > 0) == (b > 0))
  ) # introgressed-vs-not agreement
}
per_line_min <- function(a, b) min(rowMeans(a == b))

cat("\n================  nNIL caller equivalence  ================\n")
e <- concord(nh, ho)
log_info(
  "nilHMM vs Holland (File_S11): %d/%d states match (%.6f%%), %d mismatches",
  e$n - e$mism, e$n, 100 * e$pos, e$mism
)
log_info(
  "  introgressed-vs-not agreement: %.6f%% | worst per-line: %.4f%%",
  100 * e$donor, 100 * per_line_min(nh, ho)
)

# Stratify mismatches by genotype class: informative homozygous (g in {0,2})
# distinguish the states; missing (3) and het (1) are emission-degenerate (>=2
# states tied), so any Viterbi difference there is a boundary tie-break.
cat("\n--- mismatch stratification by observed genotype ---\n")
mmg <- g[nh != ho]
info_mm <- sum(mmg %in% c(0L, 2L)) # mismatches at informative homozygous markers
log_info(
  "informative homozygous markers (g=0 or 2): %d mismatches  <-- the exactness test",
  info_mm
)
log_info(
  "emission-tied positions: missing g=3: %d | het g=1: %d (states 0,2 emission-equal)",
  sum(mmg == 3L), sum(mmg == 1L)
)

cat("\n----------------  reproduction cross-check vs File S18  ----------------\n")
for (nm in c("nilHMM", "Holland")) {
  m <- if (nm == "nilHMM") nh else ho
  c18 <- concord(m, s18)
  log_info(
    "%s vs File S18: %.4f%% states, %.4f%% introgressed-vs-not (%d mismatches)",
    nm, 100 * c18$pos, 100 * c18$donor, c18$mism
  )
}

verdict <- if (e$mism == 0) "EQUIVALENT (0 mismatches)" else sprintf("%d mismatches (%.4f%%)", e$mism, 100 * (1 - e$pos))
cat(sprintf("\nVERDICT (nilHMM nnil vs Holland File_S11): %s\n", verdict))
