#!/usr/bin/env Rscript
# nNIL equivalence sweep (mirror of the RTIGER equivalence table, tab:cpp-equiv):
# at each odd-index marker-density level, decode the FULL population with both cores
# (Holland File_S11 hmmlearn vs nilHMM nnil) on the identical thinned genotypes and
# the identical per-size r, then compare the calls position-by-position. Reports, per
# size, the states compared, the Viterbi match/mismatch counts, and the mismatches
# restricted to informative homozygous markers (g in {0,2}) -- the exactness test
# (missing/het positions are emission-degenerate, so any difference there is a
# boundary tie-break, not a decoding disagreement).
#   Rscript scripts/nnil_equiv/07_equiv_sweep.R
suppressMessages({
  library(nilHMM)
  library(BEDMatrix)
})
ROOT <- "/Users/fvrodriguez/repos/zealhmm"
OUT <- file.path(ROOT, "data/nnil_equiv")
SDIR <- file.path(ROOT, "scripts/nnil_equiv")
OUTDIR <- file.path(ROOT, "results/bench")
SWEEP <- file.path(OUT, "sweep")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(SWEEP, showWarnings = FALSE, recursive = TRUE)
PYBIN <- path.expand("~/anaconda3/envs/nilhmm/bin/python")
LEVELS <- 0:6 # 64025 -> ~1000 markers, full population
log_info <- function(...) cat(sprintf("[07_equiv_sweep] %s\n", sprintf(...)))

# SPLIT: every worker reads its own pre-split thin_L<level>/ (no in-script thinning);
# materialize them once up front if absent.
if (!file.exists(file.path(OUT, "thin_L0", "geno.bed"))) {
  log_info("materializing thinned .bed per level ...")
  system2(PYBIN, file.path(SDIR, "materialize_thinned_bed.py"))
}

params <- jsonlite::fromJSON(file.path(OUT, "params.json"))
N <- length(readLines(file.path(OUT, "lines.csv"))) # lines, constant across levels

# nilHMM nnil calls + observed genotypes for one PRE-SPLIT level dir (thin_L<level>/,
# written by materialize_thinned_bed.py); loads only that dir, no in-script thinning.
nilhmm_calls <- function(dir) {
  geno <- BEDMatrix(file.path(dir, "geno.bed"))
  md <- read.csv(file.path(dir, "markers.csv"))
  M <- ncol(geno)
  chrom_idx <- split(seq_len(M), md$chrom)
  r <- 2 * 1500 / (100 * M) # per-size avg_r (File S16 formula)
  em <- emission_gt(germ = params$germ, gert = params$gert, p = params$p, mr = params$mr, nir = params$nir)
  du <- duration_geometric(rrate = r)
  priors <- list(f_1 = params$f_1, f_2 = params$f_2)
  calls <- matrix(NA_integer_, N, M)
  g <- matrix(NA_integer_, N, M) # observed genotypes (for the informative-marker strat)
  for (li in seq_len(N)) {
    g_all <- as.integer(geno[li, ])
    g_all[is.na(g_all)] <- 3L
    g[li, ] <- g_all
    for (cc in names(chrom_idx)) {
      ix <- chrom_idx[[cc]]
      gg <- g_all[ix]
      calls[li, ix] <- decode(fit(list(g = gg), em, du, priors = priors), list(g = gg))
    }
  }
  list(calls = calls, g = g, M = M)
}

rows <- list()
for (lv in LEVELS) {
  # nilHMM calls + observed genotypes from the pre-split level dir
  nc <- nilhmm_calls(file.path(OUT, sprintf("thin_L%d", lv)))
  nh <- nc$calls
  g <- nc$g
  M <- nc$M
  # Holland calls at this level, via the python worker (writes a C-order int8 binary)
  binf <- file.path(SWEEP, sprintf("holland_L%d.bin", lv))
  out <- system2(PYBIN, c(
    file.path(SDIR, "07_holland_level.py"),
    "--level", lv, "--out", binf
  ), stdout = TRUE, stderr = TRUE)
  res <- grep("^RESULT", out, value = TRUE)
  if (!length(res)) stop(sprintf("holland worker failed at level %d:\n%s", lv, paste(out, collapse = "\n")))
  if (file.info(binf)$size != N * M) {
    stop(sprintf(
      "holland worker wrote %d bytes at level %d, expected %d (N*M)",
      file.info(binf)$size, lv, N * M
    ))
  }
  ho <- matrix(readBin(binf, "integer", n = N * M, size = 1L, signed = TRUE),
    nrow = N, ncol = M, byrow = TRUE
  )
  states <- N * M
  mm <- nh != ho
  mism <- sum(mm)
  info_mm <- sum(g[mm] %in% c(0L, 2L)) # mismatches at informative homozygous markers
  rows[[length(rows) + 1L]] <- data.frame(
    level = lv, markers = M, lines = N, states = states,
    match = states - mism, mismatches = mism, informative_mismatches = info_mm,
    mismatch_rate = mism / states,
    agreement = 100 * (states - mism) / states
  )
  log_info(
    "L%d %6d markers x %d lines = %d states | match %d | mism %d | informative-mism %d (%.6f%%)",
    lv, M, N, states, states - mism, mism, info_mm, 100 * (states - mism) / states
  )
  file.remove(binf)
}
d <- do.call(rbind, rows)
write.csv(d, file.path(OUTDIR, "nnil_equiv_sweep.csv"), row.names = FALSE)
log_info("wrote results/bench/nnil_equiv_sweep.csv")
cat("\n")
print(d, row.names = FALSE)
