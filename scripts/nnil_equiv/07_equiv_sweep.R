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

geno <- BEDMatrix(file.path(OUT, "geno.bed"))
md <- read.csv(file.path(OUT, "markers.csv"))
params <- jsonlite::fromJSON(file.path(OUT, "params.json"))
N <- nrow(geno)

thin_idx <- function(level) {
  idx <- seq_len(ncol(geno))
  for (j in seq_len(level)) idx <- idx[c(TRUE, FALSE)] # odd-index thin (== python [::2])
  idx
}

# nilHMM nnil calls at one thinned size (same path as 03_nilhmm_calls / 05 worker).
nilhmm_calls <- function(idx) {
  chrom_idx <- split(seq_along(idx), md$chrom[idx])
  r <- 2 * 1500 / (100 * length(idx)) # per-size avg_r (File S16 formula)
  em <- emission_gt(germ = params$germ, gert = params$gert, p = params$p, mr = params$mr, nir = params$nir)
  du <- duration_geometric(rrate = r)
  priors <- list(f_1 = params$f_1, f_2 = params$f_2)
  calls <- matrix(NA_integer_, N, length(idx))
  for (li in seq_len(N)) {
    g_all <- as.integer(geno[li, idx])
    g_all[is.na(g_all)] <- 3L
    for (cc in names(chrom_idx)) {
      ix <- chrom_idx[[cc]]
      g <- g_all[ix]
      calls[li, ix] <- decode(fit(list(g = g), em, du, priors = priors), list(g = g))
    }
  }
  calls
}

rows <- list()
for (lv in LEVELS) {
  idx <- thin_idx(lv)
  M <- length(idx)
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
  # nilHMM calls at the same level
  nh <- nilhmm_calls(idx)
  # observed genotypes at the same markers (for the informative-marker stratification)
  g <- matrix(NA_integer_, N, M)
  for (li in seq_len(N)) {
    gi <- as.integer(geno[li, idx])
    gi[is.na(gi)] <- 3L
    g[li, ] <- gi
  }
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
