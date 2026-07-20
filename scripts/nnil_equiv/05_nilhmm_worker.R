#!/usr/bin/env Rscript
# Timing worker: nilHMM nnil at one thinned marker size (odd-index thinning, as in
# the RTIGER sweep), over the FULL 884-line population, streaming rows from the
# memory-mapped .bed. Peak RSS via parent /usr/bin/time -l.
#   Rscript 05_nilhmm_worker.R --level L
suppressMessages({
  library(nilHMM)
  library(BEDMatrix)
})
OUT <- "/Users/fvrodriguez/repos/zealhmm/data/nnil_equiv"
args <- commandArgs(trailingOnly = TRUE)
level <- as.integer(args[which(args == "--level") + 1L])

geno <- BEDMatrix(file.path(OUT, "geno.bed"))
md <- read.csv(file.path(OUT, "markers.csv"))
params <- jsonlite::fromJSON(file.path(OUT, "params.json"))
N <- nrow(geno)
idx <- seq_len(ncol(geno))
for (j in seq_len(level)) idx <- idx[c(TRUE, FALSE)] # odd-index thin
chrom_idx <- split(seq_along(idx), md$chrom[idx])
r <- 2 * 1500 / (100 * length(idx)) # per-size avg_r (File S16 formula)
em <- emission_gt(germ = params$germ, gert = params$gert, p = params$p, mr = params$mr, nir = params$nir)
du <- duration_geometric(rrate = r)
priors <- list(f_1 = params$f_1, f_2 = params$f_2)

calls <- matrix(NA_integer_, N, length(idx))
t0 <- Sys.time()
for (li in seq_len(N)) {
  g_all <- as.integer(geno[li, idx])
  g_all[is.na(g_all)] <- 3L
  for (cc in names(chrom_idx)) {
    ix <- chrom_idx[[cc]]
    g <- g_all[ix]
    calls[li, ix] <- decode(fit(list(g = g), em, du, priors = priors), list(g = g))
  }
}
dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf(
  "RESULT caller=nilhmm level=%d markers=%d lines=%d seconds=%.4f\n",
  level, length(idx), N, dt
))
