#!/usr/bin/env Rscript
# Timing worker: nilHMM nnil over the full 888-line population, streaming rows from
# the PRE-SPLIT memory-mapped .bed for one density level (thin_L<level>/, written by
# materialize_thinned_bed.py). Loads ONLY the level it benchmarks; no in-script
# thinning. Peak RSS via parent /usr/bin/time -l.
#   Rscript 05_nilhmm_worker.R --level L
suppressMessages({
  library(nilHMM)
  library(BEDMatrix)
})
OUT <- "/Users/fvrodriguez/repos/zealhmm/data/nnil_equiv"
args <- commandArgs(trailingOnly = TRUE)
level <- as.integer(args[which(args == "--level") + 1L])
DIR <- file.path(OUT, sprintf("thin_L%d", level))

geno <- BEDMatrix(file.path(DIR, "geno.bed"))
md <- read.csv(file.path(DIR, "markers.csv"))
params <- jsonlite::fromJSON(file.path(OUT, "params.json"))
N <- nrow(geno)
M <- ncol(geno)
chrom_idx <- split(seq_len(M), md$chrom) # all columns of the pre-split panel, by chrom
r <- 2 * 1500 / (100 * M) # per-size avg_r (File S16 formula)
em <- emission_gt(germ = params$germ, gert = params$gert, p = params$p, mr = params$mr, nir = params$nir)
du <- duration_geometric(rrate = r)
priors <- list(f_1 = params$f_1, f_2 = params$f_2)

calls <- matrix(NA_integer_, N, M)
t0 <- Sys.time()
for (li in seq_len(N)) {
  g_all <- as.integer(geno[li, ])
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
  level, M, N, dt
))
