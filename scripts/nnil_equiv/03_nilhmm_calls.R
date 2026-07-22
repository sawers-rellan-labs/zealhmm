#!/usr/bin/env Rscript
# Run nilHMM's nnil caller on the shared genotype table, STREAMING rows from the
# memory-mapped PLINK .bed (BEDMatrix) instead of slurping the 114 MB wide CSV.
# Each line's row is pulled from disk on demand, decoded, and discarded, so peak
# memory is ~R baseline + the output matrix -- flat in the (held) input. Uses the
# low-level engine: emission_gt (categorical, Holland's matrix verbatim) +
# duration_geometric (f-weighted transition) + priors -> fit() -> decode() (C++
# Viterbi, incumbent tie-break for GT emission). One (line, chromosome) chain at a
# time, matching Holland's per-line-per-chromosome decode.
#   Rscript scripts/nnil_equiv/03_nilhmm_calls.R
# Output: data/nnil_equiv/nilhmm_calls.rds (compact integer matrix, states {0,1,2}).

suppressMessages({
  library(nilHMM)
  library(BEDMatrix)
})
ROOT <- "/Users/fvrodriguez/repos/zealhmm"
OUT <- file.path(ROOT, "data/nnil_equiv")
log_info <- function(...) cat(sprintf("[03_nilhmm] %s\n", sprintf(...)))

geno <- BEDMatrix(file.path(OUT, "geno.bed")) # mmap; nothing loaded yet
md <- read.csv(file.path(OUT, "markers.csv")) # marker + chrom + pos, .bed column order
params <- jsonlite::fromJSON(file.path(OUT, "params.json"))
# Real line names from the sidecar (the .fam IIDs are positional/sanitised because
# line names contain spaces). Row i of the .bed == lines[i]; columns follow md.
lines <- readLines(file.path(OUT, "lines.csv"))
mk <- md$marker
stopifnot(nrow(geno) == length(lines), ncol(geno) == length(mk))
chrom_idx <- split(seq_along(mk), md$chrom)
log_info(
  "streaming %d lines x %d markers from geno.bed; r=%.3e mr=%.4f",
  length(lines), length(mk), params$r, params$mr
)

em <- emission_gt(
  germ = params$germ, gert = params$gert, p = params$p,
  mr = params$mr, nir = params$nir
)
du <- duration_geometric(rrate = params$r)
priors <- list(f_1 = params$f_1, f_2 = params$f_2)

calls <- matrix(NA_integer_, length(lines), length(mk), dimnames = list(lines, mk))
t0 <- Sys.time()
for (li in seq_along(lines)) {
  row <- geno[li, ] # one row pulled from disk (mmap)
  g_all <- as.integer(row)
  g_all[is.na(g_all)] <- 3L # NA (missing) -> 3
  for (cc in names(chrom_idx)) {
    ix <- chrom_idx[[cc]]
    g <- g_all[ix]
    calls[li, ix] <- decode(fit(list(g = g), em, du, priors = priors), list(g = g))
  }
  if (li %% 100 == 0 || li == length(lines)) {
    log_info(
      "  %d/%d lines (%.1fs)", li, length(lines),
      as.numeric(difftime(Sys.time(), t0, units = "secs"))
    )
  }
}
saveRDS(calls, file.path(OUT, "nilhmm_calls.rds"))
fr <- round(prop.table(table(factor(calls, levels = 0:2))), 4)
log_info(
  "wrote nilhmm_calls.rds %dx%d; state fractions 0=%.4f 1=%.4f 2=%.4f",
  nrow(calls), ncol(calls), fr[1], fr[2], fr[3]
)
