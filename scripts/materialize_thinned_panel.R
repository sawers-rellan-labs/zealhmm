#!/usr/bin/env Rscript
# Write odd-index-thinned copies of the shared-3 Arabidopsis panel, one directory per
# level, so the peak-RSS memory workers load ONLY the thinned working set for a given
# marker density. Without this, a worker that reads the full panel and thins in-memory
# has its peak RSS floored by the full-panel read at every size, flattening the memory
# curve; materializing the thinned set to disk (a separate, unmeasured process) makes
# each measured fit see only the data it actually uses.
#   Rscript scripts/materialize_thinned_panel.R
# Writes data/rtiger_shared3_input/thin_L{0..4}/sample{BN,Z,AU}.txt (gitignored).
ROOT <- "/Users/fvrodriguez/repos/zealhmm"
PANEL <- file.path(ROOT, "data/rtiger_shared3_input")
files <- c("sampleBN.txt", "sampleZ.txt", "sampleAU.txt")
LEVELS <- 0:4

thin <- function(n, level) {
  ix <- seq_len(n)
  for (j in seq_len(level)) ix <- ix[c(TRUE, FALSE)] # odd-index thin (matches build_obs)
  ix
}

full <- lapply(files, function(f) readLines(file.path(PANEL, f)))
names(full) <- files
n <- length(full[[1]])
stopifnot(all(vapply(full, length, 0L) == n)) # the three samples share the grid

for (level in LEVELS) {
  ix <- thin(n, level)
  outdir <- file.path(PANEL, sprintf("thin_L%d", level))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  for (f in files) writeLines(full[[f]][ix], file.path(outdir, f))
  cat(sprintf("[materialize] L%d: %d markers -> %s\n", level, length(ix), outdir))
}
