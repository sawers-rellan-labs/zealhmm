#!/usr/bin/env Rscript
# 0.1 cM ancestry-inference grid for the 118K sweeps — computed ONCE and cached.
#
# ALGORITHM. Keeping markers pairwise >=0.1 cM apart is the Maximum Distance-d
# Independent Set (MDdIS) problem: the largest vertex set whose pairwise distance is
# >= d (d=2 is ordinary MIS). On a general unit-disk graph (a 2-D point set) MDdIS is
# NP-hard -- Jena, Jallu, Das & Nandy (2018) prove NP-completeness for d>=3 and give a
# 4-factor approximation and a PTAS. That is why a generic solver over the full O(n^2)
# distance matrix is both heavy and only heuristic, and why the JLM `select_independent`
# route (which materializes that matrix) OOMs on the dense 118K panel.
#   Here the markers live on the 1-D cM axis, i.e. the disk/interval graph of an
# INTERVAL family, on which MDdIS is polynomially solvable (Jena et al. 2018, Related
# Work: interval / trapezoid / circular-arc graphs). For points on a line with a
# minimum gap d it collapses to the classic interval-scheduling greedy, which is EXACT
# and optimal: sort by cM, keep the leftmost, skip everything within d, repeat --
# O(n) time, O(n) memory, no distance matrix. So the sweep below is the exact optimum,
# not an approximation.
#
# Ref: Jena, S. K., Jallu, R. K., Das, G. K., & Nandy, S. C. (2018). "The Maximum
#   Distance-d Independent Set Problem on Unit Disk Graphs." International Workshop on
#   Frontiers in Algorithmics (FAW 2018), LNCS, Springer. (PDF: data/main.pdf)
#
# The ancestry HMMs infer segments on this thinned grid (~20x fewer markers → ~20x
# faster); segments are back-projected onto the full 118K union for the GWAS. The
# grid depends only on the cM map — not on family/coverage/caller — so it is built
# once here and read by every teonam_*_sweep_118k.R.
#
# Output: data/teonam/markers_v5_gwas118k_cm_thin01.tsv (marker, chr, pos_v5, cm)
# Run: Rscript scripts/teonam_gwas118k_thin01.R
suppressMessages(library(data.table))
setwd("/Users/fvrodriguez/repos/zealhmm")
THR <- 0.1

mc <- fread("data/teonam/markers_v5_gwas118k_cm.tsv") # marker, chr, pos_v5, cm
setorder(mc, chr, cm, pos_v5)
# per-chr left-to-right greedy: keep a marker iff it is >=THR cM past the last kept
# (ties at equal cM collapse automatically: gap 0 < THR). Optimal max independent set.
keep_idx <- mc[,
  {
    last <- -Inf
    k <- logical(.N)
    for (i in seq_len(.N)) {
      if (cm[i] - last >= THR) {
        k[i] <- TRUE
        last <- cm[i]
      }
    }
    .I[k]
  },
  by = chr
]$V1
thin <- mc[keep_idx]
setorder(thin, chr, cm)

OUT <- "data/teonam/markers_v5_gwas118k_cm_thin01.tsv"
fwrite(thin[, .(marker, chr, pos_v5, cm)], OUT, sep = "\t")
cat(sprintf(
  "0.1 cM inference grid: %d -> %d markers (%.1f%% kept) -> %s\n",
  nrow(mc), nrow(thin), 100 * nrow(thin) / nrow(mc), OUT
))
print(thin[, .(n = .N, cm_span = round(max(cm) - min(cm), 1)), by = chr][order(chr)])
