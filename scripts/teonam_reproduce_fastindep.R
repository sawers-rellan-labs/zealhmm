#!/usr/bin/env Rscript
# =============================================================================
# Reproduce the JLM marker thinning EXACTLY by running the real deterministic
# FastIndep greedy (nilHMM::select_independent, n_runs=1 -> bit-identical to the
# FastIndep CLI) on the per-chromosome cM-DISTANCE matrix, and compare the
# selected set to the 6059 markers actually in data/teonam/tassel/geno.hmp.txt.
#
# Prior diagnostic used a naive left-to-right greedy (6112) as a stand-in and
# hand-waved the 53-marker gap. Since the greedy is deterministic, that gap should
# be zero if (a) the algorithm and (b) the INPUT POOL both match. This script
# tests the full marker_info pool; if the set is not identical, it prints the
# discrepancy so we can pin down the pool the prior run actually thinned.
#
# sense="distance", threshold=0.1  -> edge iff |cm_i - cm_j| within 0.1 cM.
# Read-only. Run: Rscript agent/teonam_reproduce_fastindep.R
# =============================================================================
suppressMessages({
  library(data.table)
  devtools::load_all("/Users/fvrodriguez/repos/nilhmm", quiet = TRUE)
})
ROOT <- "/Users/fvrodriguez/repos/zealhmm"
setwd(ROOT)

THR <- 0.1

# target: the markers actually retained in the TASSEL JLM input
hmp <- fread("data/teonam/tassel/geno.hmp.txt", select = "rs#")
target <- hmp[["rs#"]]
cat(sprintf("target (geno.hmp.txt): %d markers\n", length(target)))

mc <- fread("data/teonam/marker_info_v5_cm.tsv") # candidate pool: all markers w/ cM
setnames(mc, "chr_v5", "chr")
cat(sprintf(
  "candidate pool (marker_info_v5_cm): %d markers, %d chromosomes\n\n",
  nrow(mc), uniqueN(mc$chr)
))

run_chr <- function(m) { # m: data.table(marker, cm) for one chr
  m <- m[order(cm)]
  D <- abs(outer(m$cm, m$cm, `-`)) # cM distance matrix (markers x markers)
  dimnames(D) <- list(m$marker, m$marker)
  sel <- select_independent(D,
    threshold = THR, n_runs = 1L,
    sense = "distance", max_markers = nrow(D) + 1L
  )
  as.character(sel)
}

compare <- function(sel, tag) {
  inter <- length(intersect(sel, target))
  cat(sprintf("\n[%s] selected %d markers vs %d target:\n", tag, length(sel), length(target)))
  cat(sprintf(
    "  in both: %d | selected-only: %d | target-only: %d | Jaccard: %.4f\n",
    inter, length(setdiff(sel, target)), length(setdiff(target, sel)),
    inter / length(union(sel, target))
  ))
  if (setequal(sel, target)) cat(sprintf("  ==> EXACT reproduction with the '%s' pool.\n", tag))
  setequal(sel, target)
}

# pool A: raw marker_info (every marker with a cM)
selA <- unlist(lapply(sort(unique(mc$chr)), function(ch) run_chr(mc[chr == ch, .(marker, cm)])))
compare(selA, "raw marker_info")

# pool B: unique-cM dedup per chr (the documented densify grid: mch[!duplicated(cm)])
poolB <- mc[order(chr, cm)][, .SD[!duplicated(cm)], by = chr]
selB <- unlist(lapply(sort(unique(poolB$chr)), function(ch) run_chr(poolB[chr == ch, .(marker, cm)])))
compare(selB, "unique-cM dedup")
