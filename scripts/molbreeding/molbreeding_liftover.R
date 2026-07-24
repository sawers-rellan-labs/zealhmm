#!/usr/bin/env Rscript
# Liftover the MolBreeding target sequencing marker positions to B73 NAM v5.
#
# The MolBreeding panel is on AGPv3 (GRMZM2G… gene IDs), the rest of the project
# is B73 NAM v5, and there is no direct v3→v5 chain — so lift in two hops:
#   AGPv3 ──(AGPv3_to_B73_RefGen_v4.chain)──▶ v4 ──(B73_RefGen_v4_to_…NAM-5.0.chain)──▶ v5
# Chains live in data/chain_files/ (tab-converted for rtracklayer; see its README).
#
# Input : data/GSER2026030032P01/<set>/All.GT.map  (PLINK map: chr, marker, cM, bp; v3)
# Output: data/molbreeding_45k/sites_v5_<set>.tsv
#         (marker, chr_v5, pos_v5, strand_v5, chr_v3, pos_v3) — unique 1:1 lifts on chr 1–10.
#
# Usage:  Rscript molbreeding_liftover.R [SNP mSNP ...]   (default: both)

suppressMessages({
  library(rtracklayer)
  library(GenomicRanges)
})

chain_dir <- "data/chain_files"
out_dir <- "data/molbreeding_45k"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sets <- commandArgs(trailingOnly = TRUE)
if (length(sets) == 0) sets <- c("SNP", "mSNP")

# --- one liftOver hop, keeping only unique 1:1 mappings -----------------------
lift_unique <- function(gr, chain_file, label) {
  ch <- import.chain(file.path(chain_dir, chain_file))
  lifted <- liftOver(gr, ch) # GRangesList, parallel to gr
  n <- elementNROWS(lifted)
  out <- unlist(lifted[n == 1]) # mcols (marker) carried through
  cat(sprintf(
    "  %s: %d in → %d unique (%d lost, %d multi-map)\n",
    label, length(gr), sum(n == 1), sum(n == 0), sum(n > 1)
  ))
  out
}

lift_set <- function(setname) {
  map_path <- sprintf("data/GSER2026030032P01/%s/All.GT.map", setname)
  out_path <- file.path(out_dir, sprintf("sites_v5_%s.tsv", setname))
  cat(sprintf("\n=== %s ===\n", setname))

  map <- read.table(map_path, header = FALSE, stringsAsFactors = FALSE)
  names(map) <- c("chr", "marker", "cm", "bp")
  map$chr <- as.character(map$chr)
  cat(sprintf("MolBreeding markers (v3): %d\n", nrow(map)))

  gr_v3 <- GRanges(
    seqnames = map$chr,
    ranges = IRanges(start = map$bp, width = 1),
    marker = map$marker
  )

  cat("Liftover v3 → v4 → v5:\n")
  gr_v4 <- lift_unique(gr_v3, "AGPv3_to_B73_RefGen_v4.chain", "v3→v4")
  gr_v5 <- lift_unique(gr_v4, "B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain", "v4→v5")

  res <- data.frame(
    marker = gr_v5$marker,
    chr_v5 = as.character(seqnames(gr_v5)),
    pos_v5 = start(gr_v5),
    strand_v5 = as.character(strand(gr_v5)),
    stringsAsFactors = FALSE
  )
  res <- res[res$chr_v5 %in% as.character(1:10), ]
  res <- merge(res, map[, c("marker", "chr", "bp")], by = "marker", all.x = TRUE)
  names(res)[names(res) == "chr"] <- "chr_v3"
  names(res)[names(res) == "bp"] <- "pos_v3"
  moved_chr <- res$chr_v5 != res$chr_v3
  if (any(moved_chr)) cat(sprintf("  note: %d markers changed chromosome v3→v5 (dropped)\n", sum(moved_chr)))
  res <- res[!moved_chr, ]
  res <- res[order(as.integer(res$chr_v5), res$pos_v5), ]

  write.table(res, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf(
    "Retained %d / %d markers (%.1f%%) on v5 chr1–10 → %s\n",
    nrow(res), nrow(map), 100 * nrow(res) / nrow(map), out_path
  ))
}

for (s in sets) lift_set(s)
