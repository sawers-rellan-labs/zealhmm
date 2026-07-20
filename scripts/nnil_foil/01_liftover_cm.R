#!/usr/bin/env Rscript
# Calibration foil, step 1-2: put the nNIL marker set on the maize v5 coordinate
# system and give every marker a genetic (cM) position from the TeoNAM native map.
#
#   1. Lift the 64,025 nNIL markers AGPv4 -> B73 NAM v5 (single hop) with the
#      vendored UCSC chain, keeping only unique 1:1 lifts on chr 1-10 that do not
#      change chromosome (lift_unique(), R/teonam_liftover.R). Re-sort by v5 order.
#   2. Interpolate cM at the lifted v5 bp from the TeoNAM native v5 map
#      (data/teonam/markers_v5_gwas118k_cm.tsv) via the engine's per-chromosome
#      monotone Hyman spline (nilHMM::bp_to_cm), clamped to each chr's observed
#      range -- the same native-cM convention used across the TeoNAM analyses.
#
# Both the chip-side calibration (real GBS re-keyed to this v5 order) and the
# simcross sim run on THIS marker set, so `rrate` is a per-adjacent-marker rate in
# one coordinate system on both sides of the foil.
#
# Marker IDs are renamed to v5 here and used as v5 (`S<chr>_<pos_v5>`) everywhere
# downstream; the original v4 name is kept only as provenance (`marker_v4`).
#
#   Rscript scripts/nnil_foil/01_liftover_cm.R
# Output: data/nnil_foil/markers_v5.tsv (marker, marker_v4, chr, pos_v4, pos_v5, cm),
#   v5-ordered. `marker` = v5 id; join the v4-keyed geno.bed on `marker_v4`.

suppressMessages({
  library(nilHMM)
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
})
root <- here::here()
source(file.path(root, "R/teonam_liftover.R")) # lift_unique()
source(file.path(root, "scripts/logging.R")) # log_info

CHAIN <- file.path(
  root, "data/ref/chain_files",
  "B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain"
)
MARKERS_V4 <- file.path(root, "data/nnil_equiv/markers.csv") # marker, chrom, pos (v4)
NATIVE_MAP <- file.path(root, "data/teonam/markers_v5_gwas118k_cm.tsv")
OUT <- file.path(root, "data/nnil_foil/markers_v5.tsv")

stopifnot(file.exists(CHAIN), file.exists(MARKERS_V4), file.exists(NATIVE_MAP))

# ---- 1. lift v4 -> v5 -------------------------------------------------------
m4 <- fread(MARKERS_V4, data.table = FALSE)
m4$chrom <- as.character(m4$chrom)
ok <- m4$chrom %in% as.character(1:10) & !is.na(m4$pos)
gr <- GRanges(
  seqnames = m4$chrom[ok],
  ranges = IRanges(start = as.integer(m4$pos[ok]), width = 1L),
  marker = m4$marker[ok]
)
n_in <- length(gr)
log_info("lifting %d nNIL markers AGPv4 -> v5 ...", n_in)

gr5 <- lift_unique(gr, CHAIN)
lift <- data.frame(
  marker = gr5$marker,
  chr = as.character(seqnames(gr5)),
  pos_v5 = start(gr5),
  stringsAsFactors = FALSE
)
# original v4 coords, drop chromosome-changers, keep chr 1-10
src <- data.frame(
  marker = m4$marker[ok], chr_v4 = m4$chrom[ok], pos_v4 = as.integer(m4$pos[ok]),
  stringsAsFactors = FALSE
)
lift <- merge(src, lift, by = "marker")
n_uniq <- nrow(lift)
lift <- lift[lift$chr %in% as.character(1:10) & lift$chr == lift$chr_v4, ]
n_samechr <- nrow(lift)
lift <- lift[order(as.integer(lift$chr), lift$pos_v5), ]
log_info(
  "lifted 1:1 = %d (%.1f%%); same-chromosome kept = %d (dropped %d chr-changers/off-target)",
  n_uniq, 100 * n_uniq / n_in, n_samechr, n_uniq - n_samechr
)

# ---- 2. cM from the TeoNAM native v5 map ------------------------------------
native <- fread(NATIVE_MAP, data.table = FALSE) # marker, chr, pos_v5, cm
map <- data.frame(chr = as.character(native$chr), bp = native$pos_v5, cm = native$cm)
to_cm <- bp_to_cm(map) # engine Hyman spline, clamped per chr
lift$cm <- to_cm(lift$chr, lift$pos_v5)

# how many lifted markers fall outside the native map's per-chr bp range (clamped)?
rng <- tapply(map$bp, map$chr, range)
clamp_lo <- clamp_hi <- 0L
for (ch in unique(lift$chr)) {
  r <- rng[[ch]]
  sel <- lift$chr == ch
  clamp_lo <- clamp_lo + sum(lift$pos_v5[sel] < r[1])
  clamp_hi <- clamp_hi + sum(lift$pos_v5[sel] > r[2])
}
log_info(
  "cM interpolated; %d markers clamped below / %d above the native map range",
  clamp_lo, clamp_hi
)

# rename marker IDs to v5 (S<chr>_<pos_v5>); keep the v4 name as provenance
lift$marker_v4 <- lift$marker
lift$marker <- sprintf("S%s_%d", lift$chr, lift$pos_v5)
if (anyDuplicated(lift$marker)) {
  dup <- sum(duplicated(lift$marker))
  log_info("WARNING: %d duplicate v5 ids after rename; keeping first per id", dup)
  lift <- lift[!duplicated(lift$marker), ]
}
out <- lift[, c("marker", "marker_v4", "chr", "pos_v4", "pos_v5", "cm")]
fwrite(out, OUT, sep = "\t")
log_info(
  "wrote %s : %d markers (v5 ids), cM span %.1f..%.1f (total %.0f cM)",
  OUT, nrow(out), min(out$cm), max(out$cm),
  sum(tapply(out$cm, out$chr, function(x) max(x) - min(x)))
)
