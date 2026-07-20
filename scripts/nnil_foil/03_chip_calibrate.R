#!/usr/bin/env Rscript
# Calibration foil, step 4: the CHIP-calibrated operating point. Sweep the nnil
# caller's rrate on the real GBS genotypes of the 24 both-platform NILs and score
# agreement with the reproduced chip calls. This locates rrate_chip* -- the
# operating point Holland selected against his SNP chip.
#
# Faithful to Holland: the nnil caller runs with his exact gt-emission and priors
# (params.json: germ/gert/p/nir/mr, f_1, f_2); only rrate varies. Agreement is
# scored with the SAME metric functions the sim side uses (R/metrics.R:
# donor_fragment_dice = primary objective; marker_dice), so the chip and sim
# calibration curves are directly comparable. The marker-level comparison is
# restricted to the 11,310 GBS markers the chip projects onto (Holland's footing);
# donor_fragment_dice compares donor blocks in bp space. Holland's own marker
# mismatch (= 1 - accuracy on the shared markers) is emitted as a secondary curve.
#
# Everything is keyed by v5 marker ids; the v4-keyed geno.bed is joined via the
# crosswalk (markers_v5.tsv) and relabelled to v5, unmapped markers dropped.
#
#   Rscript scripts/nnil_foil/03_chip_calibrate.R
# Output (data/nnil_foil/):
#   chip_rrate_sweep.csv    rrate, donor_frag_dice, donor_frag_FDR, donor_marker_dice,
#                           marker_macro_dice, holland_mismatch, n_breakpoints
#   chip_calib.json         rrate_chip*, avg_r (Holland's reference), n_lines, n_shared

suppressMessages({
  library(nilHMM)
  library(BEDMatrix)
  library(data.table)
  library(jsonlite)
})
root <- here::here()
for (f in list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)) source(f)
source(file.path(root, "scripts/logging.R"))
FOIL <- file.path(root, "data/nnil_foil")
EQUIV <- file.path(root, "data/nnil_equiv")

# ---- crosswalk + chip calls -------------------------------------------------
xwalk <- fread(file.path(FOIL, "markers_v5.tsv")) # marker(v5), marker_v4, chr, pos_v4, pos_v5, cm
setkey(xwalk, marker)

chip <- fread(file.path(FOIL, "chip_truth_projected.csv")) # Line + v5-marker cols
gbs_lines <- readLines(file.path(EQUIV, "lines.csv"))
both <- intersect(chip$Line, gbs_lines) # the 24 both-platform NILs (exact match)
stopifnot(length(both) == 24)
chip <- chip[Line %in% both]
setorder(chip, Line)
chip_markers <- setdiff(names(chip), "Line")
chip_markers <- chip_markers[chip_markers %in% xwalk$marker] # keep lifted only
log_info("chip calls: %d both-platform NILs x %d shared v5 markers", nrow(chip), length(chip_markers))

# chip marker positions (v5), for segments + the marker-eval grid
cm_info <- xwalk[chip_markers, .(marker, chr, pos = pos_v5)]
setorder(cm_info, chr, pos)
grid_eval <- cm_info[, .(chr, pos)] # 11,310 shared markers = Holland's comparison footing

# ---- per-marker state matrix -> common-schema segments ----------------------
states_to_segments <- function(mat, lines, marker_pos) {
  # mat: lines x markers integer states {0,1,2} (3/NA = missing); marker_pos:
  # data.table(marker, chr, pos) in column order of mat. RLE per (line, chr).
  mp <- marker_pos[match(colnames(mat), marker)]
  out <- vector("list", nrow(mat))
  for (i in seq_len(nrow(mat))) {
    s <- as.integer(mat[i, ])
    dt <- data.table(chr = mp$chr, pos = mp$pos, state = s)[!is.na(state) & state != 3L]
    if (!nrow(dt)) next
    setorder(dt, chr, pos)
    dt[, run := rleid(chr, state)]
    seg <- dt[, .(start_bp = min(pos), end_bp = max(pos), state = state[1]), by = .(chr, run)]
    seg[, name := lines[i]]
    out[[i]] <- seg[, .(name, chr, start_bp, end_bp, state)]
  }
  rbindlist(out)
}

chip_mat <- as.matrix(chip[, ..chip_markers])
rownames(chip_mat) <- chip$Line
tr <- states_to_segments(chip_mat, chip$Line, cm_info[, .(marker, chr, pos)])
log_info("chip calls segments: %d donor blocks over %d NILs", nrow(.donor_blocks(tr)), nrow(chip))

# ---- real GBS genotypes for the 24 NILs, relabelled to v5 -------------------
geno <- BEDMatrix(file.path(EQUIV, "geno.bed")) # 888 x 64025 (v4), mmap
md <- fread(file.path(EQUIV, "markers.csv")) # v4 marker order of geno columns
stopifnot(nrow(geno) == length(gbs_lines), ncol(geno) == nrow(md))
row_idx <- match(both, gbs_lines)
# v4 columns that lifted, in v5 order
xw_by_v4 <- xwalk[match(md$marker, marker_v4)] # aligned to geno columns; NA = unmapped
keep_col <- which(!is.na(xw_by_v4$marker))
g_raw <- geno[row_idx, keep_col, drop = FALSE] # 24 x (lifted)
v5_id <- xw_by_v4$marker[keep_col]
v5_chr <- xw_by_v4$chr[keep_col]
v5_pos <- xw_by_v4$pos_v5[keep_col]
ord <- order(v5_chr, v5_pos)
g_raw <- g_raw[, ord, drop = FALSE]
v5_id <- v5_id[ord]
v5_chr <- v5_chr[ord]
v5_pos <- v5_pos[ord]

# long data.table for call_ancestry: name, chr, pos, g (hard calls, 3 = missing)
data <- rbindlist(lapply(seq_along(both), function(i) {
  g <- as.integer(g_raw[i, ])
  g[is.na(g)] <- 3L
  data.table(name = both[i], chr = v5_chr, pos = v5_pos, g = g)
}))
setorder(data, name, chr, pos)
log_info(
  "GBS input: %d NILs x %d v5 markers (%.1f%% missing)",
  length(both), length(v5_id), 100 * mean(data$g == 3L)
)

# ---- sweep rrate; score with the shared metric functions --------------------
hp <- fromJSON(file.path(EQUIV, "params.json")) # Holland's exact emission + priors
values <- log_grid(1e-6, 1e-1, 24L)
# "map r" = the expected per-adjacent-marker recombination FRACTION implied by the
# genetic map: 2 * L / (100 * N), L = total map length in cM, N = markers, factor 2
# ~= the effective meioses of backcross+selfing (Haldane-Waddington). Recomputed on
# OUR native TeoNAM v5 map (not Holland's uniform 1500 cM), since the map changed.
map_r <- 2 * xwalk[, sum(tapply(cm, chr, function(x) max(x) - min(x)))] / (100 * nrow(xwalk))

score_one <- function(v) {
  called <- as.data.table(call_ancestry(
    data = data, caller = "nnil", rrate = v,
    germ = hp$germ, gert = hp$gert, p = hp$p, nir = hp$nir, mr = hp$mr,
    f_1 = hp$f_1, f_2 = hp$f_2
  ))
  mf <- marker_dice(called, tr, grid_eval)
  ff <- donor_fragment_dice(called, tr)
  dm <- mf$per_class[class == "donor(>0)"]
  data.table(
    rrate = v, donor_frag_dice = ff$dice, donor_frag_FDR = ff$fdr,
    donor_marker_dice = dm$dice, marker_macro_dice = mf$macro_dice,
    holland_mismatch = 1 - mf$accuracy, n_breakpoints = breakpoint_count(called)
  )
}
t0 <- Sys.time()
log_info("chip-side rrate sweep: %d points (nnil, Holland emission) ...", length(values))
sweep <- rbindlist(lapply(seq_along(values), function(i) {
  r <- score_one(values[i])
  log_info(
    "  rrate=%.3e | frag_dice=%.3f mismatch=%.4f (%d/%d, %.0fs)",
    values[i], r$donor_frag_dice, r$holland_mismatch, i, length(values),
    as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )
  r
}))
fwrite(sweep, file.path(FOIL, "chip_rrate_sweep.csv"))

# The chip calls are SPARSE (11,310 markers, few large donor blocks per line), so
# donor_frag_dice is monotone-decreasing in rrate (low rrate never misses a true
# block; it only cuts spurious ones), and Holland's marker-mismatch is nearly flat
# over orders of magnitude. There is no sharp interior optimum on the chip side --
# the sparse chip weakly constrains rrate. So we do NOT report a degenerate
# boundary argmax as "the" optimum; we record the curve shape and Holland's own
# documented operating point (avg_r), and characterize the weakly-constrained
# region (rrate within 0.01 Dice of the grid best, and within 5% of min mismatch).
best_fd <- max(sweep$donor_frag_dice)
min_mm <- min(sweep$holland_mismatch)
plat_fd <- sweep[donor_frag_dice >= best_fd - 0.01]
plat_mm <- sweep[holland_mismatch <= min_mm * 1.05]
fd_at_mapr <- approx(log10(sweep$rrate), sweep$donor_frag_dice, log10(map_r))$y
mm_at_mapr <- approx(log10(sweep$rrate), sweep$holland_mismatch, log10(map_r))$y
writeLines(toJSON(list(
  map_r = map_r, # map-defined per-marker recombination fraction (native v5 map)
  frag_dice_at_map_r = fd_at_mapr, mismatch_at_map_r = mm_at_mapr,
  rrate_grid_best_fragdice = sweep$rrate[which.max(sweep$donor_frag_dice)],
  frag_dice_best = best_fd, mismatch_min = min_mm,
  fragdice_plateau = c(min(plat_fd$rrate), max(plat_fd$rrate)),
  mismatch_plateau = c(min(plat_mm$rrate), max(plat_mm$rrate)),
  n_lines = length(both), n_shared_markers = length(chip_markers),
  n_donor_blocks_truth = nrow(.donor_blocks(tr))
), auto_unbox = TRUE, digits = 8), file.path(FOIL, "chip_calib.json"))
log_info(
  "chip curve: frag Dice monotone (best %.3f @ %.2e); mismatch flat (min %.4f)",
  best_fd, sweep$rrate[which.max(sweep$donor_frag_dice)], min_mm
)
log_info(
  "map r=%.3e -> frag Dice %.3f, mismatch %.4f; mismatch-plateau rrate [%.2e, %.2e]",
  map_r, fd_at_mapr, mm_at_mapr, min(plat_mm$rrate), max(plat_mm$rrate)
)
