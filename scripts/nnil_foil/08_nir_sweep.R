#!/usr/bin/env Rscript
# Calibration foil, step 8: the nir sweep at the map-defined r.
#
# Holds the recombination knob at the map-defined r = 2*L/(100*N) on the native v5
# map (the expected per-adjacent-marker recombination fraction; from 03's json) (which
# the r sweep showed is in a flat safe basin) and sweeps the emission
# non-informative rate nir, scoring the real 24-line GBS calls against the chip
# calls with THREE metrics:
#   - GBS-vs-chip calls mismatch   (Holland's ORIGINAL objective: GBS call != chip call)
#   - donor-fragment Dice
#   - donor-fragment-size KS
# Expectation (reproducing Holland's grid selection): mismatch is minimized near
# nir = 0.9. The grid is extended past 0.9 to test whether 0.9 is a true optimum
# or just where Holland's grid stopped. The founder-genotype nir ~= 0.59 is marked
# (f0 measured on the NAM-founder chip genotypes; not a "biological" prescription).
#
#   Rscript scripts/nnil_foil/08_nir_sweep.R
# Output: data/nnil_foil/nir_sweep.csv + agent/nnil_foil_nir_sweep.png

suppressMessages({
  library(nilHMM)
  library(BEDMatrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
})
root <- here::here()
for (f in list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)) source(f)
source(file.path(root, "scripts/logging.R"))
FOIL <- file.path(root, "data/nnil_foil")
EQUIV <- file.path(root, "data/nnil_equiv")

# ---- chip calls + GBS input (same load as 03_chip_calibrate.R) --------------
xw <- fread(file.path(FOIL, "markers_v5.tsv"))
setkey(xw, marker)
chip <- fread(file.path(FOIL, "chip_truth_projected.csv"))
gbs_lines <- readLines(file.path(EQUIV, "lines.csv"))
both <- intersect(chip$Line, gbs_lines)
stopifnot(length(both) == 24)
chip <- chip[Line %in% both]
setorder(chip, Line)
chip_markers <- setdiff(names(chip), "Line")
chip_markers <- chip_markers[chip_markers %in% xw$marker]
cm_info <- xw[chip_markers, .(marker, chr, pos = pos_v5)]
setorder(cm_info, chr, pos)
grid_eval <- cm_info[, .(chr, pos)]

states_to_segments <- function(mat, lines, marker_pos) {
  mp <- marker_pos[match(colnames(mat), marker)]
  out <- vector("list", nrow(mat))
  for (i in seq_len(nrow(mat))) {
    dt <- data.table(chr = mp$chr, pos = mp$pos, state = as.integer(mat[i, ]))[!is.na(state) & state != 3L]
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
tr_sizes <- donor_block_sizes(tr)

geno <- BEDMatrix(file.path(EQUIV, "geno.bed"))
md <- fread(file.path(EQUIV, "markers.csv"))
row_idx <- match(both, gbs_lines)
xw_by_v4 <- xw[match(md$marker, marker_v4)]
keep_col <- which(!is.na(xw_by_v4$marker))
g_raw <- geno[row_idx, keep_col, drop = FALSE]
v5_chr <- xw_by_v4$chr[keep_col]
v5_pos <- xw_by_v4$pos_v5[keep_col]
ord <- order(v5_chr, v5_pos)
g_raw <- g_raw[, ord, drop = FALSE]
v5_chr <- v5_chr[ord]
v5_pos <- v5_pos[ord]
data <- rbindlist(lapply(seq_along(both), function(i) {
  g <- as.integer(g_raw[i, ])
  g[is.na(g)] <- 3L
  data.table(name = both[i], chr = v5_chr, pos = v5_pos, g = g)
}))
setorder(data, name, chr, pos)

# ---- sweep nir at the map-defined r (native v5 map); other params = Holland's
hp <- fromJSON(file.path(EQUIV, "params.json"))
map_r <- fromJSON(file.path(FOIL, "chip_calib.json"))$map_r # single source (03_chip_calibrate.R)
nir_grid <- c(0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.594, 0.7, 0.8, 0.9, 0.95, 0.99)
log_info("nir sweep at map r=%.3e on %d NILs vs chip (%d shared markers)", map_r, length(both), length(chip_markers))

sweep <- rbindlist(lapply(nir_grid, function(v) {
  called <- as.data.table(call_ancestry(
    data = data, caller = "nnil", rrate = map_r,
    germ = hp$germ, gert = hp$gert, p = hp$p, nir = v, mr = hp$mr, f_1 = hp$f_1, f_2 = hp$f_2
  ))
  mf <- marker_dice(called, tr, grid_eval)
  ff <- donor_fragment_dice(called, tr)
  ks <- fragment_size_ks(donor_block_sizes(called), tr_sizes)
  data.table(
    nir = v, marker_mismatch = 1 - mf$accuracy,
    donor_frag_dice = ff$dice, frag_ks = ks, donor_marker_dice = mf$per_class[class == "donor(>0)"]$dice
  )
}))
fwrite(sweep, file.path(FOIL, "nir_sweep.csv"))
nir_mm <- sweep$nir[which.min(sweep$marker_mismatch)]
nir_fd <- sweep$nir[which.max(sweep$donor_frag_dice)]
nir_ks <- sweep$nir[which.min(sweep$frag_ks)]
log_info("best nir | mismatch-min=%.3f | fragDice-max=%.3f | KS-min=%.3f", nir_mm, nir_fd, nir_ks)
print(sweep)

# ---- fragment-size ECDF at selected nir (mirror of the r-based fragsize view) ----
# Called donor-block sizes on the 24-line GBS at the map r, for a ladder of nir,
# vs the chip-call block sizes. Feeds the notebook's Section-5 mirror figure.
nir_ecdf_grid <- c(0.1, 0.3, 0.594, 0.8, 0.9)
ecdf_nir <- rbindlist(lapply(nir_ecdf_grid, function(v) {
  called <- as.data.table(call_ancestry(
    data = data, caller = "nnil", rrate = map_r,
    germ = hp$germ, gert = hp$gert, p = hp$p, nir = v, mr = hp$mr, f_1 = hp$f_1, f_2 = hp$f_2
  ))
  data.table(size_mb = donor_block_sizes(called), series = sprintf("GBS nir %.2g", v))
}))
ecdf_nir <- rbind(ecdf_nir, data.table(size_mb = tr_sizes, series = "Chip calls"))
fwrite(ecdf_nir, file.path(FOIL, "nir_fragsize_ecdf.csv"))
log_info("wrote nir_fragsize_ecdf.csv (%d selected nir + chip calls)", length(nir_ecdf_grid))

# ---- figure: three metrics vs nir -------------------------------------------
long <- rbind(
  data.table(nir = sweep$nir, val = sweep$marker_mismatch, metric = "GBS-vs-chip calls mismatch"),
  data.table(nir = sweep$nir, val = sweep$donor_frag_dice, metric = "Donor-fragment DSC"),
  data.table(nir = sweep$nir, val = sweep$frag_ks, metric = "Fragment-size KS")
)
long[, metric := factor(metric,
  levels = c("GBS-vs-chip calls mismatch", "Donor-fragment DSC", "Fragment-size KS"),
  labels = c(
    "GBS-vs-chip calls mismatch (lower better)", "Donor-fragment DSC (higher better)",
    "Fragment-size KS (lower better)"
  )
)]
fig <- ggplot(long, aes(nir, val)) +
  geom_vline(xintercept = 0.9, linetype = "dashed", colour = "grey30", linewidth = 0.4) +
  geom_vline(xintercept = 0.594, linetype = "dotted", colour = "#009E73", linewidth = 0.5) +
  geom_line(linewidth = 0.6, colour = "#D55E00") +
  geom_point(size = 1.3, colour = "#D55E00") +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  labs(
    x = expression("non-informative rate " * italic(nir) * "   (dashed = grid-tuned 0.9, dotted = founder nir 0.59)"),
    y = NULL, title = expression("nir sweep at the map-defined " * italic(r)),
    caption = "grid-tuned nir: value selected by Holland's grid search (marker mismatch). founder nir: f0 on NAM-founder chip genotypes."
  ) +
  theme_classic(base_size = 9) +
  theme(
    text = element_text(family = "sans"), plot.title = element_text(size = 9),
    plot.caption = element_text(size = 6.5, hjust = 0)
  )
OUT <- file.path(root, "agent/nnil_foil_nir_sweep.png")
ggsave(OUT, fig, width = 190, height = 62, units = "mm", dpi = 300)
cat(sprintf("wrote %s\n", OUT))
