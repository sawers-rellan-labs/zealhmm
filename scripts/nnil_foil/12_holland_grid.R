#!/usr/bin/env Rscript
# Calibration foil, step 12: reproduce Holland's FULL grid (File_S04, 945 configs
# nir x germ x gert x p x r) on the real nNIL 24-line GBS vs chip, through the
# unified caller_grid() kernel. If our caller reproduces Holland's finding, the
# GBS-vs-chip mismatch is driven by nir (his variance decomposition: ~97.7%), with
# the optimum near his grid-tuned nir.
#
# Runs His EXACT grid levels (read from File_S04). Chunked over the emission grid:
# bounded memory (each batch's calls are scored to a mismatch row then discarded;
# only the 945-row summary is retained) + a per-batch ETA. Cores are capped so peak
# RAM stays well under 16 GB (this dataset is tiny; the cap matters more for the
# ZEAL grid). Uses the feat/viterbi-sweep-kernel nilhmm (caller_grid) via load_all.
#
#   Rscript scripts/nnil_foil/12_holland_grid.R
# Output: data/nnil_foil/holland_grid_reproduction.csv (our mismatch + his mismatchMean)
#         data/nnil_foil/holland_grid_sensitivity.csv  (variance decomposition, ours vs his)
#         agent/nnil_foil_holland_grid.png

suppressMessages({
  library(devtools)
  library(BEDMatrix)
  library(data.table)
  library(jsonlite)
  library(ggplot2)
  library(parallel)
})
devtools::load_all(path.expand("~/repos/nilhmm"), quiet = TRUE) # caller_grid (nilHMM)
root <- here::here()
for (f in list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)) source(f) # scoring helpers (R/metrics.R)
source(file.path(root, "scripts/logging.R"))

MEM_CAP_GB <- 16L
N_CORES <- min(detectCores() - 2L, 8L) # capped; per-worker footprint is tiny here
log_info("start Holland-grid reproduction | cores=%d (RAM cap %d GB; peak stays far below on this 24-line set)", N_CORES, MEM_CAP_GB)

FOIL <- file.path(root, "data/nnil_foil")
EQUIV <- file.path(root, "data/nnil_equiv")

# ---- data: identical load to scripts/nnil_foil/08_nir_sweep.R -----------------
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
hp <- fromJSON(file.path(EQUIV, "params.json"))
log_info("loaded %d NILs vs chip | %d shared markers", length(both), length(chip_markers))

# ---- Holland's exact grid, read from File_S04 --------------------------------
f4 <- fread(file.path(root, "agent/nNIL/File_S04.nNIL_gbs_vs_chip_data_HMMgridSearch.csv"))
lv <- function(col) sort(unique(f4[[col]]))
eg <- as.data.table(expand.grid(nir = lv("nir"), germ = lv("germ"), gert = lv("gert"), p = lv("p")))
eg[, mr := hp$mr]
rr <- lv("r")
n_cfg <- nrow(eg) * length(rr)
log_info("Holland grid (File_S04 levels): %d emission combos x %d r = %d configs", nrow(eg), length(rr), n_cfg)

# score one batch of emission combos: caller_grid over the r grid, marker mismatch
# vs chip per config, then DISCARD the segments (keep only the mismatch rows).
score_batch <- function(eg_batch) {
  seg <- as.data.table(caller_grid(
    data,
    caller = "nnil", emission_grid = as.data.frame(eg_batch),
    rrate = rr, f_1 = hp$f_1, f_2 = hp$f_2, threads = N_CORES
  ))
  seg[, cfg := paste(nir, germ, gert, p, rrate, sep = "_")]
  cfgs <- unique(seg[, .(nir, germ, gert, p, r = rrate, cfg)])
  out <- rbindlist(lapply(seq_len(nrow(cfgs)), function(i) {
    called <- seg[cfg == cfgs$cfg[i]]
    mf <- marker_dice(called, tr, grid_eval)
    data.table(
      nir = cfgs$nir[i], germ = cfgs$germ[i], gert = cfgs$gert[i],
      p = cfgs$p[i], r = cfgs$r[i], mismatch = 1 - mf$accuracy
    )
  }))
  rm(seg)
  gc(FALSE)
  out
}

# ---- upfront ETA probe: one emission combo x the r grid ----------------------
tp <- Sys.time()
invisible(score_batch(eg[1]))
per_cfg <- as.numeric(difftime(Sys.time(), tp, units = "secs")) / length(rr)
log_info(
  "PROBE: %.2fs/config -> projected full grid ~%.1f min (%d configs, %d cores). Starting.",
  per_cfg, per_cfg * n_cfg / 60, n_cfg, N_CORES
)

# ---- chunked run with per-batch ETA ------------------------------------------
BATCH <- 21L
batches <- split(seq_len(nrow(eg)), (seq_len(nrow(eg)) - 1L) %/% BATCH)
B <- length(batches)
t0 <- Sys.time()
res <- vector("list", B)
for (b in seq_len(B)) {
  res[[b]] <- score_batch(eg[batches[[b]]])
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  log_info(
    ">>> batch %d/%d done | elapsed %.1f min | avg %.2f min/batch | ETA ~%.1f min remaining",
    b, B, el, el / b, (el / b) * (B - b)
  )
}
score <- rbindlist(res)
log_info("scored %d configs", nrow(score))

# ---- compare to Holland's File_S04 + variance decomposition ------------------
key_cols <- c("nir", "germ", "gert", "p", "r")
score[, rk := round(r, 10)]
f4[, rk := round(r, 10)]
cmp <- merge(score, f4[, .(nir, germ, gert, p, rk, holland = mismatchMean)],
  by = c("nir", "germ", "gert", "p", "rk")
)
log_info("merged %d/%d configs with File_S04", nrow(cmp), n_cfg)
fwrite(cmp[, .(nir, germ, gert, p, r, mismatch, holland)], file.path(FOIL, "holland_grid_reproduction.csv"))

pear <- cor(cmp$mismatch, cmp$holland)
spear <- cor(cmp$mismatch, cmp$holland, method = "spearman")
our_opt <- score[which.min(mismatch)]
his_opt <- f4[which.min(mismatchMean)]
log_info("our mismatch vs Holland's: Pearson=%.3f Spearman=%.3f", pear, spear)
log_info("OUR   optimum: mismatch=%.4f at nir=%g germ=%g gert=%g p=%g r=%.2e", our_opt$mismatch, our_opt$nir, our_opt$germ, our_opt$gert, our_opt$p, our_opt$r)
log_info("HOLLAND optimum: mismatch=%.4f at nir=%g germ=%g gert=%g p=%g r=%.2e", his_opt$mismatchMean, his_opt$nir, his_opt$germ, his_opt$gert, his_opt$p, his_opt$r)

var_decomp <- function(d, y) {
  dd <- copy(d)
  dd[, (key_cols) := lapply(.SD, factor), .SDcols = key_cols]
  tab <- as.data.table(summary(aov(reformulate(key_cols, y), data = dd))[[1]], keep.rownames = "term")
  tab[, term := trimws(term)]
  tab[term %in% key_cols, .(param = term, var_share = `Sum Sq` / sum(`Sum Sq`))]
}
sens_ours <- var_decomp(score, "mismatch")[, src := "ours (caller_grid)"]
sens_his <- var_decomp(f4, "mismatchMean")[, src := "Holland File_S04"]
sens <- rbind(sens_ours, sens_his)[order(src, -var_share)]
fwrite(dcast(sens, param ~ src, value.var = "var_share"), file.path(FOIL, "holland_grid_sensitivity.csv"))
log_info("variance shares (ours | Holland):")
print(dcast(sens, param ~ src, value.var = "var_share")[order(-`ours (caller_grid)`)])

# ---- figure: marginal mean mismatch vs each parameter, ours vs Holland -------
marg <- rbindlist(lapply(key_cols, function(pp) {
  o <- score[, .(mm = mean(mismatch)), by = c(pp)][, .(level = get(pp), mm, param = pp, src = "ours (caller_grid)")]
  h <- f4[, .(mm = mean(mismatchMean)), by = c(pp)][, .(level = get(pp), mm, param = pp, src = "Holland File_S04")]
  rbind(o, h)
}))
vs <- dcast(sens, param ~ src, value.var = "var_share")
marg[, param := factor(param,
  levels = key_cols,
  labels = sprintf("%s (ours %.0f%% var)", key_cols, 100 * vs[["ours (caller_grid)"]][match(key_cols, vs$param)])
)]
fig <- ggplot(marg, aes(level, mm, colour = src)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.2) +
  facet_wrap(~param, scales = "free_x", nrow = 1) +
  scale_x_log10() +
  scale_colour_manual(values = c("ours (caller_grid)" = "#0072B2", "Holland File_S04" = "#D55E00")) +
  labs(
    x = "parameter value (log scale)", y = "mean GBS-vs-chip mismatch", colour = NULL,
    title = sprintf("Holland's 945-config grid through caller_grid on the nNIL 24-line GBS vs chip (Pearson %.2f vs File_S04)", pear)
  ) +
  theme_classic(base_size = 9) +
  theme(
    text = element_text(family = "sans"), plot.title = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom"
  )
ggsave(file.path(root, "agent/nnil_foil_holland_grid.png"), fig, width = 190, height = 60, units = "mm", dpi = 300)
log_info("wrote holland_grid_reproduction.csv + holland_grid_sensitivity.csv + agent/nnil_foil_holland_grid.png")
