# nnil sim-calibration figure: nnil non-informative-rate calibration + introgression-size recovery.
#
# Panels (2x2):
#   A  nnil calibration of nir by MARKER MISMATCH (Holland's File_S16 criterion):
#      mismatch = 1 - marker accuracy. sim = caller vs latent simcross truth;
#      chip = real GBS caller calls vs chip calls. Dotted verticals = the donors'
#      actual non-informative rate f0 (sim-realized, real founders).
#   B  Simulated GBS: sim-calibrated nnil calls vs simulated ancestry (ECDF).
#   C  Real GBS (24 both-platform NILs): chip- and sim-calibrated nnil vs chip calls.
#   D  QQ of introgression size, chip- vs sim-calibrated nnil, on the real 24 NILs.
# rtiger (count caller) has no real-data leg (no counts), so it is excluded here.
#
# Font: ONE base size (BASE) drives all panels; no per-element size overrides.
# Speed: the ~25-decode compute is CACHED to fig_nnil_sim_calibration_cache.rds. Re-renders reuse it;
#        pass NNIL_SIMCAL_RECOMPUTE=1 to recompute after the sim or parameters change.
#
#   Rscript scripts/fig_nnil_sim_calibration.R
#   NNIL_SIMCAL_RECOMPUTE=1 Rscript scripts/fig_nnil_sim_calibration.R   # force recompute

suppressMessages({
  library(data.table)
  library(here)
  library(ggplot2)
  library(ggtext)
  library(patchwork)
  library(nilHMM)
})
source(here::here("R/metrics.R"))

SIMDIR <- here::here("results/sim/nested_nnil")
DESIGN <- "BC5S4"
NIR_GRID <- c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
BASE <- 20 # single font base size for every panel
lg <- function(...) message(sprintf(...))

# ============================ compute (cached) ==============================
CACHE <- file.path(SIMDIR, "fig_nnil_sim_calibration_cache.rds")
if (file.exists(CACHE) && Sys.getenv("NNIL_SIMCAL_RECOMPUTE") == "") {
  lg("[fig] reusing cached compute (%s); NNIL_SIMCAL_RECOMPUTE=1 to recompute", basename(CACHE))
  D <- readRDS(CACHE)
} else {
  mv5 <- fread(here::here("data/nnil_foil/markers_v5.tsv"))
  v4_to_v5 <- setNames(mv5$marker, mv5$marker_v4)
  v5_to_v4 <- setNames(mv5$marker_v4, mv5$marker)

  # sim: counts -> ML hard-calls; latent truth
  sim <- readRDS(file.path(SIMDIR, sprintf("nested_nnil_%s_nnil_gbs.rds", tolower(DESIGN))))
  sim_truth <- fread(file.path(SIMDIR, sprintf("nested_nnil_%s_truth_segments.csv", tolower(DESIGN))))
  M <- nrow(sim$grid)
  sim_mk <- paste0("S", sim$grid$chr, "_", sim$grid$pos)
  sim_long <- data.table(
    name = rep(sim$names, each = M),
    chr = rep(as.integer(sim$grid$chr), times = length(sim$names)),
    pos = rep(as.integer(sim$grid$pos), times = length(sim$names)),
    n_ref = as.integer(sim$n_ref), n_alt = as.integer(sim$n_alt)
  )
  sim_long[, g := nilHMM::call_gt(n_ref, n_alt, prior = "flat", error = 0.01, return = "call")]
  sim_long[is.na(g), g := 3L]
  common <- intersect(sim_mk, mv5$marker)
  total_cM <- sum(mv5[marker %in% common, max(cm, na.rm = TRUE), by = chr]$V1)
  RRATE <- 2 * total_cM / (100 * length(common))
  sim_truth_sz <- donor_block_sizes(sim_truth)

  # real 24-line GBS (v4->v5 remap) + chip truth
  chip <- fread(here::here("data/nnil_foil/chip_truth_projected.csv"))
  geno <- fread(here::here("data/nnil_equiv/geno_recoded.csv"))
  idcol <- names(geno)[1]
  lines24 <- intersect(chip$Line, geno[[idcol]])
  mk_common <- Reduce(intersect, list(
    sim_mk, setdiff(names(chip), "Line"),
    unname(v4_to_v5[setdiff(names(geno), idcol)])
  ))
  lg(
    "[fig] %d sim NILs; chip calibration on %d both-platform NILs x %d markers; rrate %.2e",
    length(sim$names), length(lines24), length(mk_common), RRATE
  )

  chip_long <- melt(chip[Line %in% lines24, c("Line", mk_common), with = FALSE],
    id.vars = "Line", variable.name = "marker", value.name = "state"
  )
  chip_long[, `:=`(
    name = Line, chr = as.integer(sub("S(\\d+)_.*", "\\1", marker)),
    pos = as.integer(sub("S\\d+_", "", marker))
  )]
  chip_seg <- chip_long[order(name, chr, pos),
    .(start_bp = pos, end_bp = pos, state = as.integer(state)),
    by = .(name, chr)
  ]
  chip_truth_sz <- donor_block_sizes(chip_seg)

  cols_v4 <- v5_to_v4[mk_common]
  gsub24 <- geno[get(idcol) %in% lines24, c(idcol, cols_v4), with = FALSE]
  setnames(gsub24, c(idcol, cols_v4), c("name", mk_common))
  real_long <- melt(gsub24, id.vars = "name", variable.name = "marker", value.name = "g")
  real_long[, `:=`(
    chr = as.integer(sub("S(\\d+)_.*", "\\1", marker)),
    pos = as.integer(sub("S\\d+_", "", marker)), g = as.integer(g)
  )]
  real_long <- real_long[order(name, chr, pos)]

  seg_nnil <- function(long_g, nir) {
    as.data.table(call_ancestry(
      long_g[, .(name, chr, pos, g)],
      caller = "nnil", design = DESIGN, rrate = RRATE, nir = nir
    ))
  }
  grid_sim <- data.table(chr = as.integer(sim$grid$chr), pos = as.integer(sim$grid$pos))
  grid_chip <- unique(chip_seg[, .(chr, pos = start_bp)])
  mm_sim <- function(seg) 1 - marker_dice(seg, sim_truth, grid_sim)$accuracy
  mm_chip <- function(seg) 1 - marker_dice(seg, chip_seg, grid_chip)$accuracy

  # sweep nir by marker mismatch (the expensive 22-decode step). Reuse the persisted
  # CSV if present; NNIL_SIMCAL_RECOMPUTE=1 forces a fresh sweep.
  swp_csv <- file.path(SIMDIR, "fig_nnil_sim_calibration_sweep.csv")
  if (file.exists(swp_csv) && Sys.getenv("NNIL_SIMCAL_RECOMPUTE") == "") {
    lg("[fig] reusing persisted sweep (%s)", basename(swp_csv))
    sweep <- fread(swp_csv)
  } else {
    sweep <- rbindlist(list(
      rbindlist(lapply(NIR_GRID, function(v) {
        data.table(
          truth = "sim", param = v,
          mismatch = mm_sim(seg_nnil(sim_long, v))
        )
      })),
      rbindlist(lapply(NIR_GRID, function(v) {
        data.table(
          truth = "chip", param = v,
          mismatch = mm_chip(seg_nnil(real_long, v))
        )
      }))
    ))
    fwrite(sweep, swp_csv)
  }
  opt <- sweep[, .SD[which.min(mismatch)], by = truth]
  nir_sim <- opt[truth == "sim", param]
  nir_chip <- opt[truth == "chip", param]

  # caller runs at the two optima
  sz_nnil_sim <- donor_block_sizes(seg_nnil(sim_long, nir_sim)) # nnil on sim @ sim-opt
  sz_chipcal <- donor_block_sizes(seg_nnil(real_long, nir_chip)) # nnil on real GBS @ chip-opt
  sz_simcal <- donor_block_sizes(seg_nnil(real_long, nir_sim)) # nnil on real GBS @ sim-opt

  # donor non-informative rate f0 = fraction of markers where donor allele = B73/REF
  nir_sim_donor <- mean(sim$founder_f0$f0)
  fnd <- as.matrix(fread(here::here("data/nnil_foil/founders_v5.csv")), rownames = 1)
  nir_real_donor <- mean(apply(fnd, 1, function(g) mean(g[g != 3] == 0)))

  D <- list(
    sweep = sweep, opt = opt, nir_sim = nir_sim, nir_chip = nir_chip,
    sim_truth_sz = sim_truth_sz, chip_truth_sz = chip_truth_sz,
    sz_nnil_sim = sz_nnil_sim, sz_chipcal = sz_chipcal, sz_simcal = sz_simcal,
    nir_sim_donor = nir_sim_donor, nir_real_donor = nir_real_donor
  )
  saveRDS(D, CACHE)
  fwrite(opt, file.path(SIMDIR, "fig_nnil_sim_calibration_optima.csv"))
}
list2env(D, environment())
lg(
  "[fig] mismatch optima: nir_sim=%.2f nir_chip=%.2f | donor f0 sim=%.3f real=%.3f",
  nir_sim, nir_chip, nir_sim_donor, nir_real_donor
)

# ================================ plot ======================================
ANN <- BASE * 0.8 / .pt # match the legend category-label size (theme legend.text = rel(0.8) of BASE)

L_sim_anc <- sprintf("simulated ancestry\n%s truth", DESIGN)
L_nnil_on_sim <- sprintf("sim-calibrated nnil\n(nir=%.2f)", nir_sim)
L_chip_calls <- "chip calls\n(Zhong 2025)"
L_nnil_chip <- sprintf("chip-calibrated nnil\n(nir=%.2f)", nir_chip)
# colour = data domain: blue = simulation, orange = chip-calibrated, black = chip-calls truth
pal <- setNames(
  c("#0072B2", "#0072B2", "#D55E00", "black"),
  c(L_sim_anc, L_nnil_on_sim, L_nnil_chip, L_chip_calls)
)
# ground truth = dotted, callers = solid (panels B and C)
lty_pal <- setNames(
  c("dotted", "solid", "solid", "dotted"),
  c(L_sim_anc, L_nnil_on_sim, L_nnil_chip, L_chip_calls)
)

ecdf_B <- rbindlist(list(
  data.table(size_mb = sim_truth_sz, series = L_sim_anc),
  data.table(size_mb = sz_nnil_sim, series = L_nnil_on_sim)
))
ecdf_C <- rbindlist(list(
  data.table(size_mb = chip_truth_sz, series = L_chip_calls),
  data.table(size_mb = sz_chipcal, series = L_nnil_chip),
  data.table(size_mb = sz_simcal, series = L_nnil_on_sim)
))
# ground truth first so it sits at the top of each legend
ecdf_B[, series := factor(series, levels = c(L_sim_anc, L_nnil_on_sim))]
ecdf_C[, series := factor(series, levels = c(L_chip_calls, L_nnil_chip, L_nnil_on_sim))]
qs <- ppoints(200)
qq <- data.table(chip_cal = quantile(sz_chipcal, qs), sim_cal = quantile(sz_simcal, qs))

.szC <- c(chip_truth_sz, sz_chipcal, sz_simcal)
xlim_mb <- c(0.1, max(.szC[is.finite(.szC) & .szC > 0])) # window shared by B/C/D
MB_BREAKS <- c(0.1, 1, 10, 100)
MB_LABELS <- c("0.1", "1", "10", "100")

# A: calibration sweep + donor-nir reference lines (both labels right of the blue line)
p_A <- ggplot(sweep, aes(param, mismatch, colour = truth)) +
  geom_line() +
  geom_point(size = 1) +
  geom_point(
    data = opt, size = 3, shape = 21, fill = "white",
    colour = ifelse(opt$truth == "sim", "#0072B2", "#D55E00") # dot = the caller it yields in C/D
  ) +
  geom_vline(xintercept = nir_sim_donor, linetype = "dashed", colour = "#0072B2", linewidth = 0.8, alpha = 0.35) +
  geom_vline(xintercept = nir_real_donor, linetype = "dashed", colour = "black", linewidth = 0.8, alpha = 0.35) +
  annotate("text",
    x = nir_sim_donor + 0.035, y = 0.01, label = "sim donor",
    colour = "#0072B2", angle = 90, hjust = 0.5, vjust = 0.5, size = ANN, alpha = 0.5
  ) +
  annotate("text",
    x = nir_real_donor - 0.052, y = 0.01, label = "chip donor",
    colour = "black", angle = 90, hjust = 0.5, vjust = 0.5, size = ANN, alpha = 0.5
  ) +
  scale_colour_manual(
    values = c(sim = "#0072B2", chip = "#D55E00"),
    labels = c(sim = "sim GBS vs sim ancestry", chip = "real GBS vs chip calls"),
    name = "calibration"
  ) +
  scale_x_continuous(
    breaks = seq(0, 1, 0.2), minor_breaks = seq(0, 0.9, 0.1),
    guide = guide_axis(minor.ticks = TRUE)
  ) +
  labs(
    x = expression("non-informative rate " * italic(nir)),
    y = "marker mismatch rate", title = "nnil calibration of the\nnon-informative rate"
  ) +
  theme_bw(base_size = BASE) +
  theme(
    aspect.ratio = 1, legend.position = c(0.98, 0.98), legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", colour = NA),
    legend.key = element_rect(fill = "white", colour = NA)
  )

mk_ecdf <- function(d, ttl) {
  ggplot(d, aes(size_mb, colour = series, linetype = series)) +
    stat_ecdf(linewidth = 1) +
    scale_x_log10(breaks = MB_BREAKS, labels = MB_LABELS, limits = xlim_mb, oob = scales::oob_keep) +
    scale_colour_manual(values = pal, name = NULL) +
    scale_linetype_manual(values = lty_pal, name = NULL) +
    labs(x = "introgression size (Mb)", y = "ECDF", title = ttl) +
    theme_bw(base_size = BASE) +
    theme(
      aspect.ratio = 1, legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.key = element_rect(fill = "white", colour = NA)
    )
}

p_C <- mk_ecdf(ecdf_C, "Real GBS, 24 NILs:\nnnil calls vs chip calls")

# Two-sample KS on panel C: chip-calibrated vs sim-calibrated nnil introgression
# sizes (same real 24 NILs). Non-significant p => the two calibrations agree.
ks_C <- suppressWarnings(ks.test(sz_chipcal, sz_simcal))
ks_p <- ks_C$p.value
ks_txt <- if (ks_p >= 0.01) sprintf("%.2f", ks_p) else sprintf("%.0e", ks_p)
lg("[fig] panel C KS (chip-cal vs sim-cal nnil): D=%.3f p=%.3g", ks_C$statistic, ks_p)
# chip/sim coloured to match the panel-C legend; p-value symbol in italics (ggtext)
ks_rich <- sprintf(
  "calibration<br><span style='color:%s'>chip</span> vs <span style='color:%s'>sim</span><br>*p* = %s",
  pal[[L_nnil_chip]], pal[[L_nnil_on_sim]], ks_txt
)
p_C <- p_C + geom_richtext(
  data = data.frame(x = 11, y = 0.15, label = ks_rich),
  aes(x, y, label = label), inherit.aes = FALSE,
  hjust = 0, vjust = 0.5, size = ANN, lineheight = 0.9, colour = "grey20",
  fill = "white", label.color = NA, label.r = grid::unit(0, "pt"),
  label.padding = grid::unit(2, "pt")
)

p_B <- mk_ecdf(ecdf_B, "Simulated GBS:\nnnil calls vs simulated ancestry")

p_D <- ggplot(qq, aes(chip_cal, sim_cal)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey50") +
  geom_point(size = 1.3, colour = "#0072B2") +
  scale_x_log10(limits = xlim_mb, breaks = MB_BREAKS, labels = MB_LABELS, oob = scales::oob_keep) +
  scale_y_log10(limits = xlim_mb, breaks = MB_BREAKS, labels = MB_LABELS, oob = scales::oob_keep) +
  labs(
    x = "introgression size (Mb)\nfrom chip calibrated calls",
    y = "introgression size (Mb)\nfrom sim calibrated calls",
    title = "QQ plot of introgression size\nby nnil calibration method"
  ) +
  theme_bw(base_size = BASE) +
  theme(aspect.ratio = 1)

fig <- (p_A | p_B) / (p_C | p_D) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 25, face = "bold"),
    plot.tag.location = "plot", plot.tag.position = "topleft"
  )
ggsave(file.path(SIMDIR, "fig_nnil_sim_calibration.png"), fig, width = 13, height = 13, dpi = 150)
lg("[fig] wrote fig_nnil_sim_calibration.png")
