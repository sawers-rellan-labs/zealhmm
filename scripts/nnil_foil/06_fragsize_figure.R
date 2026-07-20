#!/usr/bin/env Rscript
# Calibration foil, step 6b: KS (fragment-size) view.
#
# The Dice objective (reciprocal-overlap block matching) and the KS objective
# (match of the donor-block-SIZE distribution) pick different rrate. This figure
# shows both:
#   A  KS(called block sizes, truth block sizes) vs rrate, with the KS optimum,
#      the Dice optimum, and Holland's avg_r marked.
#   B  the donor-fragment-size ECDFs the KS compares: simcross truth vs the nnil
#      calls at avg_r, the Dice optimum, and the KS optimum.
#
# Regenerates the SAME 300 calibration NILs as 04_sim_calibrate.R (identical seed
# and RNG order: for i in 1..300, .simulate_dosage() then draw_obs()), so the
# block sizes match the sweep's ks_fragsize column exactly. Reads sim_rrate_sweep.csv.
#   Rscript scripts/nnil_foil/06_fragsize_figure.R
# Output: agent/nnil_foil_fragsize.png (exploratory; promote to figures/ if kept).

suppressMessages({
  library(nilHMM)
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

DESIGN <- "BC5S2"
NIR_FOUNDER <- 0.594
N_CAL <- 300L
M_INT <- 10L
set.seed(1L) # MUST match 04_sim_calibrate.R

xw <- fread(file.path(FOIL, "markers_v5.tsv"))
markers <- data.table(chr = as.integer(xw$chr), bp = as.integer(xw$pos_v5), cm = xw$cm)[order(chr, bp)]
cmlen <- markers[, .(L = max(cm)), by = chr][order(chr)]
pd <- parse_design(DESIGN)
bpd <- .bcsft_pedigree(pd$n_bc, pd$n_self)
hp <- fromJSON(file.path(EQUIV, "params.json"))
gt_emimat <- function(germ, gert, p, mr, nir) {
  matrix(c(
    (1 - germ) * (1 - mr), p * germ * (1 - mr), (1 - p) * germ * (1 - mr), mr,
    (((1 - nir) * 0.5 * gert) + nir * (1 - germ)) * (1 - mr), (((1 - nir) * (1 - gert)) + (nir * germ * p)) * (1 - mr),
    (((1 - nir) * 0.5 * gert) + nir * germ * (1 - p)) * (1 - mr), mr,
    ((1 - nir) * germ * (1 - p) + (nir * (1 - germ))) * (1 - mr), germ * p * (1 - mr),
    ((1 - nir) * (1 - germ) + (nir * germ * (1 - p))) * (1 - mr), mr
  ), nrow = 3, byrow = TRUE)
}
emimat <- gt_emimat(hp$germ, hp$gert, hp$p, hp$mr, NIR_FOUNDER)
draw_obs <- function(s) {
  g <- integer(length(s))
  for (v in 0:2) {
    ix <- which(s == v)
    if (length(ix)) g[ix] <- sample.int(4L, length(ix), TRUE, emimat[v + 1L, ]) - 1L
  }
  g
}

# regenerate the 300 cal NILs (same RNG order as step 5's first 300)
nms <- sprintf("sim%04d", seq_len(N_CAL))
tr_l <- data_l <- vector("list", N_CAL)
for (i in seq_len(N_CAL)) {
  dosage <- .simulate_dosage(bpd$ped, cmlen, markers, m = M_INT, p = 0, nil_id = bpd$nil_id)
  tr_l[[i]] <- .truth_segments(markers, dosage, nms[i])
  data_l[[i]] <- data.table(name = nms[i], chr = markers$chr, pos = markers$bp, g = draw_obs(dosage))
}
truth <- rbindlist(tr_l)
data <- rbindlist(data_l)
setorder(data, name, chr, pos)
log_info("regenerated %d cal NILs; truth donor blocks = %d", N_CAL, nrow(.donor_blocks(truth)))

sweep <- fread(file.path(FOIL, "sim_rrate_sweep.csv"))
map_r <- fromJSON(file.path(FOIL, "chip_calib.json"))$map_r
r_dice <- sweep$rrate[which.max(sweep$donor_frag_dice)]
r_ks <- sweep$rrate[which.min(sweep$ks_fragsize)]
log_info("DSC-opt rrate=%.3e | KS-opt rrate=%.3e | map r=%.3e", r_dice, r_ks, map_r)

call_sizes <- function(v) {
  seg <- as.data.table(call_ancestry(
    data = data, caller = "nnil", rrate = v,
    germ = hp$germ, gert = hp$gert, p = hp$p, nir = hp$nir, mr = hp$mr, f_1 = hp$f_1, f_2 = hp$f_2
  ))
  donor_block_sizes(seg)
}
# the marker-mismatch-optimal r (Holland's objective) over the WIDE chip sweep:
# mismatch is monotone in r, so its optimum runs to the low edge -> under-segments.
chip_sw <- fread(file.path(FOIL, "chip_rrate_sweep.csv"))
r_mm <- chip_sw$rrate[which.min(chip_sw$holland_mismatch)]
# include an over-fragmented rrate so Panel B shows what failure looks like
r_over <- sweep$rrate[which.min(abs(sweep$rrate - 2.2e-2))]
picks <- data.table(
  rrate = c(map_r, r_mm, r_dice, r_ks, r_over),
  lab = c(
    sprintf("map r (%.1e)", map_r), sprintf("mismatch r* (%.1e)", r_mm),
    sprintf("DSC r* (%.1e)", r_dice), sprintf("KS r* (%.1e)", r_ks),
    sprintf("over-frag (%.1e)", r_over)
  )
)
log_info("marker-mismatch-optimal r (wide chip sweep) = %.2e", r_mm)
ecdf_dt <- rbindlist(lapply(seq_len(nrow(picks)), function(i) {
  data.table(size_mb = call_sizes(picks$rrate[i]), series = picks$lab[i])
}))
ecdf_dt <- rbind(ecdf_dt, data.table(size_mb = donor_block_sizes(truth), series = "Truth (simcross)"))
lev <- c("Truth (simcross)", picks$lab)
ecdf_dt[, series := factor(series, levels = lev)]
fwrite(ecdf_dt, file.path(FOIL, "fragsize_ecdf.csv")) # for the notebook

col_sim <- "#0072B2"
col_holl <- "grey30"
col_ks <- "#009E73"
base <- theme_classic(base_size = 9) + theme(
  text = element_text(family = "sans"),
  plot.tag = element_text(face = "bold"), legend.position = "top", legend.title = element_blank(),
  legend.key.height = unit(0.4, "lines")
)
lx <- scale_x_log10(breaks = 10^(-6:-1), labels = c("1e-6", "1e-5", "1e-4", "1e-3", "1e-2", "1e-1"))

col_mm <- "#CC79A7" # mismatch r* (matches its ECDF colour in panel B)
pA <- ggplot(sweep, aes(rrate, ks_fragsize)) +
  geom_vline(xintercept = r_mm, linetype = "twodash", colour = col_mm, linewidth = 0.5) +
  geom_vline(xintercept = map_r, linetype = "dashed", colour = col_holl, linewidth = 0.4) +
  geom_vline(xintercept = r_dice, linetype = "dotted", colour = col_sim, linewidth = 0.5) +
  geom_vline(xintercept = r_ks, linetype = "dotdash", colour = col_ks, linewidth = 0.5) +
  geom_line(linewidth = 0.6, colour = "grey20") +
  geom_point(size = 1, colour = "grey20") +
  lx +
  labs(
    x = expression("intermarker recombination fraction, " * italic(r)),
    y = "KS(called, truth) fragment size", tag = "A"
  ) +
  base +
  annotate("text",
    x = r_mm, y = max(sweep$ks_fragsize) * 0.9, label = "mismatch~italic(r)^'*'",
    parse = TRUE, angle = 90, hjust = 1, vjust = 1.2, size = 2.5, colour = col_mm
  ) +
  annotate("text",
    x = map_r, y = max(sweep$ks_fragsize) * 0.9, label = "map~italic(r)",
    parse = TRUE, angle = 90, hjust = 1, vjust = -0.3, size = 2.5, colour = col_holl
  ) +
  annotate("text",
    x = r_dice, y = max(sweep$ks_fragsize) * 0.6, label = "DSC~italic(r)^'*'", parse = TRUE,
    angle = 90, hjust = 1, vjust = -0.3, size = 2.5, colour = col_sim
  ) +
  annotate("text",
    x = r_ks, y = max(sweep$ks_fragsize) * 0.6, label = "KS~italic(r)^'*'", parse = TRUE,
    angle = 90, hjust = 1, vjust = 1.2, size = 2.5, colour = col_ks
  )

pB <- ggplot(ecdf_dt, aes(size_mb, colour = series)) +
  stat_ecdf(linewidth = 0.6) +
  scale_x_log10() +
  scale_colour_manual(values = setNames(
    c("black", col_holl, "#CC79A7", col_sim, col_ks, "#D55E00"), lev
  )) +
  labs(x = "donor fragment size (Mb)", y = "ECDF", tag = "B") +
  base +
  guides(colour = guide_legend(nrow = 3))

fig <- pA + pB + plot_layout(widths = c(1, 1.1))
OUT <- file.path(root, "agent/nnil_foil_fragsize.png")
ggsave(OUT, fig, width = 180, height = 82, units = "mm", dpi = 300)
cat(sprintf("wrote %s\n", OUT))
