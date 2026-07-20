#!/usr/bin/env Rscript
# Sim-only nir calibration (self-consistency check).
#
# The sim is GENERATED at nir = 0.594 (founder-genotype non-informative rate,
# injected through Holland's emission on the simcross ancestry, exactly as
# 04_sim_calibrate.R does). Here we fix r at the map value and sweep the CALLER's
# nir against the DENSE simcross truth. Question: does the caller's nir* land near
# the generation value (~0.59), rather than the ~0.9 the real chip-supervised fit
# forces? If yes, the decomposition is clean (sim -> biology; the excess to 0.9 on
# real data = data quality). This is sim-only: no chip, no GBS.
#
#   Rscript scripts/nnil_foil/09_sim_nir_sweep.R
# Output: data/nnil_foil/sim_nir_sweep.csv + agent/nnil_foil_sim_nir_sweep.png

suppressMessages({
  library(nilHMM)
  library(data.table)
  library(ggplot2)
  library(jsonlite)
})
root <- here::here()
for (f in list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)) source(f)
source(file.path(root, "scripts/logging.R"))
FOIL <- file.path(root, "data/nnil_foil")
EQUIV <- file.path(root, "data/nnil_equiv")

DESIGN <- "BC5S2"
NIR_GEN <- 0.594 # generation non-informative rate (founder f0)
N_CAL <- 300L
M_INT <- 10L
set.seed(1L) # MUST match 04_sim_calibrate.R's first N_CAL lines

# ---- grid + map r (sim-only: computed from the markers, no chip dependency) --
xw <- fread(file.path(FOIL, "markers_v5.tsv"))
markers <- data.table(chr = as.integer(xw$chr), bp = as.integer(xw$pos_v5), cm = xw$cm)[order(chr, bp)]
cmlen <- markers[, .(L = max(cm)), by = chr][order(chr)]
map_r <- 2 * markers[, sum(tapply(cm, chr, function(x) max(x) - min(x)))] / (100 * nrow(markers))
pd <- parse_design(DESIGN)
bpd <- .bcsft_pedigree(pd$n_bc, pd$n_self)
hp <- fromJSON(file.path(EQUIV, "params.json"))

# emission for GENERATION at the founder nir (Holland's other error terms)
gt_emimat <- function(germ, gert, p, mr, nir) {
  matrix(c(
    (1 - germ) * (1 - mr), p * germ * (1 - mr), (1 - p) * germ * (1 - mr), mr,
    (((1 - nir) * 0.5 * gert) + nir * (1 - germ)) * (1 - mr), (((1 - nir) * (1 - gert)) + (nir * germ * p)) * (1 - mr),
    (((1 - nir) * 0.5 * gert) + nir * germ * (1 - p)) * (1 - mr), mr,
    ((1 - nir) * germ * (1 - p) + (nir * (1 - germ))) * (1 - mr), germ * p * (1 - mr),
    ((1 - nir) * (1 - germ) + (nir * germ * (1 - p))) * (1 - mr), mr
  ), nrow = 3, byrow = TRUE)
}
emimat <- gt_emimat(hp$germ, hp$gert, hp$p, hp$mr, NIR_GEN)
draw_obs <- function(s) {
  g <- integer(length(s))
  for (v in 0:2) {
    ix <- which(s == v)
    if (length(ix)) g[ix] <- sample.int(4L, length(ix), TRUE, emimat[v + 1L, ]) - 1L
  }
  g
}

# ---- regenerate the 300 cal NILs (same RNG order as 04's first 300) ---------
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
grid_eval <- markers[, .(chr, pos = bp)]
tr_sizes <- donor_block_sizes(truth)
log_info("regenerated %d cal NILs (gen nir=%.3f); sweeping caller nir at map r=%.3e", N_CAL, NIR_GEN, map_r)

# ---- sweep the CALLER's nir at fixed map r, score vs the simcross truth ------
nir_grid <- c(0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.594, 0.7, 0.8, 0.9, 0.95, 0.99)
t0 <- Sys.time()
sweep <- rbindlist(lapply(nir_grid, function(v) {
  called <- as.data.table(call_ancestry(
    data = data, caller = "nnil", rrate = map_r,
    germ = hp$germ, gert = hp$gert, p = hp$p, nir = v, mr = hp$mr, f_1 = hp$f_1, f_2 = hp$f_2
  ))
  mf <- marker_dice(called, truth, grid_eval)
  ff <- donor_fragment_dice(called, truth)
  r <- data.table(
    nir = v, marker_mismatch = 1 - mf$accuracy,
    donor_frag_dice = ff$dice, frag_ks = fragment_size_ks(donor_block_sizes(called), tr_sizes),
    donor_marker_dice = mf$per_class[class == "donor(>0)"]$dice
  )
  log_info(
    "  caller nir=%.3f | frag_dice=%.3f mismatch=%.4f (%.0fs)", v, r$donor_frag_dice,
    r$marker_mismatch, as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )
  r
}))
fwrite(sweep, file.path(FOIL, "sim_nir_sweep.csv"))
nir_star_dice <- sweep$nir[which.max(sweep$donor_frag_dice)]
nir_star_mm <- sweep$nir[which.min(sweep$marker_mismatch)]
log_info(
  "SIM nir* | fragDice-max=%.3f | mismatch-min=%.3f  (generation nir=%.3f)",
  nir_star_dice, nir_star_mm, NIR_GEN
)

# ---- figure: 3 metrics vs caller nir, generation (0.594) + GBS-opt (0.9) marked
long <- rbind(
  data.table(nir = sweep$nir, val = sweep$marker_mismatch, metric = "Marker mismatch (lower better)"),
  data.table(nir = sweep$nir, val = sweep$donor_frag_dice, metric = "Donor-fragment DSC (higher better)"),
  data.table(nir = sweep$nir, val = sweep$frag_ks, metric = "Fragment-size KS (lower better)")
)
long[, metric := factor(metric, levels = c(
  "Marker mismatch (lower better)",
  "Donor-fragment DSC (higher better)", "Fragment-size KS (lower better)"
))]
fig <- ggplot(long, aes(nir, val)) +
  geom_vline(xintercept = NIR_GEN, linetype = "dotted", colour = "#009E73", linewidth = 0.6) +
  geom_vline(xintercept = 0.9, linetype = "dashed", colour = "grey30", linewidth = 0.4) +
  geom_line(linewidth = 0.6, colour = "#0072B2") +
  geom_point(size = 1.3, colour = "#0072B2") +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  labs(
    x = expression("caller " * italic(nir) * "   (dotted = generation 0.594, dashed = chip-supervised opt 0.9)"),
    y = NULL, title = "nir sweep on the SIMULATION (data generated at nir = 0.594), scored vs simcross truth"
  ) +
  theme_classic(base_size = 9) +
  theme(text = element_text(family = "sans"), plot.title = element_text(size = 9))
ggsave(file.path(root, "agent/nnil_foil_sim_nir_sweep.png"), fig, width = 190, height = 62, units = "mm", dpi = 300)
cat(sprintf(
  "wrote sim_nir_sweep.csv + agent/nnil_foil_sim_nir_sweep.png | nir*(DSC)=%.3f nir*(mm)=%.3f gen=%.3f\n",
  nir_star_dice, nir_star_mm, NIR_GEN
))
