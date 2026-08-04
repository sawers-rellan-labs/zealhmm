#!/usr/bin/env Rscript
# ZEAL nnil `nir` calibration on the count-from-parents 6-sib-pool sim.
#
# nnil is the CATEGORICAL caller (dispatcher: hard genotypes -> categorical), so
# the ZEAL skim counts are HARD-CALLED first via the BC2S3 DESIGN-prior MAP
# (call_gt, NOT HWE/flat), then decoded by nnil. Following S6, the free knob is a
# MINIMUM-DISTANCE estimator vs sim truth, not an MLE: we FIX rrate at the map
# value 2L/(100n) (no r sweep) and sweep only `nir`, scoring called-vs-truth by
# donor-fragment DSC (max) and fragment-size KS (min). Ported from
# scripts/nnil_foil/09_sim_nir_sweep.R onto the ZEAL Tier-A pool sim; unlike the
# nNIL case there is no scalar "generation nir" (non-informativeness emerges from
# the teosinte donor genotypes), so nir* is purely the caller value that best
# recovers the simulated ancestry.
#
#   Rscript scripts/zeal_nir_sweep.R
# Output: data/zeal/zeal_nir_sweep.csv + agent/zeal_nir_sweep.png

suppressMessages({
  library(nilHMM)
  library(data.table)
  library(ggplot2)
})
root <- here::here()
for (f in list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)) source(f)
source(file.path(root, "scripts/logging.R"))

SIM <- file.path(root, "results/sim/zeal_pool/zeal_pool_bc2s3_full.rds")
N_CAL <- 300L # calibration subset (decode ~0.02s/line; raise if wanted)
DESIGN <- "BC2S3"

sim <- readRDS(SIM)
design_prior <- breeding_prior(DESIGN)
map_r <- sim$map_r
ids <- head(sim$names, min(N_CAL, length(sim$names)))
idx <- match(ids, sim$names)
M <- nrow(sim$grid)

# ---- long count table for the calibration subset ---------------------------
skim <- data.table(
  name = rep(ids, each = M),
  chr = rep(sim$grid$chr, length(idx)),
  pos = rep(sim$grid$pos, length(idx)),
  n_ref = as.integer(sim$n_ref[, idx]),
  n_alt = as.integer(sim$n_alt[, idx])
)
# hard-call ONCE via the BC2S3 design-prior MAP (call_gt); zero-depth -> missing
skim[, g := {
  gg <- call_gt(n_ref, n_alt, prior = design_prior, error = 0.01)
  gg[is.na(gg)] <- 3L
  as.integer(gg)
}]
setorder(skim, name, chr, pos)

truth <- sim$truth[name %in% ids] # pooled-dose ancestry truth (segments)
grid_eval <- as.data.table(sim$grid)[, .(chr, pos)]
tr_sizes <- donor_block_sizes(truth)
log_info(
  "ZEAL nir sweep: %d cal lines x %d markers, hard-called; fixed map r=%.3e (design %s)",
  length(ids), M, map_r, DESIGN
)

# ---- sweep the caller's nir at fixed map r, score vs the pooled-dose truth ---
nir_grid <- c(0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.594, 0.7, 0.8, 0.9, 0.95, 0.99)
t0 <- Sys.time()
sweep <- rbindlist(lapply(nir_grid, function(v) {
  called <- as.data.table(call_ancestry(
    data = skim[, .(name, chr, pos, g)], caller = "nnil", rrate = map_r,
    design = DESIGN, germ = 0.01, gert = 0.5, p = 0, nir = v, mr = 0.1
  ))
  mf <- marker_dsc(called, truth, grid_eval)
  ff <- donor_fragment_dsc(called, truth)
  r <- data.table(
    nir = v, marker_mismatch = 1 - mf$accuracy, donor_frag_dsc = ff$dsc,
    frag_ks = fragment_size_ks(donor_block_sizes(called), tr_sizes),
    donor_marker_dsc = mf$per_class[class == "donor(>0)"]$dsc
  )
  log_info(
    "  nir=%.3f | frag_dsc=%.3f mismatch=%.4f (%.0fs)", v, r$donor_frag_dsc,
    r$marker_mismatch, as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )
  r
}))
fwrite(sweep, file.path(root, "data/zeal/zeal_nir_sweep.csv"))
nir_star_dsc <- sweep$nir[which.max(sweep$donor_frag_dsc)]
nir_star_ks <- sweep$nir[which.min(sweep$frag_ks)]
log_info("ZEAL nir* | donor-fragment-DSC max at nir=%.3f | fragment-KS min at nir=%.3f", nir_star_dsc, nir_star_ks)

# ---- figure: 3 min-distance metrics vs caller nir (every series labeled) -----
long <- rbind(
  data.table(
    nir = sweep$nir, val = sweep$marker_mismatch,
    metric = "ZEAL sim (nnil on BC2S3 6-sib pool): marker mismatch vs pooled-dose truth (lower better)"
  ),
  data.table(
    nir = sweep$nir, val = sweep$donor_frag_dsc,
    metric = "ZEAL sim (nnil on BC2S3 6-sib pool): donor-fragment DSC vs pooled-dose truth (higher better)"
  ),
  data.table(
    nir = sweep$nir, val = sweep$frag_ks,
    metric = "ZEAL sim (nnil on BC2S3 6-sib pool): fragment-size KS vs pooled-dose truth (lower better)"
  )
)
fig <- ggplot(long, aes(nir, val)) +
  geom_line(linewidth = 0.6, colour = "#0072B2") +
  geom_point(size = 1.3, colour = "#0072B2") +
  facet_wrap(~metric, scales = "free_y", nrow = 1, labeller = label_wrap_gen(38)) +
  labs(
    x = expression("caller " * italic(nir) * " (SNP50K, full set, no thinning; rrate fixed at map value)"),
    y = NULL,
    title = sprintf(
      "ZEAL nnil nir sweep on the count-from-parents BC2S3 6-sib-pool sim (Tier A het); nir*(DSC)=%.3f nir*(KS)=%.3f",
      nir_star_dsc, nir_star_ks
    )
  ) +
  theme_classic(base_size = 9) +
  theme(text = element_text(family = "sans"), plot.title = element_text(size = 8), strip.text = element_text(size = 7))
ggsave(file.path(root, "agent/zeal_nir_sweep.png"), fig, width = 190, height = 66, units = "mm", dpi = 300)
cat(sprintf(
  "wrote data/zeal/zeal_nir_sweep.csv + agent/zeal_nir_sweep.png | nir*(DSC)=%.3f nir*(KS)=%.3f\n",
  nir_star_dsc, nir_star_ks
))
