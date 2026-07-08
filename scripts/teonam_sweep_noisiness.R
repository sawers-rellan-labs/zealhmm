#!/usr/bin/env Rscript
# Noisiness of each caller's STAM Manhattan track, measured as spectral power in a
# PHYSICALLY-DEFINED noise band. Each track is a 1-D signal -log10P(bp); real ancestry
# structure lives at long wavelengths (>= introgression-fragment size), sub-resolution
# structure is meaningless, so the NOISE BAND is wavelengths in (resolution, P10_frag):
# spurious over-fragmentation shorter than ~90% of real tracts but above resolution.
#
#   * resolution R  = median adjacent-marker bp gap on the 0.1 cM inference grid
#                     (finest scale at which any caller can place a breakpoint).
#   * P10_frag      = 10th-percentile teosinte introgression-tract length (bp), from
#                     the CLEAN FSFHap ancestry mosaic (runs carrying teosinte, state>=1).
#   * HEADLINE      = FFT band-power fraction: periodogram (Hann-tapered, per chr, on a
#                     uniform bp grid) power in freq (1/P10, 1/R) / total power.
#   * wavelet check = Haar multi-level detail ENERGY in the octaves overlapping that band
#                     / total detail energy (localized, peak-robust corroboration).
# Compared at lambda=Inf (each caller's ceiling); also FFT band-power vs coverage.
# Output: results/sim/teonam/sweep_noisiness_118k.csv
# Run: Rscript scripts/teonam_sweep_noisiness.R
suppressMessages(library(data.table))
setwd("/Users/fvrodriguez/repos/zealhmm")
CALLERS <- c("control", "nnil", "rtiger", "lbimpute")

# --- 1. resolution R (bp) from the inference grid ----------------------------
thin <- fread("data/teonam/markers_v5_gwas118k_cm_thin01.tsv")
R <- thin[, as.numeric(median(diff(sort(pos_v5)))), by = chr][, median(V1)]

# --- 2. P10 introgression-tract length (bp) from the clean FSFHap mosaic ------
g <- readRDS("data/teonam/teonam_gwas118k_dosage_fsfhap.rds")
dos <- g$dos
mc <- fread("data/teonam/markers_v5_gwas118k_cm.tsv")
pos_by <- setNames(mc$pos_v5, mc$marker)
chr_by <- setNames(mc$chr, mc$marker)
ord <- order(chr_by[rownames(dos)], pos_by[rownames(dos)])
dos <- dos[ord, ]
pos <- pos_by[rownames(dos)]
chr <- chr_by[rownames(dos)]
tracts <- unlist(lapply(1:10, function(cc) {
  idx <- which(chr == cc)
  p <- pos[idx]
  unlist(lapply(seq_len(ncol(dos)), function(j) {
    teo <- as.integer(dos[idx, j] >= 1) # carries teosinte
    teo[is.na(teo)] <- 0
    r <- rle(teo)
    ends <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1
    sel <- r$values == 1 & r$lengths >= 2
    if (!any(sel)) {
      return(numeric(0))
    }
    p[ends[sel]] - p[starts[sel]]
  }))
}))
# raw P10 is sub-resolution (~34% of teo tracts are single-marker het blips in the
# mosaic); take P10 over RESOLVABLE tracts (>= R) = the smallest real introgressions.
P10 <- as.numeric(quantile(tracts[tracts >= R], 0.10))
cat(sprintf(
  "resolution R = %.0f kb (thin-grid median gap) | resolvable teosinte tract P10 = %.0f kb (median %.2f Mb; %d tracts, %.0f%% sub-R dropped)\n",
  R / 1e3, P10 / 1e3, median(tracts[tracts >= R]) / 1e6, length(tracts), 100 * mean(tracts < R)
))
cat(sprintf("NOISE BAND: wavelength (%.0f kb, %.2f Mb)  <=>  freq (%.3g, %.3g) cycles/bp\n\n", R / 1e3, P10 / 1e6, 1 / P10, 1 / R))

# --- FFT band-power fraction for one caller track at one coverage -------------
fft_band <- function(sw, cov) {
  s <- sw[coverage == cov & is.finite(P) & P > 0]
  s[, logP := -log10(P)]
  per <- lapply(1:10, function(cc) {
    z <- s[CHR == cc][order(BP)]
    if (nrow(z) < 16) {
      return(c(NA, NA))
    }
    dx <- R / 2 # Nyquist wavelength = R
    xg <- seq(min(z$BP), max(z$BP), by = dx)
    if (length(xg) < 16) {
      return(c(NA, NA))
    }
    yg <- approx(z$BP, z$logP, xout = xg, method = "constant", rule = 2)$y
    yg <- (yg - mean(yg)) * (0.5 - 0.5 * cos(2 * pi * (seq_along(yg) - 1) / (length(yg) - 1))) # Hann
    N <- length(yg)
    pw <- Mod(fft(yg))^2
    fr <- (seq_len(N) - 1) / (N * dx)
    fr <- pmin(fr, 1 / (N * dx) * (N - (seq_len(N) - 1))) # fold to [0,Nyq]
    keep <- fr > 0
    band <- fr > (1 / P10) & fr <= (1 / R)
    c(sum(pw[band]), sum(pw[keep]))
  })
  M <- do.call(rbind, per)
  M <- M[stats::complete.cases(M), , drop = FALSE]
  round(sum(M[, 1]) / sum(M[, 2]), 4)
}

# --- Haar multi-level detail energy in the band's octaves ---------------------
wav_band <- function(sw, cov) {
  s <- sw[coverage == cov & is.finite(P) & P > 0]
  s[, logP := -log10(P)]
  num <- den <- 0
  dx <- R / 2
  for (cc in 1:10) {
    z <- s[CHR == cc][order(BP)]
    if (nrow(z) < 32) next
    xg <- seq(min(z$BP), max(z$BP), by = dx)
    yg <- approx(z$BP, z$logP, xout = xg, method = "constant", rule = 2)$y
    a <- yg
    j <- 0
    while (length(a) >= 2) {
      j <- j + 1
      n <- length(a)
      m <- n - n %% 2
      a2 <- a[seq(1, m, 2)]
      b2 <- a[seq(2, m, 2)]
      d <- (b2 - a2) / sqrt(2)
      wl <- 2^j * dx # approx wavelength of level-j detail
      e <- sum(d^2)
      den <- den + e
      if (wl > R & wl < P10) num <- num + e
      a <- (a2 + b2) / sqrt(2)
    }
  }
  round(num / den, 4)
}

sweeps <- lapply(CALLERS, function(c) fread(sprintf("results/sim/teonam/stam_gwas_%s_118k_sweep.csv", c)))
names(sweeps) <- CALLERS
covs <- sort(unique(sweeps[[1]]$coverage))

head_tab <- rbindlist(lapply(CALLERS, function(c) {
  data.table(caller = c, fft_band = fft_band(sweeps[[c]], Inf), wav_band = wav_band(sweeps[[c]], Inf))
}))[order(fft_band)]
cat("=== NOISE-BAND power at lambda=Inf (fraction in the (R, P10) band), sorted ===\n")
print(head_tab)

grid_tab <- rbindlist(lapply(CALLERS, function(c) {
  data.table(caller = c, t(sapply(covs, function(v) fft_band(sweeps[[c]], v))))
}))
setnames(grid_tab, c("caller", paste0("l", covs)))
cat("\n=== FFT band-power fraction vs coverage ===\n")
print(grid_tab)

fwrite(head_tab, "results/sim/teonam/sweep_noisiness_118k.csv")
fwrite(grid_tab, "results/sim/teonam/sweep_noisiness_vscov_118k.csv")
fwrite(
  data.table(R_kb = round(R / 1e3), P10_kb = round(P10 / 1e3), median_frag_Mb = round(median(tracts[tracts >= R]) / 1e6, 2)),
  "results/sim/teonam/sweep_noisiness_band_118k.csv"
)
cat("\nwrote sweep_noisiness_118k.csv + _vscov_118k.csv + _band_118k.csv\n")
