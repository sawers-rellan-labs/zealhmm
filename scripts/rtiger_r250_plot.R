#!/usr/bin/env Rscript
# Figures for the RTIGER benchmark at the operating rigidity r=250. The nilHMM
# C++/Rcpp core is measured at all five sizes; the upstream-original Julia is
# measured at the three small sizes and PROJECTED (power-law fit) to the two large
# sizes (open symbols, dashed projection). Two figures, mirroring the r=2 pair:
#   rtiger_marker_scaling_r250.png     (2-panel: per-iter time + peak memory)
#   rtiger_equivalence_r250.png        (per-iter time + EM parameter-delta 1:1,
#                                       at the largest MEASURED size)
#   Rscript scripts/rtiger_r250_plot.R
ROOT <- "/Users/fvrodriguez/repos/zealhmm"
B <- file.path(ROOT, "results/bench")
FIGDIR <- file.path(ROOT, "nilhmm-paper/figures")
speed <- read.csv(file.path(B, "rtiger_marker_scaling_r250.csv"))
trace <- read.csv(file.path(B, "rtiger_conv_trace_r250.csv"))
cppm <- read.csv(file.path(B, "rtiger_memory_markers_r250.csv"))
julm <- read.csv(file.path(B, "rtiger_julia_memory_markers_r250.csv"))
COL_ORIG <- "firebrick"
COL_CPP <- "steelblue"

sc <- speed[speed$core == "cpp", ]
sc <- sc[order(sc$mps), ]
so <- speed[speed$core == "orig", ]
so <- so[order(so$mps), ]
allmps <- sort(unique(sc$mps))
proj_mps <- setdiff(allmps, so$mps) # the large sizes the original is projected to
xt <- c(10000, 20000, 50000, 100000)
xl <- formatC(xt, big.mark = ",", format = "d")
et <- function(d, y) coef(stats::lm(log10(d[[y]]) ~ log10(d$mps)))[2]
em <- function(d, y = "peak_rss_mib") coef(stats::lm(log10(d[[y]]) ~ log10(d$markers)))[2]

png(file.path(FIGDIR, "rtiger_marker_scaling_r250.png"), width = 1500, height = 640, res = 130)
par(mfrow = c(1, 2), mar = c(4.6, 5.2, 3.2, 1))

# ---- panel A: per-iteration wall time vs markers (r=250) ----
fo <- stats::lm(log10(per_iter) ~ log10(mps), so)
fc <- stats::lm(log10(per_iter) ~ log10(mps), sc)
o_proj <- data.frame(mps = proj_mps, per_iter = 10^predict(fo, data.frame(mps = proj_mps)))
yr <- range(c(sc$per_iter, so$per_iter, o_proj$per_iter))
plot(NA,
  log = "xy", xlim = range(allmps), ylim = yr, xaxt = "n", yaxt = "n",
  xlab = "markers per sample", ylab = "per-iteration wall time (s)",
  main = "RTIGER time scaling (shared panel, r=250)"
)
grid(col = "grey92")
axis(1, xt, xl)
yt <- c(0.01, 0.1, 1, 10, 100, 1000)
axis(2, yt, c("0.01", "0.1", "1", "10", "100", "1000"), las = 1)
xx <- 10^seq(log10(min(allmps)), log10(max(allmps)), length = 50)
lines(xx, 10^predict(fo, data.frame(mps = xx)), col = COL_ORIG, lty = 2, lwd = 2)
lines(xx, 10^predict(fc, data.frame(mps = xx)), col = COL_CPP, lty = 2, lwd = 2)
points(so$mps, so$per_iter, col = COL_ORIG, pch = 17, cex = 1.3)
points(sc$mps, sc$per_iter, col = COL_CPP, pch = 19, cex = 1.3)
legend("topleft", bty = "n", cex = 0.9, legend = c(
  sprintf("original Julia ~ markers^%.2f", et(so, "per_iter")),
  sprintf("nilHMM C++/Rcpp ~ markers^%.2f", et(sc, "per_iter"))
), col = c(COL_ORIG, COL_CPP), pch = c(17, 19), lty = 2, lwd = 2)

# ---- panel B: peak RSS vs markers (r=250) ----
cppm <- cppm[order(cppm$markers), ]
julm <- julm[order(julm$markers), ]
fom <- stats::lm(log10(peak_rss_mib) ~ log10(markers), julm)
jul_proj_mps <- setdiff(cppm$markers, julm$markers)
jul_proj <- 10^predict(fom, data.frame(markers = jul_proj_mps))
plot(NA,
  log = "xy", xlim = range(cppm$markers), ylim = range(c(cppm$peak_rss_mib, julm$peak_rss_mib, jul_proj)),
  xaxt = "n", xlab = "markers per sample", ylab = "peak RSS (MiB)",
  main = "RTIGER memory scaling (shared panel, r=250)"
)
grid(col = "grey92")
axis(1, xt, xl)
xxm <- 10^seq(log10(min(cppm$markers)), log10(max(cppm$markers)), length = 50)
lines(xxm, 10^predict(fom, data.frame(markers = xxm)), col = COL_ORIG, lty = 2, lwd = 2)
fcm <- stats::lm(log10(peak_rss_mib) ~ log10(markers), cppm)
lines(xxm, 10^predict(fcm, data.frame(markers = xxm)), col = COL_CPP, lty = 2, lwd = 2)
points(julm$markers, julm$peak_rss_mib, col = COL_ORIG, pch = 17, cex = 1.3)
points(cppm$markers, cppm$peak_rss_mib, col = COL_CPP, pch = 19, cex = 1.3)
legend("left", bty = "n", cex = 0.9, legend = c(
  sprintf("original Julia (retains) ~ markers^%.2f", em(julm)),
  sprintf("nilHMM C++/Rcpp (streams) ~ markers^%.2f", em(cppm))
), col = c(COL_ORIG, COL_CPP), pch = c(17, 19), lty = 2, lwd = 2)
dev.off()
cat("wrote rtiger_marker_scaling_r250.png\n")

# ---- equivalence figure at the largest MEASURED size ----
meas <- max(so$level == min(so$level)) # placeholder
lvl <- min(so$level) # smallest level index == largest markers among measured original
tc <- trace[trace$core == "cpp" & trace$level == lvl, ]
to <- trace[trace$core == "orig" & trace$level == lvl, ]
mps_eq <- unique(tc$mps)
if (nrow(tc) >= 2 && nrow(to) >= 2) {
  n <- min(nrow(tc), nrow(to))
  dc <- tc$delta[seq_len(n)]
  do <- to$delta[seq_len(n)]
  drel <- max(abs(do - dc) / pmax(abs(dc), 1e-12))
  dt <- c(1e-3, 0.01, 0.1, 1, 10, 100, 1000)
  dl <- c("0.001", "0.01", "0.1", "1", "10", "100", "1000")
  png(file.path(FIGDIR, "rtiger_equivalence_r250.png"), width = 1400, height = 640, res = 130)
  par(mfrow = c(1, 2), mar = c(4.4, 5.0, 3.6, 1))
  plot(seq_len(nrow(tc)), tc$per_iter,
    log = "y", type = "b", pch = 19, col = COL_CPP,
    ylim = range(c(to$per_iter, tc$per_iter)), xaxt = "n", yaxt = "n", xlab = "EM iteration",
    ylab = "per-iteration wall time (s)",
    main = sprintf("Per-iteration time -- %s markers/sample, r=250\norig ~%.0f s vs C++ ~%.2f s per iter", format(mps_eq, big.mark = ","), mean(to$per_iter), mean(tc$per_iter))
  )
  axis(1, at = seq_len(max(nrow(tc), nrow(to)))) # EM iterations are integers
  axis(2, dt, dl, las = 1)
  points(seq_len(nrow(to)), to$per_iter, type = "b", pch = 17, col = COL_ORIG, cex = 1.2)
  legend("right", c(sprintf("nilHMM C++/Rcpp (%d it)", nrow(tc)), sprintf("original Julia (%d it)", nrow(to))),
    col = c(COL_CPP, COL_ORIG), pch = c(19, 17), bty = "n", cex = 0.85
  )
  rng <- range(c(do, dc))
  cols <- rev(hcl.colors(n, "viridis"))
  plot(dc, do,
    log = "xy", pch = 19, col = cols, cex = 1.3, xlim = rng, ylim = rng, xaxt = "n", yaxt = "n",
    xlab = "nilHMM C++/Rcpp EM parameter delta", ylab = "original Julia EM parameter delta",
    main = sprintf("EM parameter-convergence delta (per iteration), r=250\n1:1 => identical (max rel %.0e)", drel)
  )
  axis(1, dt, dl)
  axis(2, dt, dl, las = 1)
  abline(0, 1, col = "grey50", lty = 2, lwd = 2)
  legend("topleft", c("early iterations", "late iterations"), pch = 19, col = c(cols[1], cols[n]), bty = "n", cex = 0.85)
  dev.off()
  cat(sprintf("wrote rtiger_equivalence_r250.png (level %d, %d markers, max rel delta %.1e)\n", lvl, mps_eq, drel))
} else {
  cat("skipped equivalence figure (need cpp+orig delta trace at a shared level)\n")
}
