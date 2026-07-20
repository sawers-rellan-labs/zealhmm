#!/usr/bin/env Rscript
# RTIGER scaling figure (mirrors nnil_marker_scaling.png): two panels, wall time and peak
# memory vs markers per sample on the Arabidopsis shared-3 panel, original
# array-retaining Julia vs nilHMM C++/Rcpp. Both cores measured on this machine.
# Reads:
#   results/bench/rtiger_marker_scaling.csv          (per-iteration time; core=cpp/orig)
#   results/bench/rtiger_memory_markers.csv          (C++ peak RSS)
#   results/bench/rtiger_julia_memory_markers.csv    (original Julia peak RSS)
#   Rscript scripts/rtiger_memory_plot.R
ROOT <- "/Users/fvrodriguez/repos/zealhmm"
tm <- read.csv(file.path(ROOT, "results/bench/rtiger_marker_scaling.csv"))
cppm <- read.csv(file.path(ROOT, "results/bench/rtiger_memory_markers.csv"))
julm <- read.csv(file.path(ROOT, "results/bench/rtiger_julia_memory_markers.csv"))
FIG <- file.path(ROOT, "nilhmm-paper/figures/rtiger_marker_scaling.png")

cols <- c(cpp = "steelblue", orig = "firebrick")
pch <- c(cpp = 19, orig = 17)
xt <- c(10000, 20000, 50000, 100000)
xl <- formatC(xt, big.mark = ",", format = "d")
axfmt <- function(side, at) axis(side, at, formatC(at, big.mark = ",", format = "d"), las = 1)

png(FIG, width = 1500, height = 640, res = 130)
par(mfrow = c(1, 2), mar = c(4.6, 5.2, 3.2, 1))

# ---- panel A: per-iteration wall time vs markers ----
tc <- tm[tm$core == "cpp", ]
to <- tm[tm$core == "orig", ]
et <- function(d) coef(stats::lm(log10(per_iter) ~ log10(mps), d))[2]
plot(NA,
  log = "xy", xlim = range(tm$mps), ylim = range(tm$per_iter), xaxt = "n", yaxt = "n",
  xlab = "markers per sample", ylab = "per-iteration wall time (s)",
  main = "RTIGER time scaling (shared panel)"
)
grid(col = "grey92")
axis(1, xt, xl)
yt <- c(0.01, 0.1, 1, 10, 100)
axis(2, yt, c("0.01", "0.1", "1", "10", "100"), las = 1)
for (nm in c("orig", "cpp")) {
  d <- tm[tm$core == nm, ]
  d <- d[order(d$mps), ]
  f <- stats::lm(log10(per_iter) ~ log10(mps), d)
  xx <- 10^seq(log10(min(d$mps)), log10(max(d$mps)), length = 50)
  lines(xx, 10^predict(f, data.frame(mps = xx)), col = cols[nm], lty = 2, lwd = 2)
  points(d$mps, d$per_iter, col = cols[nm], pch = pch[nm], cex = 1.3)
}
legend("topleft",
  bty = "n",
  legend = c(
    sprintf("original Julia ~ markers^%.2f", et(to)),
    sprintf("nilHMM C++/Rcpp ~ markers^%.2f", et(tc))
  ),
  col = cols[c("orig", "cpp")], pch = pch[c("orig", "cpp")], lty = 2, lwd = 2
)

# ---- panel B: peak RSS vs markers ----
em <- function(d) coef(stats::lm(log10(peak_rss_mib) ~ log10(markers), d))[2]
plot(NA,
  log = "xy", xlim = range(cppm$markers), ylim = range(c(cppm$peak_rss_mib, julm$peak_rss_mib)),
  xaxt = "n", xlab = "markers per sample", ylab = "peak RSS (MiB)",
  main = "RTIGER memory scaling (shared panel)"
)
grid(col = "grey92")
axis(1, xt, xl)
for (dd in list(list(d = julm, nm = "orig"), list(d = cppm, nm = "cpp"))) {
  d <- dd$d[order(dd$d$markers), ]
  nm <- dd$nm
  f <- stats::lm(log10(peak_rss_mib) ~ log10(markers), d)
  xx <- 10^seq(log10(min(d$markers)), log10(max(d$markers)), length = 50)
  lines(xx, 10^predict(f, data.frame(markers = xx)), col = cols[nm], lty = 2, lwd = 2)
  points(d$markers, d$peak_rss_mib, col = cols[nm], pch = pch[nm], cex = 1.3)
}
legend("left",
  bty = "n",
  legend = c(
    sprintf("original Julia (retains) ~ markers^%.2f", em(julm)),
    sprintf("nilHMM C++/Rcpp (streams) ~ markers^%.2f", em(cppm))
  ),
  col = cols[c("orig", "cpp")], pch = pch[c("orig", "cpp")], lty = 2, lwd = 2
)
dev.off()
cat("wrote", FIG, "\n")
