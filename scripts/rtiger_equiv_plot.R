#!/usr/bin/env Rscript
# Re-plot the RTIGER equivalence figure (rtiger_equivalence_110K_cpp.png) from the
# preserved convergence trace, WITHOUT re-fitting. Mirrors the plotting block of
# scripts/bench_rtiger_cpp_vs_julia.R, but labels the right panel correctly: the
# per-iteration delta is the EM PARAMETER-convergence delta, max |Delta(alpha,beta)|
# over the BetaBinomial emission parameters between successive EM iterations (RTIGER's
# own convergence criterion, compared against eps), NOT a log-likelihood delta.
#   Rscript scripts/rtiger_equiv_plot.R
ROOT <- "/Users/fvrodriguez/repos/zealhmm"
conv_trace <- read.csv(file.path(ROOT, "results/bench/rtiger_conv_trace.csv"))
FIGDIR <- file.path(ROOT, "nilhmm-paper/figures")
COL_ORIG <- "firebrick"
COL_CPP <- "steelblue"

tr0c <- conv_trace[conv_trace$core == "cpp" & conv_trace$level == 0L, ]
tr0o <- conv_trace[conv_trace$core == "orig" & conv_trace$level == 0L, ]
stopifnot(nrow(tr0c) >= 2, nrow(tr0o) >= 2)
piC <- tr0c$per_iter
piO <- tr0o$per_iter
n <- min(nrow(tr0c), nrow(tr0o))
dc <- tr0c$delta[seq_len(n)]
do <- tr0o$delta[seq_len(n)]
drel <- max(abs(do - dc) / pmax(abs(dc), 1e-12))
dt <- c(1e-3, 0.01, 0.1, 1, 10, 100, 1000)
dl <- c("0.001", "0.01", "0.1", "1", "10", "100", "1000")

png(file.path(FIGDIR, "rtiger_equivalence_110K_cpp.png"), width = 1400, height = 640, res = 130)
par(mfrow = c(1, 2), mar = c(4.4, 5.0, 3.6, 1))
# left: per-iteration wall time
plot(seq_along(piC), piC,
  log = "y", type = "b", pch = 19, col = COL_CPP,
  ylim = range(c(piO, piC)), yaxt = "n", xlab = "EM iteration",
  ylab = "per-iteration wall time (s)",
  main = sprintf(
    "Per-iteration time — 109,703 markers/sample\norig ~%.0f s vs C++ ~%.1f s per iter",
    mean(piO), mean(piC)
  )
)
axis(2, dt, dl, las = 1)
points(seq_along(piO), piO, type = "b", pch = 17, col = COL_ORIG, cex = 1.2)
legend("right", c(
  sprintf("nilHMM C++/Rcpp (%d it)", nrow(tr0c)),
  sprintf("original Julia (%d it)", nrow(tr0o))
),
col = c(COL_CPP, COL_ORIG), pch = c(19, 17), bty = "n", cex = 0.85
)
# right: EM parameter-convergence delta, 1:1
rng <- range(c(do, dc))
cols <- rev(hcl.colors(n, "viridis"))
plot(dc, do,
  log = "xy", pch = 19, col = cols, cex = 1.3, xlim = rng, ylim = rng,
  xaxt = "n", yaxt = "n",
  xlab = "nilHMM C++/Rcpp EM parameter delta",
  ylab = "original Julia EM parameter delta",
  main = sprintf("EM parameter-convergence delta (per iteration)\n1:1 => identical (max rel %.0e)", drel)
)
axis(1, dt, dl)
axis(2, dt, dl, las = 1)
abline(0, 1, col = "grey50", lty = 2, lwd = 2)
legend("topleft", c("early iterations", "late iterations"),
  pch = 19, col = c(cols[1], cols[n]), bty = "n", cex = 0.85
)
dev.off()
cat(sprintf("wrote rtiger_equivalence_110K_cpp.png (max rel delta %.1e)\n", drel))
