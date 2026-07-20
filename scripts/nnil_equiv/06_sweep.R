#!/usr/bin/env Rscript
# nNIL time + memory scaling with marker density, comparable to the RTIGER §1
# marker sweep: odd-index thinning of the 64,025-marker panel to seven sizes, over
# the full 888-line population, both callers (Holland hmmlearn vs nilHMM nnil) on
# the same memory-mapped .bed, each under `/usr/bin/time -l` for peak RSS.
# Two-panel figure: wall time vs markers (log-log, power-law fits + exponents) and
# peak RSS vs markers.
#   Rscript scripts/nnil_equiv/06_sweep.R
suppressMessages(library(nilHMM)) # (paths only)
ROOT <- "/Users/fvrodriguez/repos/zealhmm"
SDIR <- file.path(ROOT, "scripts/nnil_equiv")
OUTDIR <- file.path(ROOT, "results/bench")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
FIGDIR <- file.path(ROOT, "nilhmm-paper/figures")
PYBIN <- path.expand("~/anaconda3/envs/nilhmm/bin/python")
LEVELS <- 0:6 # 64025 -> ~1000 markers, full 888 lines
log_info <- function(...) cat(sprintf("[06_sweep] %s\n", sprintf(...)))

run1 <- function(caller, level) {
  cmd <- if (caller == "holland") {
    c("-l", PYBIN, file.path(SDIR, "05_holland_worker.py"), "--level", level)
  } else {
    c("-l", "Rscript", file.path(SDIR, "05_nilhmm_worker.R"), "--level", level)
  }
  out <- suppressWarnings(system2("/usr/bin/time", cmd, stdout = TRUE, stderr = TRUE))
  res <- grep("^RESULT", out, value = TRUE)
  rss <- grep("maximum resident set size", out, value = TRUE)
  g <- function(f, s) as.numeric(sub(paste0(".*", f, "=([0-9.]+).*"), "\\1", s))
  data.frame(
    caller = caller, level = level, markers = g("markers", res),
    lines = g("lines", res), seconds = g("seconds", res),
    rss_gib = as.numeric(sub("^\\s*([0-9]+).*", "\\1", rss)) / 1024^3
  )
}

log_info("marker sweep, full population, levels %s x {nilhmm, holland} ...", paste(LEVELS, collapse = ","))
rows <- list()
for (lv in LEVELS) {
  for (ca in c("nilhmm", "holland")) {
    r <- run1(ca, lv)
    rows[[length(rows) + 1L]] <- r
    log_info("  %-7s L%d  %6.0f markers  %6.2fs  %.2f GiB", ca, lv, r$markers, r$seconds, r$rss_gib)
  }
}
d <- do.call(rbind, rows)
NLINES <- d$lines[1]
write.csv(d, file.path(OUTDIR, "nnil_scaling.csv"), row.names = FALSE)
log_info("wrote nnil_scaling.csv")

# ---- two-panel figure ----
cols <- c(nilhmm = "steelblue", holland = "firebrick")
pch <- c(nilhmm = 19, holland = 17)
expo <- function(ca) coef(stats::lm(log10(seconds) ~ log10(markers), d[d$caller == ca, ]))[2]
png(file.path(FIGDIR, "nnil_marker_scaling.png"), width = 1500, height = 640, res = 130)
par(mfrow = c(1, 2), mar = c(4.6, 5.2, 3.2, 1))

# panel A: wall time vs markers, log-log, power-law fits (RTIGER-style)
plot(NA,
  log = "xy", xlim = range(d$markers), ylim = range(d$seconds),
  xlab = "markers per line", ylab = "wall time (s)",
  main = sprintf("nNIL decode time (%d lines)", NLINES)
)
grid(col = "grey92")
for (ca in c("holland", "nilhmm")) {
  s <- d[d$caller == ca, ]
  s <- s[order(s$markers), ]
  f <- stats::lm(log10(seconds) ~ log10(markers), s)
  xx <- 10^seq(log10(min(s$markers)), log10(max(s$markers)), length = 50)
  lines(xx, 10^predict(f, data.frame(markers = xx)), col = cols[ca], lty = 2, lwd = 2)
  points(s$markers, s$seconds, col = cols[ca], pch = pch[ca], cex = 1.3)
}
legend("topleft",
  bty = "n",
  legend = c(
    sprintf("Holland (hmmlearn) ~ markers^%.2f", expo("holland")),
    sprintf("nilHMM nnil ~ markers^%.2f", expo("nilhmm"))
  ),
  col = cols[c("holland", "nilhmm")], pch = pch[c("holland", "nilhmm")], lty = 2, lwd = 2
)

# panel B: peak RSS vs markers (log-log, power-law fit; memory has a runtime floor
# plus the held matrix, so the exponent is an effective slope, not a pure power law)
expm <- function(ca) coef(stats::lm(log10(rss_gib) ~ log10(markers), d[d$caller == ca, ]))[2]
plot(NA,
  log = "xy", xlim = range(d$markers), ylim = range(d$rss_gib),
  xlab = "markers per line", ylab = "peak RSS (GiB)",
  main = sprintf("nNIL peak memory (%d lines)", NLINES)
)
grid(col = "grey92")
for (ca in c("holland", "nilhmm")) {
  s <- d[d$caller == ca, ]
  s <- s[order(s$markers), ]
  f <- stats::lm(log10(rss_gib) ~ log10(markers), s)
  xx <- 10^seq(log10(min(s$markers)), log10(max(s$markers)), length = 50)
  lines(xx, 10^predict(f, data.frame(markers = xx)), col = cols[ca], lty = 2, lwd = 2)
  points(s$markers, s$rss_gib, col = cols[ca], pch = pch[ca], cex = 1.3)
}
legend("topleft",
  bty = "n",
  legend = c(
    sprintf("Holland (hmmlearn) ~ markers^%.2f", expm("holland")),
    sprintf("nilHMM nnil ~ markers^%.2f", expm("nilhmm"))
  ),
  col = cols[c("holland", "nilhmm")], pch = pch[c("holland", "nilhmm")], lty = 2, lwd = 2
)
dev.off()
log_info("wrote figures/nnil_marker_scaling.png")
