#!/usr/bin/env Rscript
# =============================================================================
# nilHMM C++/Rcpp RTIGER core vs upstream-original Julia RTIGER
# =============================================================================
# Reproduces, for the C++/Rcpp core, the two figures the fork produced for its
# optimized-Julia core (docs/optimization.md: marker_scaling.png,
# equivalence_110K.png). The comparison is now UPSTREAM-ORIGINAL Julia (red) vs
# nilHMM C++/Rcpp (blue), on the shared-3 Arabidopsis Col x Ler panel.
#
# The original-Julia baseline is NOT re-run (its 110k convergence fit is ~3.4 h
# upstream); it is read from the preserved fork dumps in data/bench_ref/. Only
# the C++ core is measured live here.
#
# Method (mirrors the fork's shared-panel marker sweep):
#   * one shared panel = 109,703 loci covered in all three BN/Z/AU samples,
#     thinned by odd index to 6,857 / 13,713 / 27,426 / 54,852 / 109,703.
#   * throughput: fixed-iteration fit (eps=0.01, r=2, max_iter=6) per size,
#     per-iteration mean wall time -> marker-scaling figure.
#   * equivalence: fit to convergence from the IDENTICAL deterministic init the
#     original used (alpha=[20,20,1], beta=[1,20,20], A=0.1+10I normalised,
#     uniform pi), decode (Viterbi + border postprocessing), and compare the
#     paths position-by-position against the stored original-Julia oracle.
#   * per-iteration delta trace (cold-restart reconstruction: fit(max_iter=m)
#     for m=1..K; delta_m = max|Delta(alpha,beta)| between successive fits, the
#     same metric Julia logs) -> equivalence-trajectory figure.
#
# Outputs (results/bench/, gitignored) + figures (nilhmm-paper/figures/):
#   rtiger_marker_scaling.csv, rtiger_conv_trace.csv, rtiger_equivalence.csv
#   rtiger_marker_scaling_cpp.png, rtiger_equivalence_110K_cpp.png
#
# Usage:  Rscript scripts/bench_rtiger_cpp_vs_julia.R [--smoke] [--generate]
#   --smoke     only L4+L3, no L0 delta trace (fast sanity path)
#   --generate  full sweep (default)
# Persist a run:  ... 2>&1 | tee agent/bench_rtiger_<ts>_<pid>.log
# =============================================================================

suppressMessages(library(nilHMM))
ROOT <- tryCatch(here::here(), error = function(e) getwd())
if (!file.exists(file.path(ROOT, "scripts/logging.R"))) {
  ROOT <- "/Users/fvrodriguez/repos/zealhmm"
}
source(file.path(ROOT, "scripts/logging.R"))

args <- commandArgs(trailingOnly = TRUE)
SMOKE <- "--smoke" %in% args
REPS <- 3L # throughput timing reps (median)

PANELDIR <- file.path(ROOT, "data/rtiger_shared3_input")
REFDIR <- file.path(ROOT, "data/bench_ref")
OUTDIR <- file.path(ROOT, "results/bench")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
FIGDIR <- file.path(ROOT, "nilhmm-paper/figures")
files <- c(S1 = "sampleBN.txt", S2 = "sampleZ.txt", S3 = "sampleAU.txt")
chrs <- paste0("Chr", 1:5)
RIG <- 2L
EPS <- 0.01
INIT_A <- c(20, 20, 1)
INIT_B <- c(1, 20, 20) # == original-Julia deterministic init
LEVELS <- if (SMOKE) c(4L, 3L) else 0:4 # 0 = full 109,703; 4 = 6,857
if (SMOKE) { # never clobber canonical deliverables from a smoke run
  OUTDIR <- file.path(OUTDIR, "smoke")
  dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
  FIGDIR <- OUTDIR
}

log_info(
  "nilHMM C++/Rcpp vs original-Julia RTIGER benchmark  (%s)",
  if (SMOKE) "SMOKE" else "GENERATE"
)
log_info("panel: %s   ref(orig): %s   out: %s", PANELDIR, REFDIR, OUTDIR)

# ---- machine + toolchain record (written with every full/generate run) -----
write_env <- function() {
  sh <- function(cmd) {
    tryCatch(paste(system(cmd, intern = TRUE), collapse = " "),
      error = function(e) NA_character_
    )
  }
  pv <- function(p) tryCatch(format(packageVersion(p)), error = function(e) NA_character_)
  cpu <- sh("sysctl -n machdep.cpu.brand_string")
  cores <- sh("sysctl -n hw.logicalcpu")
  mem <- tryCatch(paste0(round(as.numeric(sh("sysctl -n hw.memsize")) / 1024^3), " GB"),
    error = function(e) NA_character_
  )
  lines <- c(
    sprintf("cpu\t%s (%s cores)", cpu, cores),
    sprintf("memory\t%s", mem),
    sprintf("os\t%s %s (%s)", Sys.info()[["sysname"]], Sys.info()[["release"]], Sys.info()[["machine"]]),
    sprintf("compiler\t%s %s", sh("R CMD config CXX"), sh("R CMD config CXXFLAGS")),
    sprintf("R\t%s", getRversion()),
    sprintf("julia\t%s", sub("julia version ", "", sh("julia --version"))),
    sprintf("nilHMM\t%s", pv("nilHMM")),
    sprintf("Rcpp\t%s", pv("Rcpp")),
    sprintf("RcppParallel\t%s", pv("RcppParallel")),
    "threads\t1"
  )
  writeLines(lines, file.path(OUTDIR, "rtiger_env.txt"))
  log_info(
    "env: %s, %s, %s %s, R %s, clang -O2 arm64, threads=1",
    cpu, mem, Sys.info()[["sysname"]], Sys.info()[["release"]], getRversion()
  )
}

# ---- panel -----------------------------------------------------------------
raw <- lapply(files, function(f) {
  read.table(file.path(PANELDIR, f),
    sep = "\t",
    col.names = c("chr", "pos", "refA", "refC", "altA", "altF"), stringsAsFactors = FALSE
  )
})
thin <- function(n, level) {
  ix <- seq_len(n)
  for (j in seq_len(level)) ix <- ix[c(TRUE, FALSE)]
  ix
}
mps_of <- function(level) length(thin(nrow(raw[[1]]), level))
build_obs <- function(level) {
  N <- nrow(raw[[1]])
  ix <- thin(N, level)
  obs <- list()
  for (s in names(files)) {
    d <- raw[[s]]
    obs[[s]] <- list()
    for (c in chrs) {
      rows <- ix[d$chr[ix] == c]
      if (!length(rows)) next
      k <- as.integer(round(d$refC[rows]))
      n <- k + as.integer(round(d$altF[rows]))
      obs[[s]][[c]] <- list(k = k, n = n)
    }
  }
  obs
}

# ---- original-Julia reference parsers --------------------------------------
parse_orig_log <- function(level) { # per-iteration delta + elapsed
  p <- file.path(REFDIR, "orig_conv", sprintf("log_orig_conv_L%d.log", level))
  if (!file.exists(p)) {
    return(NULL)
  }
  ln <- grep("^iter ", readLines(p), value = TRUE)
  if (!length(ln)) {
    return(NULL)
  }
  g <- function(f) as.numeric(sub(paste0(".*", f, "=([0-9.eE+-]+).*"), "\\1", ln))
  el <- g("elapsed")
  data.frame(
    iter = seq_along(ln), delta = g("delta"), elapsed = el,
    per_iter = diff(c(0, el))
  )
}
parse_orig_dump <- function(level) { # converged params + Viterbi oracle
  p <- file.path(REFDIR, "orig_conv", sprintf("conv_orig_L%d.txt", level))
  if (!file.exists(p)) {
    return(NULL)
  }
  L <- readLines(p)
  num <- function(key) {
    as.numeric(strsplit(gsub(
      ".*=\\[|\\].*", "",
      grep(paste0("^", key, "="), L, value = TRUE)[1]
    ), ",")[[1]])
  }
  meta <- L[1]
  iters <- as.integer(sub(".*iters=([0-9]+).*", "\\1", meta))
  vit <- L[grepl("^vit\\[", L)]
  nm <- sub("=.*", "", vit)
  paths <- lapply(
    sub("^[^=]*=", "", vit),
    function(s) as.integer(strsplit(gsub("[^0-9,]", "", s), ",")[[1]])
  )
  names(paths) <- nm
  list(
    iters = iters, alpha = num("alpha"), beta = num("beta"), pi = num("pi"),
    transition = num("transition"), vit = paths
  )
}
parse_orig_sweep <- function() { # fixed-iter throughput points
  p <- file.path(REFDIR, "30_panel_sweep_orig.txt")
  ln <- grep("^SWEEP", readLines(p), value = TRUE)
  g <- function(f) as.numeric(sub(paste0(".*", f, "=([0-9.eE+-]+).*"), "\\1", ln))
  data.frame(
    core = "orig", level = g("level"), mps = g("markers_per_sample"),
    iters = g("iters"), runtime = g("runtime"), per_iter = g("per_iter")
  )
}

fit_cpp <- function(obs, max_iter) {
  nilHMM:::.rtiger_fit(obs,
    r = RIG, nstates = 3L, eps = EPS, max_iter = as.integer(max_iter),
    init_alpha = INIT_A, init_beta = INIT_B
  )
}

# =============================================================================
# 0. MEMORY SWEEP  (C++ peak RSS vs sample count, fixed markers/sample)
# =============================================================================
# Two modes on this same script so the sweep is one tracked artifact:
#   --memworker N MK : build N cycled samples at MK markers each, fit once, quit
#                      (peak RSS of THIS process is measured by the parent).
#   --memory         : driver -- for each N run the worker under `/usr/bin/time -l`,
#                      parse peak RSS, write results/bench/rtiger_memory_sweep.csv.
SCRIPT <- file.path(ROOT, "scripts/bench_rtiger_cpp_vs_julia.R")
FIXED_MK <- 50000L
MEM_N <- c(2L, 4L, 8L, 16L, 30L, 60L)
build_obs_cycled <- function(N, markers) { # cycle the 3 real samples to N
  Ntot <- nrow(raw[[1]])
  keep <- unique(round(seq(1, Ntot, length.out = markers)))
  obs <- list()
  for (i in seq_len(N)) {
    d <- raw[[((i - 1L) %% 3L) + 1L]]
    nm <- sprintf("S%03d", i)
    obs[[nm]] <- list()
    for (c in chrs) {
      rows <- keep[d$chr[keep] == c]
      if (!length(rows)) next
      k <- as.integer(round(d$refC[rows]))
      n <- k + as.integer(round(d$altF[rows]))
      obs[[nm]][[c]] <- list(k = k, n = n)
    }
  }
  obs
}
if ("--memworker" %in% args) {
  ai <- which(args == "--memworker")
  N <- as.integer(args[ai + 1])
  mk <- as.integer(args[ai + 2])
  invisible(fit_cpp(build_obs_cycled(N, mk), 6L))
  cat(sprintf("MEMWORKER_DONE N=%d markers=%d\n", N, mk))
  quit(save = "no")
}
if ("--memory" %in% args) {
  log_info("--- C++ peak-RSS sweep: N cycled samples at %d markers each ---", FIXED_MK)
  rows <- list()
  for (N in MEM_N) {
    out <- system2("/usr/bin/time", c("-l", "Rscript", SCRIPT, "--memworker", N, FIXED_MK),
      stdout = TRUE, stderr = TRUE
    )
    rl <- grep("maximum resident set size", out, value = TRUE)
    by <- as.numeric(sub("^\\s*([0-9]+).*", "\\1", rl))
    mib <- by / 1024^2
    rows[[length(rows) + 1L]] <- data.frame(
      n_samples = N, markers = FIXED_MK,
      peak_rss_mib = round(mib, 1)
    )
    log_info("  N=%2d  %d mk/sample  peak RSS = %.1f MiB", N, FIXED_MK, mib)
  }
  mem <- do.call(rbind, rows)
  write.csv(mem, file.path(OUTDIR, "rtiger_memory_sweep.csv"), row.names = FALSE)
  lf <- stats::lm(peak_rss_mib ~ n_samples, mem)
  log_info(
    "wrote rtiger_memory_sweep.csv | peak RSS ~ %.0f MiB + %.2f MiB/sample (R^2=%.3f)",
    coef(lf)[1], coef(lf)[2], summary(lf)$r.squared
  )
  quit(save = "no")
}

# =============================================================================
# 1. THROUGHPUT  (fixed-iteration, matches the original sweep protocol)
# =============================================================================
write_env()
log_info("--- throughput sweep (C++), %d rep(s), fixed max_iter=6 ---", REPS)
thr <- list()
for (level in LEVELS) {
  obs <- build_obs(level)
  mps <- mps_of(level)
  reps <- vapply(seq_len(REPS), function(i) {
    system.time(f <- fit_cpp(obs, 6L))[["elapsed"]]
  }, numeric(1))
  fit <- fit_cpp(obs, 6L)
  it <- fit$iterations
  rt <- stats::median(reps)
  thr[[length(thr) + 1L]] <- data.frame(
    core = "cpp", level = level, mps = mps,
    iters = it, runtime = rt, per_iter = rt / it
  )
  log_info("  L%d  %6d mk  iters=%d  runtime=%.3fs  per_iter=%.4fs", level, mps, it, rt, rt / it)
}
cpp_thr <- do.call(rbind, thr)
orig_thr <- parse_orig_sweep()
scaling <- rbind(orig_thr[orig_thr$level %in% LEVELS, ], cpp_thr)
scaling <- scaling[order(scaling$core, scaling$mps), ]
write.csv(scaling, file.path(OUTDIR, "rtiger_marker_scaling.csv"), row.names = FALSE)
log_info("wrote rtiger_marker_scaling.csv (%d rows)", nrow(scaling))

# =============================================================================
# 2. CONVERGENCE: delta trace (cold-restart) + equivalence vs original oracle
# =============================================================================
log_info("--- convergence + equivalence (C++ vs original-Julia oracle) ---")
trace_rows <- list()
equiv_rows <- list()
t_all <- Sys.time()
for (li in seq_along(LEVELS)) {
  level <- LEVELS[li]
  obs <- build_obs(level)
  mps <- mps_of(level)
  orig_dump <- parse_orig_dump(level)
  orig_log <- parse_orig_log(level)

  # ---- C++ delta trace via cold restarts: fit(max_iter=m) for m=1..K --------
  # (skip the expensive full-panel trace in SMOKE)
  do_trace <- !(SMOKE && level == 0L)
  p_prev <- list(alpha = INIT_A, beta = INIT_B)
  cum <- 0
  last_fit <- NULL
  Kcpp <- NA_integer_
  if (do_trace) {
    t0 <- Sys.time()
    for (m in 1:50) {
      el <- system.time(fm <- fit_cpp(obs, m))[["elapsed"]]
      if (fm$iterations < m) {
        Kcpp <- fm$iterations
        break
      } # converged before m
      dlt <- max(abs(c(fm$alpha, fm$beta) - c(p_prev$alpha, p_prev$beta)))
      trace_rows[[length(trace_rows) + 1L]] <- data.frame(
        core = "cpp", level = level, mps = mps, iter = m, delta = dlt, per_iter = el - cum
      )
      cum <- el
      p_prev <- fm
      last_fit <- fm
      Kcpp <- m
    }
    log_info(
      "  L%d trace: K=%d iters  (%.1fs)", level, Kcpp,
      as.numeric(difftime(Sys.time(), t0, units = "secs"))
    )
  } else {
    last_fit <- fit_cpp(obs, 50L)
    Kcpp <- last_fit$iterations
  }

  # original per-iteration trace (from the preserved log)
  if (!is.null(orig_log)) {
    for (r in seq_len(nrow(orig_log))) {
      trace_rows[[length(trace_rows) + 1L]] <- data.frame(
        core = "orig", level = level, mps = mps, iter = orig_log$iter[r],
        delta = orig_log$delta[r], per_iter = orig_log$per_iter[r]
      )
    }
  }

  # ---- equivalence: decode C++ + compare vs original oracle -----------------
  paths <- nilHMM:::.rtiger_decode(obs, last_fit, r = RIG, postprocess = TRUE)
  tot <- 0L
  mis <- 0L
  for (s in names(obs)) {
    for (c in names(obs[[s]])) {
      key <- sprintf("vit[%s][%s]", s, c)
      a <- paths[[s]][[c]]
      b <- orig_dump$vit[[key]]
      if (is.null(b)) next
      n <- min(length(a), length(b))
      tot <- tot + n
      mis <- mis + sum(a[seq_len(n)] != b[seq_len(n)])
    }
  }
  pdiff <- max(abs(c(last_fit$alpha, last_fit$beta) - c(orig_dump$alpha, orig_dump$beta)))
  equiv_rows[[length(equiv_rows) + 1L]] <- data.frame(
    level = level, mps = mps, iters_cpp = Kcpp, iters_orig = orig_dump$iters,
    viterbi_total = tot, mismatches = mis, param_maxdiff = pdiff
  )
  log_info(
    "  L%d equiv: iters cpp/orig=%d/%d  Viterbi %d/%d match  mismatch=%d  paramdiff=%.2e",
    level, Kcpp, orig_dump$iters, tot - mis, tot, mis, pdiff
  )
}
conv_trace <- do.call(rbind, trace_rows)
equiv <- do.call(rbind, equiv_rows)
equiv <- equiv[order(-equiv$mps), ]
write.csv(conv_trace, file.path(OUTDIR, "rtiger_conv_trace.csv"), row.names = FALSE)
write.csv(equiv, file.path(OUTDIR, "rtiger_equivalence.csv"), row.names = FALSE)
log_info(
  "wrote rtiger_conv_trace.csv (%d rows) + rtiger_equivalence.csv (%d rows)",
  nrow(conv_trace), nrow(equiv)
)
log_info("convergence block: %.1f min", as.numeric(difftime(Sys.time(), t_all, units = "mins")))

# =============================================================================
# 3. FIGURES  (base R, mirroring docs/optimization.md)
# =============================================================================
COL_ORIG <- "firebrick"
COL_CPP <- "steelblue"
expo <- function(d) {
  f <- stats::lm(log10(per_iter) ~ log10(mps), d)
  coef(f)[2]
}

# ---- (a) marker scaling: original vs C++ -----------------------------------
op <- scaling[scaling$core == "orig", ]
op <- op[order(op$mps), ]
cp <- scaling[scaling$core == "cpp", ]
cp <- cp[order(cp$mps), ]
if (nrow(cp) >= 2 && nrow(op) >= 2) {
  bo <- expo(op)
  bc <- expo(cp)
  png(file.path(FIGDIR, "rtiger_marker_scaling_cpp.png"), width = 1100, height = 780, res = 130)
  par(mar = c(4.6, 5.4, 3.2, 1))
  plot(NA,
    log = "xy", xlim = c(5000, 115000), ylim = c(1e-3, 1000), xaxt = "n", yaxt = "n",
    xlab = "markers per sample", ylab = "per-iteration wall time (s)",
    main = "RTIGER scaling, shared panel: original (Julia) vs nilHMM C++/Rcpp"
  )
  yt <- c(1e-3, 1e-2, 0.1, 1, 10, 100, 1000)
  xt <- c(5000, 10000, 20000, 50000, 100000)
  abline(h = yt, v = xt, col = "grey92")
  axis(2, yt, c("0.001", "0.01", "0.1", "1", "10", "100", "1000"), las = 1)
  axis(1, xt, c("5000", "10000", "20000", "50000", "100000"))
  xx <- 10^seq(log10(5000), log10(115000), length = 100)
  fo <- stats::lm(log10(per_iter) ~ log10(mps), op)
  fc <- stats::lm(log10(per_iter) ~ log10(mps), cp)
  lines(xx, 10^predict(fo, data.frame(mps = xx)), col = COL_ORIG, lty = 2, lwd = 2)
  lines(xx, 10^predict(fc, data.frame(mps = xx)), col = COL_CPP, lty = 2, lwd = 2)
  points(op$mps, op$per_iter, pch = 17, col = COL_ORIG, cex = 1.3)
  points(cp$mps, cp$per_iter, pch = 19, col = COL_CPP, cex = 1.3)
  legend("topleft",
    bty = "n", cex = 0.9,
    legend = c(
      sprintf("original (Julia)  ~ markers^%.2f", bo),
      sprintf("nilHMM C++/Rcpp    ~ markers^%.2f", bc)
    ),
    col = c(COL_ORIG, COL_CPP), pch = c(17, 19), lty = 2, lwd = 2
  )
  dev.off()
  log_info("wrote rtiger_marker_scaling_cpp.png  (orig ^%.2f, cpp ^%.2f)", bo, bc)
}

# ---- (b) equivalence at the full panel (L0): time + trajectory -------------
tr0c <- conv_trace[conv_trace$core == "cpp" & conv_trace$level == 0L, ]
tr0o <- conv_trace[conv_trace$core == "orig" & conv_trace$level == 0L, ]
if (nrow(tr0c) >= 2 && nrow(tr0o) >= 2) {
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
  # right: convergence trajectory 1:1
  rng <- range(c(do, dc))
  cols <- rev(hcl.colors(n, "viridis"))
  plot(dc, do,
    log = "xy", pch = 19, col = cols, cex = 1.3, xlim = rng, ylim = rng,
    xaxt = "n", yaxt = "n", xlab = "nilHMM C++/Rcpp delta", ylab = "original Julia delta",
    main = sprintf("Convergence trajectory (per-iter delta)\n1:1 => identical (max rel %.0e)", drel)
  )
  axis(1, dt, dl)
  axis(2, dt, dl, las = 1)
  abline(0, 1, col = "grey50", lty = 2, lwd = 2)
  legend("topleft", c("early iterations", "late iterations"),
    pch = 19,
    col = c(cols[1], cols[n]), bty = "n", cex = 0.85
  )
  dev.off()
  log_info("wrote rtiger_equivalence_110K_cpp.png  (max rel delta %.1e)", drel)
} else {
  log_info("skipped equivalence_110K figure (need L0 trace; run without --smoke)")
}

log_info("DONE.")
