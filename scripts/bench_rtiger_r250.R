#!/usr/bin/env Rscript
# RTIGER benchmark at the DATASET OPERATING RIGIDITY r=250 (the r=2 companion is
# bench_rtiger_cpp_vs_julia.R). The nilHMM C++/Rcpp core is measured at all five
# thinned sizes; the upstream-original Julia is measured at the small sizes only
# (levels 4,3,2 via scripts/rtiger_julia_conv_worker.jl -> results/bench/orig_conv_r250/)
# and its large sizes are PROJECTED from the power-law fit (done in the plotting
# script). Because the rigidity is collapsed (window product, not s*r enumeration),
# the C++ per-iteration cost is ~r-independent, so r=250 is cheap for it.
#   Rscript scripts/bench_rtiger_r250.R
suppressMessages(library(nilHMM))
ROOT <- "/Users/fvrodriguez/repos/zealhmm"
source(file.path(ROOT, "scripts/logging.R"))
PANELDIR <- file.path(ROOT, "data/rtiger_shared3_input")
ORIGDIR <- file.path(ROOT, "results/bench/orig_conv_r250")
OUTDIR <- file.path(ROOT, "results/bench")
files <- c(S1 = "sampleBN.txt", S2 = "sampleZ.txt", S3 = "sampleAU.txt")
chrs <- paste0("Chr", 1:5)
RIG <- 250L
EPS <- 0.01
INIT_A <- c(20, 20, 1)
INIT_B <- c(1, 20, 20)
REPS <- 3L
LEVELS <- 0:4
ORIG_LEVELS <- c(4L, 3L, 2L) # the small sizes we actually run the original at

# Each benchmark reads its PRE-SPLIT panel (thin_L<level>/, written by
# scripts/materialize_thinned_panel.R) and does NO in-script thinning: it loads only
# the dataset it is benchmarking, never the full panel.
build_obs_from <- function(dir) {
  raw <- lapply(files, function(f) {
    read.table(file.path(dir, f),
      sep = "\t", col.names = c("chr", "pos", "refA", "refC", "altA", "altF"),
      stringsAsFactors = FALSE
    )
  })
  obs <- list()
  for (s in names(files)) {
    d <- raw[[s]]
    obs[[s]] <- list()
    for (cc in chrs) {
      rows <- which(d$chr == cc)
      if (!length(rows)) next
      k <- as.integer(round(d$refC[rows]))
      n <- k + as.integer(round(d$altF[rows]))
      obs[[s]][[cc]] <- list(k = k, n = n)
    }
  }
  obs
}
fit_cpp <- function(obs, max_iter) {
  nilHMM:::.rtiger_fit(obs,
    r = RIG, nstates = 3L, eps = EPS, max_iter = as.integer(max_iter),
    init_alpha = INIT_A, init_beta = INIT_B
  )
}

# ---- original r=250 reference parsers (results/bench/orig_conv_r250/) --------
parse_orig_dump <- function(level) {
  p <- file.path(ORIGDIR, sprintf("conv_orig_L%d.txt", level))
  if (!file.exists(p)) {
    return(NULL)
  }
  L <- readLines(p)
  num <- function(key) {
    as.numeric(strsplit(gsub(
      ".*=\\[|\\].*", "", grep(paste0("^", key, "="), L, value = TRUE)[1]
    ), ",")[[1]])
  }
  meta <- L[1]
  vit <- list()
  for (ln in grep("^vit\\[", L, value = TRUE)) {
    key <- sub("=.*", "", ln)
    vit[[key]] <- as.integer(strsplit(gsub(".*=\\[|\\].*", "", ln), ",")[[1]])
  }
  list(
    iters = as.integer(sub(".*iters=([0-9]+).*", "\\1", meta)),
    runtime = as.numeric(sub(".*runtime=([0-9.]+).*", "\\1", meta)),
    alpha = num("alpha"), beta = num("beta"), vit = vit
  )
}
parse_orig_log <- function(level) {
  p <- file.path(ORIGDIR, sprintf("log_orig_conv_L%d.log", level))
  if (!file.exists(p)) {
    return(NULL)
  }
  ln <- grep("^iter ", readLines(p), value = TRUE)
  if (!length(ln)) {
    return(NULL)
  }
  g <- function(f) as.numeric(sub(paste0(".*", f, "=([0-9.eE+-]+).*"), "\\1", ln))
  el <- g("elapsed")
  data.frame(iter = seq_along(ln), delta = g("delta"), per_iter = diff(c(0, el)))
}

trace_rows <- list()
speed_rows <- list()
equiv_rows <- list()

for (level in LEVELS) {
  obs <- build_obs_from(file.path(PANELDIR, sprintf("thin_L%d", level)))
  mps <- sum(vapply(obs[[1]], function(x) length(x$k), 0L))

  # C++ per-iteration delta trace via cold restarts (fit(max_iter=m), m=1..K)
  p_prev <- list(alpha = INIT_A, beta = INIT_B)
  cum <- 0
  last_fit <- NULL
  Kcpp <- NA_integer_
  for (m in 1:50) {
    el <- system.time(fm <- fit_cpp(obs, m))[["elapsed"]]
    if (fm$iterations < m) {
      Kcpp <- fm$iterations
      break
    }
    dlt <- max(abs(c(fm$alpha, fm$beta) - c(p_prev$alpha, p_prev$beta)))
    trace_rows[[length(trace_rows) + 1L]] <- data.frame(
      core = "cpp", level = level, mps = mps, iter = m, delta = dlt, per_iter = el - cum
    )
    cum <- el
    p_prev <- fm
    last_fit <- fm
    Kcpp <- m
  }
  if (is.null(last_fit)) last_fit <- fit_cpp(obs, 50L)

  # C++ full-fit throughput: median of REPS runs to convergence
  tt <- vapply(seq_len(REPS), function(i) system.time(fit_cpp(obs, 50L))[["elapsed"]], numeric(1))
  cpp_total <- stats::median(tt)
  speed_rows[[length(speed_rows) + 1L]] <- data.frame(
    core = "cpp", level = level, mps = mps, iters = Kcpp,
    total = cpp_total, per_iter = cpp_total / Kcpp
  )
  log_info("cpp L%d %d markers: %d iters, %.3fs total (%.4fs/iter)", level, mps, Kcpp, cpp_total, cpp_total / Kcpp)

  # original r=250 (small sizes only): throughput + per-iter trace + equivalence
  od <- parse_orig_dump(level)
  ol <- parse_orig_log(level)
  if (!is.null(od)) {
    speed_rows[[length(speed_rows) + 1L]] <- data.frame(
      core = "orig", level = level, mps = mps, iters = od$iters,
      total = od$runtime, per_iter = od$runtime / od$iters
    )
    paths <- nilHMM:::.rtiger_decode(obs, last_fit, r = RIG, postprocess = TRUE)
    tot <- 0L
    mis <- 0L
    for (s in names(obs)) {
      for (cc in names(obs[[s]])) {
        b <- od$vit[[sprintf("vit[%s][%s]", s, cc)]]
        if (is.null(b)) next
        a <- paths[[s]][[cc]]
        n <- min(length(a), length(b))
        tot <- tot + n
        mis <- mis + sum(a[seq_len(n)] != b[seq_len(n)])
      }
    }
    pdiff <- max(abs(c(last_fit$alpha, last_fit$beta) - c(od$alpha, od$beta)))
    equiv_rows[[length(equiv_rows) + 1L]] <- data.frame(
      level = level, mps = mps, iters_cpp = Kcpp, iters_orig = od$iters,
      viterbi_total = tot, mismatches = mis, param_maxdiff = pdiff
    )
    log_info(
      "  L%d equiv (r=250): iters cpp/orig=%d/%d  Viterbi %d/%d match  mism=%d  paramdiff=%.2e",
      level, Kcpp, od$iters, tot - mis, tot, mis, pdiff
    )
  }
  if (!is.null(ol)) {
    ol$core <- "orig"
    ol$level <- level
    ol$mps <- mps
    trace_rows[[length(trace_rows) + 1L]] <- ol[, c("core", "level", "mps", "iter", "delta", "per_iter")]
  }
}

conv_trace <- do.call(rbind, trace_rows)
speed <- do.call(rbind, speed_rows)
equiv <- do.call(rbind, equiv_rows)
write.csv(conv_trace, file.path(OUTDIR, "rtiger_conv_trace_r250.csv"), row.names = FALSE)
write.csv(speed, file.path(OUTDIR, "rtiger_marker_scaling_r250.csv"), row.names = FALSE)
write.csv(equiv, file.path(OUTDIR, "rtiger_equiv_r250.csv"), row.names = FALSE)
log_info("wrote rtiger_{conv_trace,marker_scaling,equiv}_r250.csv")
cat("\n== speed ==\n")
print(speed, row.names = FALSE)
cat("\n== equivalence ==\n")
print(equiv, row.names = FALSE)
