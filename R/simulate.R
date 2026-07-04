# simcross NIL truth generator + per-source degradation.
#
# Ported from the zealtiger pipeline (R/02_simulate.R, R/05_make_rtiger_input.R)
# into base/data.table (zealhmm has no tidyverse dep). Emits, for one design and
# source: the donor-mosaic TRUTH (common-schema segments), a marker grid, and
# degraded per-sample allele counts in the `read_counts()` 6-col layout.
#
# DESIGN (decided): BC2S2, n = 1500, maize consensus cM map (~1783 cM; the
# `maize_map_v5_clean.rds` staged under data/ref/), Stahl interference m = 10,
# p = 0. Consumed by analysis/02-simulation-calibration.qmd.
#
# Source-safe: top level defines functions only (no library()/stopifnot); the
# analysis notes source() all of R/ without side effects. Uses data.table + the
# kbroman/simcross package (Suggests).

SIM_DEFAULT_SEED <- 1L

# --- pedigree ---------------------------------------------------------------
# Founders 1 = recurrent (B73), 2 = donor (teosinte). Two backcrosses to 1, then
# selfing: id5 = BC2, id6 = BC2S1, id7 = BC2S2, id8 = BC2S3. A base data.frame
# is required (simcross::check_pedigree indexes it as `ped$col`).
.nil_pedigree <- function() {
  data.frame(
    id = 1:8,
    mom = c(0, 0, 1, 1, 1, 5, 6, 7),
    dad = c(0, 0, 2, 3, 4, 5, 6, 7),
    sex = c(0, 1, 0, 1, 0, 1, 0, 1),
    gen = c(0, 0, 1, 2, 3, 4, 5, 6)
  )
}
.nil_id_for <- function(design) {
  switch(design,
    BC2S2 = "7",
    BC2S3 = "8"
  )
}

# --- consensus map ----------------------------------------------------------
#' Load the staged consensus map (marker -> chr, cm, bp)
#' @export
load_consensus_map <- function(path = here::here("data/ref/maize_map_v5_clean.rds")) {
  if (!file.exists(path)) {
    stop("consensus map not staged: ", path, " (stage maize_map_v5_clean.rds)")
  }
  m <- data.table::as.data.table(readRDS(path))
  m <- m[, .(chr = as.integer(chr), cm = as.numeric(cm), bp = as.numeric(bp))]
  m[order(chr, cm)]
}

# per-chromosome cM length (the simcross `L` vector), ordered by chr
.chr_cm_lengths <- function(map) map[, .(L = max(cm)), by = chr][order(chr)]

# monotone bp -> cM interpolator for one chromosome (Hyman spline, clamped)
.bp_to_cm_fun <- function(chr_map) {
  d <- chr_map[, .(cm = mean(cm)), by = bp][order(bp)]
  f <- stats::splinefun(d$bp, d$cm, method = "hyman")
  rng <- range(d$bp)
  function(bp) f(pmin(pmax(bp, rng[1]), rng[2]))
}

#' Build a physical marker grid: n_markers spread by span, cM per marker
#' @export
build_marker_grid <- function(map, n_markers = 2500L, chr_prefix = "chr") {
  spans <- map[, .(min_bp = min(bp), max_bp = max(bp), span = max(bp) - min(bp)), by = chr]
  spans[, quota := pmax(2L, as.integer(round(n_markers * span / sum(span))))]
  data.table::rbindlist(lapply(seq_len(nrow(spans)), function(i) {
    ch <- spans$chr[i]
    bp <- round(seq(spans$min_bp[i], spans$max_bp[i], length.out = spans$quota[i]))
    f <- .bp_to_cm_fun(map[chr == ch])
    data.table::data.table(
      chr = as.integer(ch), chr_label = paste0(chr_prefix, ch),
      bp = as.integer(bp), cm = f(bp)
    )
  }))[order(chr, bp)]
}

# --- per-NIL true dosage (0/1/2) at markers ---------------------------------
.hap_allele_at <- function(hap, cm) {
  idx <- findInterval(cm, hap$locations, left.open = TRUE) + 1L
  idx[idx > length(hap$alleles)] <- length(hap$alleles)
  idx[idx < 1L] <- 1L
  hap$alleles[idx]
}

# one realization of the pedigree -> donor dosage (0/1/2) at each marker
.simulate_dosage <- function(ped, cmlen, markers, m, p, donor_allele = 2L, nil_id = "7") {
  sim <- simcross::sim_from_pedigree(ped, L = cmlen$L, m = m, p = p)
  dosage <- integer(nrow(markers))
  for (ci in seq_len(nrow(cmlen))) {
    ch <- cmlen$chr[ci]
    rows <- which(markers$chr == ch)
    if (!length(rows)) next
    nil <- sim[[ci]][[nil_id]]
    cm <- markers$cm[rows]
    dosage[rows] <- (.hap_allele_at(nil$mat, cm) == donor_allele) +
      (.hap_allele_at(nil$pat, cm) == donor_allele)
  }
  dosage
}

# RLE the per-marker dosage into common-schema truth segments (one sample)
.truth_segments <- function(markers, dosage, name) {
  dt <- data.table::data.table(
    chr = markers$chr, bp = markers$bp,
    state = as.integer(dosage)
  )[order(chr, bp)]
  dt[,
    {
      r <- rle(state)
      e <- cumsum(r$lengths)
      s <- c(1L, utils::head(e, -1L) + 1L)
      list(start_bp = bp[s], end_bp = bp[e], state = r$values)
    },
    by = chr
  ][, .(source = "sim", donor = NA_character_, name = name, chr, start_bp, end_bp, state)]
}

# --- degrade: draw REF/ALT counts under a source's coverage/error regime -----
# missing(lambda) = pi_floor + (1-pi_floor) e^{-k lambda}; present-site depth is
# 1 + Poisson(cond_mean - 1) so present markers carry >= 1 read (the low-lambda
# regime where a het often shows a single read and looks homozygous).
.source_regime <- function(source) {
  switch(source,
    skim   = list(lambda_mean = 0.43, shape = 8, pi_floor = 0, k_decay = 0.8, error = 0.01),
    brb    = list(lambda_mean = 1.0, shape = 2, pi_floor = 0.30, k_decay = 0.8, error = 0.02),
    target = list(lambda_mean = 20, shape = 8, pi_floor = 0, k_decay = 1.0, error = 0.002)
  )
}
.draw_counts <- function(dosage, lambda, pi_floor, k_decay, error) {
  n <- length(dosage)
  present_prob <- (1 - pi_floor) * (1 - exp(-k_decay * lambda))
  cond_mean <- lambda / present_prob
  present <- stats::runif(n) < present_prob
  depth <- integer(n)
  depth[present] <- 1L + stats::rpois(sum(present), max(cond_mean - 1, 0))
  p_alt <- c(0, 0.5, 1)[dosage + 1L] # donor allele = ALT
  p_eff <- p_alt * (1 - error) + (1 - p_alt) * error
  alt <- stats::rbinom(n, depth, p_eff)
  list(ref = as.integer(depth - alt), alt = as.integer(alt))
}

#' Simulate NILs for one design + degrade to one source
#'
#' Writes, under `outdir`:
#'   <source>.rds                 a compact bundle (see below), NOT per-sample files
#'   <design>_truth_segments.csv  (name, chr, start_bp, end_bp, state) — TRUTH
#'
#' The counts are stored as ONE `.rds` -- `list(grid = data.frame(chr, pos),
#' n_ref = <M x n integer matrix>, n_alt = <M x n integer matrix>, names, meta)` --
#' instead of `n` per-sample TSVs. The marker grid is identical across samples, so
#' storing it once (rather than repeating it `n` times) plus the sparse integer
#' count matrices is ~14x smaller and ~95x faster to load (380 MB / 38 s of gz
#' files -> ~27 MB / 0.4 s). Reconstruct the caller-ready long table with
#' [sim_counts()]; open the bundle with [load_sim()].
#'
#' @param design "BC2S2" (default) or "BC2S3".
#' @param source "skim" (default), "brb", or "target".
#' @param n Number of NILs (default 1500).
#' @param m,p Stahl interference (default m = 10, p = 0).
#' @param seed RNG seed (pinned; same seed => same truth across source calls).
#' @param n_markers Marker-grid size.
#' @param outdir,map_path Output dir and staged consensus map.
#' @param compress `saveRDS` compression ("gzip" default; "bzip2" is ~40% smaller
#'   but slower to read).
#' @return (invisibly) list(rds, truth, grid).
#' @export
simulate_source <- function(design = c("BC2S2", "BC2S3"),
                            source = c("skim", "brb", "target"),
                            n = 1500L, m = 10L, p = 0, seed = SIM_DEFAULT_SEED,
                            n_markers = 50000L,
                            outdir = here::here("results/sim"),
                            map_path = here::here("data/ref/maize_map_v5_clean.rds"),
                            compress = "gzip") {
  design <- match.arg(design)
  source <- match.arg(source)
  if (!requireNamespace("simcross", quietly = TRUE)) {
    stop("simulate_source() needs the 'simcross' package (kbroman/simcross).")
  }
  set.seed(seed)
  map <- load_consensus_map(map_path)
  cmlen <- .chr_cm_lengths(map)
  markers <- build_marker_grid(map, n_markers)
  ped <- .nil_pedigree()
  nid <- .nil_id_for(design)
  reg <- .source_regime(source)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  M <- nrow(markers)
  nms <- sprintf("sim%04d", seq_len(n))
  n_ref <- matrix(0L, M, n, dimnames = list(NULL, nms))
  n_alt <- matrix(0L, M, n, dimnames = list(NULL, nms))
  truth <- vector("list", n)
  for (i in seq_len(n)) {
    dosage <- .simulate_dosage(ped, cmlen, markers, m, p, nil_id = nid)
    truth[[i]] <- .truth_segments(markers, dosage, nms[i])
    lambda <- max(0.01, stats::rgamma(1, shape = reg$shape, scale = reg$lambda_mean / reg$shape))
    ac <- .draw_counts(dosage, lambda, reg$pi_floor, reg$k_decay, reg$error)
    n_ref[, i] <- as.integer(ac$ref)
    n_alt[, i] <- as.integer(ac$alt)
  }
  truth <- data.table::rbindlist(truth)
  grid <- data.frame(chr = as.integer(markers$chr), pos = as.integer(markers$bp))

  sim <- list(
    grid = grid, n_ref = n_ref, n_alt = n_alt, names = nms,
    source = source, design = design, n_markers = M, seed = seed,
    regime = reg
  )
  rds <- file.path(outdir, paste0(source, ".rds"))
  saveRDS(sim, rds, compress = compress)
  data.table::fwrite(truth, file.path(outdir, sprintf("%s_truth_segments.csv", tolower(design))))
  invisible(list(rds = rds, truth = truth, grid = grid))
}

#' Load a simulated-counts bundle written by [simulate_source()]
#' @param path Path to the `<source>.rds`.
#' @return The bundle list: `grid, n_ref, n_alt, names, source, design, ...`.
#' @export
load_sim <- function(path) {
  sim <- readRDS(path)
  if (!all(c("grid", "n_ref", "n_alt", "names") %in% names(sim))) {
    stop("load_sim(): '", path, "' is not a simulate_source() bundle")
  }
  sim
}

#' Reconstruct the caller-ready long count table from a sim bundle
#'
#' Expands the stored (grid + count matrices) into the `read_counts()`-style long
#' table `(name, chr, pos, n_ref, n_alt)` that the callers/metrics consume. Pass
#' `samples` to materialize only a subset (the grid is shared, so a calibration
#' subset costs only its own columns).
#'
#' @param sim A bundle from [load_sim()], or a path to the `.rds`.
#' @param samples Optional character vector of sample names (default: all).
#' @return data.frame `(name, chr, pos, n_ref, n_alt)`.
#' @export
sim_counts <- function(sim, samples = NULL) {
  if (is.character(sim)) sim <- load_sim(sim)
  idx <- if (is.null(samples)) seq_along(sim$names) else match(samples, sim$names)
  if (anyNA(idx)) {
    stop(
      "sim_counts(): unknown sample(s): ",
      paste(samples[is.na(idx)], collapse = ", ")
    )
  }
  M <- nrow(sim$grid)
  N <- length(idx)
  data.frame(
    name = rep(sim$names[idx], each = M),
    chr = rep(sim$grid$chr, N),
    pos = rep(sim$grid$pos, N),
    n_ref = as.integer(sim$n_ref[, idx, drop = FALSE]),
    n_alt = as.integer(sim$n_alt[, idx, drop = FALSE]),
    stringsAsFactors = FALSE
  )
}

# CLI: Rscript -e 'source("R/simulate.R"); simulate_source("BC2S2","skim")'
