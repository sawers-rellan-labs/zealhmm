# simcross NIL truth generator + per-source degradation.
#
# Ported from the zealtiger pipeline (R/02_simulate.R, R/05_make_rtiger_input.R)
# into base/data.table (zealhmm has no tidyverse dep). Emits, for one design and
# source: the donor-mosaic TRUTH (common-schema segments), a marker grid, and
# degraded per-sample allele counts in the `read_counts()` 6-col layout.
#
# DESIGN (decided): BC2S2, n = 1500, maize consensus cM map (~1783 cM; the
# `maize_map_v5_clean.rds` staged under data/ref/), Stahl interference m = 10,
# p = 0. Consumed by analysis/simulation-calibration.qmd.
#
# Source-safe: top level defines functions only (no library()/stopifnot); the
# analysis notes source() all of R/ without side effects. Uses data.table + the
# kbroman/simcross package (Suggests).

SIM_DEFAULT_SEED <- 1L

# --- pedigree ---------------------------------------------------------------
# Founders 1 = recurrent (B73), 2 = donor (teosinte). F1 = 1 x 2; then `n_bc`
# backcrosses to the recurrent (1); then `n_self` generations of selfing. Builds
# the pedigree for ANY BC_n S_m design (BC2S2 = the old hardcoded 8-row table;
# BC1S4 = TeoNAM). For anything else simcross can simulate (RIL, AIL, MAGIC),
# pass a ready `pedigree` + `nil_id` to simulate_source() directly. simcross
# selfing is sex-agnostic (mom == dad); backcrosses use recurrent(sex 0) x
# prev(sex 1). A base data.frame is required (simcross::check_pedigree indexes it
# as `ped$col`).
.bcsft_pedigree <- function(n_bc, n_self) {
  id <- c(1L, 2L)
  mom <- c(0L, 0L)
  dad <- c(0L, 0L)
  gen <- c(0L, 0L)
  add <- function(mm, dd, g) { # append one progeny, return its id
    i <- length(id) + 1L
    id[i] <<- i
    mom[i] <<- mm
    dad[i] <<- dd
    gen[i] <<- g
    i
  }
  cur <- add(1L, 2L, 1L) # F1
  for (k in seq_len(n_bc)) cur <- add(1L, cur, 1L + k) # BCk = recurrent x prev
  g <- 1L + n_bc
  for (j in seq_len(n_self)) {
    g <- g + 1L
    cur <- add(cur, cur, g)
  } # self x n_self
  # sex by role: anyone ever used as a `mom` is female (0), else male (1). sex is
  # cosmetic for our use -- sim_from_pedigree() takes the maternal gamete from
  # `mom` and paternal from `dad` regardless of sex, and backcrosses always put
  # the recurrent (1) as mom, so dosage is correct. (simcross::check_pedigree
  # cannot validate selfing at all: a selfed individual is both mom and dad, so it
  # trips either "moms male" or "dads female" -- which is why the original never
  # called it. sim_from_pedigree is the workhorse and handles mom == dad selfing.)
  sex <- ifelse(id %in% mom[mom > 0L], 0L, 1L)
  list(
    ped = data.frame(id = id, mom = mom, dad = dad, sex = sex, gen = gen),
    nil_id = as.character(cur)
  )
}

# "BC1S4" -> list(n_bc = 1, n_self = 4)
.parse_design <- function(design) {
  mm <- regmatches(design, regexec("^BC(\\d+)S(\\d+)$", design))[[1]]
  if (!length(mm)) {
    stop(".parse_design(): design must be 'BC<n>S<m>' (e.g. BC1S4, BC2S2): ", design)
  }
  list(n_bc = as.integer(mm[2]), n_self = as.integer(mm[3]))
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
  # span as double: on bp-scale maize maps `n_markers * span` and `sum(span)`
  # overflow 32-bit integer arithmetic (span ~3e8, sum ~2e9 > .Machine$integer.max),
  # yielding NA quotas. Double arithmetic keeps the quota computation exact.
  spans <- map[, .(min_bp = min(bp), max_bp = max(bp), span = as.numeric(max(bp) - min(bp))), by = chr]
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
#' @param design Breeding design as `"BC<n>S<m>"` (e.g. `"BC2S2"`, `"BC1S4"`):
#'   F1, then n backcrosses to the recurrent, then m selfs. Ignored if `pedigree`
#'   is supplied.
#' @param pedigree,nil_id Optional escape hatch for any simcross-simulable design:
#'   a ready pedigree data.frame (`id, mom, dad, sex, gen`) and the id of the
#'   individual to sample. When given, `design` is used only to name the output.
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
simulate_source <- function(design = "BC2S2",
                            source = c("skim", "brb", "target"),
                            n = 1500L, m = 10L, p = 0, seed = SIM_DEFAULT_SEED,
                            n_markers = 50000L,
                            pedigree = NULL, nil_id = NULL,
                            outdir = here::here("results/sim"),
                            map_path = here::here("data/ref/maize_map_v5_clean.rds"),
                            compress = "gzip") {
  source <- match.arg(source)
  if (!requireNamespace("simcross", quietly = TRUE)) {
    stop("simulate_source() needs the 'simcross' package (kbroman/simcross).")
  }
  set.seed(seed)
  map <- load_consensus_map(map_path)
  cmlen <- .chr_cm_lengths(map)
  markers <- build_marker_grid(map, n_markers)
  if (!is.null(pedigree)) { # custom-pedigree escape hatch
    if (is.null(nil_id)) stop("simulate_source(): supply `nil_id` with a custom `pedigree`.")
    ped <- pedigree
    nid <- as.character(nil_id)
  } else {
    pd <- .parse_design(design)
    bp <- .bcsft_pedigree(pd$n_bc, pd$n_self)
    ped <- bp$ped
    nid <- bp$nil_id
  }
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
