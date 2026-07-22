# Nested-NIL simulator with REAL founder alleles (the count-from-parents model).
#
# This is the fidelity upgrade over R/simulate.R's `simulate_source()`, which
# hard-wires `donor = ALT everywhere` (p_alt = {0,0.5,1}[ancestry]) and so assumes
# every marker is informative (nir = 0). Here the observed allelic dose is the
# ANCESTRY dose filtered through the actual founder genotype:
#
#     observed_alt_dose = ancestry_dose * founder_allele,   founder_allele in {0,1}
#
# so wherever a NAM founder carries the B73/REF allele (founder = 0/0), the site is
# NON-INFORMATIVE regardless of ancestry -- the donor mosaic is invisible there.
# Non-informativeness (the nNIL `nir`/f0 ~= 0.59) thus EMERGES per-donor and
# spatially from the 24 real founders, instead of a scalar knob. Truth = the
# ancestry mosaic; the caller's job is to recover it despite the masked markers.
#
# Breeding design defaults to BC5S2 = the genotyped nNIL generation (Zhong 2025;
# f_1 = 1/128). Nested: each of the 24 founders is its own B73 x donor family,
# simulated apart, then unioned.
#
# Reuses R/simulate.R helpers (.bcsft_pedigree, .simulate_dosage, .draw_counts,
# .truth_segments, .source_regime). Founder input from
# scripts/nnil_foil/10_founder_genotypes.py (data/nnil_foil/founders_v5.csv).
#
#   Rscript scripts/simulate_nested_nils.R            # smoke (small n)
#   Rscript scripts/simulate_nested_nils.R --generate # full run

suppressMessages({
  library(data.table)
  library(here)
  library(jsonlite)
})
source(here::here("R/simulate.R")) # defines the reused helpers (side-effect free)

#' Simulate nested NILs with real founder alleles
#'
#' @param founders_path founders_v5.csv (founder x v5-marker genotype {0,1,2,3}).
#' @param markers_path   markers_v5.tsv (marker, chr, pos_v5, cm) for the grid/cM.
#' @param design         breeding design string; default "BC5S2" (genotyped gen).
#' @param n_per_family   NILs simulated per founder family.
#' @param source         degradation regime ("skim"/"brb"/"target"); the nNIL-GBS
#'                        regime fit is step 2, so this stands in for now.
#' @param het_policy     founder het (0/1) handling: "fix_per_family" (draw the
#'                        donor haplotype allele once per family, default),
#'                        "noninformative" (treat as REF), or "drop" (mask).
#' @param m,p            Stahl interference (default m = 10, p = 0).
#' @param seed           RNG seed.
#' @param outdir         output directory.
#' @return (invisibly) list(rds, truth, grid, founder_f0).
simulate_nested_nils <- function(founders_path = here::here("data/nnil_foil/founders_v5.csv"),
                                 markers_path = here::here("data/nnil_foil/markers_v5.tsv"),
                                 design = "BC5S4", # Zhong 2025 Methods: "backcrossed to B73 for five generations (BC5) and ... self-pollinated for four generations (BC5F4)" -> 4 selfings
                                 n_per_family = 40L,
                                 source = c("nnil_gbs", "skim", "brb", "target"),
                                 regime_path = here::here("data/nnil_foil/nnil_gbs_regime.json"),
                                 het_policy = c("fix_per_family", "noninformative", "drop"),
                                 m = 10L, p = 0, seed = 1L,
                                 outdir = here::here("results/sim/nested_nnil")) {
  source <- match.arg(source)
  het_policy <- match.arg(het_policy)
  if (!requireNamespace("simcross", quietly = TRUE)) {
    stop("simulate_nested_nils() needs the 'simcross' package (kbroman/simcross).")
  }
  set.seed(seed)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # ---- founders + marker grid (co-registered on the v5 marker names) --------
  fmat <- as.matrix(data.table::fread(founders_path), rownames = 1) # founder x marker {0,1,2,3}
  mv5 <- data.table::fread(markers_path) # marker, marker_v4, chr, pos_v4, pos_v5, cm
  mk <- data.table::data.table(marker = colnames(fmat))
  mk <- merge(mk, mv5[, .(marker, chr = as.integer(chr), bp = as.integer(pos_v5), cm = as.numeric(cm))],
    by = "marker", sort = FALSE
  )
  # keep only markers with a cM (grid the sim runs on) and reorder both together
  mk <- mk[!is.na(cm) & !is.na(chr)][order(chr, bp)]
  fmat <- fmat[, mk$marker, drop = FALSE]
  markers <- data.table::data.table(chr = mk$chr, bp = mk$bp, cm = mk$cm)
  cmlen <- markers[, .(L = max(cm)), by = chr][order(chr)]
  M <- nrow(markers)
  message(sprintf(
    "[nested] %d founders x %d markers on %d chr; design %s, %d/family, source %s, het=%s",
    nrow(fmat), M, nrow(cmlen), design, n_per_family, source, het_policy
  ))

  # ---- pedigree for the design ----------------------------------------------
  pd <- nilHMM::parse_design(design)
  bp <- .bcsft_pedigree(pd$n_bc, pd$n_self)
  ped <- bp$ped
  nid <- bp$nil_id
  # nnil_gbs = the mr-fitted regime (scripts/nnil_foil/11_fit_gbs_regime.R); the
  # others are R/simulate.R's hand-set skim/brb/target.
  reg <- if (source == "nnil_gbs") {
    j <- jsonlite::fromJSON(regime_path)
    list(
      lambda_mean = j$lambda_mean, shape = j$shape, pi_floor = j$pi_floor,
      k_decay = j$k_decay, error = j$error
    )
  } else {
    .source_regime(source)
  }

  # ---- per-founder allele vector (donor haplotype ALT-carrying, {0,1,NA}) ----
  # g = 0 -> 0 (REF = non-informative); g = 2 -> 1 (ALT = informative);
  # g = 3 (founder missing) -> 0 (treat as non-informative: the founder allele is
  # unknown, but a missing FOUNDER call is not a missing NIL genotype, so the site
  # still emits reads rather than dropping out); g = 1 (het) -> per het_policy.
  # Only het under het_policy = "drop" yields NA (a genuinely masked marker).
  founder_allele_vec <- function(g) {
    fa <- rep(0, length(g)) # default non-informative (covers g == 0 and g == 3)
    fa[g == 2L] <- 1
    if (het_policy == "noninformative") {
      fa[g == 1L] <- 0
    } else if (het_policy == "fix_per_family") {
      h <- which(g == 1L)
      if (length(h)) fa[h] <- stats::rbinom(length(h), 1L, 0.5)
    } else { # "drop": mask het markers
      fa[g == 1L] <- NA_real_
    }
    fa
  }

  fnames <- rownames(fmat)
  n_total <- nrow(fmat) * n_per_family
  n_ref <- matrix(0L, M, n_total)
  n_alt <- matrix(0L, M, n_total)
  smp_names <- character(n_total)
  smp_founder <- character(n_total)
  truth <- vector("list", n_total)
  f0_obs <- numeric(nrow(fmat)) # realized non-informative rate per family

  col <- 0L
  for (fi in seq_along(fnames)) {
    fa <- founder_allele_vec(as.integer(fmat[fi, ]))
    f0_obs[fi] <- mean(fa == 0, na.rm = TRUE)
    informative <- !is.na(fa) # markers where the donor mosaic can show
    for (k in seq_len(n_per_family)) {
      col <- col + 1L
      a <- .simulate_dosage(ped, cmlen, markers, m, p, nil_id = nid) # ancestry 0/1/2
      truth[[col]] <- .truth_segments(markers, a, sprintf("%s_%04d", fnames[fi], k))
      obs <- a * ifelse(informative, fa, 0) # observed ALT dose {0,1,2}; masked -> 0
      lambda <- max(0.01, stats::rgamma(1, shape = reg$shape, scale = reg$lambda_mean / reg$shape))
      ac <- .draw_counts(obs, lambda, reg$pi_floor, reg$k_decay, reg$error)
      ref <- as.integer(ac$ref)
      alt <- as.integer(ac$alt)
      ref[!informative] <- 0L # masked markers carry no reads (founder het-drop/missing)
      alt[!informative] <- 0L
      n_ref[, col] <- ref
      n_alt[, col] <- alt
      smp_names[col] <- sprintf("%s_%04d", fnames[fi], k)
      smp_founder[col] <- fnames[fi]
    }
  }
  colnames(n_ref) <- colnames(n_alt) <- smp_names
  truth <- data.table::rbindlist(truth)
  grid <- data.frame(chr = markers$chr, pos = markers$bp)

  sim <- list(
    grid = grid, n_ref = n_ref, n_alt = n_alt, names = smp_names,
    founder = smp_founder, source = source, design = design, n_markers = M,
    seed = seed, het_policy = het_policy, regime = reg,
    founder_f0 = data.frame(founder = fnames, f0 = f0_obs)
  )
  rds <- file.path(outdir, sprintf("nested_nnil_%s_%s.rds", tolower(design), source))
  saveRDS(sim, rds, compress = "gzip")
  data.table::fwrite(truth, file.path(outdir, sprintf("nested_nnil_%s_truth_segments.csv", tolower(design))))
  sim_mr <- mean((n_ref + n_alt) == 0L) # depth-0 = missing call
  message(sprintf(
    "[nested] wrote %s : %d markers x %d NILs (%d families); realized f0 mean %.3f [%.3f, %.3f]",
    basename(rds), M, n_total, nrow(fmat), mean(f0_obs), min(f0_obs), max(f0_obs)
  ))
  message(sprintf(
    "[nested] simulated missing rate = %.4f  (real nNIL GBS mr = 0.0928)", sim_mr
  ))
  invisible(list(rds = rds, truth = truth, grid = grid, founder_f0 = sim$founder_f0))
}

# --- CLI --------------------------------------------------------------------
if (sys.nframe() == 0L) {
  smoke <- !("--generate" %in% commandArgs(trailingOnly = TRUE))
  if (smoke) {
    message("[nested] SMOKE run (5/family). Pass --generate for the full run.")
    simulate_nested_nils(n_per_family = 5L)
  } else {
    simulate_nested_nils(n_per_family = 40L)
  }
}
