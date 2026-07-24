#!/usr/bin/env Rscript
# ZEAL simple BC2S3 single-individual simulator (calibration baseline; NO pooling).
#
# Simplification of scripts/simulate_zeal_pool.R: each line is ONE BC2S3 individual
# (integer donor-ancestry dose 0/1/2), not a 6-sib pool. Only the pooling is removed;
# the HET DONOR SIMULATION is kept.
#
# HET DONOR SIMULATION: teosinte donors are ~20% heterozygous, so the inbred
# `observed = ancestry_dose * founder_allele` shortcut is not exact. All donor ancestry
# in a line descends from ONE founding donor->F1 gamete, so at each heterozygous donor
# site that gamete carried ALT or REF, decided once per line by a Bernoulli(0.5): ALT ->
# the site is informative (behaves like a hom-ALT donor), REF -> silent (non-informative).
# Averaged over het sites this gives the authentic per-donor non-informative rate ~0.68
# (f_REF + 0.5 f_HET). The recurrent parent is B73 (inbred REF); the donor is a teosinte
# reference individual drawn with replacement as a taxon-matched stand-in for the real
# accession (we lack the accession's own genotypes; snp50k_genotype_identifiability.qmd s2.6).
#
# Truth = the individual's latent ancestry mosaic (integer dose); non-informativeness is
# observation masking, not truth. Skim counts drawn under the SNP50K coverage regime; the
# marker grid is the full SNP50K native-cM set (NO thinning).
#
#   Rscript scripts/simulate_zeal_nil.R            # SMOKE (5/family, ~2 families/taxon)
#   Rscript scripts/simulate_zeal_nil.R --generate # full cohort

suppressMessages({
  library(data.table)
  library(here)
  library(nilHMM) # breeding_prior, parse_design, call_gt, call_ancestry
})
source(here::here("R/simulate.R")) # .bcsft_pedigree, .simulate_dosage, .truth_segments, .hap_allele_at
source(here::here("scripts/logging.R")) # log_info

# ---- SNP50K teosinte reference-panel genotypes -> {0,1,2,3} matrix (cached) ----
# (shared with scripts/simulate_zeal_pool.R)
load_ref_panel_gt <- function(vcf = here::here("data/zeal/reference_panel/bzea_50K_ref_panel.vcf.gz"),
                              markers,
                              cache = here::here("data/zeal/reference_panel/ref_panel_gt.rds")) {
  if (file.exists(cache)) {
    g <- readRDS(cache)
    if (nrow(g$mat) == nrow(markers)) {
      return(g)
    }
  }
  tmp <- tempfile(fileext = ".tsv")
  system2("bcftools", c("query", "-f", "'%CHROM\\t%POS[\\t%GT]\\n'", vcf), stdout = tmp)
  samples <- system2("bcftools", c("query", "-l", vcf), stdout = TRUE)
  raw <- fread(tmp, header = FALSE)
  setnames(raw, c("chr", "pos", samples))
  raw[, chr := as.integer(sub("^chr", "", chr))]
  recode <- function(x) {
    out <- rep(3L, length(x))
    out[x %in% c("0/0", "0|0")] <- 0L
    out[x %in% c("0/1", "1/0", "0|1", "1|0")] <- 1L
    out[x %in% c("1/1", "1|1")] <- 2L
    out
  }
  gmat <- as.matrix(raw[, lapply(.SD, recode), .SDcols = samples])
  idx <- match(paste(markers$chr, markers$bp, sep = "_"), paste(raw$chr, raw$pos, sep = "_"))
  out_mat <- matrix(3L, nrow(markers), length(samples), dimnames = list(NULL, samples))
  ok <- !is.na(idx)
  out_mat[ok, ] <- gmat[idx[ok], ]
  g <- list(mat = out_mat, samples = samples, n_absent = sum(!ok))
  saveRDS(g, cache, compress = "gzip")
  g
}

# ---- SNP50K skim coverage regime (shared with simulate_zeal_pool.R) ------------
fit_snp50k_regime <- function(path = here::here("data/missing_data/snp50k_per_sample.tsv"),
                              pi_floor = 0.01, k_decay = 1.0, error = 0.01) {
  d <- fread(path)
  lam <- d$lambda[is.finite(d$lambda) & d$lambda > 0]
  list(
    lambda_mean = mean(lam), shape = mean(lam)^2 / stats::var(lam),
    pi_floor = pi_floor, k_decay = k_decay, error = error,
    missing_mean = mean(d$missing_obs, na.rm = TRUE)
  )
}

# ---- skim counts from a continuous observed ALT fraction ----------------------
.draw_counts_pfrac <- function(p_alt, lambda, pi_floor, k_decay, error) {
  n <- length(p_alt)
  present_prob <- (1 - pi_floor) * (1 - exp(-k_decay * lambda))
  cond_mean <- lambda / present_prob
  present <- stats::runif(n) < present_prob
  depth <- integer(n)
  depth[present] <- 1L + stats::rpois(sum(present), max(cond_mean - 1, 0))
  p_eff <- p_alt * (1 - error) + (1 - p_alt) * error
  alt <- stats::rbinom(n, depth, p_eff)
  list(ref = as.integer(depth - alt), alt = as.integer(alt))
}

#' Simulate the ZEAL BC2S3 single-individual cohort (het donor simulation, no pool)
simulate_zeal_nil <- function(design = "BC2S3",
                              n_per_family = 5L, n_families = NULL,
                              taxa = c("Zx", "Zv", "Zd", "Zl", "Zh"),
                              m = 10L, p = 0, seed = 1L,
                              markers_path = here::here("data/zeal/markers_snp50k_cm.tsv"),
                              cohort_meta = here::here("data/zeal/reference_panel/bzea_50K_cohort_ref_metadata.csv"),
                              ref_ids = here::here("data/zeal/reference_panel/reference_ids_by_taxon.csv"),
                              outdir = here::here("results/sim/zeal_nil"), tag = "smoke") {
  if (!requireNamespace("simcross", quietly = TRUE)) stop("needs the 'simcross' package.")
  set.seed(seed)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  t0 <- Sys.time()

  # marker grid + cM (full SNP50K, no thinning); map-derived r
  mk <- fread(markers_path)
  markers <- data.table(chr = as.integer(mk$chr), bp = as.integer(mk$pos), cm = as.numeric(mk$cm))[order(chr, bp)]
  cmlen <- markers[, .(L = max(cm)), by = chr][order(chr)]
  M <- nrow(markers)
  map_r <- 2 * markers[, sum(tapply(cm, chr, function(x) max(x) - min(x)))] / (100 * M)

  # families = cohort donor accessions x taxon-matched teosinte reference stand-in
  cm_meta <- fread(cohort_meta)
  fams <- unique(cm_meta[
    is_cohort == TRUE & !is.na(donor_accession) & maizegdb_prefix %in% taxa,
    .(family = donor_accession, taxon = maizegdb_prefix)
  ])
  refs <- fread(ref_ids)[maizegdb_prefix %in% taxa, .(sample, taxon = maizegdb_prefix)]
  if (!is.null(n_families)) {
    per_tax <- max(1L, ceiling(n_families / length(taxa)))
    fams <- fams[, head(.SD, per_tax), by = taxon][seq_len(min(.N, n_families))]
  }
  fams[, stand_in := {
    pool <- refs[taxon == .BY$taxon]$sample
    if (!length(pool)) stop("no reference stand-in for taxon ", .BY$taxon)
    sample(pool, .N, replace = TRUE)
  }, by = taxon]

  gp <- load_ref_panel_gt(markers = markers)
  reg <- fit_snp50k_regime()
  pd <- parse_design(design)
  ped <- .bcsft_pedigree(pd$n_bc, pd$n_self) # SINGLE BC2S3 individual (no pool)
  log_info(
    "grid: %d SNP50K markers; map r=%.3e | %d families x %d taxa | %d lines/family (single BC2S3 individual)",
    M, map_r, nrow(fams), uniqueN(fams$taxon), n_per_family
  )

  n_lines <- nrow(fams) * n_per_family
  n_ref <- matrix(0L, M, n_lines)
  n_alt <- matrix(0L, M, n_lines)
  dose_mat <- matrix(0L, M, n_lines)
  gtrue_mat <- matrix(0L, M, n_lines) # level-2 true genotype g in {0,1,2} = dose * founder_allele
  line_names <- character(n_lines)
  line_family <- character(n_lines)
  truth_l <- vector("list", n_lines)
  col <- 0L
  for (fi in seq_len(nrow(fams))) {
    fa <- fams[fi]
    g_vec <- gp$mat[, fa$stand_in] # teosinte stand-in donor genotype {0,1,2,3}
    for (k in seq_len(n_per_family)) {
      col <- col + 1L
      nm <- sprintf("%s_%03d", fa$family, k)
      dose <- .simulate_dosage(ped$ped, cmlen, markers, m, p, nil_id = ped$nil_id) # M-vector {0,1,2}
      # HET DONOR SIMULATION: founding-gamete allele, resolved once per line
      founder_allele <- integer(M)
      founder_allele[g_vec == 2L] <- 1L
      het <- which(g_vec == 1L)
      if (length(het)) founder_allele[het] <- stats::rbinom(length(het), 1L, 0.5)
      p_alt <- (dose / 2) * founder_allele # observed ALT fraction (0 where non-informative)
      gtrue_mat[, col] <- as.integer(dose * founder_allele) # true genotype g in {0,1,2} (no RNG consumed)
      lambda <- max(0.01, stats::rgamma(1, shape = reg$shape, scale = reg$lambda_mean / reg$shape))
      ac <- .draw_counts_pfrac(p_alt, lambda, reg$pi_floor, reg$k_decay, reg$error)
      n_ref[, col] <- ac$ref
      n_alt[, col] <- ac$alt
      dose_mat[, col] <- dose
      line_names[col] <- nm
      line_family[col] <- fa$family
      truth_l[[col]] <- .truth_segments(markers, dose, nm) # latent ancestry mosaic
    }
  }
  colnames(n_ref) <- colnames(n_alt) <- colnames(dose_mat) <- colnames(gtrue_mat) <- line_names

  sim <- list(
    grid = data.frame(chr = markers$chr, pos = markers$bp), n_ref = n_ref, n_alt = n_alt,
    dose = dose_mat, g_true = gtrue_mat, names = line_names, family = line_family, truth = rbindlist(truth_l),
    design = design, n_markers = M, seed = seed, regime = reg, map_r = map_r,
    donor_model = "het_donor_simulation", pooling = FALSE
  )
  rds <- file.path(outdir, sprintf("zeal_nil_%s_%s.rds", tolower(design), tag))
  saveRDS(sim, rds, compress = "gzip")
  fwrite(fams, file.path(outdir, sprintf("zeal_nil_%s_%s_families.csv", tolower(design), tag)))

  sim_mr <- mean((n_ref + n_alt) == 0L)
  log_info(
    "wrote %s : %d markers x %d lines (%d families); missing=%.3f (target ~%.3f); donor sites=%.3f; mean dose=%.3f",
    basename(rds), M, n_lines, nrow(fams), sim_mr, reg$missing_mean, mean(dose_mat > 0), mean(dose_mat)
  )
  log_info("total %.1fs", as.numeric(difftime(Sys.time(), t0, units = "secs")))
  invisible(list(rds = rds, sim = sim, families = fams, markers = markers, map_r = map_r))
}

if (sys.nframe() == 0L) {
  full <- "--generate" %in% commandArgs(trailingOnly = TRUE)
  if (full) {
    simulate_zeal_nil(n_per_family = 17L, n_families = NULL, tag = "full")
  } else {
    log_info("[zeal_nil] SMOKE (5/family, ~2 families/taxon). Pass --generate for the full cohort.")
    simulate_zeal_nil(n_per_family = 5L, n_families = 10L, tag = "smoke")
  }
}
