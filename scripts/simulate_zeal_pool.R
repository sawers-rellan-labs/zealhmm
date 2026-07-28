# ZEAL count-from-parents BC2S3 6-sib-pool simulator (Phase-1 calibration input).
#
# Fidelity upgrade of R/simulate.R's inbred-donor model for the ZEAL population.
# The recurrent parent is B73 (inbred, REF); the donor is a TEOSINTE reference
# individual drawn (with replacement) as a TAXON-MATCHED stand-in for the real
# donor accession (we lack the accession's own genotypes -- see
# snp50k_genotype_identifiability.qmd s2.6). Teosinte donors are ~20% heterozygous
# per site, so the inbred `observed = ancestry_dose * founder_allele` shortcut is
# not exact at het sites.
#
# HET INJECTION -- TIER A (per-site Bernoulli resolution; decided 2026-07-22):
#   All donor ancestry in a line traces to ONE founding donor->F1 gamete, so at a
#   donor-het site that gamete carried ALT or REF (decided once per founding
#   lineage) and every donor-ancestry haplotype in the line carries that same
#   allele. We therefore resolve each het site once per LINE by a Bernoulli(0.5)
#   draw: ~half the donor-het sites become informative (behave like a hom-ALT
#   donor), ~half go silent (non-informative). Averaged over het sites this gives
#   E[obs] = 0.5 * ancestry_dose -- the "F1 looks like a BC1 at het sites" halving.
#   Tier A reproduces this MARGINAL but not the donor-side het-block LD; the LD
#   upgrade (Tier B) needs phased donor haplotypes (the schnable2023 lineage is
#   unphased top to bottom, so Tier B would run Beagle on the SNP50K panel first).
#
# 6-PLANT POOL (emission detail, not a design change): each line = 6 BC2S3 sibs
# sharing the BC2 backbone (one simcross realization, 6 selfing leaves), ancestry
# dosages averaged to a pooled dose, THEN skim counts drawn. Truth = the ancestry
# mosaic (rounded pooled dose); non-informativeness is observation masking, not truth.
#
# STORAGE: the interval master (per-sib ancestry segments, bgzip+tabix) + the
# families/founder-draw table + the RNG seed are the grid-independent, lossless
# record; the raster (.rds of pooled dose + REF/ALT counts) is a regenerable cache.
#
#   Rscript scripts/simulate_zeal_pool.R            # SMOKE (n_per_family=5, few families)
#   Rscript scripts/simulate_zeal_pool.R --generate # full cohort (82 families)

suppressMessages({
  library(data.table)
  library(here)
  library(nilHMM) # call_gt, breeding_prior, parse_design, call_ancestry
})
source(here::here("R/simulate.R")) # .hap_allele_at, .truth_segments (+ others)
source(here::here("scripts/logging.R")) # log_info

# --- BC2S3 6-sib pool pedigree (6 selfing leaves off a shared BC backbone) ----
.bc2s3_pool_pedigree <- function(n_bc = 2L, n_self = 3L, n_sib = 6L) {
  id <- c(1L, 2L)
  mom <- c(0L, 0L)
  dad <- c(0L, 0L)
  gen <- c(0L, 0L)
  add <- function(mm, dd, g) {
    i <- length(id) + 1L
    id[i] <<- i
    mom[i] <<- mm
    dad[i] <<- dd
    gen[i] <<- g
    i
  }
  cur <- add(1L, 2L, 1L) # F1
  for (k in seq_len(n_bc)) cur <- add(1L, cur, 1L + k) # BCk = recurrent x prev
  backbone <- cur
  g0 <- 1L + n_bc
  nil_ids <- integer(n_sib)
  for (s in seq_len(n_sib)) { # 6 independent selfing chains off the shared backbone
    c2 <- backbone
    g <- g0
    for (j in seq_len(n_self)) {
      g <- g + 1L
      c2 <- add(c2, c2, g)
    }
    nil_ids[s] <- c2
  }
  sex <- ifelse(id %in% mom[mom > 0L], 0L, 1L)
  list(
    ped = data.frame(id = id, mom = mom, dad = dad, sex = sex, gen = gen),
    nil_id = as.character(nil_ids)
  )
}

# one simcross realization -> M x n_sib donor-ancestry dosage (shared backbone)
.simulate_pool_dosage <- function(ped, cmlen, markers, m, p, nil_ids, donor_allele = 2L) {
  sim <- simcross::sim_from_pedigree(ped, L = cmlen$L, m = m, p = p)
  M <- nrow(markers)
  out <- matrix(0L, M, length(nil_ids))
  for (ci in seq_len(nrow(cmlen))) {
    ch <- cmlen$chr[ci]
    rows <- which(markers$chr == ch)
    if (!length(rows)) next
    cm <- markers$cm[rows]
    for (s in seq_along(nil_ids)) {
      nil <- sim[[ci]][[nil_ids[s]]]
      out[rows, s] <- (.hap_allele_at(nil$mat, cm) == donor_allele) +
        (.hap_allele_at(nil$pat, cm) == donor_allele)
    }
  }
  out
}

# skim counts from a CONTINUOUS pooled ALT fraction (generalizes R/simulate.R's
# .draw_counts, which maps integer dosage -> {0,.5,1}).
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

# --- SNP50K teosinte reference-panel genotypes -> {0,1,2,3} matrix (cached) ----
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
  st <- system2("bcftools", c("query", "-f", "'%CHROM\\t%POS[\\t%GT]\\n'", vcf), stdout = tmp)
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
  key_v <- paste(raw$chr, raw$pos, sep = "_")
  key_m <- paste(markers$chr, markers$bp, sep = "_")
  idx <- match(key_m, key_v)
  out_mat <- matrix(3L, nrow(markers), length(samples), dimnames = list(NULL, samples))
  ok <- !is.na(idx)
  out_mat[ok, ] <- gmat[idx[ok], ]
  g <- list(mat = out_mat, samples = samples, n_absent = sum(!ok))
  saveRDS(g, cache, compress = "gzip")
  g
}

# fit the SNP50K skim coverage regime (per-sample lambda already tabulated) ------
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

#' Simulate the ZEAL 6-sib-pool cohort (BC2S3, Tier A het injection)
#' @return (invisibly) list(rds, intervals, families, grid).
simulate_zeal_pool <- function(design = "BC2S3",
                               n_per_family = 5L,
                               n_families = NULL, # NULL = all cohort accessions
                               n_sib = 6L,
                               taxa = c("Zx", "Zv", "Zd", "Zl", "Zh"),
                               m = 10L, p = 0, seed = 1L,
                               markers_path = here::here("data/zeal/markers_snp50k_cm.tsv"),
                               cohort_meta = here::here("data/zeal/reference_panel/bzea_50K_cohort_ref_metadata.csv"),
                               ref_ids = here::here("data/zeal/reference_panel/reference_ids_by_taxon.csv"),
                               outdir = here::here("results/sim/zeal_pool"),
                               tag = "smoke") {
  if (!requireNamespace("simcross", quietly = TRUE)) {
    stop("simulate_zeal_pool() needs the 'simcross' package (kbroman/simcross).")
  }
  set.seed(seed)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  t0 <- Sys.time()

  # ---- SNP50K marker grid + cM (the calibration grid; full 50K, NO thinning) --
  mk <- fread(markers_path)
  markers <- data.table(chr = as.integer(mk$chr), bp = as.integer(mk$pos), cm = as.numeric(mk$cm))[order(chr, bp)]
  cmlen <- markers[, .(L = max(cm)), by = chr][order(chr)]
  M <- nrow(markers)
  map_r <- 2 * markers[, sum(tapply(cm, chr, function(x) max(x) - min(x)))] / (100 * M)
  log_info("grid: %d SNP50K markers on %d chr; map r = %.3e (design %s)", M, nrow(cmlen), map_r, design)

  # ---- families = cohort donor accessions x taxon-matched reference stand-in --
  cm_meta <- fread(cohort_meta)
  fams <- unique(cm_meta[
    is_cohort == TRUE & !is.na(donor_accession) & maizegdb_prefix %in% taxa,
    .(family = donor_accession, taxon = maizegdb_prefix)
  ])
  refs <- fread(ref_ids)[maizegdb_prefix %in% taxa, .(sample, taxon = maizegdb_prefix)]
  if (!is.null(n_families)) { # smoke: up to ~2 families per taxon, covering all taxa
    per_tax <- max(1L, ceiling(n_families / length(taxa)))
    fams <- fams[, head(.SD, per_tax), by = taxon][seq_len(min(.N, n_families))]
  }
  # draw one stand-in ref individual per family (WITH replacement within taxon)
  fams[, stand_in := {
    pool <- refs[taxon == .BY$taxon]$sample
    if (!length(pool)) stop("no reference stand-in for taxon ", .BY$taxon)
    sample(pool, .N, replace = TRUE)
  }, by = taxon]
  log_info("families: %d accessions across %d taxa; %d lines/family x %d sibs", nrow(fams), uniqueN(fams$taxon), n_per_family, n_sib)

  # ---- inputs: reference-panel genotypes + skim regime + design prior ---------
  gp <- load_ref_panel_gt(markers = markers)
  reg <- fit_snp50k_regime()
  design_prior <- breeding_prior(design)
  ped <- .bc2s3_pool_pedigree(parse_design(design)$n_bc, parse_design(design)$n_self, n_sib)
  log_info(
    "ref panel: %d markers x %d samples (%d absent->missing); regime lambda_mean=%.3f shape=%.2f",
    M, length(gp$samples), gp$n_absent, reg$lambda_mean, reg$shape
  )

  # ---- simulate every line: 6-sib pool -> pooled dose -> skim counts ----------
  n_lines <- nrow(fams) * n_per_family
  n_ref <- matrix(0L, M, n_lines)
  n_alt <- matrix(0L, M, n_lines)
  pool_dose <- matrix(0, M, n_lines)
  line_names <- character(n_lines)
  line_family <- character(n_lines)
  truth_l <- vector("list", n_lines) # pooled-dose truth (rounded) per line
  ivl_l <- vector("list", n_lines) # per-sib ancestry segments (interval master)
  col <- 0L
  for (fi in seq_len(nrow(fams))) {
    fa <- fams[fi]
    g_vec <- gp$mat[, fa$stand_in] # stand-in donor genotype {0,1,2,3}
    for (k in seq_len(n_per_family)) {
      col <- col + 1L
      nm <- sprintf("%s_%03d", fa$family, k)
      dmat <- .simulate_pool_dosage(ped$ped, cmlen, markers, m, p, ped$nil_id) # M x n_sib
      # per-LINE het resolution (Tier A): het site informative w.p. 0.5, fixed for the line
      founder_allele <- integer(M)
      founder_allele[g_vec == 2L] <- 1L
      het <- which(g_vec == 1L)
      if (length(het)) founder_allele[het] <- stats::rbinom(length(het), 1L, 0.5)
      pd <- rowMeans(dmat) # pooled ancestry dose in [0,2]
      p_alt <- (pd / 2) * founder_allele # observed pooled ALT fraction
      lambda <- max(0.01, stats::rgamma(1, shape = reg$shape, scale = reg$lambda_mean / reg$shape))
      ac <- .draw_counts_pfrac(p_alt, lambda, reg$pi_floor, reg$k_decay, reg$error)
      n_ref[, col] <- ac$ref
      n_alt[, col] <- ac$alt
      pool_dose[, col] <- pd
      line_names[col] <- nm
      line_family[col] <- fa$family
      truth_l[[col]] <- .truth_segments(markers, as.integer(round(pd)), nm)
      # interval master: per-sib ancestry dosage segments (grid = SNP50K positions)
      ivl_l[[col]] <- rbindlist(lapply(seq_len(n_sib), function(s) {
        seg <- .truth_segments(markers, dmat[, s], nm)
        seg[, .(chrom = chr, start = start_bp, end = end_bp, line = nm, sib = s, state, donor = fa$taxon)]
      }))
    }
  }
  colnames(n_ref) <- colnames(n_alt) <- colnames(pool_dose) <- line_names

  # ---- write raster bundle (.rds, regenerable cache) --------------------------
  sim <- list(
    grid = data.frame(chr = markers$chr, pos = markers$bp), n_ref = n_ref, n_alt = n_alt,
    pool_dose = pool_dose, names = line_names, family = line_family,
    truth = rbindlist(truth_l), design = design, n_markers = M, n_sib = n_sib,
    seed = seed, regime = reg, map_r = map_r, het_injection = "tierA_bernoulli"
  )
  rds <- file.path(outdir, sprintf("zeal_pool_%s_%s.rds", tolower(design), tag))
  saveRDS(sim, rds, compress = "gzip")

  # ---- write interval master (bgzip + tabix BED-like TSV) + families table ----
  ivl <- rbindlist(ivl_l)[order(chrom, start)]
  ivl_tsv <- file.path(outdir, sprintf("zeal_pool_%s_%s_intervals.tsv", tolower(design), tag))
  fwrite(ivl, ivl_tsv, sep = "\t")
  bg <- suppressWarnings(system2("bgzip", c("-f", ivl_tsv)))
  if (bg == 0L) system2("tabix", c("-f", "-s", "1", "-b", "2", "-e", "3", "-S", "1", paste0(ivl_tsv, ".gz")))
  fwrite(fams, file.path(outdir, sprintf("zeal_pool_%s_%s_families.csv", tolower(design), tag)))

  sim_mr <- mean((n_ref + n_alt) == 0L)
  log_info(
    "wrote %s : %d markers x %d lines (%d families); missing=%.3f (target ~%.3f), pooled-donor sites=%.3f",
    basename(rds), M, n_lines, nrow(fams), sim_mr, reg$missing_mean, mean(pool_dose > 0)
  )
  log_info(
    "interval master: %d segments (bgzip=%s); total %.1fs", nrow(ivl),
    if (bg == 0L) "yes" else "no (bgzip missing)", as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )

  invisible(list(rds = rds, sim = sim, families = fams, markers = markers, design_prior = design_prior, map_r = map_r))
}

# --- CLI -------------------------------------------------------------------
if (sys.nframe() == 0L) {
  full <- "--generate" %in% commandArgs(trailingOnly = TRUE)
  if (full) {
    res <- simulate_zeal_pool(n_per_family = 17L, n_families = NULL, tag = "full")
  } else {
    log_info("[zeal_pool] SMOKE (5/family, ~2 families/taxon). Pass --generate for the full cohort.")
    res <- simulate_zeal_pool(n_per_family = 5L, n_families = 10L, tag = "smoke")

    # time ONE nnil decode (hard-call via BC2S3 design prior, then nnil) to size the sweep
    tryCatch(
      {
        sim <- res$sim
        one <- sim$names[1]
        idx <- match(one, sim$names)
        dat <- data.table(
          name = one, chr = sim$grid$chr, pos = sim$grid$pos,
          n_ref = sim$n_ref[, idx], n_alt = sim$n_alt[, idx]
        )
        dat[, g := {
          gg <- call_gt(n_ref, n_alt, prior = res$design_prior, error = 0.01)
          gg[is.na(gg)] <- 3L
          as.integer(gg)
        }]
        log_info("[sanity] mean pooled donor dose = %.3f (BC2S3 expectation ~0.25 = 2 x 12.5%%)", mean(sim$pool_dose))
        td <- Sys.time()
        called <- call_ancestry(
          data = dat[, .(name, chr, pos, g)], caller = "nnil", rrate = sim$map_r,
          design = sim$design, germ = 0.01, gert = 0.5, p = 0, nir = 0.3, mr = 0.1
        )
        dt <- as.numeric(difftime(Sys.time(), td, units = "secs"))
        log_info(
          "[decode] one nnil decode = %.2fs -> full nir sweep (11 pts x %d lines) ~ %.1f min on 10 cores",
          dt, length(sim$names), dt * 11 * length(sim$names) / 10 / 60
        )
      },
      error = function(e) log_info("[decode] timing skipped: %s", conditionMessage(e))
    )
  }
}
