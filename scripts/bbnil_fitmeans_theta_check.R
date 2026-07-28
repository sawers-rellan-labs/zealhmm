#!/usr/bin/env Rscript
# Does bbnil's EM-fit emission land at the compressed geometry {ε, d/2, d}?
# ---------------------------------------------------------------------------
# snp50k_genotype_identifiability.qmd §2.4/§3 predicts the count-emission state means are NOT
# {0, ½, 1} but {ε, d/2, d} with the transmitted donor dose d = f2 + ½f1 ≈ 0.28–0.35 per taxon,
# so a fixed c(err,0.5,1-err) emission is mis-specified and fit_means is the correct fix (not err/
# conc, which cannot move HET off 0.5). This diagnostic closes the loop empirically: it runs the
# ACTUAL EM (nilHMM:::.em_fit_means, the exact code call_ancestry uses) on the real skim counts,
# pulls the fitted theta = (REF, HET, HOM), and checks it against the panel d per taxon.
#
# It also reports the marginal log-likelihood BEFORE (fixed theta0 = c(err,0.5,1-err)) vs AFTER the
# fit, via an exact 3-state forward pass -- the EM-likelihood trace, for inspection.
#
# No sweep, no figure: a table. rrate is pinned to the map value; conc = 20; design = BC2S3.
#   Rscript scripts/bbnil_fitmeans_theta_check.R

suppressMessages({
  library(devtools)
  library(data.table)
  library(here)
})
devtools::load_all(path.expand("~/repos/nilhmm"), quiet = TRUE)
source(here::here("scripts/logging.R"))

OUT <- here::here("results/sim/zeal_nil")
R_MOLB <- 1.67e-3 # map-derived rrate (same as the molb leg of the calibration figures)
ERR <- 0.01
CONC <- 20

# Panel d per taxon (identifiability §2.3) and the read-based d (§2.5), for reference.
PANEL <- data.table(
  taxon = c("Zv", "Zx", "Zl", "Zd", "Zn", "Zh"),
  panel_d = c(0.282, 0.330, 0.333, 0.346, 0.313, 0.37),
  read_HET = c(0.157, 0.187, 0.164, 0.198, NA, 0.339), # §2.5 HET alt-frac (~d/2)
  read_HOM = c(0.342, 0.492, 0.363, 0.461, NA, 0.894) # §2.5 HOM alt-frac (~d)
)

# ---- exact 3-state (geometric, n_sub = 1) forward log-likelihood ----------------------------
logsumexp <- function(v) {
  m <- max(v)
  if (!is.finite(m)) {
    return(-Inf)
  }
  m + log(sum(exp(v - m)))
}
chr_loglik <- function(em, log_start, log_trans) { # em: T x 3 emission loglik
  Tn <- nrow(em)
  a <- log_start + em[1, ]
  for (t in 2:Tn) {
    a <- em[t, ] + vapply(1:3, function(s) logsumexp(a + log_trans[, s]), numeric(1))
  }
  logsumexp(a)
}
sample_loglik <- function(obs_list, theta, log_start, log_trans) {
  sum(vapply(obs_list, function(o) {
    em <- nilHMM:::count_emission_loglik_cpp(as.integer(o$n), as.integer(o$a), theta, CONC)
    chr_loglik(em, log_start, log_trans)
  }, numeric(1)))
}

# ---- faithfully reconstruct the engine's fit path (mirrors call_states one_sample) ----------
spec <- caller_spec("bbnil", rrate = R_MOLB, err = ERR, conc = CONC, fit_means = TRUE)
priors <- nilHMM:::.state_freqs("BC2S3", NULL, NULL)
td <- nilHMM:::.duration_transition(spec$duration, priors) # n_sub = 1 (geometric)
theta0 <- nilHMM:::.emission_theta(spec$emission) # c(err, 0.5, 1 - err)
stopifnot(td$n_sub == 1L)
log_info(
  "[theta] theta0 (fixed, mis-specified) = {%.3f, %.3f, %.3f}; rrate=%.2e conc=%g design=BC2S3",
  theta0[1], theta0[2], theta0[3], R_MOLB, CONC
)

# ---- real skim samples + taxon from the calibration pairing ---------------------------------
pair <- fread(here::here("data/zeal/correspondence/calibration_pairing.csv"))[in_calibration == TRUE]
pair[, taxon := sub("\\..*", "", true_pedigree)]

fit_one <- function(s) {
  cf <- fread(here::here("data/zeal/skim/counts", paste0(s, ".tsv")),
    header = FALSE,
    col.names = c("contig", "pos", "rb", "rc", "ab", "ac")
  )
  d <- data.table(
    chr = as.integer(sub("chr", "", cf$contig)), pos = cf$pos,
    n = as.integer(cf$rc + cf$ac), a = as.integer(cf$ac)
  )
  setorder(d, chr, pos)
  obs_list <- lapply(split(d, d$chr), function(dc) list(chr = dc$chr[1], pos = dc$pos, n = dc$n, a = dc$a))
  theta_fit <- nilHMM:::.em_fit_means(obs_list, spec$emission, td, theta0)
  ll0 <- sample_loglik(obs_list, theta0, td$log_start, td$log_trans)
  ll1 <- sample_loglik(obs_list, theta_fit, td$log_start, td$log_trans)
  data.table(
    sample = s,
    theta_REF = theta_fit[1], theta_HET = theta_fit[2], theta_HOM = theta_fit[3],
    ll_fixed = ll0, ll_fit = ll1, ll_gain = ll1 - ll0
  )
}

res <- rbindlist(lapply(pair$test_sample, function(s) {
  r <- fit_one(s)
  log_info(
    "[theta] %-14s fit {%.3f, %.3f, %.3f}  LL %.0f -> %.0f (+%.0f)",
    s, r$theta_REF, r$theta_HET, r$theta_HOM, r$ll_fixed, r$ll_fit, r$ll_gain
  )
  r
}))
res <- merge(res, pair[, .(sample = test_sample, taxon)], by = "sample")

# ---- per-taxon summary vs the panel prediction ----------------------------------------------
by_tax <- res[, .(
  n = .N,
  fit_REF = median(theta_REF), fit_HET = median(theta_HET), fit_HOM = median(theta_HOM),
  ll_gain = median(ll_gain)
), by = taxon]
by_tax <- merge(by_tax, PANEL, by = "taxon", all.x = TRUE)
by_tax[, `:=`(exp_HET = panel_d / 2, exp_HOM = panel_d)]
setcolorder(by_tax, c(
  "taxon", "n", "fit_REF", "fit_HET", "exp_HET", "read_HET",
  "fit_HOM", "exp_HOM", "read_HOM", "panel_d", "ll_gain"
))

cat("\n=== per-sample fitted theta (real skim, bbnil EM, rrate=map, conc=20) ===\n")
print(res[order(taxon, sample)], digits = 3)
cat("\n=== per-taxon fitted theta vs panel/read-based d (identifiability §2.3/§2.5) ===\n")
cat("fit_* = median EM-fit theta; exp_* = panel prediction {d/2, d}; read_* = §2.5 read-based\n\n")
print(by_tax[order(taxon)], digits = 3)

fwrite(res, file.path(OUT, "bbnil_fitmeans_theta_check_persample.csv"))
fwrite(by_tax, file.path(OUT, "bbnil_fitmeans_theta_check_bytaxon.csv"))
log_info("[theta] wrote per-sample + per-taxon CSVs to %s", OUT)
