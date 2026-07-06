#!/usr/bin/env Rscript
# =============================================================================
# STAM GWAS degradation sweep: GL+HWE interpolation CONTROL vs coverage (TeoNAM)
# Plan: agent/teonam-control-sweep-plan.md  (sibling of the RTIGER sweep,
#       scripts/teonam_rtiger_sweep.R, which this mirrors)
# -----------------------------------------------------------------------------
# The no-HMM control for the coverage-degradation comparison. The shared failure
# mode across ancestry callers is HET EXCESS, so the control is the naive caller
# that is MOST het-prone: a per-site genotype-likelihood call with an HWE posterior
# (a single ALT read -> HET via the 2p(1-p) prior term). Closed-form on the
# simulated (n_ref, n_alt); NO external caller (GATK/bcftools) is run.
#
# For each TeoNAM family x coverage lambda:
#   1. simulate low-coverage reads from the real (truth) 0/1/2 genotypes
#      [R/simulate.R::.draw_counts(pi_floor=0, k_decay=1, error=0.01)]
#   2. GL+HWE-call genotypes [nilHMM::call_gl(prior="hwe", af=<per-marker truth
#      teosinte AF>, error=0.01)] -> 0/1/2, NA at zero depth (het-excess control)
#   3. PER-RIL step-interpolate each RIL's called (non-NA) markers onto the union
#      cM grid [nilHMM::interpolate_genotype(mode="step")]. Unlike RTIGER (which
#      decodes EVERY marker -> a complete rectangular block), call_gl leaves NA at
#      uncovered markers, so covered markers vary per RIL and interpolation is
#      per-RIL. Cheap: no HMM.
# Then assemble the union matrix (~47,750 markers x 1,237 lines) at each lambda and
# run the same GWAS scan as the baseline (STAM ~ Family + marker, 1-df F).
#
# Writes results/sim/teonam/stam_gwas_control_lambda<L>.csv (per lambda), the
# combined long table stam_gwas_control_sweep.csv (adds `coverage`; includes the
# lambda=Inf baseline rows copied from data/teonam/stam_gwas_scan_interpolated.csv), and
# stam_control_het_fraction.csv (per-lambda HET-call fraction + call rate).
#
# KEY DECISIONS (recorded per plan):
#  - Caller: GL + HWE posterior (het-excess control). error = 0.01. Per-marker af =
#    teosinte (ALT) allele frequency from the TRUTH genotypes = rowMeans(D)/2 over
#    that family's RILs (population parameter, per family per marker).
#  - Coverage grid EXACTLY {0.1, 0.2, 0.5, 1, 5, 10, 20} (+ lambda=Inf baseline).
#  - Read model: pi_floor=0, k_decay=1, error=0.01; 1 replicate per (family,lambda),
#    RNG seed = 1000 + 100*family_index + lambda_index (same scheme as the RTIGER
#    sweep, so both sweeps degrade the identical truth mosaics).
#  - Zero-depth markers -> NA (call_gl); "covered" = non-NA. No min_cov floor:
#    interpolation fills uncovered union markers from each RIL's flanking calls.
#
# Run:  Rscript scripts/teonam_control_sweep.R --generate
# =============================================================================

suppressMessages({
  library(data.table)
  library(parallel)
  library(readxl)
  devtools::load_all("/Users/fvrodriguez/repos/nilhmm", quiet = TRUE)
})

ROOT <- "/Users/fvrodriguez/repos/zealhmm"
setwd(ROOT)
source(file.path(ROOT, "R/simulate.R")) # .draw_counts()
OUTDIR <- file.path(ROOT, "results/sim/teonam")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

if (!("--generate" %in% commandArgs(TRUE))) {
  message("teonam_control_sweep.R: pass --generate to (re)compute the sweep CSVs.")
  quit(save = "no", status = 0)
}

LAMBDAS <- c(0.1, 0.2, 0.5, 1, 5, 10, 20)
ERROR <- 0.01 # GL per-read error (matches read model)
THREADS <- max(1L, detectCores() - 2L)
READ_PARS <- list(pi_floor = 0, k_decay = 1, error = 0.01)

# --- marker map + union target grid (identical to the RTIGER sweep) ----------
mc <- fread(file.path(ROOT, "data/teonam/marker_info_v5_cm.tsv"))
setnames(mc, "chr_v5", "chr")
cm_by <- setNames(mc$cm, mc$marker)
pos_by <- setNames(mc$pos_v5, mc$marker)

setorder(mc, chr, cm)
u <- mc[, .SD[!duplicated(cm)], by = chr] # union grid: cm-dedup per chr
setorder(u, chr, cm)
target_df <- data.frame(chr = as.integer(u$chr), cm = as.numeric(u$cm))
union_markers <- u$marker
union_pos <- as.integer(pos_by[union_markers])
union_chr <- as.integer(u$chr)
message(sprintf("union grid: %d markers x 10 chromosomes", nrow(u)))

# --- per-family: genotypes, robust keys, cm-dedup marker table, truth dosage --
fams <- c(
  TIL01 = "W22TIL01_genotype.csv", TIL03 = "W22TIL03_genotype.csv",
  TIL11 = "W22TIL11_genotype.csv", TIL14 = "W22TIL14_genotype.csv",
  TIL25 = "W22TIL25_genotype.csv"
)

load_family <- function(fam) {
  g <- fread(file.path(ROOT, "data/teonam", fams[fam]))
  g <- g[!duplicated(g[[1]])]
  keys <- paste0(fam, sub("^.*Line_", "", g[[1]])) # e.g. TIL01A001
  mk <- intersect(names(g)[-(1:3)], mc$marker)
  mt <- data.table(
    marker = mk,
    chr = as.integer(mc$chr[match(mk, mc$marker)]),
    pos = as.integer(pos_by[mk]),
    cm = as.numeric(cm_by[mk])
  )
  setorder(mt, chr, cm)
  mt <- mt[, .SD[!duplicated(cm)], by = chr] # strictly increasing cm/chr
  setorder(mt, chr, cm)
  D <- t(as.matrix(g[, mt$marker, with = FALSE])) # markers(mt order) x RILs
  storage.mode(D) <- "double"
  colnames(D) <- keys
  af <- rowMeans(D) / 2 # per-marker truth teosinte AF
  list(mt = mt, D = D, keys = keys, af = af)
}
message("loading families ...")
fam_data <- lapply(names(fams), load_family)
names(fam_data) <- names(fams)
for (f in names(fams)) {
  message(sprintf(
    "  %s: %d markers x %d RILs (mean truth teosinte AF = %.3f)", f,
    nrow(fam_data[[f]]$mt), length(fam_data[[f]]$keys),
    mean(fam_data[[f]]$af)
  ))
}

# --- one (family,lambda) cell: simulate -> GL+HWE -> per-RIL interpolate ------
recover_block <- function(fam, li) {
  lambda <- LAMBDAS[li]
  fi <- match(fam, names(fams))
  fd <- fam_data[[fam]]
  mt <- fd$mt
  D <- fd$D
  keys <- fd$keys
  af <- fd$af
  M <- nrow(mt)
  N <- length(keys)

  set.seed(1000L + 100L * fi + li) # per-cell seed (as RTIGER sweep)
  ac <- .draw_counts(as.vector(D),
    lambda = lambda,
    pi_floor = READ_PARS$pi_floor, k_decay = READ_PARS$k_decay,
    error = READ_PARS$error
  )
  n_ref <- matrix(as.integer(ac$ref), M, N)
  n_alt <- matrix(as.integer(ac$alt), M, N)

  # GL+HWE call: 0/1/2, NA at zero depth. af recycled column-wise over M x N.
  calls <- call_gl(n_ref, n_alt, prior = "hwe", af = af, error = ERROR) # M x N int, NA=uncovered

  covered <- !is.na(calls)
  n_called <- sum(covered)
  n_het <- sum(calls == 1L, na.rm = TRUE)

  # per-RIL step-interpolation of the RIL's covered calls onto the union grid.
  block <- matrix(NA_real_, nrow(u), N, dimnames = list(union_markers, keys))
  for (j in seq_len(N)) {
    obs <- which(covered[, j]) # this RIL's called markers (mt order)
    if (!length(obs)) next # (never at these lambda; guarded)
    geno <- matrix(as.double(calls[obs, j]), ncol = 1L)
    obs_df <- data.frame(chr = mt$chr[obs], cm = mt$cm[obs])
    # a target chr with no observed markers for this RIL would error; restrict the
    # target to chromosomes this RIL covers, leave the rest NA (essentially never).
    ok_chr <- unique(obs_df$chr)
    tsel <- target_df$chr %in% ok_chr
    block[tsel, j] <- interpolate_genotype(geno, obs_df, target_df[tsel, , drop = FALSE],
      mode = "step"
    )[, 1L]
  }
  list(
    lambda = lambda, fam = fam, block = block,
    n_called = n_called, n_het = n_het, n_cells = M * N
  )
}

grid <- expand.grid(
  fam = names(fams), li = seq_along(LAMBDAS),
  stringsAsFactors = FALSE
)
message(sprintf(
  "GL+HWE control sweep: %d (family,lambda) cells, %d threads ...",
  nrow(grid), THREADS
))
t0 <- Sys.time()
cells <- mclapply(seq_len(nrow(grid)), function(i) {
  recover_block(grid$fam[i], grid$li[i])
}, mc.cores = THREADS)
bad <- vapply(cells, function(x) inherits(x, "try-error") || is.null(x), logical(1))
if (any(bad)) {
  stop(
    "call_gl cell(s) failed: ",
    paste(which(bad), collapse = ", "), " -> ", cells[[which(bad)[1]]]
  )
}
message(sprintf(
  "  control sweep done in %.1f min",
  as.numeric(Sys.time() - t0, units = "mins")
))

# --- phenotype (STAM) --------------------------------------------------------
ph <- as.data.frame(read_excel(
  file.path(ROOT, "data/teonam/9250682/TeoNAM_1257RILs_22traits_phenotype_data.xlsx")
))
names(ph)[1] <- "line"
stam_by <- setNames(ph$STAM, ph$line)

# --- GWAS scan: STAM ~ Family + marker (additive, 1 df) ----------------------
gwas_scan <- function(G) {
  y <- suppressWarnings(as.numeric(stam_by[colnames(G)]))
  fam <- factor(substr(colnames(G), 1, 5))
  ok <- !is.na(y)
  y <- y[ok]
  fam <- fam[ok]
  Gm <- G[, ok, drop = FALSE]
  Xr <- model.matrix(~fam)
  n <- length(y)
  RSS0 <- sum(lm.fit(Xr, y)$residuals^2)
  scan1 <- function(i) {
    g <- Gm[i, ]
    if (anyNA(g) || sd(g) == 0) {
      return(NA_real_)
    }
    fit <- lm.fit(cbind(Xr, g), y)
    RSS1 <- sum(fit$residuals^2)
    df2 <- n - fit$rank
    pf(((RSS0 - RSS1) / 1) / (RSS1 / df2), 1, df2, lower.tail = FALSE)
  }
  P <- unlist(mclapply(seq_len(nrow(Gm)), scan1, mc.cores = THREADS))
  data.table(SNP = rownames(G), CHR = union_chr, BP = union_pos, P = P)
}

# --- tb1 peak helper (canonical STAM QTL) ------------------------------------
tb1 <- data.table(chr = 1L, start = 272330564L) # tb1 (Zm00001eb054440)
tb1_peak <- function(scan) {
  w <- scan[CHR == tb1$chr & abs(BP - tb1$start) <= 5e5 & is.finite(P) & P > 0]
  if (!nrow(w)) {
    return(NA_real_)
  }
  round(max(-log10(w$P)), 2)
}

# --- per-lambda: cbind family blocks, scan, write; accumulate het stats -------
sweep_list <- vector("list", length(LAMBDAS))
het_list <- vector("list", length(LAMBDAS))
for (li in seq_along(LAMBDAS)) {
  lambda <- LAMBDAS[li]
  idx <- which(grid$li == li)
  blocks <- lapply(idx, function(i) cells[[i]]$block) # one per family
  G <- do.call(cbind, blocks)
  rownames(G) <- union_markers
  n_na <- sum(is.na(G))
  if (n_na) {
    message(sprintf(
      "  lambda=%-4g : %d NA union cells (empty-chr RILs) imputed by per-marker mean",
      lambda, n_na
    ))
  }
  if (n_na) { # defensive: essentially never
    rm <- rowMeans(G, na.rm = TRUE)
    na_ij <- which(is.na(G), arr.ind = TRUE)
    G[na_ij] <- rm[na_ij[, 1]]
  }
  scan <- gwas_scan(G)[order(CHR, BP)]
  fwrite(scan, file.path(OUTDIR, sprintf("stam_gwas_control_lambda%s.csv", lambda)))
  scan[, coverage := lambda]
  sweep_list[[li]] <- scan

  n_called <- sum(vapply(idx, function(i) cells[[i]]$n_called, numeric(1)))
  n_het <- sum(vapply(idx, function(i) cells[[i]]$n_het, numeric(1)))
  n_cells <- sum(vapply(idx, function(i) cells[[i]]$n_cells, numeric(1)))
  het_list[[li]] <- data.table(
    coverage = lambda,
    het_frac = n_het / n_called, # HET / covered calls
    call_rate = n_called / n_cells, # covered / all cells
    n_het = n_het, n_called = n_called
  )
  message(sprintf(
    "  lambda=%-4g : %d markers, het-frac = %.3f, call-rate = %.3f, tb1 peak -log10P = %s, global max = %.1f",
    lambda, nrow(scan), n_het / n_called, n_called / n_cells,
    tb1_peak(scan), max(-log10(scan[is.finite(P) & P > 0, P]))
  ))
}

# --- combined sweep + lambda=Inf baseline ------------------------------------
baseline <- fread(file.path(ROOT, "data/teonam/stam_gwas_scan_interpolated.csv"))
baseline[, coverage := Inf]
message(sprintf(
  "  lambda=Inf  : %d markers, tb1 peak -log10P = %s (baseline, panel C)",
  nrow(baseline), tb1_peak(baseline)
))

sweep <- rbindlist(c(sweep_list, list(baseline)), use.names = TRUE)
fwrite(sweep, file.path(OUTDIR, "stam_gwas_control_sweep.csv"))
message(sprintf(
  "wrote %s (%d rows, %d coverage levels)",
  file.path(OUTDIR, "stam_gwas_control_sweep.csv"),
  nrow(sweep), uniqueN(sweep$coverage)
))

het <- rbindlist(het_list)
fwrite(het, file.path(OUTDIR, "stam_control_het_fraction.csv"))
message("wrote stam_control_het_fraction.csv:")
print(het)
