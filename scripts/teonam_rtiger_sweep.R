#!/usr/bin/env Rscript
# =============================================================================
# STAM GWAS degradation sweep: RTIGER ancestry vs sequencing coverage (TeoNAM)
# Plan: agent/teonam-rtiger-degradation-plan.md
# -----------------------------------------------------------------------------
# For each TeoNAM family x coverage lambda:
#   1. simulate low-coverage reads from the real (truth) 0/1/2 genotypes
#      [R/simulate.R::.draw_counts(pi_floor=0, k_decay=1, error=0.01)]
#   2. RTIGER-call ancestry [call_ancestry(caller="rtiger", rigidity=2)]
#   3. step-interpolate the recovered per-family block onto the union cM grid
#      [nilHMM::interpolate_genotype(mode="step")]
# Then assemble the union matrix (~47,750 markers x 1,237 lines) at each lambda
# and run the same GWAS scan as the baseline (STAM ~ Family + marker, 1-df F).
#
# Writes results/sim/teonam/stam_gwas_rtiger_lambda<L>.csv (per lambda) and the
# combined long table stam_gwas_rtiger_sweep.csv (adds `coverage`; includes the
# lambda=Inf baseline rows copied from data/teonam/stam_gwas_scan_interpolated.csv).
#
# KEY DECISIONS (recorded per plan S7):
#  - Coverage grid: {0.1, 0.2, 0.5, 1, 5, 10, 20} (+ lambda=Inf baseline = panel C).
#  - RTIGER rigidity = 2  (rigidity_star from results/sim/calib_params.csv).
#  - Design priors: RTIGER fits its OWN start frequencies -- call_states() returns
#    from the caller=="rtiger" branch BEFORE any design/f_1/f_2 is resolved, so
#    the design priors are IGNORED by this caller. (design_priors() has no "BC1S4"
#    entry anyway; TeoNAM is BC1S4 ~15% teo-hom / ~8% het, which would map to
#    f_2~=0.15, f_1~=0.08 had a prior-consuming caller been used.)
#  - Read model: pi_floor=0, k_decay=1, error=0.01; 1 replicate per (family,lambda),
#    RNG seed = 1000 + 100*family_index + lambda_index (deterministic per cell).
#  - min_cov = 0L for the RTIGER call: decode EVERY family marker (uncovered
#    markers carry a flat emission and are filled by the rigidity HMM from
#    neighbours). This yields a COMPLETE rectangular per-family block for step-
#    interpolation and avoids the 2*rigidity per-(sample,chr) coverage floor
#    aborting the batch at low lambda.
#
# Batch RTIGER: ONE call_ancestry per (family,lambda) (joint EM across all family
# RILs), not one fit per RIL. The 5 families x 7 lambda = 35 batch calls are
# fanned out with parallel::mclapply.
#
# Run:  Rscript scripts/teonam_rtiger_sweep.R --generate
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
  message("teonam_rtiger_sweep.R: pass --generate to (re)compute the sweep CSVs.")
  quit(save = "no", status = 0)
}

LAMBDAS <- c(0.1, 0.2, 0.5, 1, 5, 10, 20)
RIGIDITY <- 2L # rigidity_star (calib_params.csv)
THREADS <- max(1L, detectCores() - 2L)
READ_PARS <- list(pi_floor = 0, k_decay = 1, error = 0.01)

# --- marker map + union target grid -----------------------------------------
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
  list(mt = mt, D = D, keys = keys)
}
message("loading families ...")
fam_data <- lapply(names(fams), load_family)
names(fam_data) <- names(fams)
for (f in names(fams)) {
  message(sprintf(
    "  %s: %d markers x %d RILs", f,
    nrow(fam_data[[f]]$mt), length(fam_data[[f]]$keys)
  ))
}

# --- one (family,lambda) cell: simulate -> RTIGER -> interpolate to union -----
recover_block <- function(fam, li) {
  lambda <- LAMBDAS[li]
  fi <- match(fam, names(fams))
  fd <- fam_data[[fam]]
  mt <- fd$mt
  D <- fd$D
  keys <- fd$keys
  M <- nrow(mt)
  N <- length(keys)

  set.seed(1000L + 100L * fi + li) # per-cell seed
  ac <- .draw_counts(as.vector(D),
    lambda = lambda,
    pi_floor = READ_PARS$pi_floor, k_decay = READ_PARS$k_decay,
    error = READ_PARS$error
  )
  n_ref <- as.integer(ac$ref)
  n_alt <- as.integer(ac$alt)

  long <- data.table(
    name = rep(keys, each = M),
    chr = rep(mt$chr, N),
    pos = rep(mt$pos, N),
    n_ref = n_ref, n_alt = n_alt
  )

  st <- call_states(long,
    caller = "rtiger", rigidity = RIGIDITY,
    min_cov = 0L, threads = 1L
  ) # RTIGER ignores priors
  W <- dcast(as.data.table(st), chr + pos ~ name, value.var = "state")
  W <- W[mt[, .(chr, pos)], on = c("chr", "pos")] # align rows to mt (cm order)
  R <- as.matrix(W[, keys, with = FALSE])
  storage.mode(R) <- "double"

  block <- interpolate_genotype(R, data.frame(chr = mt$chr, cm = mt$cm),
    target_df,
    mode = "step"
  )
  colnames(block) <- keys
  list(lambda = lambda, fam = fam, block = block)
}

grid <- expand.grid(
  fam = names(fams), li = seq_along(LAMBDAS),
  stringsAsFactors = FALSE
)
message(sprintf(
  "RTIGER sweep: %d (family,lambda) batch calls, %d threads ...",
  nrow(grid), THREADS
))
t0 <- Sys.time()
cells <- mclapply(seq_len(nrow(grid)), function(i) {
  recover_block(grid$fam[i], grid$li[i])
}, mc.cores = THREADS)
bad <- vapply(cells, function(x) inherits(x, "try-error") || is.null(x), logical(1))
if (any(bad)) {
  stop(
    "RTIGER failed for cell(s): ",
    paste(which(bad), collapse = ", "), " -> ", cells[[which(bad)[1]]]
  )
}
message(sprintf(
  "  RTIGER sweep done in %.1f min",
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
    if (sd(g) == 0) {
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

# --- per-lambda: cbind family blocks, scan, write ----------------------------
tb1 <- data.table(chr = 1L, start = 272330564L) # tb1 (Zm00001eb054440)
tb1_peak <- function(scan) {
  w <- scan[CHR == tb1$chr & abs(BP - tb1$start) <= 5e5 & is.finite(P) & P > 0]
  if (!nrow(w)) {
    return(NA_real_)
  }
  round(max(-log10(w$P)), 2)
}

sweep_list <- vector("list", length(LAMBDAS))
for (li in seq_along(LAMBDAS)) {
  lambda <- LAMBDAS[li]
  idx <- which(grid$li == li)
  blocks <- lapply(idx, function(i) cells[[i]]$block) # one per family
  G <- do.call(cbind, blocks)
  rownames(G) <- union_markers
  scan <- gwas_scan(G)[order(CHR, BP)]
  fwrite(scan, file.path(OUTDIR, sprintf("stam_gwas_rtiger_lambda%s.csv", lambda)))
  scan[, coverage := lambda]
  sweep_list[[li]] <- scan
  message(sprintf(
    "  lambda=%-4g : %d markers, tb1 peak -log10P = %s, global max = %.1f",
    lambda, nrow(scan), tb1_peak(scan),
    max(-log10(scan[is.finite(P) & P > 0, P]))
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
fwrite(sweep, file.path(OUTDIR, "stam_gwas_rtiger_sweep.csv"))
message(sprintf(
  "wrote %s (%d rows, %d coverage levels)",
  file.path(OUTDIR, "stam_gwas_rtiger_sweep.csv"),
  nrow(sweep), uniqueN(sweep$coverage)
))
