#!/usr/bin/env Rscript
# =============================================================================
# STAM GWAS degradation sweep: LB-Impute ancestry vs sequencing coverage (TeoNAM)
# Sibling of scripts/teonam_rtiger_sweep.R — same simulate -> call -> interpolate
# -> GWAS pipeline, swapping the caller to LB-Impute (Fragoso et al. 2014, G3;
# nilHMM native port: coverage-aware emission + distance-dependent transition with
# a double-recombination penalty, full-chromosome Viterbi).
# -----------------------------------------------------------------------------
# For each TeoNAM family x coverage lambda:
#   1. simulate low-coverage reads from the real (truth) 0/1/2 genotypes
#      [R/simulate.R::.draw_counts(pi_floor=0, k_decay=1, error=0.01)]
#   2. LB-Impute-call ancestry  [call_ancestry(caller="lbimpute", unit="cm",
#      recombdist=recombdist_star, genotypeerr=0.05, drp=drp) — recombdist + drp
#      CALIBRATED, read from calib_params.csv]
#   3. step-interpolate the recovered per-family block onto the union grid
#      [nilHMM::interpolate_genotype(mode="step")]
# Then assemble the union matrix at each lambda and run the same GWAS scan as the
# baseline (STAM ~ Family + marker, 1-df F).
#
# GRID: the FULL 51,004-marker GWAS set (teonam_map_v5_gwas) — matches the 51K
# rebuild of the interpolated GWAS/MLM. interpolate_genotype accepts the ~3,300
# duplicate-cM markers as a target directly (they come out as identical terraced
# rows), so no unique-cM dedup of the TARGET (the per-family OBS stays cm-dedup'd,
# since obs$cm must be strictly increasing). Baseline (lambda=Inf) = the 51,004
# interpolated OLS scan, stam_gwas_scan_interpolated.csv.
#
# LB-IMPUTE PARAMS (map-aware; owner decision 2026-07-06):
#  - unit = "cm": transition decays over v5 consensus cM, so local recombination-
#    rate variation (maize centromeric suppression) is captured — consistent with
#    the cM step-interpolation and RTIGER/nNIL's map-based operation.
#  - recombdist + drp: CALIBRATED (donor-fragment-Dice optimal on the BC1S4 sim,
#    scripts/02_calibrate.R) and read from calib_params.csv (drp=TRUE for RILs).
#  - err = 0.01 (read/allele error, = the read-sim error), genotypeerr = 0.05
#    (LB-Impute default), min_reads = 0L
#    (decode every family marker; zero-coverage markers carry a flat emission and
#    are filled by the distance transition -> complete rectangular block).
#  - Flat start (no design prior): LB-Impute has no state-frequency prior in its
#    transition, so the start distribution is left flat.
#  - Read model + seeds identical to the RTIGER sweep (deterministic per cell:
#    seed = 1000 + 100*family_index + lambda_index).
#
# Writes results/sim/teonam/stam_gwas_lbimpute_lambda<L>.csv (per lambda) and the
# combined long table stam_gwas_lbimpute_sweep.csv (adds `coverage`; Inf baseline).
#
# Run:  Rscript scripts/teonam_lbimpute_sweep.R --generate
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
source(file.path(ROOT, "scripts/map_tools.R")) # DEFAULT_TEONAM_MAP (native est.map)
OUTDIR <- file.path(ROOT, "results/sim/teonam")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

if (!("--generate" %in% commandArgs(TRUE))) {
  message("teonam_lbimpute_sweep.R: pass --generate to (re)compute the sweep CSVs.")
  quit(save = "no", status = 0)
}

LAMBDAS <- c(0.1, 0.2, 0.5, 1, 5, 10, 20)
# recombdist + drp are CALIBRATED (donor-fragment-Dice optimal on the BC1S4 sim,
# scripts/02_calibrate.R) and READ from calib_params.csv — do not hardcode.
cp <- fread(file.path(ROOT, "results/sim/calib_params.csv"))
RECOMBDIST <- as.numeric(cp$value[cp$key == "recombdist_star"]) # cM (unit-aware transition scale)
if (!isTRUE(is.finite(RECOMBDIST))) {
  stop("recombdist_star not in results/sim/calib_params.csv — run scripts/02_calibrate.R first")
}
DRP <- isTRUE(toupper(cp$value[cp$key == "lbimpute_drp"]) == "TRUE") # double-recomb penalty (RIL)
GENOERR <- 0.05 # LB-Impute genotypeerr default
THREADS <- max(1L, detectCores() - 2L)
READ_PARS <- list(pi_floor = 0, k_decay = 1, error = 0.01)
message(sprintf("LB-Impute: recombdist_star = %.4g cM, drp = %s (calib_params.csv)", RECOMBDIST, DRP))

# --- marker map + native GWAS union target grid ------------------------------
# GWAS grid = the FULL genotype set (51,004 markers): every genotyped marker with
# a v5 position, scanned at every coverage level (a blind test). Roster + v5
# physical positions come from data/teonam/markers_v5.tsv (the map-neutral v2->v5
# liftover). cM is taken ENTIRELY from the NATIVE TeoNAM est.map: native cM for the
# markers it placed, and for those it did not (quirky/non-Mendelian/unplaced) a cM
# interpolated ON THE NATIVE MAP via its monotone bp->cM Marey spline
# (nilHMM::bp_to_cm, Hyman, clamped) fit per chr to the placed markers.
mc <- fread(file.path(ROOT, "data/teonam/markers_v5.tsv")) # roster + v5 bp (liftover)
setnames(mc, "chr_v5", "chr")
nat_cm <- fread(file.path(ROOT, DEFAULT_TEONAM_MAP)) # native est.map: cM for placed markers
mc[, cm := nat_cm$cm[match(marker, nat_cm$marker)]] # native cM; NA where est.map didn't place it
# place est.map-unplaced markers on the NATIVE cM scale via a per-chr Marey spline
# fit on the placed markers (bp_to_cm splits by chr internally; a chr with <2
# placed markers has no spline and its unplaced markers stay NA).
placed <- mc[!is.na(cm), .(chr, bp = pos_v5, cm)]
fit_chr <- placed[, .N, by = chr][N >= 2L, chr]
to_cm <- bp_to_cm(placed[chr %in% fit_chr])
mc[is.na(cm) & chr %in% fit_chr, cm := to_cm(chr, pos_v5)]
cm_by <- setNames(mc$cm, mc$marker)
pos_by <- setNames(mc$pos_v5, mc$marker)

gcols <- names(fread(file.path(ROOT, "data/teonam/TeoNAM_genotype_clean.csv"), nrows = 0))[-(1:3)]
GWAS_MK <- intersect(gcols, mc$marker) # FULL 51,004 GWAS set
message(sprintf(
  "native cM coordinate: %d placed + %d native-Marey-spline = %d full GWAS markers",
  sum(GWAS_MK %in% nat_cm$marker), sum(!(GWAS_MK %in% nat_cm$marker)), length(GWAS_MK)
))
u <- mc[marker %in% GWAS_MK] # FULL GWAS set — duplicate cM kept
setorder(u, chr, cm)
target_df <- data.frame(chr = as.integer(u$chr), cm = as.numeric(u$cm))
union_markers <- u$marker
union_pos <- as.integer(pos_by[union_markers])
union_chr <- as.integer(u$chr)
message(sprintf("union grid (51K GWAS set): %d markers x 10 chromosomes", nrow(u)))

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
  # sort by (chr, pos): LB-Impute decodes runs in bp order and the cM transition
  # coordinate must be non-decreasing along it. cm-dedup keeps obs$cm strictly
  # increasing for the downstream interpolate_genotype.
  setorder(mt, chr, pos)
  mt <- mt[, .SD[!duplicated(cm)], by = chr]
  setorder(mt, chr, pos)
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

# --- one (family,lambda) cell: simulate -> LB-Impute -> interpolate to union --
recover_block <- function(fam, li) {
  lambda <- LAMBDAS[li]
  fi <- match(fam, names(fams))
  fd <- fam_data[[fam]]
  mt <- fd$mt
  D <- fd$D
  keys <- fd$keys
  M <- nrow(mt)
  N <- length(keys)

  set.seed(1000L + 100L * fi + li) # per-cell seed (identical to the RTIGER sweep)
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
    cm = rep(mt$cm, N),
    n_ref = n_ref, n_alt = n_alt
  )

  st <- call_states(long,
    caller = "lbimpute", unit = "cm", recombdist = RECOMBDIST,
    err = READ_PARS$error, genotypeerr = GENOERR, drp = DRP,
    min_reads = 0L, threads = 1L
  )
  W <- dcast(as.data.table(st), chr + pos ~ name, value.var = "state")
  W <- W[mt[, .(chr, pos)], on = c("chr", "pos")] # align rows to mt (bp order)
  R <- as.matrix(W[, keys, with = FALSE])
  storage.mode(R) <- "double"

  # interpolate on cm; obs must be cm-sorted & strictly increasing
  o <- order(mt$chr, mt$cm)
  block <- interpolate_genotype(R[o, , drop = FALSE],
    data.frame(chr = mt$chr[o], cm = mt$cm[o]),
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
  "LB-Impute sweep: %d (family,lambda) batch calls, %d threads ...",
  nrow(grid), THREADS
))
t0 <- Sys.time()
cells <- mclapply(seq_len(nrow(grid)), function(i) {
  recover_block(grid$fam[i], grid$li[i])
}, mc.cores = THREADS)
bad <- vapply(cells, function(x) inherits(x, "try-error") || is.null(x), logical(1))
if (any(bad)) {
  stop(
    "LB-Impute failed for cell(s): ",
    paste(which(bad), collapse = ", "), " -> ", cells[[which(bad)[1]]]
  )
}
message(sprintf(
  "  LB-Impute sweep done in %.1f min",
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
  if ("--save-matrix" %in% commandArgs(TRUE)) { # opt-in: emit the union matrix for downstream MLM/LOCO
    dir.create(file.path(OUTDIR, "cache"), recursive = TRUE, showWarnings = FALSE)
    saveRDS(
      list(G = G, markers = union_markers, chr = union_chr, bp = union_pos),
      file.path(OUTDIR, "cache", sprintf("geno_lb_native%s.rds", lambda))
    )
  }
  scan <- gwas_scan(G)[order(CHR, BP)]
  fwrite(scan, file.path(OUTDIR, sprintf("stam_gwas_lbimpute_native_lambda%s.csv", lambda)))
  scan[, coverage := lambda]
  sweep_list[[li]] <- scan
  message(sprintf(
    "  lambda=%-4g : %d markers, tb1 peak -log10P = %s, global max = %.1f",
    lambda, nrow(scan), tb1_peak(scan),
    max(-log10(scan[is.finite(P) & P > 0, P]))
  ))
}

# --- combined sweep + lambda=Inf baseline (51K interpolated OLS) --------------
baseline <- fread(file.path(ROOT, "data/teonam/stam_gwas_scan_interpolated.csv"))
baseline[, coverage := Inf]
message(sprintf(
  "  lambda=Inf  : %d markers, tb1 peak -log10P = %s (baseline)",
  nrow(baseline), tb1_peak(baseline)
))

sweep <- rbindlist(c(sweep_list, list(baseline)), use.names = TRUE)
fwrite(sweep, file.path(OUTDIR, "stam_gwas_lbimpute_native_sweep.csv"))
message(sprintf(
  "wrote %s (%d rows, %d coverage levels)",
  file.path(OUTDIR, "stam_gwas_lbimpute_native_sweep.csv"),
  nrow(sweep), uniqueN(sweep$coverage)
))
