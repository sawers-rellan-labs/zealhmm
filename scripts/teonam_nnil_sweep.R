#!/usr/bin/env Rscript
# =============================================================================
# STAM GWAS degradation sweep: nNIL ancestry vs sequencing coverage (TeoNAM)
# Plan: agent/teonam-control-sweep-plan.md (Part 3) -- the middle panel of the
#       control / nNIL / RTIGER composite. A trivial variant of the RTIGER sweep
#       (scripts/teonam_rtiger_sweep.R): swap caller = "rtiger" for caller = "nnil".
# -----------------------------------------------------------------------------
# For each TeoNAM family x coverage lambda:
#   1. simulate low-coverage reads from the real (truth) 0/1/2 genotypes
#      [R/simulate.R::.draw_counts(pi_floor=0, k_decay=1, error=0.01)]
#   2. nNIL-call ancestry [call_ancestry(caller="nnil", rrate=rrate_star, err=0.01)]
#      -- count emission + geometric duration (Holland's nNIL).
#   3. step-interpolate the recovered per-family block onto the union cM grid
#      [nilHMM::interpolate_genotype(mode="step")]
# Then assemble the union matrix (51,004 markers x 1,237 lines; full GWAS set) at each lambda
# and run the same GWAS scan as the baseline (STAM ~ Family + marker, 1-df F).
#
# Writes results/sim/teonam/stam_gwas_nnil_lambda<L>.csv (per lambda) and the
# combined long table stam_gwas_nnil_sweep.csv (adds `coverage`; includes the
# lambda=Inf baseline rows copied from data/teonam/stam_gwas_scan_interpolated.csv).
#
# KEY DECISIONS:
#  - Coverage grid: {0.1, 0.2, 0.5, 1, 5, 10, 20} (+ lambda=Inf baseline = panel C).
#  - rrate = rrate_star, READ from results/sim/calib_params.csv (do not invent).
#  - Design priors: nNIL (geometric duration) CONSUMES f_1/f_2 (unlike RTIGER,
#    which fits its own start freqs). We use the Chen 2019 OBSERVED frequencies,
#    not the Mendelian breeding_prior("BC1S4") expectation: real TeoNAM shows het
#    excess / teosinte-hom deficit (~8% het / ~15% teosinte-hom), so f_1 = 0.08,
#    f_2 = 0.15.
#  - Read model: pi_floor=0, k_decay=1, error=0.01; 1 replicate per (family,lambda),
#    RNG seed = 1000 + 100*family_index + lambda_index (identical scheme to the
#    RTIGER + control sweeps, so all three degrade the same truth mosaics).
#  - min_reads = 0L: decode EVERY family marker -> a COMPLETE rectangular per-family
#    block for step-interpolation (uncovered markers carry a flat count emission and
#    are filled by the geometric HMM from neighbours), matching the RTIGER sweep.
#
# Run:  Rscript scripts/teonam_nnil_sweep.R --generate
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
  message("teonam_nnil_sweep.R: pass --generate to (re)compute the sweep CSVs.")
  quit(save = "no", status = 0)
}

LAMBDAS <- c(0.1, 0.2, 0.5, 1, 5, 10, 20)
cp <- fread(file.path(ROOT, "results/sim/calib_params.csv"))
RRATE <- as.numeric(cp$value[cp$key == "rrate_star"]) # calibrated nNIL rrate
if (!is.finite(RRATE)) stop("rrate_star not found in results/sim/calib_params.csv")
F1 <- 0.08
F2 <- 0.15 # TeoNAM BC1S4 (Chen 2019 obs freqs)
THREADS <- max(1L, detectCores() - 2L)
READ_PARS <- list(pi_floor = 0, k_decay = 1, error = 0.01)
message(sprintf("nNIL: rrate_star = %.5g (calib_params.csv), f_1 = %.2f, f_2 = %.2f", RRATE, F1, F2))

# --- marker map + union target grid (identical to the other sweeps) ----------
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
# fit on the placed markers (bp_to_cm splits by chr internally; a chr with
# <2 placed markers has no spline and its unplaced markers stay NA).
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
u <- mc[marker %in% GWAS_MK] # FULL GWAS set (51K) — duplicate cM kept as terraced target rows
setorder(u, chr, cm)
target_df <- data.frame(chr = as.integer(u$chr), cm = as.numeric(u$cm))
union_markers <- u$marker
union_pos <- as.integer(pos_by[union_markers])
union_chr <- as.integer(u$chr)
message(sprintf("union grid (51K GWAS set): %d markers x 10 chromosomes", nrow(u)))

# --- per-family: genotypes, robust keys, cm-dedup marker table ---------------
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

# --- one (family,lambda) cell: simulate -> nNIL -> interpolate to union -------
recover_block <- function(fam, li) {
  lambda <- LAMBDAS[li]
  fi <- match(fam, names(fams))
  fd <- fam_data[[fam]]
  mt <- fd$mt
  D <- fd$D
  keys <- fd$keys
  M <- nrow(mt)
  N <- length(keys)

  set.seed(1000L + 100L * fi + li) # per-cell seed (as RTIGER sweep)
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
    caller = "nnil", f_1 = F1, f_2 = F2,
    rrate = RRATE, err = READ_PARS$error,
    min_reads = 0L, threads = 1L
  )
  W <- dcast(as.data.table(st), chr + pos ~ name, value.var = "state")
  W <- W[mt[, .(chr, pos)], on = c("chr", "pos")] # align rows to mt (cm order)
  R <- as.matrix(W[, keys, with = FALSE])
  storage.mode(R) <- "double"

  block <- interpolate_genotype(R, data.frame(chr = mt$chr, cm = mt$cm),
    target_df,
    mode = "step"
  )
  colnames(block) <- keys
  n_het <- sum(R == 1L)
  n_cells <- length(R)
  list(lambda = lambda, fam = fam, block = block, n_het = n_het, n_cells = n_cells)
}

grid <- expand.grid(
  fam = names(fams), li = seq_along(LAMBDAS),
  stringsAsFactors = FALSE
)
message(sprintf(
  "nNIL sweep: %d (family,lambda) batch calls, %d threads ...",
  nrow(grid), THREADS
))
t0 <- Sys.time()
cells <- mclapply(seq_len(nrow(grid)), function(i) {
  recover_block(grid$fam[i], grid$li[i])
}, mc.cores = THREADS)
bad <- vapply(cells, function(x) inherits(x, "try-error") || is.null(x), logical(1))
if (any(bad)) {
  stop(
    "nNIL failed for cell(s): ",
    paste(which(bad), collapse = ", "), " -> ", cells[[which(bad)[1]]]
  )
}
message(sprintf(
  "  nNIL sweep done in %.1f min",
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

# --- tb1 peak helper ---------------------------------------------------------
tb1 <- data.table(chr = 1L, start = 272330564L) # tb1 (Zm00001eb054440)
tb1_peak <- function(scan) {
  w <- scan[CHR == tb1$chr & abs(BP - tb1$start) <= 5e5 & is.finite(P) & P > 0]
  if (!nrow(w)) {
    return(NA_real_)
  }
  round(max(-log10(w$P)), 2)
}

# --- per-lambda: cbind family blocks, scan, write ----------------------------
sweep_list <- vector("list", length(LAMBDAS))
het_list <- vector("list", length(LAMBDAS))
for (li in seq_along(LAMBDAS)) {
  lambda <- LAMBDAS[li]
  idx <- which(grid$li == li)
  blocks <- lapply(idx, function(i) cells[[i]]$block) # one per family
  G <- do.call(cbind, blocks)
  rownames(G) <- union_markers
  scan <- gwas_scan(G)[order(CHR, BP)]
  fwrite(scan, file.path(OUTDIR, sprintf("stam_gwas_nnil_native_lambda%s.csv", lambda)))
  scan[, coverage := lambda]
  sweep_list[[li]] <- scan

  n_het <- sum(vapply(idx, function(i) cells[[i]]$n_het, numeric(1)))
  n_cells <- sum(vapply(idx, function(i) cells[[i]]$n_cells, numeric(1)))
  het_list[[li]] <- data.table(
    coverage = lambda, het_frac = n_het / n_cells,
    n_het = n_het, n_cells = n_cells
  )
  message(sprintf(
    "  lambda=%-4g : %d markers, het-frac = %.3f, tb1 peak -log10P = %s, global max = %.1f",
    lambda, nrow(scan), n_het / n_cells, tb1_peak(scan),
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
fwrite(sweep, file.path(OUTDIR, "stam_gwas_nnil_native_sweep.csv"))
message(sprintf(
  "wrote %s (%d rows, %d coverage levels)",
  file.path(OUTDIR, "stam_gwas_nnil_native_sweep.csv"),
  nrow(sweep), uniqueN(sweep$coverage)
))

het <- rbindlist(het_list)
fwrite(het, file.path(OUTDIR, "stam_nnil_native_het_fraction.csv"))
message("wrote stam_nnil_native_het_fraction.csv:")
print(het)
