#!/usr/bin/env Rscript
# =============================================================================
# STAM GWAS degradation sweep on the AUTHENTIC 118K panel (TeoNAM, RTIGER)
# simulate reads -> RTIGER ancestry-segment inference -> back-project -> GWAS.
#   * TRUTH = the AUTHENTIC 118K per-SNP genotypes (teonam_gwas118k_dosage_polar.rds,
#     W22<->teo polarized) -- the REAL genotypes with real per-SNP structure, NOT an
#     ancestry-smoothed mosaic. Reads are Poisson-sampled from these at each SNP.
#   * The smoothing / clean line geometry comes from the caller's ancestry inference
#     + back-projection (block output), NOT from pre-smoothing the truth.
#   * GRID = the full 118,514-marker v5 set on the native cM (markers_v5_gwas118k_cm.tsv).
#   * BASELINE (lambda=Inf) = this SAME caller at perfect coverage (deterministic
#     decisive counts) -> its own ancestry ceiling. Per-caller, produced in-pipeline;
#     no separate baseline file.
#
# Run:  Rscript scripts/teonam_rtiger_sweep_118k.R --generate      # full 35-cell grid
#       Rscript scripts/teonam_rtiger_sweep_118k.R --smoke         # 1 family x 1 lambda, timed
# =============================================================================
suppressMessages({
  library(data.table)
  library(parallel)
  library(readxl)
  devtools::load_all("/Users/fvrodriguez/repos/nilhmm", quiet = TRUE)
})
ROOT <- "/Users/fvrodriguez/repos/zealhmm"
setwd(ROOT)
source(file.path(ROOT, "scripts/logging.R"))
t0 <- Sys.time()
source(file.path(ROOT, "R/simulate.R")) # .draw_counts()
source(file.path(ROOT, "scripts/map_tools.R"))
source(file.path(ROOT, "scripts/emmax_qk.R")) # emmax_qk_scan (MLM Q+K, Chen Fig-4C)
OUTDIR <- file.path(ROOT, "results/sim/teonam")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

ARGS <- commandArgs(TRUE)
SMOKE <- "--smoke" %in% ARGS
if (!SMOKE && !("--generate" %in% ARGS)) {
  log_info("pass --generate (full grid) or --smoke (1 cell, timed).")
  quit(save = "no", status = 0)
}

LAMBDAS <- if (SMOKE) c(1, Inf) else c(0.1, 0.2, 0.5, 1, 5, 10, 20, Inf) # Inf = perfect-coverage ceiling
# PER-COVERAGE rigidity from the simulation calibration (scripts/teonam_rtiger_calib_bycov.R):
# BC1S4 ground truth on this 118K grid, rigidity chosen by MIN donor-fragment FDR (err-long),
# min_reads=1. rigidity* rises with coverage (16 -> 128) -- deep coverage exposes more
# per-marker detail, so a longer minimum run is needed to stay stable. Because it was tuned
# with min_reads=1 (rigidity counts COVERED markers), the sweep below also uses min_reads=1
# and fills uncovered markers back to the full grid by within-chromosome carry-forward.
calib <- fread(file.path(ROOT, "results/sim/teonam/rtiger_calib_bycov.csv"))
RIG_BY_COV <- setNames(as.integer(calib$rigidity), as.character(calib$coverage))
THREADS <- max(1L, detectCores() - 2L)
READ_PARS <- list(pi_floor = 0, k_decay = 1, error = 0.01)

# --- 118K cM grid (native est.map + Marey spline; built by teonam_gwas118k_cm_grid.R) ---
mc <- fread(file.path(ROOT, "data/teonam/markers_v5_gwas118k_cm.tsv")) # marker, chr, pos_v5, cm
setnames(mc, "pos_v5", "pos")
setorder(mc, chr, cm)
# per-chr strictly-increasing cM (dedup) -> the shared marker table for all families
mt_all <- mc[, .SD[!duplicated(cm)], by = chr]
setorder(mt_all, chr, cm)
# union target = full 118K set (duplicate cM kept as terraced target rows)
u <- copy(mc)
setorder(u, chr, cm)
target_df <- data.frame(chr = as.integer(u$chr), cm = as.numeric(u$cm))
union_markers <- u$marker
union_pos <- as.integer(u$pos)
union_chr <- as.integer(u$chr)

# --- ancestry-inference grid: the FULL 118K marker set, per chromosome (NO thinning,
# NO back-projection). RTIGER runs the HMM directly on every marker; the recovered
# states ARE the union genotypes for the GWAS. Slower than the 0.1 cM thin grid but
# removes the thin->union interpolation and the grid-density rigidity coupling.
mt_thin <- copy(mc)[, .(marker, chr, pos, cm)] # full union = inference grid
stopifnot(all(is.finite(RIG_BY_COV)))
log_info(
  "118K grid: RTIGER on the FULL %d markers/genome, per chromosome, min_reads=1; per-coverage rigidity = %s",
  nrow(mt_thin), paste(sprintf("%s:%d", names(RIG_BY_COV), RIG_BY_COV), collapse = " ")
)

# --- dense polarized 118K truth, split by family ------------------------------
g118 <- readRDS(file.path(ROOT, "data/teonam/teonam_gwas118k_dosage_polar.rds"))
dos <- g118$dos # AUTHENTIC per-SNP genotypes, markers x lines, 0/1/2 (0=W22,2=teo), ~2.7% NA
FAMS <- c("TIL01", "TIL03", "TIL11", "TIL14", "TIL25")

load_family <- function(fam) {
  keys <- colnames(dos)[substr(colnames(dos), 1, 5) == fam]
  D <- dos[mt_thin$marker, keys, drop = FALSE] # truth at the 0.1 cM INFERENCE grid x RILs
  storage.mode(D) <- "double"
  # dense truth, but the imputed HapMap still carries ~2.7% NA -> fill per RIL by
  # the family modal call so .draw_counts has a complete truth (same spirit as the
  # 51K family blocks being complete within family).
  if (anyNA(D)) {
    rm <- apply(D, 1, function(z) {
      z <- z[!is.na(z)]
      if (!length(z)) 0 else as.numeric(names(which.max(table(z))))
    })
    for (j in seq_len(ncol(D))) {
      na <- is.na(D[, j])
      if (any(na)) D[na, j] <- rm[na]
    }
  }
  list(mt = mt_thin, D = D, keys = keys)
}
log_info("loading families (authentic per-SNP truth on the 0.1 cM inference grid) ...")
fam_data <- lapply(FAMS, load_family)
names(fam_data) <- FAMS
for (f in FAMS) log_info("  %s: %d markers x %d RILs", f, nrow(fam_data[[f]]$mt), length(fam_data[[f]]$keys))

# Memory-safe recovery: fit the joint per-family RTIGER emission ONCE, then decode
# each chromosome in a separate low-memory worker (fit_rtiger + rtiger_fit=). RTIGER's
# HMM is per-chromosome, so this is identical to a whole-family call (verified,
# agent/verify_fit_reuse.R) at ~1/10 the peak memory. Families run SEQUENTIALLY (the
# fit is the memory-heavy step); chromosomes decode in parallel.
DECODE_CORES <- max(1L, min(detectCores() - 1L, 8L)) # per-chromosome workers (each ~1 GB)
CHRS <- sort(unique(mt_thin$chr))

# covered-marker reads (n>=1) for one family at coverage li -- same read model + seeds
# as before; filtering to covered markers gives the fit and the per-chr decode the
# identical min_reads=1 support.
build_reads <- function(fam, li) {
  lambda <- LAMBDAS[li]
  fi <- match(fam, FAMS)
  fd <- fam_data[[fam]]
  D <- fd$D
  keys <- fd$keys
  M <- nrow(fd$mt)
  N <- length(keys)
  set.seed(1000L + 100L * fi + li)
  if (is.infinite(lambda)) {
    DINF <- 100L
    p_alt <- c(0, 0.5, 1)[as.vector(D) + 1L]
    p_eff <- p_alt * (1 - READ_PARS$error) + (1 - p_alt) * READ_PARS$error
    n_alt <- as.integer(round(DINF * p_eff))
    n_ref <- DINF - n_alt
  } else {
    ac <- .draw_counts(as.vector(D),
      lambda = lambda,
      pi_floor = READ_PARS$pi_floor, k_decay = READ_PARS$k_decay, error = READ_PARS$error
    )
    n_ref <- as.integer(ac$ref)
    n_alt <- as.integer(ac$alt)
  }
  long <- data.table(
    name = rep(keys, each = M), chr = rep(fd$mt$chr, N),
    pos = rep(fd$mt$pos, N), n_ref = n_ref, n_alt = n_alt
  )
  long[n_ref + n_alt > 0L] # covered markers only (shared min_reads=1 support: fit + decode)
}

# one family -> full-grid recovered genotype block (markers x family RILs), integer.
recover_family <- function(fam, li) {
  keys <- fam_data[[fam]]$keys
  rig <- RIG_BY_COV[[if (is.infinite(LAMBDAS[li])) "Inf" else as.character(LAMBDAS[li])]]
  reads <- build_reads(fam, li)
  fit <- fit_rtiger(reads, rig, threads = 1L, seed = 1L) # joint emission, once per family
  blocks <- mclapply(CHRS, function(ch) {
    st <- call_states(reads[chr == ch],
      caller = "rtiger", rigidity = rig, rtiger_fit = fit, min_reads = 1L, threads = 1L
    )
    mtc <- mt_thin[chr == ch]
    W <- dcast(as.data.table(st), chr + pos ~ name, value.var = "state")
    W <- W[mtc[, .(chr, pos)], on = c("chr", "pos")] # full chr grid, NA at uncovered markers
    b <- as.matrix(W[, keys, with = FALSE])
    apply(b, 2L, function(v) nafill(nafill(v, "locf"), "nocb")) # carry-forward fill within chr
  }, mc.cores = DECODE_CORES)
  if (any(vapply(blocks, function(x) inherits(x, "try-error") || is.null(x), logical(1)))) {
    stop(sprintf("recover_family(%s, lambda=%s): a chromosome decode failed", fam, LAMBDAS[li]))
  }
  block <- do.call(rbind, blocks) # CHRS-ordered -> union (chr, cM) order
  storage.mode(block) <- "integer"
  block
}

# --- phenotype + GWAS scan (STAM ~ Family + marker, 1 df) --------------------
ph <- as.data.frame(read_excel(file.path(ROOT, "data/teonam/9250682/TeoNAM_1257RILs_22traits_phenotype_data.xlsx")))
names(ph)[1] <- "line"
stam_by <- setNames(ph$STAM, ph$line)
gwas_scan <- function(G) {
  y <- suppressWarnings(as.numeric(stam_by[colnames(G)]))
  fam <- factor(substr(colnames(G), 1, 5))
  ok <- !is.na(y)
  y <- y[ok]
  fam <- droplevels(fam[ok])
  Gm <- G[, ok, drop = FALSE]
  Xr <- if (nlevels(fam) > 1) model.matrix(~fam) else matrix(1, length(y), 1)
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

tb1 <- data.table(chr = 1L, start = 272330564L) # tb1 (Zm00001eb054440), v5
tb1_peak <- function(scan) {
  w <- scan[CHR == tb1$chr & abs(BP - tb1$start) <= 5e5 & is.finite(P) & P > 0]
  if (!nrow(w)) NA_real_ else round(max(-log10(w$P)), 2)
}

FAM_USE <- if (SMOKE) FAMS[1] else FAMS
log_info(
  "RTIGER-118K sweep: %d coverages x %d families; fit-once-per-family, per-chromosome decode (%d cores)",
  length(LAMBDAS), length(FAM_USE), DECODE_CORES
)
t0 <- Sys.time()
sweep_list <- mlm_list <- vector("list", length(LAMBDAS))
for (li in seq_along(LAMBDAS)) {
  lambda <- LAMBDAS[li]
  covlab <- if (is.infinite(lambda)) "Inf" else as.character(lambda)
  tl <- Sys.time()
  G <- do.call(cbind, lapply(FAM_USE, function(fam) recover_family(fam, li))) # families sequential
  rownames(G) <- union_markers
  scan <- gwas_scan(G)[order(CHR, BP)] # OLS (Family + marker)
  fwrite(scan, file.path(OUTDIR, sprintf("stam_gwas_rtiger_118k_lambda%s.csv", lambda)))
  scan[, coverage := lambda]
  sweep_list[[li]] <- scan
  null_li <- readRDS(file.path(ROOT, sprintf("data/teonam/mlm_null_118k_l%s.rds", lambda))) # coverage-matched Family+K (GL-dosage of downsampled reads)
  mlm <- emmax_qk_scan(G, null_li, union_chr, union_pos)[order(CHR, BP)] # MLM (Family+K)
  mlm[, coverage := lambda]
  mlm_list[[li]] <- mlm
  rm(G)
  invisible(gc())
  log_info(
    "  [%d/%d] lambda=%-4s (rig=%d, %.1f min): OLS tb1 %s / MLM tb1 %s (OLS max %.1f)",
    li, length(LAMBDAS), covlab, RIG_BY_COV[[covlab]], as.numeric(difftime(Sys.time(), tl, units = "mins")),
    tb1_peak(scan), tb1_peak(mlm), max(-log10(scan[is.finite(P) & P > 0, P]))
  )
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  log_info(">>> %d/%d done | elapsed %.1f min | avg %.1f min | ETA ~%.1f min remaining", li, length(LAMBDAS), el, el / li, (el / li) * (length(LAMBDAS) - li))
}

if (SMOKE) {
  log_info("smoke ok.")
  quit(save = "no", status = 0)
}

fwrite(rbindlist(sweep_list, use.names = TRUE), file.path(OUTDIR, "stam_gwas_rtiger_118k_sweep.csv"))
fwrite(rbindlist(mlm_list, use.names = TRUE), file.path(OUTDIR, "stam_gwas_rtiger_118k_mlm_sweep.csv"))
log_info("%s", paste0("wrote OLS + MLM(Q+K) sweeps, ", uniqueN(rbindlist(sweep_list)$coverage), " coverage levels"))
