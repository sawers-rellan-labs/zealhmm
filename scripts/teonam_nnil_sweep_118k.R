#!/usr/bin/env Rscript
# =============================================================================
# STAM GWAS degradation sweep on the AUTHENTIC 118K panel (TeoNAM, nNIL)
# simulate reads -> nNIL ancestry inference -> assemble -> GWAS. Memory-safe
# rewrite mirroring scripts/teonam_rtiger_sweep_118k.R (shared design there):
#   * TRUTH = authentic 118K per-SNP genotypes (teonam_gwas118k_dosage_polar.rds).
#   * GRID  = the full 118,514-marker v5 set on native cM; nNIL runs on EVERY marker
#     per chromosome (NO 0.1 cM thinning, NO back-projection -- the recovered states
#     ARE the union genotypes, same as the RTIGER 118K sweep).
#   * nNIL has NO EM (count emission fixed from err + BC1S4 priors, geometric
#     duration), so -- unlike RTIGER -- there is NO fit-once step: decode each
#     (family, chromosome) directly. min_reads=1 (covered markers) + within-chr
#     carry-forward fill to the full grid, matching RTIGER's assembly.
#   * rrate: nNIL is the MAP-BLIND baseline and its rrate is NOT worth per-coverage
#     tuning -- calibration (teonam_nnil_calib_bycov.R) showed donor_frag_FDR is a
#     flat low-rrate plateau (recall ~1 throughout), so a SINGLE fixed rrate=3.3e-5
#     (the all-coverage plateau floor) is used for every coverage. See memory
#     nnil-rrate-flat-plateau.
#
# Run:  Rscript scripts/teonam_nnil_sweep_118k.R --generate   # full grid
#       Rscript scripts/teonam_nnil_sweep_118k.R --smoke      # 1 family x {1, Inf}, timed
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
source(file.path(ROOT, "R/simulate.R")) # .draw_counts()
for (f in list.files(file.path(ROOT, "R"), "\\.R$", full.names = TRUE)) source(f) # single_locus_expectation()
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
# FIXED plateau rrate (see header + memory nnil-rrate-flat-plateau): min-FDR is flat
# across the low-rrate plateau with recall ~1, so per-coverage tuning is noise; use the
# all-coverage plateau floor for every coverage.
RRATE <- 3.3e-5
# BC1S4 start/transition priors -- SAME source as the rrate calibration
# (single_locus_expectation(1,4)); consistency with the tuning run, not the pre-calib
# sweep's Chen observed freqs. f_1 = HET, f_2 = ALT (homozygous donor).
EXP <- single_locus_expectation(1L, 4L)
F1 <- as.numeric(EXP["HET"])
F2 <- as.numeric(EXP["ALT"])
THREADS <- max(1L, detectCores() - 2L)
READ_PARS <- list(pi_floor = 0, k_decay = 1, error = 0.01)
log_info("nNIL-118K: FIXED rrate = %.5g (plateau floor), f_1 = %.3f, f_2 = %.3f", RRATE, F1, F2)

# --- 118K cM grid (native est.map + Marey spline) ----------------------------
mc <- fread(file.path(ROOT, "data/teonam/markers_v5_gwas118k_cm.tsv")) # marker, chr, pos_v5, cm
setnames(mc, "pos_v5", "pos")
setorder(mc, chr, cm)
u <- copy(mc)
setorder(u, chr, cm)
union_markers <- u$marker
union_pos <- as.integer(u$pos)
union_chr <- as.integer(u$chr)
mt_thin <- copy(mc)[, .(marker, chr, pos, cm)] # full union = inference grid (no thinning)
CHRS <- sort(unique(mt_thin$chr))
log_info(
  "118K grid: nNIL on the FULL %d markers/genome, per chromosome, min_reads=1 + carry-forward fill",
  nrow(mt_thin)
)

# --- dense polarized 118K truth, split by family (shared with rtiger_118k) ----
g118 <- readRDS(file.path(ROOT, "data/teonam/teonam_gwas118k_dosage_polar.rds"))
dos <- g118$dos # AUTHENTIC per-SNP genotypes, markers x lines, 0/1/2 (0=W22,2=teo), ~2.7% NA
FAMS <- c("TIL01", "TIL03", "TIL11", "TIL14", "TIL25")

load_family <- function(fam) {
  keys <- colnames(dos)[substr(colnames(dos), 1, 5) == fam]
  D <- dos[mt_thin$marker, keys, drop = FALSE] # truth at the full 118K grid x RILs
  storage.mode(D) <- "double"
  if (anyNA(D)) { # imputed HapMap still ~2.7% NA -> fill per RIL by the family modal call
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
log_info("loading families (authentic per-SNP truth on the full 118K grid) ...")
fam_data <- lapply(FAMS, load_family)
names(fam_data) <- FAMS
for (f in FAMS) log_info("  %s: %d markers x %d RILs", f, nrow(fam_data[[f]]$mt), length(fam_data[[f]]$keys))

DECODE_CORES <- max(1L, min(detectCores() - 1L, 8L)) # per-chromosome workers (each ~1 GB)

# covered-marker reads (n>=1) for one family at coverage li -- same read model + seeds
# as the RTIGER sweep; filtering to covered markers gives the per-chr decode the
# min_reads=1 support the rrate was calibrated at.
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
  long[n_ref + n_alt > 0L] # covered markers only (min_reads=1 support)
}

# one family -> full-grid recovered genotype block (markers x family RILs), integer.
# Decode per chromosome in low-memory workers; families run sequentially.
recover_family <- function(fam, li) {
  keys <- fam_data[[fam]]$keys
  reads <- build_reads(fam, li)
  blocks <- mclapply(CHRS, function(ch) {
    st <- call_states(reads[chr == ch],
      caller = "nnil", f_1 = F1, f_2 = F2, rrate = RRATE,
      err = READ_PARS$error, min_reads = 1L, threads = 1L
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
  "nNIL-118K sweep: %d coverages x %d families; per-chromosome decode (%d cores), families sequential",
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
  fwrite(scan, file.path(OUTDIR, sprintf("stam_gwas_nnil_118k_lambda%s.csv", lambda)))
  scan[, coverage := lambda]
  sweep_list[[li]] <- scan
  null_li <- readRDS(file.path(ROOT, sprintf("data/teonam/mlm_null_118k_l%s.rds", lambda))) # coverage-matched Family+K
  mlm <- emmax_qk_scan(G, null_li, union_chr, union_pos)[order(CHR, BP)] # MLM (Family+K)
  mlm[, coverage := lambda]
  mlm_list[[li]] <- mlm
  rm(G)
  invisible(gc())
  log_info(
    "  [%d/%d] lambda=%-4s (%.1f min): OLS tb1 %s / MLM tb1 %s (OLS max %.1f)",
    li, length(LAMBDAS), covlab, as.numeric(difftime(Sys.time(), tl, units = "mins")),
    tb1_peak(scan), tb1_peak(mlm), max(-log10(scan[is.finite(P) & P > 0, P]))
  )
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  log_info(">>> %d/%d done | elapsed %.1f min | avg %.1f min | ETA ~%.1f min remaining", li, length(LAMBDAS), el, el / li, (el / li) * (length(LAMBDAS) - li))
}

if (SMOKE) {
  log_info("smoke ok.")
  quit(save = "no", status = 0)
}

fwrite(rbindlist(sweep_list, use.names = TRUE), file.path(OUTDIR, "stam_gwas_nnil_118k_sweep.csv"))
fwrite(rbindlist(mlm_list, use.names = TRUE), file.path(OUTDIR, "stam_gwas_nnil_118k_mlm_sweep.csv"))
log_info("%s", paste0("wrote OLS + MLM(Q+K) sweeps, ", uniqueN(rbindlist(sweep_list)$coverage), " coverage levels"))
