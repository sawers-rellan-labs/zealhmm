#!/usr/bin/env Rscript
# =============================================================================
# STAM GWAS degradation sweep on the AUTHENTIC 118K panel (TeoNAM, binhmm)
# simulate reads -> binhmm (1 Mb bins) ancestry -> back-project -> GWAS. Fifth
# caller of the coverage-degradation composite, alongside control/nNIL/RTIGER/
# LB-Impute (shared design in scripts/teonam_rtiger_sweep_118k.R).
#   * TRUTH = authentic 118K per-SNP genotypes (teonam_gwas118k_dosage_polar.rds).
#   * binhmm bins the COVERAGE-SUBSAMPLED simulated reads into 1 Mb windows
#     (per-bin ALT_FREQ), then genotypes each bin with the default `gauss` backend
#     (anchored 3-state Gaussian-emission HMM). NO tuning: bin_size = 1 Mb fixed
#     (there is no per-coverage calibration -- unlike rtiger/nnil/lbimpute).
#   * binhmm's gauss backend anchors REF from the sample's WHOLE-GENOME bin set, so
#     it decodes per FAMILY (all chromosomes together), not per-chromosome. It is
#     light (~2,100 bins/sample vs 118K markers), so no per-chr memory split is
#     needed; families run sequentially.
#   * ASSEMBLE: each 118K marker inherits its 1 Mb bin's state (bin = pos %/% 1e6);
#     bins with no covered markers (low coverage) are filled by within-chromosome
#     carry-forward, matching the other callers' full-grid assembly.
#   * min_reads is a no-op for binhmm (it bins alt/total; zero-read markers add
#     nothing to a bin's ALT_FREQ). Reads are filtered to covered markers (lighter).
#
# Run:  Rscript scripts/teonam_binhmm_sweep_118k.R --generate   # full grid
#       Rscript scripts/teonam_binhmm_sweep_118k.R --smoke      # 1 family x {1, Inf}, timed
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
source(file.path(ROOT, "scripts/emmax_qk.R")) # emmax_qk_scan (MLM Family+K)
OUTDIR <- file.path(ROOT, "results/sim/teonam")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

ARGS <- commandArgs(TRUE)
SMOKE <- "--smoke" %in% ARGS
if (!SMOKE && !("--generate" %in% ARGS)) {
  log_info("pass --generate (full grid) or --smoke (1 cell, timed).")
  quit(save = "no", status = 0)
}

LAMBDAS <- if (SMOKE) c(1, Inf) else c(0.1, 0.2, 0.5, 1, 5, 10, 20, Inf) # Inf = perfect-coverage ceiling
BIN_SIZE <- 1e6 # 1 Mb bins, FIXED (no tuning) -- binhmm's segmentation scale
# BC1S4 start priors (start distribution for binhmm's sticky transition) -- same
# source as the nNIL/LB-Impute calibrations: single_locus_expectation(1,4).
EXP <- single_locus_expectation(1L, 4L)
F1 <- as.numeric(EXP["HET"])
F2 <- as.numeric(EXP["ALT"])
THREADS <- max(1L, detectCores() - 2L)
READ_PARS <- list(pi_floor = 0, k_decay = 1, error = 0.01)
log_info("binhmm-118K: bin_size = %.0f bp (1 Mb, gauss backend), f_1 = %.3f, f_2 = %.3f", BIN_SIZE, F1, F2)

# --- 118K cM grid (native est.map + Marey spline) ----------------------------
mc <- fread(file.path(ROOT, "data/teonam/markers_v5_gwas118k_cm.tsv")) # marker, chr, pos_v5, cm
setnames(mc, "pos_v5", "pos")
setorder(mc, chr, cm)
u <- copy(mc)
setorder(u, chr, cm)
union_markers <- u$marker
union_pos <- as.integer(u$pos)
union_chr <- as.integer(u$chr)
mt_thin <- copy(mc)[, .(marker, chr, pos, cm)] # full union = back-projection target
CHRS <- sort(unique(mt_thin$chr))
log_info("118K grid: binhmm on %d markers -> 1 Mb bins -> back-project to the full grid", nrow(mt_thin))

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

# covered-marker reads (n>=1) for one family at coverage li -- same read model + seeds
# as the other sweeps. binhmm bins alt/total, so covered markers carry all the signal
# (zero-read markers add nothing to a bin's ALT_FREQ); min_reads is a no-op.
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
  long[n_ref + n_alt > 0L] # covered markers only
}

# one family -> full-grid recovered genotype block (markers x family RILs), integer.
# binhmm decodes the WHOLE family (gauss REF anchoring is per-sample, all chromosomes),
# then each 118K marker inherits its 1 Mb bin's state; empty bins filled by within-chr
# carry-forward. Families run sequentially (binhmm is light; the read table is the cost).
recover_family <- function(fam, li) {
  keys <- fam_data[[fam]]$keys
  reads <- build_reads(fam, li)
  st <- as.data.table(call_states(reads,
    caller = "binhmm", bin_size = BIN_SIZE, f_1 = F1, f_2 = F2
  )) # per-(sample, bin): name, chr, pos(=bin index), start_bp, end_bp, state
  Bw <- dcast(
    st[, .(name, chr = as.integer(chr), bin = as.integer(pos), state = as.integer(state))],
    chr + bin ~ name,
    value.var = "state"
  ) # bins x RILs (NA where a sample has no informative markers in a bin)
  miss <- setdiff(keys, names(Bw)) # a RIL with no bins at all (extreme low coverage)
  if (length(miss)) Bw[, (miss) := NA_integer_]
  mk <- mt_thin[, .(chr = as.integer(chr), bin = as.integer(pos %/% BIN_SIZE))]
  idx <- Bw[mk, on = c("chr", "bin"), which = TRUE] # bin row per union marker (NA if bin absent)
  block <- as.matrix(Bw[idx, ..keys]) # markers x RILs (union order), NA at missing bins
  storage.mode(block) <- "double"
  for (ch in CHRS) { # carry-forward fill within chromosome (union order = chr,cm asc = bp asc)
    r <- which(mt_thin$chr == ch)
    block[r, ] <- apply(block[r, , drop = FALSE], 2L, function(v) nafill(nafill(v, "locf"), "nocb"))
  }
  nna <- sum(is.na(block))
  if (nna) { # residual NA only for an all-empty RIL/chr -> neutral REF (0), logged
    log_warn("  recover_family(%s, lambda=%s): %d residual NA cells -> REF(0)", fam, LAMBDAS[li], nna)
    block[is.na(block)] <- 0
  }
  storage.mode(block) <- "integer"
  block
}

# --- phenotype + GWAS scan (STAM ~ Family + marker, 1 df) --------------------
ph <- as.data.frame(read_excel(file.path(ROOT, "data/teonam/9250682/TeoNAM_1257RILs_22traits_phenotype_data.xlsx")))
names(ph)[1] <- "line"
TRAIT <- toupper(Sys.getenv("TRAIT", "STAM"))
TTAG <- tolower(TRAIT) # phenotype col; STAM default, e.g. DTA
NTAG <- if (TRAIT == "STAM") "" else paste0("_", TTAG) # mlm-null trait tag ("" keeps STAM paths)
if (!TRAIT %in% names(ph)) stop("TRAIT '", TRAIT, "' is not a phenotype column")
stam_by <- setNames(ph[[TRAIT]], ph$line)
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
  "binhmm-118K sweep: %d coverages x %d families; whole-family bin decode, families sequential",
  length(LAMBDAS), length(FAM_USE)
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
  fwrite(scan, file.path(OUTDIR, sprintf("%s_gwas_binhmm_118k_lambda%s.csv", TTAG, lambda)))
  scan[, coverage := lambda]
  sweep_list[[li]] <- scan
  null_li <- readRDS(file.path(ROOT, sprintf("data/teonam/mlm_null%s_118k_l%s.rds", NTAG, lambda))) # coverage-matched Family+K
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

fwrite(rbindlist(sweep_list, use.names = TRUE), file.path(OUTDIR, sprintf("%s_gwas_binhmm_118k_sweep.csv", TTAG)))
fwrite(rbindlist(mlm_list, use.names = TRUE), file.path(OUTDIR, sprintf("%s_gwas_binhmm_118k_mlm_sweep.csv", TTAG)))
log_info("%s", paste0("wrote OLS + MLM(Family+K) sweeps, ", uniqueN(rbindlist(sweep_list)$coverage), " coverage levels"))
