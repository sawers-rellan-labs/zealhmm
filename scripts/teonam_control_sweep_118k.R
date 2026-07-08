#!/usr/bin/env Rscript
# =============================================================================
# STAM GWAS degradation sweep on the AUTHENTIC 118K panel (TeoNAM, GL+HWE control)
# 118K variant of scripts/teonam_control_sweep.R -- the no-HMM het-excess control:
# per-site GL call with an HWE posterior (a single ALT read -> HET via 2p(1-p)).
# TRUTH/GRID/BASELINE are the dense polarized 118K panel (see
# scripts/teonam_rtiger_sweep_118k.R for the shared design).
#   caller: call_gl(n_ref, n_alt, prior="hwe", af=<per-marker truth teosinte AF>,
#           error=0.01) -> 0/1/2, NA at zero depth. Covered markers vary per RIL,
#           so interpolation is PER-RIL (no HMM).
#
# Run:  Rscript scripts/teonam_control_sweep_118k.R --generate   # full 35-cell grid
#       Rscript scripts/teonam_control_sweep_118k.R --smoke      # 1 family x 1 lambda
# =============================================================================
suppressMessages({
  library(data.table)
  library(parallel)
  library(readxl)
  devtools::load_all("/Users/fvrodriguez/repos/nilhmm", quiet = TRUE)
})
ROOT <- "/Users/fvrodriguez/repos/zealhmm"
setwd(ROOT)
source(file.path(ROOT, "R/simulate.R"))
source(file.path(ROOT, "scripts/map_tools.R"))
OUTDIR <- file.path(ROOT, "results/sim/teonam")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

ARGS <- commandArgs(TRUE)
SMOKE <- "--smoke" %in% ARGS
if (!SMOKE && !("--generate" %in% ARGS)) {
  message("pass --generate (full grid) or --smoke (1 cell).")
  quit(save = "no", status = 0)
}

LAMBDAS <- if (SMOKE) c(1, Inf) else c(0.1, 0.2, 0.5, 1, 5, 10, 20, Inf) # Inf = perfect coverage
ERROR <- 0.01 # GL per-read error (matches read model)
THREADS <- max(1L, detectCores() - 2L)
READ_PARS <- list(pi_floor = 0, k_decay = 1, error = 0.01)

# --- 118K cM grid + dense polarized truth (shared with rtiger_118k) ----------
mc <- fread(file.path(ROOT, "data/teonam/markers_v5_gwas118k_cm.tsv"))
setnames(mc, "pos_v5", "pos")
setorder(mc, chr, cm)
mt_all <- mc[, .SD[!duplicated(cm)], by = chr]
setorder(mt_all, chr, cm)
u <- copy(mc)
setorder(u, chr, cm)
target_df <- data.frame(chr = as.integer(u$chr), cm = as.numeric(u$cm))
union_markers <- u$marker
union_pos <- as.integer(u$pos)
union_chr <- as.integer(u$chr)
mt_thin <- fread(file.path(ROOT, "data/teonam/markers_v5_gwas118k_cm_thin01.tsv")) # cached 0.1 cM inference grid
setnames(mt_thin, "pos_v5", "pos")
message(sprintf("118K grid: %d union markers (back-projection target) | inference grid %d markers @0.1 cM", nrow(u), nrow(mt_thin)))

g118 <- readRDS(file.path(ROOT, "data/teonam/teonam_gwas118k_dosage_polar.rds")) # AUTHENTIC per-SNP genotypes
dos <- g118$dos
FAMS <- c("TIL01", "TIL03", "TIL11", "TIL14", "TIL25")

load_family <- function(fam) {
  keys <- colnames(dos)[substr(colnames(dos), 1, 5) == fam]
  D <- dos[mt_thin$marker, keys, drop = FALSE] # truth at the 0.1 cM inference grid
  storage.mode(D) <- "double"
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
  af <- rowMeans(D) / 2 # per-marker truth teosinte AF (polarized: 2 = teosinte)
  list(mt = mt_thin, D = D, keys = keys, af = af)
}
message("loading families (authentic per-SNP truth on the 0.1 cM inference grid) ...")
fam_data <- lapply(FAMS, load_family)
names(fam_data) <- FAMS
for (f in FAMS) message(sprintf("  %s: %d markers x %d RILs (mean truth teosinte AF = %.3f)", f, nrow(fam_data[[f]]$mt), length(fam_data[[f]]$keys), mean(fam_data[[f]]$af)))

recover_block <- function(fam, li) {
  lambda <- LAMBDAS[li]
  fi <- match(fam, FAMS)
  fd <- fam_data[[fam]]
  mt <- fd$mt
  D <- fd$D
  keys <- fd$keys
  af <- fd$af
  M <- nrow(mt)
  N <- length(keys)
  set.seed(1000L + 100L * fi + li)
  if (is.infinite(lambda)) {
    DINF <- 100L
    p_alt <- c(0, 0.5, 1)[as.vector(D) + 1L]
    p_eff <- p_alt * (1 - READ_PARS$error) + (1 - p_alt) * READ_PARS$error
    av <- as.integer(round(DINF * p_eff))
    n_ref <- matrix(DINF - av, M, N)
    n_alt <- matrix(av, M, N)
  } else {
    ac <- .draw_counts(as.vector(D),
      lambda = lambda,
      pi_floor = READ_PARS$pi_floor, k_decay = READ_PARS$k_decay, error = READ_PARS$error
    )
    n_ref <- matrix(as.integer(ac$ref), M, N)
    n_alt <- matrix(as.integer(ac$alt), M, N)
  }
  # GL+HWE call: 0/1/2, NA at zero depth; af recycled column-wise over M x N.
  calls <- call_gl(n_ref, n_alt, prior = "hwe", af = af, error = ERROR)
  covered <- !is.na(calls)
  # per-RIL step-interpolation of the RIL's covered calls onto the union grid.
  block <- matrix(NA_real_, nrow(u), N, dimnames = list(union_markers, keys))
  for (j in seq_len(N)) {
    obs <- which(covered[, j])
    if (!length(obs)) next
    geno <- matrix(as.double(calls[obs, j]), ncol = 1L)
    obs_df <- data.frame(chr = mt$chr[obs], cm = mt$cm[obs])
    ok_chr <- unique(obs_df$chr)
    tsel <- target_df$chr %in% ok_chr
    block[tsel, j] <- interpolate_genotype(geno, obs_df, target_df[tsel, , drop = FALSE], mode = "step")[, 1L]
  }
  n_het <- sum(calls == 1L, na.rm = TRUE) # het CALLS (pre-interpolation)
  n_called <- sum(covered)
  list(lambda = lambda, fam = fam, block = block, n_het = n_het, n_called = n_called, n_cells = M * N)
}

grid <- expand.grid(fam = if (SMOKE) FAMS[1] else FAMS, li = seq_along(LAMBDAS), stringsAsFactors = FALSE)
message(sprintf("GL+HWE-118K control sweep: %d cells, %d threads ...", nrow(grid), THREADS))
t0 <- Sys.time()
cells <- mclapply(seq_len(nrow(grid)), function(i) recover_block(grid$fam[i], grid$li[i]), mc.cores = THREADS)
bad <- vapply(cells, function(x) inherits(x, "try-error") || is.null(x), logical(1))
if (any(bad)) stop("cell(s) failed: ", paste(which(bad), collapse = ", "), " -> ", cells[[which(bad)[1]]])
message(sprintf("  recover done in %.1f min", as.numeric(Sys.time() - t0, units = "mins")))

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
    ok2 <- !is.na(g)
    if (sum(ok2) < 20 || sd(g[ok2]) == 0) {
      return(NA_real_)
    }
    Xi <- cbind(Xr[ok2, , drop = FALSE], g[ok2])
    fit <- lm.fit(Xi, y[ok2])
    RSS1 <- sum(fit$residuals^2)
    RSS0i <- sum(lm.fit(Xr[ok2, , drop = FALSE], y[ok2])$residuals^2)
    df2 <- sum(ok2) - fit$rank
    if (df2 <= 0 || RSS1 <= 0) {
      return(NA_real_)
    }
    pf(((RSS0i - RSS1) / 1) / (RSS1 / df2), 1, df2, lower.tail = FALSE)
  }
  P <- unlist(mclapply(seq_len(nrow(Gm)), scan1, mc.cores = THREADS))
  data.table(SNP = rownames(G), CHR = union_chr, BP = union_pos, P = P)
}

tb1 <- data.table(chr = 1L, start = 272330564L)
tb1_peak <- function(scan) {
  w <- scan[CHR == tb1$chr & abs(BP - tb1$start) <= 5e5 & is.finite(P) & P > 0]
  if (!nrow(w)) NA_real_ else round(max(-log10(w$P)), 2)
}

sweep_list <- vector("list", length(LAMBDAS))
het_list <- vector("list", length(LAMBDAS))
for (li in seq_along(LAMBDAS)) {
  lambda <- LAMBDAS[li]
  idx <- which(grid$li == li)
  G <- do.call(cbind, lapply(idx, function(i) cells[[i]]$block))
  rownames(G) <- union_markers
  scan <- gwas_scan(G)[order(CHR, BP)]
  fwrite(scan, file.path(OUTDIR, sprintf("stam_gwas_control_118k_lambda%s.csv", lambda)))
  scan[, coverage := lambda]
  sweep_list[[li]] <- scan
  nh <- sum(vapply(idx, function(i) cells[[i]]$n_het, numeric(1)))
  nca <- sum(vapply(idx, function(i) cells[[i]]$n_called, numeric(1)))
  ncl <- sum(vapply(idx, function(i) cells[[i]]$n_cells, numeric(1)))
  het_list[[li]] <- data.table(
    coverage = lambda, het_frac = nh / nca, call_rate = nca / ncl, n_het = nh, n_called = nca
  )
  message(sprintf(
    "  lambda=%-4g : %d markers, tb1 peak -log10P = %s, global max = %.1f (het/called %.3f)",
    lambda, nrow(scan), tb1_peak(scan), max(-log10(scan[is.finite(P) & P > 0, P]), na.rm = TRUE), nh / nca
  ))
}

if (SMOKE) {
  message("smoke ok.")
  quit(save = "no", status = 0)
}

fwrite(rbindlist(het_list), file.path(OUTDIR, "stam_control_het_fraction_118k.csv"))
message("wrote stam_control_het_fraction_118k.csv")

sweep <- rbindlist(sweep_list, use.names = TRUE) # lambda=Inf ceiling already in sweep_list
fwrite(sweep, file.path(OUTDIR, "stam_gwas_control_118k_sweep.csv"))
message(sprintf("wrote %s (%d rows, %d coverage levels)", file.path(OUTDIR, "stam_gwas_control_118k_sweep.csv"), nrow(sweep), uniqueN(sweep$coverage)))
