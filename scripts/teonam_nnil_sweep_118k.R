#!/usr/bin/env Rscript
# =============================================================================
# STAM GWAS degradation sweep on the AUTHENTIC 118K panel (TeoNAM, nNIL)
# 118K variant of scripts/teonam_nnil_sweep.R -- identical simulate -> call ->
# interpolate -> GWAS pipeline, but TRUTH/GRID/BASELINE are the dense polarized
# 118K panel (see scripts/teonam_rtiger_sweep_118k.R for the shared design).
#   caller: call_states(caller="nnil", f_1, f_2, rrate=rrate_star) -- count
#           emission + geometric duration (Holland's nNIL). min_cov=0L.
#
# Run:  Rscript scripts/teonam_nnil_sweep_118k.R --generate   # full 35-cell grid
#       Rscript scripts/teonam_nnil_sweep_118k.R --smoke      # 1 family x 1 lambda
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

LAMBDAS <- if (SMOKE) c(1) else c(0.1, 0.2, 0.5, 1, 5, 10, 20)
cp <- fread(file.path(ROOT, "results/sim/calib_params.csv"))
RRATE <- as.numeric(cp$value[cp$key == "rrate_star"]) # calibrated nNIL rrate
if (!is.finite(RRATE)) stop("rrate_star not in results/sim/calib_params.csv")
F1 <- 0.08
F2 <- 0.15 # TeoNAM BC1S4 (Chen 2019 obs freqs)
THREADS <- max(1L, detectCores() - 2L)
READ_PARS <- list(pi_floor = 0, k_decay = 1, error = 0.01)
message(sprintf("nNIL-118K: rrate_star = %.5g, f_1 = %.2f, f_2 = %.2f", RRATE, F1, F2))

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
message(sprintf("118K grid: %d union markers | %d strictly-increasing-cM target rows", nrow(u), nrow(mt_all)))

g118 <- readRDS(file.path(ROOT, "data/teonam/teonam_gwas118k_dosage_polar.rds"))
dos <- g118$dos
FAMS <- c("TIL01", "TIL03", "TIL11", "TIL14", "TIL25")

load_family <- function(fam) {
  keys <- colnames(dos)[substr(colnames(dos), 1, 5) == fam]
  D <- dos[mt_all$marker, keys, drop = FALSE]
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
  list(mt = mt_all, D = D, keys = keys)
}
message("loading families (dense 118K truth) ...")
fam_data <- lapply(FAMS, load_family)
names(fam_data) <- FAMS
for (f in FAMS) message(sprintf("  %s: %d markers x %d RILs", f, nrow(fam_data[[f]]$mt), length(fam_data[[f]]$keys)))

recover_block <- function(fam, li) {
  lambda <- LAMBDAS[li]
  fi <- match(fam, FAMS)
  fd <- fam_data[[fam]]
  mt <- fd$mt
  D <- fd$D
  keys <- fd$keys
  M <- nrow(mt)
  N <- length(keys)
  set.seed(1000L + 100L * fi + li)
  ac <- .draw_counts(as.vector(D),
    lambda = lambda,
    pi_floor = READ_PARS$pi_floor, k_decay = READ_PARS$k_decay, error = READ_PARS$error
  )
  long <- data.table(
    name = rep(keys, each = M), chr = rep(mt$chr, N),
    pos = rep(mt$pos, N), n_ref = as.integer(ac$ref), n_alt = as.integer(ac$alt)
  )
  st <- call_states(long,
    caller = "nnil", f_1 = F1, f_2 = F2, rrate = RRATE,
    err = READ_PARS$error, min_cov = 0L, threads = 1L
  )
  W <- dcast(as.data.table(st), chr + pos ~ name, value.var = "state")
  W <- W[mt[, .(chr, pos)], on = c("chr", "pos")]
  R <- as.matrix(W[, keys, with = FALSE])
  storage.mode(R) <- "double"
  block <- interpolate_genotype(R, data.frame(chr = mt$chr, cm = mt$cm), target_df, mode = "step")
  colnames(block) <- keys
  list(lambda = lambda, fam = fam, block = block)
}

grid <- expand.grid(fam = if (SMOKE) FAMS[1] else FAMS, li = seq_along(LAMBDAS), stringsAsFactors = FALSE)
message(sprintf("nNIL-118K sweep: %d cells, %d threads ...", nrow(grid), THREADS))
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

tb1 <- data.table(chr = 1L, start = 272330564L)
tb1_peak <- function(scan) {
  w <- scan[CHR == tb1$chr & abs(BP - tb1$start) <= 5e5 & is.finite(P) & P > 0]
  if (!nrow(w)) NA_real_ else round(max(-log10(w$P)), 2)
}

sweep_list <- vector("list", length(LAMBDAS))
for (li in seq_along(LAMBDAS)) {
  lambda <- LAMBDAS[li]
  idx <- which(grid$li == li)
  G <- do.call(cbind, lapply(idx, function(i) cells[[i]]$block))
  rownames(G) <- union_markers
  scan <- gwas_scan(G)[order(CHR, BP)]
  fwrite(scan, file.path(OUTDIR, sprintf("stam_gwas_nnil_118k_lambda%s.csv", lambda)))
  scan[, coverage := lambda]
  sweep_list[[li]] <- scan
  message(sprintf(
    "  lambda=%-4g : %d markers, tb1 peak -log10P = %s, global max = %.1f",
    lambda, nrow(scan), tb1_peak(scan), max(-log10(scan[is.finite(P) & P > 0, P]))
  ))
}

if (SMOKE) {
  message("smoke ok.")
  quit(save = "no", status = 0)
}

baseline <- fread(file.path(ROOT, "data/teonam/stam_gwas_scan_118k.csv"))[, .(SNP, CHR, BP, P)]
baseline[, coverage := Inf]
message(sprintf("  lambda=Inf  : %d markers, tb1 peak -log10P = %s (authentic 118K baseline)", nrow(baseline), tb1_peak(baseline)))
sweep <- rbindlist(c(sweep_list, list(baseline)), use.names = TRUE)
fwrite(sweep, file.path(OUTDIR, "stam_gwas_nnil_118k_sweep.csv"))
message(sprintf("wrote %s (%d rows, %d coverage levels)", file.path(OUTDIR, "stam_gwas_nnil_118k_sweep.csv"), nrow(sweep), uniqueN(sweep$coverage)))
