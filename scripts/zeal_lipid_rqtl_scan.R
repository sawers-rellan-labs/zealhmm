#!/usr/bin/env Rscript
# =============================================================================
# ZEAL lipid QTL scan — JOINT multi-trait R/qtl (bcsft) HK interval mapping on the
# rtiger ancestry mosaic. Reads the wide lipid phenotype matrix
# (data/zeal/pheno_lipid_all.csv, from scripts/zeal_lipid_pheno.R) and emits ONE
# significant-QTL summary across all 133 lipid traits.
#
# Map/genotype: TeoNAM JLM 0.1-cM grid (markers S<chr>_<pos_v5>), ZEAL mosaic projected
# on (nearest block), pre-scan HMM error correction; taxon covariate; per-trait 1000-perm
# threshold. Peaks/CIs via scripts/detect_peaks.R.
#
# PARALLELISM (MODE env): the 1000-perm null is the bottleneck; R/qtl's forked n.cluster
# segfaults here (macOS Accelerate BLAS is not fork-safe). The per-permutation cost is
# dominated by the marker loop, NOT the trait count (measured: 14 traits cost the same
# ~5 s/perm as 133), so we split the PERMUTATIONS across SEPARATE PROCESSES (each scans
# ALL traits for a slice of the perms, with a DISTINCT seed), then pool the null.
#   MODE=prep   build the cross ONCE (projection + read.cross + calc.genoprob), compute
#               the observed scan (all traits), save chunks/cross_gp.rds + meta.rds +
#               zeal_lipid_scanone_taxon.csv + the cross csv (notebook reads it). Run first.
#   MODE=chunk  CHUNK="i/N": load cross_gp.rds, seed=base+i, run NPERM/N perms over ALL
#               traits, save chunks/perms_i.rds.
#   MODE=merge  rbind the chunk perm matrices -> full NPERM-row null; per-trait summary.
#   MODE=full   (default) single-process path (fallback).
# Orchestrate: bash scripts/zeal_lipid_scan_parallel.sh
# =============================================================================
suppressMessages({
  library(qtl)
  library(data.table)
  library(here)
})
source(here("scripts/logging.R"))
source(here("scripts/detect_peaks.R"))
SEED <- 1234567890L

NPERM <- as.integer(Sys.getenv("NPERM", "1000"))
ERR_PROB <- as.numeric(Sys.getenv("ERR_PROB", "0.01"))
ALPHA <- as.numeric(Sys.getenv("ALPHA", "0.05"))
MODEL <- Sys.getenv("MODEL", "taxon")
MODE <- Sys.getenv("MODE", "full")
CHUNK <- Sys.getenv("CHUNK", "")
TAG <- Sys.getenv("TAG", "taxon") # output namespace; "na4" for the 4-heavy-trait NA rerun
OUT <- here("results/sim/zeal/rqtl/lipid")
CH <- file.path(OUT, if (TAG == "taxon") "chunks" else paste0("chunks_", TAG))
dir.create(CH, recursive = TRUE, showWarnings = FALSE)
csv <- file.path(OUT, if (TAG == "taxon") "zeal_lipid_cross.csv" else sprintf("zeal_lipid_cross_%s.csv", TAG))
sc_out <- file.path(OUT, sprintf("zeal_lipid_scanone_%s.csv", TAG))
pm_out <- file.path(OUT, sprintf("zeal_lipid_perms_%s.rds", TAG))
summ_out <- file.path(OUT, if (TAG == "taxon") "zeal_lipid_qtl_summary.csv" else sprintf("zeal_lipid_qtl_summary_%s.csv", TAG))

ph <- fread(Sys.getenv("PHENO_FILE", here("data/zeal/pheno_lipid_all.csv")))
meta_cols <- c("Genotype", "nil_id", "taxon", "Material")
TRAITS <- setdiff(names(ph), meta_cols)
TO <- Sys.getenv("TRAITS_ONLY", "") # comma-sep subset (default: all)
if (nzchar(TO)) TRAITS <- intersect(TRAITS, trimws(strsplit(TO, ",", fixed = TRUE)[[1]]))
trait_class <- fifelse(
  grepl("^[A-Za-z]+_[0-9]+_[0-9]+$", TRAITS), "species",
  fifelse(grepl("^[A-Za-z]+$", TRAITS), "sum", "ratio")
)
names(trait_class) <- TRAITS

build_summary <- function(one_df, perms_mat, mk, per_trait_n, trait_class, traits) {
  rbindlist(lapply(traits, function(tr) {
    one_t <- data.frame(chr = one_df$chr, pos = one_df$pos, LOD = one_df[[tr]], row.names = one_df$marker)
    class(one_t) <- c("scanone", "data.frame")
    perms_t <- matrix(perms_mat[, tr], ncol = 1, dimnames = list(NULL, tr))
    class(perms_t) <- "scanoneperm"
    thr <- as.numeric(summary(perms_t, alpha = ALPHA))[1]
    pk <- tryCatch(get_peak_table(one_t, perms_t, mk, alpha = ALPHA), error = function(e) NULL)
    if (is.null(pk) || nrow(pk) == 0) {
      return(NULL)
    }
    ref <- tryCatch(refine_peaks(pk, one_t, mk), error = function(e) NULL)
    tab <- if (!is.null(ref) && nrow(ref) > 0) as.data.table(ref) else as.data.table(pk)
    tab[, `:=`(
      trait = tr, class = unname(trait_class[tr]), thr5 = thr,
      peak_mb = mk$pos[match(marker, mk$marker)] / 1e6,
      ci_left_mb = ci_left / 1e6, ci_right_mb = ci_right / 1e6, n = per_trait_n[[tr]]
    )]
    tab
  }), fill = TRUE)
}
write_summary <- function(summ, ntr) {
  if (nrow(summ) == 0) {
    log_info("no significant QTL at alpha=%.2f", ALPHA)
    return(invisible())
  }
  setnames(summ, "pos", "peak_cm", skip_absent = TRUE)
  cols <- c(
    "trait", "class", "chr", "peak_cm", "peak_mb", "lod",
    if ("thresh" %in% names(summ)) "thresh" else "thr5",
    "ci.low", "ci.high", "ci_left_mb", "ci_right_mb", "width_mb", "marker", "n"
  )
  summ <- summ[order(-lod), intersect(cols, names(summ)), with = FALSE]
  fwrite(summ, summ_out)
  log_info("SIGNIFICANT QTL (alpha=%.2f): %d peaks across %d/%d traits", ALPHA, nrow(summ), uniqueN(summ$trait), ntr)
  log_info("by class: %s", paste(names(table(summ$class)), table(summ$class), sep = "=", collapse = " "))
  cat("\n===== top 25 lipid QTL (by LOD) =====\n")
  print(head(summ, 25))
}

# =============================================================================
# MERGE: observed scan (from prep) + pooled perms (rbind chunks) -> summary
# =============================================================================
if (MODE == "merge") {
  meta <- readRDS(file.path(CH, "meta.rds"))
  one_full <- fread(sc_out)
  pm_files <- sort(Sys.glob(file.path(CH, "perms_*.rds")))
  perms_full <- do.call(rbind, lapply(pm_files, function(f) {
    m <- readRDS(f)
    class(m) <- NULL
    as.matrix(m)
  }))
  class(perms_full) <- "scanoneperm"
  saveRDS(perms_full, pm_out)
  TR <- setdiff(names(one_full), c("marker", "chr", "pos"))
  log_info("merge: %d perm chunks -> pooled null %d x %d; %d traits", length(pm_files), nrow(perms_full), ncol(perms_full), length(TR))
  write_summary(build_summary(one_full, perms_full, meta$mk, meta$per_trait_n, meta$trait_class, TR), length(TR))
  log_info("=== MERGE DONE ===")
  quit(save = "no")
}

# =============================================================================
# CHUNK: load pre-built genoprob cross, run this chunk's slice of the perms
# =============================================================================
if (MODE == "chunk") {
  parts <- as.integer(strsplit(CHUNK, "/", fixed = TRUE)[[1]])
  i <- parts[1]
  N <- parts[2]
  set.seed(SEED + i) # DISTINCT seed per chunk -> independent permutations to pool
  cross <- readRDS(file.path(CH, "cross_gp.rds"))
  meta <- readRDS(file.path(CH, "meta.rds"))
  addcov <- if (MODEL == "taxon") meta$covar else NULL
  nper <- NPERM %/% N + as.integer(i <= (NPERM %% N)) # split the perms, distribute remainder
  pcol <- match(TRAITS, names(cross$pheno))
  stopifnot(!anyNA(pcol))
  BATCH <- as.integer(Sys.getenv("PERM_BATCH", "25")) # perms per progress step
  out <- file.path(CH, sprintf("perms_%02d.rds", i))
  done <- 0L
  acc <- NULL
  t0 <- proc.time()
  while (done < nper) {
    b <- min(BATCH, nper - done)
    pm <- scanone(cross, pheno.col = pcol, method = "hk", addcovar = addcov, n.perm = b, n.cluster = 1L, verbose = FALSE)
    acc <- if (is.null(acc)) unclass(pm) else rbind(acc, unclass(pm))
    done <- done + b
    x <- acc
    class(x) <- "scanoneperm"
    saveRDS(x, out) # incremental checkpoint after every batch
    el <- (proc.time() - t0)[["elapsed"]]
    log_info("chunk %d/%d: %d/%d perms | %.0fs | %.1fs/perm | ETA %.0fs", i, N, done, nper, el, el / done, (el / done) * (nper - done))
  }
  quit(save = "no")
}

# =============================================================================
# prep / full: projection + read.cross + calc.genoprob (the expensive setup, ONCE)
# =============================================================================
set.seed(SEED)
mo <- readRDS(here("data/zeal/zeal_rtiger_mosaic.rds"))
mk <- as.data.table(mo$markers)[, .(marker, chr = as.integer(chr), pos = as.integer(pos))]
st <- mo$state
lines <- as.data.table(mo$lines)
jlm_rs <- fread(here("data/teonam/tassel/geno.hmp.txt"), select = 1L)[[1]]
teonam_map <- fread(here("data/teonam/teonam_v5_native.tsv"))
jlm <- teonam_map[match(jlm_rs, marker), .(marker, chr = chr_v5, pos = pos_v5, cm)][!is.na(pos) & !is.na(cm)][order(chr, pos)]
jlm[, marker := sprintf("S%d_%d", chr, pos)]
stopifnot(!anyDuplicated(jlm$marker))
nearest_idx <- function(target, sorted) {
  j <- findInterval(target, sorted, all.inside = TRUE)
  ifelse(abs(target - sorted[j]) <= abs(sorted[j + 1L] - target), j, j + 1L)
}
proj <- matrix(NA_integer_, nrow = nrow(jlm), ncol = ncol(st))
for (ch in sort(unique(jlm$chr))) {
  zi <- which(mk$chr == ch)
  zi <- zi[order(mk$pos[zi])]
  ji <- which(jlm$chr == ch)
  proj[ji, ] <- st[zi[nearest_idx(jlm$pos[ji], mk$pos[zi])], , drop = FALSE]
}
mk <- jlm
st <- proj
G <- t(st)
Gc <- matrix(c("A", "H", "B")[G + 1L], nrow = nrow(G))
colnames(Gc) <- mk$marker

yidx <- match(lines$pedigree, ph$Genotype)
Y <- as.matrix(ph[yidx, ..TRAITS])
per_trait_n <- colSums(is.finite(Y))
ss <- fread(here("data/zeal/samplesheet_3way.csv"))[, .(pedigree, taxon)]
taxv <- ss[match(lines$pedigree, pedigree), taxon]
taxv[taxv == "" | is.na(taxv)] <- NA
tax <- factor(taxv)
lv <- levels(tax)
covar <- vapply(lv[-1], function(l) as.integer(tax == l), integer(length(tax)))
colnames(covar) <- lv[-1]
addcov <- if (MODEL == "taxon") covar else NULL

np <- length(TRAITS)
top <- rbind(
  c("id", TRAITS, mk$marker),
  c(rep("", np + 1L), as.character(mk$chr)),
  c(rep("", np + 1L), sprintf("%.6f", mk$cm))
)
Yc <- ifelse(is.na(Y), "", formatC(Y, format = "f", digits = 6))
fwrite(as.data.table(rbind(top, cbind(lines$pedigree, Yc, Gc))), csv, col.names = FALSE, quote = FALSE)
cross <- read.cross(format = "csv", file = csv, genotypes = c("A", "H", "B"), alleles = c("A", "B"), estimate.map = FALSE, BC.gen = 2, F.gen = 3)
cross <- jittermap(cross)
cross <- calc.genoprob(cross, step = 1, error.prob = ERR_PROB, map.function = "haldane")
pcol <- match(TRAITS, names(cross$pheno))
one <- scanone(cross, pheno.col = pcol, method = "hk", addcovar = addcov) # observed, once
fwrite(data.table(marker = rownames(one), one), sc_out)

if (MODE == "prep") {
  saveRDS(cross, file.path(CH, "cross_gp.rds"))
  saveRDS(list(mk = mk, per_trait_n = per_trait_n, trait_class = trait_class, covar = covar), file.path(CH, "meta.rds"))
  log_info("prep: cross+genoprob+observed scan saved; %d ind x %d markers x %d traits", nrow(G), nrow(mk), np)
  quit(save = "no")
}

# ---- MODE=full: single-process perms + summary ------------------------------
t0 <- proc.time()
perms <- scanone(cross, pheno.col = pcol, method = "hk", addcovar = addcov, n.perm = NPERM, n.cluster = 1L, verbose = FALSE)
saveRDS(perms, pm_out)
log_info("full: %d perms, %.1f s", NPERM, (proc.time() - t0)[["elapsed"]])
write_summary(build_summary(data.table(marker = rownames(one), one), perms, mk, per_trait_n, trait_class, TRAITS), np)
log_info("=== DONE ===")
