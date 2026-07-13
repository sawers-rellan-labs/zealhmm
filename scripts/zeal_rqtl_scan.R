#!/usr/bin/env Rscript
# =============================================================================
# ZEAL trait scan — JOINT multi-trait R/qtl (bcsft) interval mapping on the RTIGER
# ancestry mosaic. Every trait shares ONE cross (identical genotype = the projected
# mosaic; only the phenotype columns differ), so calc.genoprob, scanone, and the
# 1000x permutation null are each computed ONCE over the whole panel — the permutation
# sample is drawn once and reused across traits (the map_single_marker.R design),
# instead of re-running the null per trait. Per-trait peak tables are then sliced out
# and written with the same zeal_<trait>_* filenames the notebooks read.
#
# Marker set/map as before: TeoNAM JLM 0.1-cM grid, markers v5-renamed S<chr>_<pos_v5>;
# pre-scan HMM error correction (calc.genoprob error.prob); per-taxon additive scan on
# the HMM expected dosage. Peaks via detect_peaks.R (get_peak_table / refine_peaks).
#
# Env: TRAITS (comma-sep; default the continuous SpATS-BLUE panel), NPERM (1000),
#      NCORES, ERR_PROB (0.01), INTRO_MIN (0.05), MIN_CLASS (3)
# Out: results/sim/zeal/rqtl/  shared: zeal_rqtl_{scanone_<tag>.csv, perms_<tag>.rds};
#      per-trait (what the notebooks read): zeal_<trait>_{scanone,perms,peaks,peaks_ci,
#      peaks_refined}_<tag>, _cross.rds, _scanone_by_taxon.csv, _peaks_refined_by_taxon.csv
# Run: Rscript scripts/zeal_rqtl_scan.R    (or TRAITS="DTA,PH" Rscript scripts/zeal_rqtl_scan.R)
# =============================================================================
suppressMessages({
  library(qtl)
  library(data.table)
  library(here)
})
source(here("scripts/logging.R"))
source(here("scripts/detect_peaks.R")) # get_peak_table / refine_peaks (airmine multi-peak)
set.seed(1234567890)

TRAITS <- trimws(strsplit(Sys.getenv("TRAITS", "DTA,DTS,PH,EH,EN,LAE,NBR,Prolif,SPAD"), ",")[[1]])
NPERM <- as.integer(Sys.getenv("NPERM", "1000"))
INTRO_MIN <- as.numeric(Sys.getenv("INTRO_MIN", "0.05")) # per-family min fraction introgressed (HET or ALT)
MIN_CLASS <- as.integer(Sys.getenv("MIN_CLASS", "3")) # per-family min samples in any present genotype class
ERR_PROB <- as.numeric(Sys.getenv("ERR_PROB", "0.01")) # HMM genotyping-error rate (pre-scan error correction)
NCORES <- as.integer(Sys.getenv("NCORES", as.character(max(1L, parallel::detectCores() - 2L))))
OUT <- here("results/sim/zeal/rqtl")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
log_info("JOINT scan over %d traits: %s | NPERM=%d NCORES=%d", length(TRAITS), paste(TRAITS, collapse = ","), NPERM, NCORES)

# ---- inputs -----------------------------------------------------------------
mo <- readRDS(here("data/zeal/zeal_rtiger_mosaic.rds"))
mk <- as.data.table(mo$markers)[, .(marker, chr = as.integer(chr), pos = as.integer(pos))]
st <- mo$state
lines <- as.data.table(mo$lines)

# ---- marker set + map = the TeoNAM JLM 0.1cM set (v5-renamed markers) --------
jlm_rs <- fread(here("data/teonam/tassel/geno.hmp.txt"), select = 1L)[[1]]
teonam_map <- fread(here("data/teonam/teonam_v5_native.tsv"))
jlm <- teonam_map[match(jlm_rs, marker), .(marker, chr = chr_v5, pos = pos_v5, cm)][
  !is.na(pos) & !is.na(cm)
][order(chr, pos)]
jlm[, marker_v2 := marker]
jlm[, marker := sprintf("S%d_%d", chr, pos)] # self-documenting v5 id
stopifnot(!anyDuplicated(jlm$marker))
log_info("JLM marker set: %d markers, %d chr, TeoNAM cM range %.1f-%.1f", nrow(jlm), uniqueN(jlm$chr), min(jlm$cm), max(jlm$cm))

# ---- project the ZEAL ancestry mosaic onto the JLM positions (nearest block) ----
nearest_idx <- function(target, sorted) {
  j <- findInterval(target, sorted, all.inside = TRUE)
  ifelse(abs(target - sorted[j]) <= abs(sorted[j + 1L] - target), j, j + 1L)
}
proj <- matrix(NA_integer_, nrow = nrow(jlm), ncol = ncol(st))
for (ch in sort(unique(jlm$chr))) {
  zi <- which(mk$chr == ch)
  zi <- zi[order(mk$pos[zi])]
  zpos <- mk$pos[zi]
  ji <- which(jlm$chr == ch)
  proj[ji, ] <- st[zi[nearest_idx(jlm$pos[ji], zpos)], , drop = FALSE]
}
mk <- jlm
st <- proj
log_info("projected ZEAL mosaic onto %d JLM sites (%d lines); NA=%d", nrow(mk), ncol(st), sum(is.na(st)))

# ---- recode 0/1/2 -> A/H/B, individuals x markers ---------------------------
G <- t(st)
log_info("genotype fractions: A=%.3f H=%.3f B=%.3f", mean(G == 0), mean(G == 1), mean(G == 2))
Gc <- matrix(c("A", "H", "B")[G + 1L], nrow = nrow(G))
colnames(Gc) <- mk$marker

# ---- multi-trait phenotype matrix aligned to the mosaic lines ---------------
# Most traits are SpATS BLUEs (pheno_<t>_blue.csv); the SPAD-date traits are direct
# per-line values and the binary stem/kinki traits empirical logits (as in the GWAS
# drivers). The <TRAIT>_mean column is read case-insensitively.
PHENO_BY <- c(SPAD20DAS = "direct", SPAD36DAS = "direct", STPI = "elogit", STPU = "elogit", KINKI = "elogit")
phtype <- function(tr) {
  p <- unname(PHENO_BY[toupper(tr)])
  if (is.na(p)) "blue" else p
}
Y <- vapply(TRAITS, function(tr) {
  ph <- fread(here(sprintf("data/zeal/pheno_%s_%s.csv", tolower(tr), phtype(tr))))
  mcol <- names(ph)[tolower(names(ph)) == tolower(paste0(tr, "_mean"))][1]
  stopifnot(!is.na(mcol))
  setNames(ph[[mcol]], ph$Genotype)[lines$pedigree]
}, numeric(nrow(lines)))
colnames(Y) <- TRAITS
for (tr in TRAITS) log_info("  %-6s matched %d/%d lines", tr, sum(is.finite(Y[, tr])), nrow(Y))

# ---- taxon covariate (5-family factor), aligned to lines --------------------
ss <- fread(here("data/zeal/samplesheet_3way.csv"))[, .(pedigree, taxon)]
taxv <- ss[match(lines$pedigree, pedigree), taxon]
taxv[taxv == "" | is.na(taxv)] <- NA
tax <- factor(taxv)
lv <- levels(tax)
covar <- vapply(lv[-1], function(l) as.integer(tax == l), integer(length(tax)))
colnames(covar) <- lv[-1]
log_info("taxon: %d lines assigned (%s)", sum(!is.na(tax)), paste(lv, collapse = "/"))

# ---- write R/qtl csv with ALL trait phenotype columns -----------------------
np <- length(TRAITS)
csv <- file.path(OUT, "zeal_rqtl_cross.csv")
top <- rbind(
  c("id", TRAITS, mk$marker),
  c(rep("", np + 1L), as.character(mk$chr)),
  c(rep("", np + 1L), sprintf("%.6f", mk$cm))
)
Yc <- ifelse(is.na(Y), "", formatC(Y, format = "f", digits = 4))
bodym <- cbind(lines$pedigree, Yc, Gc)
fwrite(as.data.table(rbind(top, bodym)), csv, col.names = FALSE, quote = FALSE)
log_info("wrote joint cross csv: %d ind x %d markers x %d traits", nrow(G), nrow(mk), np)

# ---- build cross, one error-corrected genoprob ------------------------------
cross <- read.cross(
  format = "csv", file = csv, genotypes = c("A", "H", "B"), alleles = c("A", "B"),
  estimate.map = FALSE, BC.gen = 2, F.gen = 3
)
cross <- jittermap(cross)
log_info("cross: %s | %d ind, %d markers, %d chr", paste(class(cross), collapse = "/"), nind(cross), totmar(cross), nchr(cross))
cross <- calc.genoprob(cross, step = 1, error.prob = ERR_PROB, map.function = "haldane") # pre-scan HMM error correction

# light cross (no bulky genoprob) for the notebooks — written to each per-trait path.
cross_light <- cross
for (i in seq_along(cross_light$geno)) cross_light$geno[[i]]$prob <- NULL
for (tr in TRAITS) saveRDS(cross_light, file.path(OUT, sprintf("zeal_%s_cross.rds", tolower(tr))))

pcol <- match(TRAITS, names(cross$pheno))
stopifnot(!anyNA(pcol))

# ---- JOINT scanone + ONE shared permutation null per model, sliced per trait --
run_joint <- function(tag, addcov) {
  one <- scanone(cross, pheno.col = pcol, method = "hk", addcovar = addcov)
  fwrite(data.table(marker = rownames(one), one), file.path(OUT, sprintf("zeal_rqtl_scanone_%s.csv", tag)))
  log_info("[%s] joint scanone done (%d traits); starting %d SHARED permutations on %d cores ...", tag, np, NPERM, NCORES)
  t0 <- proc.time()
  perms <- scanone(cross,
    pheno.col = pcol, method = "hk", addcovar = addcov,
    n.perm = NPERM, n.cluster = NCORES, verbose = FALSE
  )
  el <- (proc.time() - t0)[["elapsed"]]
  saveRDS(perms, file.path(OUT, sprintf("zeal_rqtl_perms_%s.rds", tag)))
  log_info("[%s] SHARED PERM null: %.1f s for all %d traits (one sample, reused)", tag, el, np)
  for (tr in TRAITS) {
    ttag <- tolower(tr)
    one_t <- one[, c("chr", "pos", tr)]
    class(one_t) <- c("scanone", "data.frame")
    perms_t <- matrix(perms[, tr], ncol = 1, dimnames = list(NULL, tr))
    class(perms_t) <- "scanoneperm"
    fwrite(
      data.table(marker = rownames(one_t), chr = one_t$chr, pos = one_t$pos, lod = one_t[[3]]),
      file.path(OUT, sprintf("zeal_%s_scanone_%s.csv", ttag, tag))
    )
    saveRDS(perms_t, file.path(OUT, sprintf("zeal_%s_perms_%s.rds", ttag, tag)))
    th <- summary(perms_t, alpha = c(0.05, 0.10))
    peaks <- summary(one_t, perms = perms_t, alpha = 0.05, format = "tabByChr", pvalues = TRUE)
    fwrite(as.data.table(peaks, keep.rownames = "locus"), file.path(OUT, sprintf("zeal_%s_peaks_%s.csv", ttag, tag)))
    ci_tab <- get_peak_table(one_t, perms_t, mk, alpha = 0.05)
    ref_tab <- refine_peaks(ci_tab, one_t, mk)
    fwrite(if (is.null(ci_tab)) data.table() else as.data.table(ci_tab), file.path(OUT, sprintf("zeal_%s_peaks_ci_%s.csv", ttag, tag)))
    fwrite(if (is.null(ref_tab)) data.table() else as.data.table(ref_tab), file.path(OUT, sprintf("zeal_%s_peaks_refined_%s.csv", ttag, tag)))
    log_info(
      "  [%s/%s] 5%%=%.2f maxLOD=%.2f | CI-peaks=%d refined=%d", tag, tr, th[1],
      max(one_t[[3]], na.rm = TRUE), if (is.null(ci_tab)) 0L else nrow(ci_tab), if (is.null(ref_tab)) 0L else nrow(ref_tab)
    )
  }
  invisible(el)
}
log_info("=== joint scan: no-covariate ===")
run_joint("nocovar", NULL)
log_info("=== joint scan: taxon covariate (%d families) ===", ncol(covar))
run_joint("taxon", covar)

# ---- per-taxon ADDITIVE scan (HMM expected dosage), per trait ---------------
# st_hard (hard calls, QC only) and st_dose (expected dosage, the LOD predictor) are
# genotype-only, so computed ONCE and reused across traits.
st_hard <- t(pull.geno(cross)[, mk$marker, drop = FALSE]) - 1L
st_dose <- do.call(rbind, lapply(names(cross$geno), function(ch) {
  pr <- cross$geno[[ch]]$prob
  t(pr[, , 2] * 1 + pr[, , 3] * 2)
}))[mk$marker, , drop = FALSE]
log_info("=== per-taxon ADDITIVE scans (frac_intro>%.2f & class>=%d), %d traits ===", INTRO_MIN, MIN_CLASS, np)

for (tr in TRAITS) {
  ttag <- tolower(tr)
  yv_all <- Y[, tr]
  by_taxon <- rbindlist(lapply(levels(tax), function(t) {
    idx <- which(tax == t & is.finite(yv_all))
    if (length(idx) < 20) {
      return(NULL)
    }
    yv <- yv_all[idx]
    if (sd(yv) == 0) {
      return(NULL)
    }
    Gh <- st_hard[, idx, drop = FALSE]
    nAA <- rowSums(Gh == 0L, na.rm = TRUE)
    nAB <- rowSums(Gh == 1L, na.rm = TRUE)
    nBB <- rowSums(Gh == 2L, na.rm = TRUE)
    ntot <- nAA + nAB + nBB
    frac_intro <- (nAB + nBB) / ntot
    BIG <- .Machine$integer.max
    min_class <- pmin(fifelse(nAA > 0L, nAA, BIG), fifelse(nAB > 0L, nAB, BIG), fifelse(nBB > 0L, nBB, BIG))
    keep <- which(frac_intro > INTRO_MIN & min_class >= MIN_CLASS)
    if (length(keep) < 10) {
      return(NULL)
    }
    n <- length(yv)
    Xk <- t(st_dose[keep, idx, drop = FALSE]) # HMM expected dosage
    lod <- as.vector(-(n / 2) * log10(1 - cor(yv, Xk)^2))
    maxperm <- replicate(NPERM, max(-(n / 2) * log10(1 - cor(sample(yv), Xk)^2), na.rm = TRUE))
    thr5 <- as.numeric(quantile(maxperm, 0.95, na.rm = TRUE))
    data.table(taxon = t, marker = mk$marker[keep], chr = mk$chr[keep], pos = mk$cm[keep], lod = lod, n = n, thr5 = thr5)
  }))
  fwrite(by_taxon, file.path(OUT, sprintf("zeal_%s_scanone_by_taxon.csv", ttag)))
  by_taxon_ref <- data.table()
  if (nrow(by_taxon) > 0) {
    by_taxon_ref <- rbindlist(lapply(split(by_taxon, by_taxon$taxon), function(d) {
      thr <- d$thr5[1]
      if (max(d$lod, na.rm = TRUE) <= thr) {
        return(NULL)
      }
      rbindlist(lapply(sort(unique(d$chr)), function(ch) {
        dc <- d[chr == ch]
        if (nrow(dc) < 3 || max(dc$lod, na.rm = TRUE) <= thr) {
          return(NULL)
        }
        sp <- find_subpeaks(dc$pos, dc$lod, thresh = thr)
        if (is.null(sp)) {
          return(NULL)
        }
        sp <- sp[sp$lod > thr, , drop = FALSE]
        if (nrow(sp) == 0) {
          return(NULL)
        }
        pkm <- Map(.nearest_marker, list(mk), ch, sp$pos_cm)
        cibp <- t(mapply(function(lo, hi) .ci_bp(mk, ch, lo, hi), sp$ci_low_cm, sp$ci_high_cm))
        data.table(
          taxon = d$taxon[1], chr = ch, pos = sp$pos_cm, lod = sp$lod, thr5 = thr,
          ci.low = sp$ci_low_cm, ci.high = sp$ci_high_cm,
          marker = vapply(pkm, `[[`, "", "marker"),
          ci_left = cibp[, 1], ci_right = cibp[, 2], width_mb = (cibp[, 2] - cibp[, 1]) / 1e6
        )
      }))
    }), fill = TRUE)
  }
  fwrite(by_taxon_ref, file.path(OUT, sprintf("zeal_%s_peaks_refined_by_taxon.csv", ttag)))
  log_info("  per-taxon [%s]: %d taxa, %d rows -> %d refined QTL", tr, uniqueN(by_taxon$taxon), nrow(by_taxon), nrow(by_taxon_ref))
}

log_info("=== DONE (joint %d-trait scan) ===", np)
