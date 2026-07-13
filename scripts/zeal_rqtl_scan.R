#!/usr/bin/env Rscript
# =============================================================================
# ZEAL trait scan — classic R/qtl (bcsft) interval mapping on the RTIGER ancestry
# mosaic. Trait selected by env TRAIT (DTA | DTS | PH | ...); reads the matching
# data/zeal/pheno_<trait>_blue.csv (column <TRAIT>_mean) and writes zeal_<trait>_*.
#
# Design: ZEAL/BZea is a BC2S3 (B73 x 5 teosinte taxa). Genotype = RTIGER mosaic
# states 0/1/2 recoded B73=A / het=H / teosinte=B.
#
# Genetic map: we do NOT estimate a map from ZEAL genotypes (excess heterozygotes
# ~9% H vs ~2% hom-teo make a ZEAL est.map unreliable). ALL distances here are the
# **TeoNAM** genetic distances (cM), NOT native to ZEAL. Specifically, the marker
# set + map IS the TeoNAM JLM set (9,063 markers, 0.1 cM, TeoNAM cM from
# teonam_jlm_build.R) — the SAME grid + TeoNAM map as the TeoNAM JLM reproduction,
# so the two analyses are directly comparable. ZEAL SNP50K and TeoNAM markers are
# disjoint platforms, so the ZEAL RTIGER-mosaic ancestry is PROJECTED onto the JLM
# positions (nearest-block); no cM is computed from ZEAL.
#
# scanone (Haley-Knott) + N permutations (single subsampling reused across the
# scan; TIMED). Threshold from permutations, not a naive Bonferroni/BH bar.
#
# Peak detection: beyond R/qtl's one-peak-per-chromosome summary(), we run the
# airmine multi-peak search (scripts/detect_peaks.R) — lodint 1.5-LOD confidence
# intervals (get_peak_table) plus Akima-interpolation + pracma::findpeaks +
# interval-graph deconvolution (refine_peaks) to recover multiple linked QTL on a
# single chromosome. Applied to the whole-genome scans AND the per-taxon scans.
#
# Env: TRAIT (default DTA), NPERM (default 1000), NCORES (default detectCores()-2)
# Out: results/sim/zeal/rqtl/  (zeal_<trait>_ cross rds, scanone csv, perms rds,
#      peaks csv, peaks_ci csv [lodint CIs], peaks_refined csv [multi-peak], timing)
# Run: TRAIT=DTS Rscript scripts/zeal_rqtl_scan.R
# =============================================================================

suppressMessages({
  library(qtl)
  library(data.table)
  library(here)
})
source(here("scripts/logging.R"))
source(here("scripts/detect_peaks.R")) # get_peak_table / refine_peaks (airmine multi-peak)
set.seed(1234567890)

TRAIT <- toupper(Sys.getenv("TRAIT", "DTA")) # DTA | DTS | PH | ...
trait_lc <- tolower(TRAIT)
PFX <- sprintf("zeal_%s", trait_lc) # output filename prefix
NPERM <- as.integer(Sys.getenv("NPERM", "1000"))
INTRO_MIN <- as.numeric(Sys.getenv("INTRO_MIN", "0.05")) # per-family min fraction introgressed (HET or ALT)
MIN_CLASS <- as.integer(Sys.getenv("MIN_CLASS", "3")) # per-family min samples in any present genotype class
ERR_PROB <- as.numeric(Sys.getenv("ERR_PROB", "0.01")) # HMM genotyping-error rate (pre-scan error correction)
NCORES <- as.integer(Sys.getenv("NCORES", as.character(max(1L, parallel::detectCores() - 2L))))
OUT <- here("results/sim/zeal/rqtl")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- inputs ---------------------------------------------------------------
mo <- readRDS(here("data/zeal/zeal_rtiger_mosaic.rds"))
mk <- as.data.table(mo$markers)[, .(marker, chr = as.integer(chr), pos = as.integer(pos))]
st <- mo$state
lines <- as.data.table(mo$lines)
ph <- fread(here(sprintf("data/zeal/pheno_%s_blue.csv", trait_lc)))
ph <- ph[, .(pedigree = Genotype, y = get(paste0(TRAIT, "_mean")))]

# ---- marker set + map = the TeoNAM JLM 0.1cM set (9,063 markers), TeoNAM cM ----
# ALL of this ZEAL R/qtl work uses the SAME marker grid + genetic map as the TeoNAM
# JLM reproduction (data/teonam/tassel/geno.hmp.txt built by teonam_jlm_build.R),
# NOT a fresh thinning of ZEAL's own SNP50K. The cM are TeoNAM's (not native to
# ZEAL). ZEAL SNP50K and TeoNAM markers are disjoint platforms, so the JLM positions
# + TeoNAM cM define the grid and the ZEAL genotypes are projected onto them (below).
jlm_rs <- fread(here("data/teonam/tassel/geno.hmp.txt"), select = 1L)[[1]]
teonam_map <- fread(here("data/teonam/teonam_v5_native.tsv")) # TeoNAM est.map (cm = TeoNAM cM)
jlm <- teonam_map[match(jlm_rs, marker), .(marker, chr = chr_v5, pos = pos_v5, cm)][
  !is.na(pos) & !is.na(cm)
][order(chr, pos)]
# Rename markers to their v5 coordinates: S<chr_v5>_<pos_v5>. The original TeoNAM
# rs# (e.g. S1_128373) embeds the *v2* bp — confusing here, where every position is
# v5. So the ZEAL marker id IS its v5 position: self-documenting in the scanone /
# peaks / refine_peaks output, no lookup needed to read a peak's Mb. marker_v2 keeps
# the TeoNAM id for provenance / cross-referencing the TeoNAM JLM run.
jlm[, marker_v2 := marker]
jlm[, marker := sprintf("S%d_%d", chr, pos)]
stopifnot(!anyDuplicated(jlm$marker))
log_info(
  "JLM marker set: %d markers, %d chr, TeoNAM cM range %.1f-%.1f",
  nrow(jlm), uniqueN(jlm$chr), min(jlm$cm), max(jlm$cm)
)

# ---- project the ZEAL ancestry mosaic onto the JLM positions ----
# The mosaic is piecewise-constant ancestry (RTIGER blocks); each JLM site takes the
# state of the nearest ZEAL marker on the same chromosome (= its covering block).
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
mk <- jlm # analysis markers = JLM set (marker, chr, pos, cm)
st <- proj # projected ancestry states, 9063 x lines
log_info("projected ZEAL mosaic onto %d JLM sites (%d lines); NA=%d", nrow(mk), ncol(st), sum(is.na(st)))

# ---- recode genotypes 0/1/2 -> A/H/B, individuals x markers ----------------
G <- t(st) # lines x markers
fA <- mean(G == 0)
fH <- mean(G == 1)
fB <- mean(G == 2)
log_info("genotype fractions: A/B73=%.3f  H/het=%.3f  B/teo=%.3f  (H:B ratio=%.1f)", fA, fH, fB, fH / fB)
Gc <- matrix(c("A", "H", "B")[G + 1L], nrow = nrow(G))
colnames(Gc) <- mk$marker

# ---- phenotype aligned to mosaic lines -------------------------------------
# The <TRAIT>_mean BLUE is already per-taxa fence-cleaned at the raw-plot level
# upstream (zeal_spats_blues.R: per field x taxa 3xIQR fence before SpATS), so no
# further phenotype cleaning here.
y <- ph[match(lines$pedigree, pedigree), y]
log_info("%s matched for %d of %d lines (BLUEs fence-cleaned upstream)", TRAIT, sum(!is.na(y)), length(y))

# ---- taxon covariate (5-family factor: Zd/Zh/Zl/Zv/Zx), aligned to lines ----
ss <- fread(here("data/zeal/samplesheet_3way.csv"))[, .(pedigree, taxon)]
taxv <- ss[match(lines$pedigree, pedigree), taxon]
taxv[taxv == "" | is.na(taxv)] <- NA
tax <- factor(taxv)
log_info(
  "taxon: %d lines assigned (%s); %d without taxon (checks/B73)",
  sum(!is.na(tax)), paste(levels(tax), collapse = "/"), sum(is.na(tax))
)
# indicator design matrix (k-1 cols); NA rows propagate -> scanone drops them
lv <- levels(tax)
covar <- vapply(lv[-1], function(l) as.integer(tax == l), integer(length(tax)))
colnames(covar) <- lv[-1]

# ---- write R/qtl csv (row1 header, row2 chr, row3 cM, then individuals) -----
csv <- file.path(OUT, sprintf("%s_cross.csv", PFX))
top <- rbind(
  c("id", TRAIT, mk$marker),
  c("", "", as.character(mk$chr)),
  c("", "", sprintf("%.6f", mk$cm))
)
bodym <- cbind(lines$pedigree, ifelse(is.na(y), "", sprintf("%.4f", y)), Gc)
fwrite(as.data.table(rbind(top, bodym)), csv, col.names = FALSE, quote = FALSE)
log_info("wrote cross csv: %s (%d ind x %d markers)", csv, nrow(G), nrow(mk))

# ---- build bcsft cross, interpolated map (no est.map) ----------------------
cross <- read.cross(
  format = "csv", file = csv, genotypes = c("A", "H", "B"), alleles = c("A", "B"),
  estimate.map = FALSE, BC.gen = 2, F.gen = 3
)
cross <- jittermap(cross)
log_info(
  "cross: %s | %d ind, %d markers, %d chr", paste(class(cross), collapse = "/"),
  nind(cross), totmar(cross), nchr(cross)
)
# ---- pre-scan genotype error correction via the HMM (calc.genoprob) ---------
# Genotypes are RTIGER ancestry projected onto the JLM grid; the hard calls carry
# het/hom-teosinte miscalls that punch artificial LOD dives -- splitting peaks and
# clipping their CIs. We correct this at the GENOTYPE level, before the scan, instead
# of smoothing peaks afterwards. calc.genoprob runs the R/qtl HMM with a realistic
# ERR_PROB (not the old 1e-4, which trusted every call and did no correction) and
# returns genotype PROBABILITIES that down-weight calls improbable given the flanking
# markers + recombination. The whole-genome Haley-Knott scan regresses on these probs
# directly; the per-taxon additive scan uses the HMM EXPECTED dosage derived from them
# (st_dose, below). (Hard deletion -- cleanGeno / calc.errorlod -- was tried: cleanGeno
# supports only 2-genotype crosses, and errorlod deletion barely moved the dive because
# it is a gradual het/hom cline, not isolated flips. The soft probabilities fix it.)
cross <- calc.genoprob(cross, step = 1, error.prob = ERR_PROB, map.function = "haldane")
saveRDS(cross, file.path(OUT, sprintf("%s_cross.rds", PFX)))

# ---- scanone + TIMED permutations: no-covariate vs taxon-covariate ---------
run_scan <- function(tag, addcov) {
  one <- scanone(cross, pheno.col = TRAIT, method = "hk", addcovar = addcov)
  fwrite(data.table(marker = rownames(one), one), file.path(OUT, sprintf("%s_scanone_%s.csv", PFX, tag)))
  log_info("[%s] starting %d permutations on %d cores ...", tag, NPERM, NCORES)
  t0 <- proc.time()
  perms <- scanone(cross,
    pheno.col = TRAIT, method = "hk", addcovar = addcov,
    n.perm = NPERM, n.cluster = NCORES, verbose = FALSE
  )
  el <- (proc.time() - t0)[["elapsed"]]
  saveRDS(perms, file.path(OUT, sprintf("%s_perms_%s.rds", PFX, tag)))
  th <- summary(perms, alpha = c(0.05, 0.10))
  peaks <- summary(one, perms = perms, alpha = 0.05, format = "tabByChr", pvalues = TRUE)
  fwrite(as.data.table(peaks, keep.rownames = "locus"), file.path(OUT, sprintf("%s_peaks_%s.csv", PFX, tag)))
  # airmine peak detection: lodint CIs (one-per-chr) + multi-peak refinement
  ci_tab <- get_peak_table(one, perms, mk, alpha = 0.05)
  ref_tab <- refine_peaks(ci_tab, one, mk)
  fwrite(
    if (is.null(ci_tab)) data.table() else as.data.table(ci_tab),
    file.path(OUT, sprintf("%s_peaks_ci_%s.csv", PFX, tag))
  )
  fwrite(
    if (is.null(ref_tab)) data.table() else as.data.table(ref_tab),
    file.path(OUT, sprintf("%s_peaks_refined_%s.csv", PFX, tag))
  )
  log_info(
    "[%s] PERM TIMING: %.1f s (%.4f s/perm) | 5%%=%.2f 10%%=%.2f | maxLOD=%.2f | nQTL=%d | CI-peaks=%d refined=%d",
    tag, el, el / NPERM, th[1], th[2], max(one$lod, na.rm = TRUE),
    if (is.data.frame(peaks)) nrow(peaks) else 0L,
    if (is.null(ci_tab)) 0L else nrow(ci_tab),
    if (is.null(ref_tab)) 0L else nrow(ref_tab)
  )
  list(tag = tag, elapsed = el, thr = th, peaks = peaks, refined = ref_tab, maxlod = max(one$lod, na.rm = TRUE))
}

log_info("=== %s scan: no-covariate ===", TRAIT)
r0 <- run_scan("nocovar", NULL)
log_info("=== %s scan: taxon covariate (%d families) ===", TRAIT, ncol(covar))
r1 <- run_scan("taxon", covar)

# Per-taxon predictors, both markers x lines in mk order:
#  - st_hard: hard calls 0/1/2 (A/H/B), used ONLY for the QC class-count filters.
#  - st_dose: HMM EXPECTED teosinte dosage E[g] = P(het) + 2 P(hom-teo) from the
#    error-aware genoprobs, used for the additive LOD -- this is the pre-scan error
#    correction (soft probabilities absorb the het/hom miscalls that split peaks).
st_hard <- t(pull.geno(cross)[, mk$marker, drop = FALSE]) - 1L
st_dose <- do.call(rbind, lapply(names(cross$geno), function(ch) {
  pr <- cross$geno[[ch]]$prob # ind x positions x 3 (A/H/B)
  t(pr[, , 2] * 1 + pr[, , 3] * 2) # positions x ind: expected teosinte dosage
}))[mk$marker, , drop = FALSE] # real markers only (drops step=1 pseudomarkers), mk order

# ---- per-taxon ADDITIVE-ONLY scan (family-resolved; cf. Andosol bcfst_by_donor.R) --
# Each taxon's NILs are scanned separately to show WHICH families carry each QTL.
# Model = ADDITIVE, regressed on the HMM EXPECTED teosinte dosage E[g] = P(het) +
# 2*P(hom-teo) (continuous 0-2, from calc.genoprob) -- NOT the hard 0/1/2 calls. The
# soft dosage carries the same pre-scan error correction as the whole-genome HK scan;
# calc.genoprob CORRECTS (soft-weights) the projected calls, it does not re-genotype.
# LOD is the 1-df regression of the trait on that dosage. The HARD calls (st_hard) are
# used ONLY for the two per-family QC filters, both required:
#   (a) fraction introgressed (HET or ALT) > INTRO_MIN -- HET or ALT because the
#       RTIGER caller calls het over hom-ALT, so ALT-only undercounts introgression;
#   (b) >= MIN_CLASS samples in EVERY present hard genotype class -- so no single 1-2
#       sample class (e.g. a lone B/B with an extreme value) can drive a peak.
log_info("=== per-taxon ADDITIVE-only scans (frac_intro>%.2f & class>=%d) ===", INTRO_MIN, MIN_CLASS)
by_taxon <- rbindlist(lapply(levels(tax), function(t) {
  idx <- which(tax == t & !is.na(y))
  if (length(idx) < 20) {
    log_info("  %s: skipped (n=%d)", t, length(idx))
    return(NULL)
  }
  yv <- y[idx]
  if (sd(yv) == 0) {
    log_info("  %s: skipped (no %s variation)", t, TRAIT)
    return(NULL)
  }
  Gh <- st_hard[, idx, drop = FALSE] # hard calls, for the QC class-count filters only
  nAA <- rowSums(Gh == 0L, na.rm = TRUE)
  nAB <- rowSums(Gh == 1L, na.rm = TRUE)
  nBB <- rowSums(Gh == 2L, na.rm = TRUE)
  ntot <- nAA + nAB + nBB
  frac_intro <- (nAB + nBB) / ntot
  BIG <- .Machine$integer.max
  min_class <- pmin(
    fifelse(nAA > 0L, nAA, BIG), fifelse(nAB > 0L, nAB, BIG), fifelse(nBB > 0L, nBB, BIG)
  )
  keep <- which(frac_intro > INTRO_MIN & min_class >= MIN_CLASS)
  if (length(keep) < 10) {
    log_info("  %s: skipped (only %d markers pass filters)", t, length(keep))
    return(NULL)
  }
  n <- length(yv)
  Xk <- t(st_dose[keep, idx, drop = FALSE]) # lines x kept markers: HMM expected dosage
  # additive 1-df LOD via the correlation identity LOD = -(n/2) log10(1 - r^2)
  lod <- as.vector(-(n / 2) * log10(1 - cor(yv, Xk)^2))
  # per-taxon permutation 5% threshold: max additive LOD over markers, NPERM shuffles of DTA
  maxperm <- replicate(NPERM, max(-(n / 2) * log10(1 - cor(sample(yv), Xk)^2), na.rm = TRUE))
  thr5 <- as.numeric(quantile(maxperm, 0.95, na.rm = TRUE))
  log_info(
    "  %s: n=%d | kept %d markers | maxLOD=%.2f | perm 5%% thr=%.2f (%s)",
    t, n, length(keep), max(lod, na.rm = TRUE), thr5, ifelse(max(lod, na.rm = TRUE) > thr5, "SIG", "ns")
  )
  data.table(taxon = t, marker = mk$marker[keep], chr = mk$chr[keep], pos = mk$cm[keep], lod = lod, n = n, thr5 = thr5)
}))
fwrite(by_taxon, file.path(OUT, sprintf("%s_scanone_by_taxon.csv", PFX)))
log_info("per-taxon (additive): %d taxa scanned, %d rows", uniqueN(by_taxon$taxon), nrow(by_taxon))

# ---- per-taxon multi-peak refinement (airmine find_subpeaks on each family) ----
# Each taxon x chr LOD profile is deconvolved into non-overlapping QTL (Akima +
# findpeaks + 1.5-LOD-drop interval graph), thresholded at that taxon's perm 5%.
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
        ci_left = cibp[, 1], ci_right = cibp[, 2],
        width_mb = (cibp[, 2] - cibp[, 1]) / 1e6
      )
    }))
  }), fill = TRUE)
  fwrite(by_taxon_ref, file.path(OUT, sprintf("%s_peaks_refined_by_taxon.csv", PFX)))
  log_info(
    "per-taxon refined: %d QTL across %d taxa",
    nrow(by_taxon_ref), if (nrow(by_taxon_ref) > 0) uniqueN(by_taxon_ref$taxon) else 0L
  )
}

log_info("=== DONE ===")
log_info("no-covariate : 5%%=%.2f  nQTL=%d  (%.1f s)", r0$thr[1], nrow(r0$peaks), r0$elapsed)
log_info("taxon-covar  : 5%%=%.2f  nQTL=%d  (%.1f s)", r1$thr[1], nrow(r1$peaks), r1$elapsed)
cat("\n--- taxon-covariate significant peaks (alpha=0.05) ---\n")
print(r1$peaks)
