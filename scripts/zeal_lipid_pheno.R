#!/usr/bin/env Rscript
# =============================================================================
# ZEAL lipid phenotypes — ID crosswalk + scale choice for QTL mapping.
#
# SCALE (verified 2026-09-03 against the ORIGINAL High-PC1 QTL data-prep script,
# PTxB73_lipid_QTL_GxE, B73xPT_base_phenotype.csv builder): the original QTL
# phenotypes are RAW — raw species, raw rowSums class sums, and raw ratios
# (num/den per plant, then averaged). The only log10 in that script is on a
# throwaway matrix `fullm` used SOLELY to select interesting ratios via a Gaussian
# graphical model (GGMselect); it is never a mapped phenotype. Fig_3.Rmd likewise
# z-scores raw concentrations and log10's ONLY the PCs/LPCs ratio (and only for the
# Fig 3F overlay plot); it explicitly comments out a log10(signal) option. So the
# scale-faithful reproduction is TRANSFORM=none (raw); the TRANSFORM=log path below
# is our own count+log variant, NOT the High-PC1 method — keep for comparison only.
#
# Input : data/zeal/BZeaLipid_AdjustedMeans_wSums.csv — per-Material NORMALIZED MS
#         SIGNAL (peak-area of the ion current, HPLC-MS), averaged over ~3 reps. These
#         are NOT BLUEs; there is no spatial / mixed-model correction. 133 traits =
#         115 species + 10 class sums + 8 ratios. A ~4% fraction of species values is
#         <=0 (below-detection: normalization pushed the signal at/under the noise floor).
# IDs   : Material (origin-name pedigree) -> line_id -> nil_id via the verified GDL
#         crosswalks (agent/gdl_flowering/out/{id_crosswalk_variant_to_lineid,
#         register_bc2s3}.csv). Scan joins genotypes on `pedigree` (= line_id), so
#         Genotype := line_id; nil_id attached as the stable serial label.
# Prep  TRANSFORM=none (DEFAULT-for-reproduction, raw scale, matches the original):
#         all 133 provided columns used AS-IS (no count scaling, no log, no below-
#         detection recoding); negatives/near-zeros kept. Sums/ratios are the file's
#         plant-level-then-averaged columns, NOT recomputed from species.
#       TRANSFORM=log (our count+log variant, NOT the High-PC1 method; comparison only,
#         decisions 2026-08-25):
#         species: below-detection (<=0) -> NA or 0 (BELOW_DET); scale to integer counts
#         K=1e9; log10(count+1). class sums: same. ratios: log10(ratio) only.
# Out   : data/zeal/pheno_lipid_all.csv (Genotype[=line_id], nil_id, Material, 133 traits)
#         data/zeal/lipid_material_nilid_crosswalk.csv (traceability)
# Scale : TRANSFORM=log (default) applies the High-PC1 log recipe above; TRANSFORM=none
#         emits the RAW adjusted means/sums/ratios as-is (same crosswalk + line_id
#         aggregation, no log) to data/zeal/pheno_lipid_raw.csv for a non-log-scale scan.
# Run   : Rscript scripts/zeal_lipid_pheno.R                 # log scale (default)
#         TRANSFORM=none Rscript scripts/zeal_lipid_pheno.R  # raw scale
# =============================================================================
suppressMessages({
  library(here)
  library(data.table)
})
source(here("scripts/logging.R"))

K <- as.numeric(Sys.getenv("LIPID_K", "1e9")) # integer-count scaling constant (decision D1)
# TRANSFORM: "log" (default) = High-PC1 log pipeline below (QTL-mapping scale).
#            "none"          = raw analysis: use the provided columns AS-IS (no count
#                              scaling, no log), skipping below-detection recoding, so the
#                              scan runs on the untransformed adjusted means/sums/ratios.
TRANSFORM <- Sys.getenv("TRANSFORM", "log")
stopifnot(TRANSFORM %in% c("log", "none"))

# ---- read + classify traits --------------------------------------------------
d <- fread(here("data/zeal/BZeaLipid_AdjustedMeans_wSums.csv"))
setnames(d, 1, "Material")
traits <- setdiff(names(d), "Material")
is_species <- grepl("^[A-Za-z]+_[0-9]+_[0-9]+$", traits) # CLASS_carbons_desat
is_sum <- grepl("^[A-Za-z]+$", traits) # bare class token
is_ratio <- grepl("^[A-Za-z]+_[A-Za-z]+$", traits) # classA_classB
species <- traits[is_species]
sums <- traits[is_sum]
ratios <- traits[is_ratio]
stopifnot(length(species) + length(sums) + length(ratios) == length(traits)) # partition covers all
log_info("traits: %d species + %d sums + %d ratios = %d", length(species), length(sums), length(ratios), length(traits))

if (TRANSFORM == "none") {
  # ---- raw scale: provided columns as-is, no count scaling, no log ------------
  X <- as.matrix(d[, ..traits])
  X[!is.finite(X)] <- NA
  colnames(X) <- traits
  log_info("TRANSFORM=none: raw scale, %d traits used as-is (no below-detection recoding, no log)", length(traits))
} else {
  # ---- species: below-detection handling, integer counts (xK), log10(count+1) ---
  # BELOW_DET env: "zero" (default) sets <=0 -> 0 (complete data => R/qtl's FAST batched
  # scan path, ~5 s/perm) ; "na" sets <=0 -> NA (statistically cleaner but ~13x slower,
  # since scanone re-runs the marker loop per distinct missing-pattern). The 0-spike from
  # "zero" only materially affects the few heavy-censoring species (DG_34_6, PC_32_4,
  # LPC_18_1, PC_32_1); use "na" (optionally a HEAVY-trait-only subset) to check those.
  BELOW_DET <- Sys.getenv("BELOW_DET", "zero")
  S <- as.matrix(d[, ..species])
  S[!is.finite(S)] <- NA
  belowdet <- sum(S <= 0, na.rm = TRUE)
  S[S <= 0] <- if (BELOW_DET == "na") NA else 0
  cnt <- round(S * K)
  log_info(
    "species: %d below-detection (<=0) -> %s (%.1f%%); K=%.0e -> counts median=%.0f max=%.0f",
    belowdet, toupper(BELOW_DET), 100 * belowdet / sum(is.finite(as.matrix(d[, ..species]))), K,
    median(cnt[cnt > 0], na.rm = TRUE), max(cnt, na.rm = TRUE)
  )
  species_log <- log10(cnt + 1) # NA (if BELOW_DET=na) stays NA
  colnames(species_log) <- species

  # ---- class sums: USE PROVIDED columns; below-detection->NA, counts, log10(count+1) ----
  # Do NOT recompute from species — the file's sums are plant-level then averaged.
  SM <- as.matrix(d[, ..sums])
  SM[!is.finite(SM)] <- NA
  bd_sm <- sum(SM <= 0, na.rm = TRUE)
  SM[SM <= 0] <- NA
  sums_log <- log10(round(SM * K) + 1)
  colnames(sums_log) <- sums
  log_info("class sums: PROVIDED columns (not recomputed); %d below-detection (<=0) -> NA", bd_sm)

  # ---- ratios: USE PROVIDED columns (per-plant ratio averaged); log10(ratio) only ----
  # The mean of per-plant ratios has far lower variance than a ratio of means -> sharper
  # peaks (per-plant/spike-in normalization; cf. hpc1 in the PNAS paper). No count scaling.
  R <- as.matrix(d[, ..ratios])
  R[!is.finite(R)] <- NA
  bd_r <- sum(R <= 0, na.rm = TRUE)
  R[R <= 0] <- NA # a ratio is strictly positive; <=0 is an artifact -> NA
  ratio_log <- log10(R)
  colnames(ratio_log) <- ratios
  log_info("ratios: PROVIDED per-plant-averaged columns; log10(ratio); %d <=0 -> NA", bd_r)

  # ---- assemble transformed trait matrix in original column order --------------
  X <- cbind(species_log, sums_log, ratio_log)[, traits, drop = FALSE]
} # end TRANSFORM == "log"

out <- data.table(Material = d$Material, X)

# ---- Material -> line_id -> nil_id -------------------------------------------
cw <- fread(here("agent/gdl_flowering/out/id_crosswalk_variant_to_lineid.csv"))
reg <- fread(here("agent/gdl_flowering/out/register_bc2s3.csv"))
out[, line_id := cw$line_id[match(Material, cw$variant)]]
out[, nil_id := reg$nil_id[match(line_id, reg$pedigree)]]
out[, taxon := reg$taxon[match(line_id, reg$pedigree)]]
n_res <- out[!is.na(nil_id), .N]
log_info("Material -> nil_id resolved: %d / %d (%.1f%%)", n_res, nrow(out), 100 * n_res / nrow(out))
unresolved <- out[is.na(nil_id), Material]
if (length(unresolved)) log_info("unresolved (kept, e.g. checks): %s", paste(head(unresolved, 6), collapse = ", "))

# crosswalk table for traceability
fwrite(
  out[, .(Material, line_id, nil_id, taxon, resolved = !is.na(nil_id))],
  here("data/zeal/lipid_material_nilid_crosswalk.csv")
)

# ---- aggregate to the scan key (line_id; unresolved keep Material) -----------
out[, Genotype := fifelse(is.na(line_id), Material, line_id)]
ncollapse <- nrow(out) - uniqueN(out$Genotype)
if (ncollapse > 0) log_info("%d Materials collapse onto a shared line_id (mean over the analysis-scale values)", ncollapse)
agg <- out[, c(
  lapply(.SD, function(v) mean(v, na.rm = TRUE)),
  .(nil_id = nil_id[1], taxon = taxon[1], Material = paste(Material, collapse = ";"))
), by = Genotype, .SDcols = traits]
for (tr in traits) agg[[tr]][is.nan(agg[[tr]])] <- NA # all-NA groups -> NA not NaN

setcolorder(agg, c("Genotype", "nil_id", "taxon", "Material", traits))
default_out <- if (TRANSFORM == "none") "data/zeal/pheno_lipid_raw.csv" else "data/zeal/pheno_lipid_all.csv"
outfile <- Sys.getenv("OUTFILE", here(default_out))
fwrite(agg, outfile)
log_info(
  "wrote %s: %d genotypes x %d traits (TRANSFORM=%s%s)", outfile, nrow(agg), length(traits),
  TRANSFORM, if (TRANSFORM == "log") sprintf(", BELOW_DET=%s", BELOW_DET) else ""
)
