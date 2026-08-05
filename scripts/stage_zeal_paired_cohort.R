#!/usr/bin/env Rscript
# Stage the ZEAL skim-x-BrB PAIRED COHORT for the two-panel fragment-size cM validation figure.
#
# Cohort = NILs with BOTH SNP50K skim and BRB-seq data (skim_brbseq_correspondence.csv,
# in_skim & in_brbseq) = 335, minus the exclusions in reference/zeal_paired_cohort_exclusions.csv
# = 332 NILs.
#
# THE NIL ID IS THE PEDIGREE STRING. The PN#_SID# prefixes are plate/well coordinates LOCAL to
# each platform: of 238 prefix strings shared between the skim and BrB namespaces, 238 (100%)
# name a DIFFERENT NIL in each. 36 of those collisions fall inside this cohort. So every staged
# file is NAMED by pedigree, and the in-file sample identifier (`SAMPLE` in bins, `sample` in
# pangene) is REWRITTEN to the pedigree. Downstream code may then take `name` from the basename
# without ambiguity. Never join skim to BrB on a prefix.
#
# Sources (all on the rsstu automount; run with the sandbox DISABLED or the files read as missing):
#   skim counts  $BZ/50K/results/allelic_counts/<skim_prefix>_allele_counts.tsv   (49,002 sites)
#   skim bins    $BZ/ancestry/<skim_prefix>_bin_genotypes.tsv
#   BrB pangene  $CAS/results/<species>/pangene/<brbseq_prefix>.pangene_counts.tsv
# NOT data/skimsweep/skim/counts_50k/: those zealtiger-derived copies are a DIFFERENT frame
# (51,991 sites) whose ALT column is a synthetic placeholder. Never mix the two frames.
#
# Resumable: a file already present is skipped, so re-running after an interrupt is cheap.
#   Rscript scripts/stage_zeal_paired_cohort.R
#   ZEAL_STAGE_FORCE=1 Rscript scripts/stage_zeal_paired_cohort.R   # re-copy everything

suppressMessages({
  library(data.table)
  library(here)
})
source(here::here("scripts/logging.R"))

BZ <- "/Volumes/rsstu/users/r/rrellan/BZea/bzeaseq"
CAS <- "/Volumes/rsstu/users/r/rrellan/tlaloc/cassini"
CORR <- here::here("data/zeal/correspondence")
OUT <- here::here("data/zeal/paired_cohort")
FORCE <- Sys.getenv("ZEAL_STAGE_FORCE") != ""

for (d in c("skim/counts_50k", "skim/bins", "brb/pangene")) {
  dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
}

# canonical pedigree: drop the trailing ".B" library suffix
canon <- function(x) sub("\\.B$", "", trimws(x))

# ---- 1. build the cohort, keyed on pedigree ---------------------------------
corr <- fread(file.path(CORR, "skim_brbseq_correspondence.csv"))
coh <- corr[in_skim == TRUE & in_brbseq == TRUE,
  .(brbseq_prefix = trimws(unlist(strsplit(brbseq_prefix, ";")))),
  by = .(pedigree = canon(pedigree), skim_prefix = trimws(skim_prefix))
]
log_info("[stage] paired cohort before exclusions: %d NILs, %d BrB libraries", uniqueN(coh$pedigree), nrow(coh))

excl <- fread(here::here("reference/zeal_paired_cohort_exclusions.csv"))
excl[, pedigree := canon(pedigree)]
log_info("[stage] exclusions (%s):", basename("reference/zeal_paired_cohort_exclusions.csv"))
print(excl[, .(pedigree, reason_code, decided_by)])
coh <- coh[!pedigree %in% excl$pedigree]
log_info("[stage] cohort after exclusions: %d NILs, %d BrB libraries", uniqueN(coh$pedigree), nrow(coh))

if (uniqueN(coh$pedigree) != nrow(coh)) {
  log_error("[stage] a pedigree still maps to >1 BrB library after exclusions; refusing to guess a collapse rule")
  print(coh[pedigree %in% coh[, .N, by = pedigree][N > 1]$pedigree])
  stop("multi-library pedigrees remain")
}

# ---- 2. integrity: does each namespace agree with the correspondence pedigree? ----
skped <- fread(file.path(CORR, "skim_sample_pedigree.csv"))[, .(skim_prefix = sample, ped_skim = canon(pedigree))]
brbm <- fread(file.path(CORR, "brbseq_metadata_master.csv"))[, .(brbseq_prefix = sample_id, ped_brb = canon(Genotype))]
chk <- merge(merge(coh, skped, by = "skim_prefix", all.x = TRUE), brbm, by = "brbseq_prefix", all.x = TRUE)
bad <- chk[(!is.na(ped_skim) & ped_skim != pedigree) | (!is.na(ped_brb) & ped_brb != pedigree)]
if (nrow(bad)) {
  log_error("[stage] %d cohort lines where a namespace disagrees with the correspondence pedigree:", nrow(bad))
  print(head(bad[, .(pedigree, skim_prefix, ped_skim, brbseq_prefix, ped_brb)], 20))
  stop("pedigree disagreement across namespaces; resolve before staging")
}
log_info("[stage] integrity OK: skim and BrB namespaces both agree with the correspondence pedigree on all %d lines", nrow(chk))

# ---- 3. locate sources ------------------------------------------------------
coh[, taxon := sub("\\..*$", "", pedigree)]
log_info("[stage] cohort by taxon: %s", paste(sprintf("%s=%d", coh[, .N, by = taxon][order(taxon)]$taxon, coh[, .N, by = taxon][order(taxon)]$N), collapse = " "))

# pangene: index every species dir by prefix (a library lives under its donor taxon's dir)
pg_idx <- rbindlist(lapply(list.dirs(file.path(CAS, "results"), recursive = FALSE), function(d) {
  fs <- list.files(file.path(d, "pangene"), pattern = "\\.pangene_counts\\.tsv$", full.names = TRUE)
  if (!length(fs)) NULL else data.table(brbseq_prefix = sub("\\.pangene_counts\\.tsv$", "", basename(fs)), species = basename(d), src_pangene = fs)
}))
log_info("[stage] cassini pangene index: %d files, %d unique prefixes", nrow(pg_idx), uniqueN(pg_idx$brbseq_prefix))
coh <- merge(coh, pg_idx, by = "brbseq_prefix", all.x = TRUE)

# skim counts: the BAM read tally (GATK CollectAllelicCounts), extracted from the aggregated
# 2.5 GB table by scripts/extract_gatk_counts_50k.sh into a shared prefix-keyed cache, already
# filtered to the 49,002-site biallelic panel. NOT 50K/results/allelic_counts/, which is `AD` read
# out of the joint-called cohort.vcf.gz (mediated by bcftools call) and disagrees with the tally on
# ~21% of covered sites. One counter everywhere. See DATA.md "The SNP50K panel".
GATK_CACHE <- here::here("data/zeal/gatk_counts_50k")
dir.create(GATK_CACHE, recursive = TRUE, showWarnings = FALSE)
lst <- tempfile()
writeLines(sort(unique(coh$skim_prefix)), lst)
rc <- system2("bash", c(here::here("scripts/extract_gatk_counts_50k.sh"), lst, GATK_CACHE))
if (rc != 0) stop("extract_gatk_counts_50k.sh failed (rc=", rc, ")")
coh[, src_counts := file.path(GATK_CACHE, sprintf("%s.tsv", skim_prefix))]
coh[, src_bins := file.path(BZ, "ancestry", sprintf("%s_bin_genotypes.tsv", skim_prefix))]
coh[, `:=`(
  has_counts = file.exists(src_counts), has_bins = file.exists(src_bins),
  has_pangene = !is.na(src_pangene) & file.exists(fifelse(is.na(src_pangene), "", src_pangene))
)]
log_info(
  "[stage] source availability: counts %d/%d | bins %d/%d | pangene %d/%d",
  sum(coh$has_counts), nrow(coh), sum(coh$has_bins), nrow(coh), sum(coh$has_pangene), nrow(coh)
)
if (coh[!(has_counts & has_bins & has_pangene), .N]) {
  log_warn("[stage] %d lines missing at least one input; they are dropped from the staged cohort:", coh[!(has_counts & has_bins & has_pangene), .N])
  print(coh[!(has_counts & has_bins & has_pangene), .(pedigree, skim_prefix, brbseq_prefix, has_counts, has_bins, has_pangene)])
}
coh <- coh[has_counts & has_bins & has_pangene]
setorder(coh, pedigree)
log_info("[stage] staging %d complete lines", nrow(coh))

# ---- 4. stage, renaming the in-file sample id to the pedigree ---------------
coh[, `:=`(
  dst_counts = file.path(OUT, "skim/counts_50k", sprintf("%s.tsv", pedigree)),
  dst_bins = file.path(OUT, "skim/bins", sprintf("%s.tsv", pedigree)),
  dst_pangene = file.path(OUT, "brb/pangene", sprintf("%s.pangene_counts.tsv", pedigree))
)]

t0 <- Sys.time()
n_copied <- 0L
for (i in seq_len(nrow(coh))) {
  r <- coh[i]
  # counts: no in-file sample id (name comes from the basename), so a plain copy is correct
  if (FORCE || !file.exists(r$dst_counts)) {
    file.copy(r$src_counts, r$dst_counts, overwrite = TRUE)
    n_copied <- n_copied + 1L
  }
  # bins: rewrite SAMPLE -> pedigree
  if (FORCE || !file.exists(r$dst_bins)) {
    b <- fread(r$src_bins)
    b[, SAMPLE := r$pedigree]
    fwrite(b, r$dst_bins, sep = "\t")
    n_copied <- n_copied + 1L
  }
  # pangene: rewrite sample -> pedigree
  if (FORCE || !file.exists(r$dst_pangene)) {
    p <- fread(r$src_pangene)
    p[, sample := r$pedigree]
    fwrite(p, r$dst_pangene, sep = "\t")
    n_copied <- n_copied + 1L
  }
  if (i %% 25 == 0 || i == nrow(coh)) {
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    log_info(
      ">>> %d/%d lines | %d files written | elapsed %.1f min | ETA ~%.1f min",
      i, nrow(coh), n_copied, el, (el / i) * (nrow(coh) - i)
    )
  }
}

# ---- 5. manifest ------------------------------------------------------------
man <- coh[, .(
  pedigree, taxon, species, skim_prefix, brbseq_prefix,
  counts_50k = basename(dst_counts), bins = basename(dst_bins), pangene = basename(dst_pangene)
)]
fwrite(man, file.path(OUT, "cohort_manifest.csv"))
log_info("[stage] wrote %s (%d lines)", file.path(OUT, "cohort_manifest.csv"), nrow(man))
sz <- sum(file.size(c(coh$dst_counts, coh$dst_bins, coh$dst_pangene)), na.rm = TRUE)
log_info("[stage] staged %d NILs, %.2f GB on disk", nrow(man), sz / 1e9)
log_info("[stage] taxon x species crosswalk:")
print(man[, .N, by = .(taxon, species)][order(taxon)])
