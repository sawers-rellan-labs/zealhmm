#!/usr/bin/env Rscript
# =============================================================================
# Phase B — assemble the shareable BZea genotype release.
#
# 50K "variations" (SNP50K panel, B73 v5):
#   <caller>_mosaic  ANCESTRY inference (call_ancestry): rtiger, nnil, binhmm, lbimpute
#                    -> marker x line 0/1/2 matrix -> our VCF -> PLINK + tidy TSV + rds
#   gphwe            the GENOTYPE: authoritative cohort VCF bzea_50K_cohort.vcf.gz
#                    (bcftools mpileup | call -mv, HWE-prior MAP) shipped VERBATIM + PLINK
#                    from it + the gwas_nil subset rds. NO single-sample GL is shipped.
# 250K: the *previous* 2-state RTIGER introgression segments (Nirwan / inv4m paper).
#
# plus shared markers/ + lines/ tables, the 250K segments, a MANIFEST and README.
#
# IMPORTANT (see TERMINOLOGY.md): a mosaic is ANCESTRY, not a genotype. The PLINK/VCF
# 0/1/2 for a mosaic is the ANCESTRY dosage (0=B73, 1=het, 2=teosinte), encoded on the
# SNP's ref/alt alleles for tooling compatibility — it does NOT report the true allele
# at invariant sites. Only gphwe is an actual genotype. This is documented in the README.
#
# Requires plink2 on PATH (override with env PLINK2). Bulk output is gitignored (/release/).
# =============================================================================
suppressMessages({
  library(here)
  library(data.table)
})
source(here("scripts/logging.R"))

OUT <- here("release/bzea_genotypes")
PLINK2 <- Sys.getenv("PLINK2", "plink2")
for (d in c("snp50k", "markers", "lines", "250k")) {
  dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
}

# --- shared marker metadata: marker, chr, pos, ref, alt, cM (B73 v5) ----------
D <- readRDS(here("data/zeal/zeal_snp50k_dosage.rds"))
mkref <- as.data.table(D$markers)[, .(marker, chr = as.integer(chr), pos = as.integer(pos), ref, alt)]
cm <- fread(here("data/zeal/markers_snp50k_cm.tsv")) # marker, cm
mkref[, cm := cm$cm[match(marker, cm$marker)]]
setkey(mkref, marker)
fwrite(mkref, file.path(OUT, "markers", "snp50k_markers.tsv"), sep = "\t")
log_info("markers: %d SNP50K sites (B73 v5) with ref/alt/cM", nrow(mkref))

# --- shared line metadata -----------------------------------------------------
ss <- fread(here("data/zeal/samplesheet_3way.csv"))
linemeta <- unique(ss[!is.na(skim_id), .(pedigree, taxon, donor_accession, taxa_code, skim_id)])
fwrite(linemeta, file.path(OUT, "lines", "snp50k_lines.tsv"), sep = "\t")

# --- per-object exporter: 0/1/2 matrix -> TSV.gz + VCF -> PLINK bed ------------
export_50k <- function(obj, layer) {
  M <- readRDS(here(sprintf("data/zeal/zeal_%s.rds", obj)))
  state <- M$state
  mk <- data.table(marker = rownames(state))
  # de-duplicate line columns (keep first) so sample IDs are unique
  dup <- duplicated(colnames(state))
  if (any(dup)) {
    log_warn("%s: dropping %d duplicate line columns", obj, sum(dup))
    state <- state[, !dup, drop = FALSE]
  }
  # join marker coords/alleles; keep only markers with ref/alt (all SNP50K should)
  mk <- mkref[mk, on = "marker"]
  keep <- !is.na(mk$chr) & !is.na(mk$ref) & !is.na(mk$alt)
  if (any(!keep)) {
    log_warn("%s: dropping %d markers lacking ref/alt/coords", obj, sum(!keep))
    state <- state[keep, , drop = FALSE]
    mk <- mk[keep]
  }
  ord <- order(mk$chr, mk$pos)
  state <- state[ord, , drop = FALSE]
  mk <- mk[ord]
  base <- sprintf("bzea_snp50k_%s", obj)

  # 1) tidy TSV.gz (marker column + one column per line)
  tsv <- data.table(marker = mk$marker)
  tsv <- cbind(tsv, as.data.table(state))
  tsvf <- file.path(OUT, "snp50k", sprintf("%s_012.tsv.gz", base))
  fwrite(tsv, tsvf, sep = "\t")

  # 2) VCF (GT hardcall: 0->0/0, 1->0/1, 2->1/1, NA->./.) -> plink2 --make-bed
  gtmap <- c("0" = "0/0", "1" = "0/1", "2" = "1/1")
  GT <- matrix(gtmap[as.character(state)], nrow = nrow(state))
  GT[is.na(state)] <- "./."
  vcf <- data.table(
    `#CHROM` = paste0("chr", mk$chr), POS = mk$pos, ID = mk$marker, REF = mk$ref, ALT = mk$alt,
    QUAL = ".", FILTER = ".", INFO = ".", FORMAT = "GT"
  )
  vcf <- cbind(vcf, as.data.table(GT))
  setnames(vcf, c(names(vcf)[1:9], colnames(state)))
  vcff <- file.path(OUT, "snp50k", sprintf("%s.vcf", base))
  hdr <- c(
    "##fileformat=VCFv4.2",
    sprintf('##FILTER=<ID=PASS,Description="All filters passed">'),
    sprintf("##contig=<ID=chr%d>", sort(unique(mk$chr))),
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="0/1/2 = REF/HET/ALT (ancestry dosage for _mosaic; genotype for _gt)">'
  )
  writeLines(hdr, vcff)
  fwrite(vcf, vcff, sep = "\t", append = TRUE, col.names = TRUE, quote = FALSE)
  outp <- file.path(OUT, "snp50k", base)
  rc <- system2(PLINK2, c(
    "--vcf", vcff, "--make-bed", "--out", outp,
    "--allow-extra-chr", "--double-id", "--set-all-var-ids", "@:#"
  ), stdout = FALSE, stderr = FALSE)
  file.remove(vcff)
  if (rc != 0 || !file.exists(paste0(outp, ".bed"))) stop("plink2 failed for ", obj)

  # 3) copy the native RDS
  file.copy(here(sprintf("data/zeal/zeal_%s.rds", obj)),
    file.path(OUT, "snp50k", sprintf("%s.rds", base)),
    overwrite = TRUE
  )
  comp <- prop.table(table(factor(state, levels = 0:2)))
  log_info(
    "%s [%s]: %d markers x %d lines | 0=%.1f%% 1=%.1f%% 2=%.1f%% NA=%.1f%% -> bed+tsv+rds",
    obj, layer, nrow(state), ncol(state),
    100 * comp["0"], 100 * comp["1"], 100 * comp["2"], 100 * mean(is.na(state))
  )
}

for (o in c("rtiger_mosaic", "nnil_mosaic", "binhmm_mosaic", "lbimpute_mosaic")) export_50k(o, "ancestry mosaic")

# gpHWE GENOTYPE: ship the AUTHORITATIVE bcftools cohort VCF verbatim (no reconstruction) +
# PLINK derived from it. Genotypes are `bcftools mpileup | call -mv` (HWE-prior MAP). The
# gwas_nil analysis subset (zeal_gphwe_gt.rds) rides along. NO single-sample GL is shipped.
export_gphwe <- function() {
  src <- here("data/zeal/bzea_50K_cohort.vcf.gz")
  base <- "bzea_snp50k_gphwe"
  vcf <- file.path(OUT, "snp50k", paste0(base, ".vcf.gz"))
  file.copy(src, vcf, overwrite = TRUE)
  if (file.exists(paste0(src, ".csi"))) file.copy(paste0(src, ".csi"), paste0(vcf, ".csi"), overwrite = TRUE)
  outp <- file.path(OUT, "snp50k", base)
  rc <- system2(PLINK2, c(
    "--vcf", vcf, "--make-bed", "--out", outp,
    "--allow-extra-chr", "--double-id", "--set-all-var-ids", "@:#"
  ), stdout = FALSE, stderr = FALSE)
  if (rc != 0 || !file.exists(paste0(outp, ".bed"))) stop("plink2 failed for gphwe")
  unlink(paste0(outp, ".log"))
  file.copy(here("data/zeal/zeal_gphwe_gt.rds"),
    file.path(OUT, "snp50k", paste0(base, "_gwasnil.rds")),
    overwrite = TRUE
  )
  n <- length(readLines(paste0(outp, ".fam")))
  log_info("gphwe [genotype]: authoritative bcftools VCF (%d samples) + PLINK + gwasnil.rds", n)
}
export_gphwe()

# --- 250K (previous): 2-state RTIGER introgression segments -> long TSV.gz -----
k250 <- readRDS(here("data/zeal/rtiger_250K_calls_introfinder.rds"))
seg <- rbindlist(lapply(names(k250), function(nm) {
  d <- as.data.table(k250[[nm]])
  data.table(
    pedigree = sub("\\.B$", "", nm),
    chr = as.integer(sub("^chr", "", d[[1]])), start = as.integer(d[[2]]),
    end = as.integer(d[[3]]), state = d[[4]]
  )
}))[order(pedigree, chr, start)]
fwrite(seg, file.path(OUT, "250k", "bzea_250k_rtiger_introgression_segments.tsv.gz"), sep = "\t")
file.copy(here("data/zeal/rtiger_250K_calls_introfinder.rds"),
  file.path(OUT, "250k", "bzea_250k_rtiger_introgression.rds"),
  overwrite = TRUE
)
log_info("250K (previous): %d lines, %d segments (2-state B73/Introgression)", uniqueN(seg$pedigree), nrow(seg))

# --- README (tracked template) + MANIFEST (file, bytes, sha256) ----------------
unlink(list.files(file.path(OUT, "snp50k"), pattern = "\\.log$", full.names = TRUE)) # drop plink2 logs
file.copy(here("scripts/release_README.md"), file.path(OUT, "README.md"), overwrite = TRUE)
files <- list.files(OUT, recursive = TRUE, full.names = TRUE)
files <- sort(files[!grepl("/MANIFEST\\.tsv$", files)])
sha256 <- vapply(files, function(f) sub(" .*", "", system2("shasum", c("-a", "256", shQuote(f)), stdout = TRUE)), "")
man <- data.table(file = sub(paste0(OUT, "/"), "", files), bytes = file.size(files), sha256 = sha256)
fwrite(man, file.path(OUT, "MANIFEST.tsv"), sep = "\t")
log_info("wrote README + MANIFEST (%d files) | release assembled under %s", nrow(man), OUT)
