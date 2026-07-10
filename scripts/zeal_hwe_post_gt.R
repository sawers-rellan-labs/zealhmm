#!/usr/bin/env Rscript
# =============================================================================
# ZEAL/BZea — the HWE-posterior GENOTYPE object, extracted from the AUTHORITATIVE cohort VCF.
#
# bzea_50K_cohort.vcf.gz was produced by `bcftools mpileup -f <B73 v5> -R <sites> |
# bcftools call -mv` (HWE-prior MAP genotypes = "HWE-posterior"), the cohort set Fausto sent to
# Jim Holland. This is the REAL genotype the project uses/shares — NOT a call_gt
# reconstruction, and NOT single-sample argmax-GL. REF = B73 allele, so GT dosage
# 0/1/2 = B73 / het / teosinte, matching the ancestry-mosaic state convention.
#
# Restricts to the gwas_nil panel; renames columns skim_id -> pedigree.
# Output: data/zeal/zeal_hwe_post_gt.rds  list(markers, state[marker x line 0/1/2], lines)
# =============================================================================
suppressMessages({
  library(here)
  library(data.table)
})
source(here("scripts/logging.R"))
VCF <- here("data/zeal/bzea_50K_cohort.vcf.gz")
BCFTOOLS <- Sys.getenv("BCFTOOLS", "bcftools")

ss <- fread(here("data/zeal/samplesheet_3way.csv"))[gwas_nil == TRUE & !is.na(skim_id)]
vsamp <- system2(BCFTOOLS, c("query", "-l", shQuote(VCF)), stdout = TRUE)
nil <- ss[skim_id %in% vsamp][!duplicated(skim_id)]
sf <- tempfile()
writeLines(nil$skim_id, sf)
log_info("gwas_nil in VCF: %d / %d lines", nrow(nil), uniqueN(ss$skim_id))

# extract GT (subset to the panel via -S; sample order follows the file)
tf <- tempfile(fileext = ".tsv")
system2(BCFTOOLS, c(
  "query", "-S", shQuote(sf), "-f", shQuote("%CHROM\\t%POS[\\t%GT]\\n"), shQuote(VCF)
), stdout = tf)
G <- fread(tf, header = FALSE)
setnames(G, c("chrom", "pos", nil$skim_id))
G[, marker := sprintf("S%d_%d", as.integer(sub("^chr", "", chrom)), pos)]

# restrict to the SNP50K panel roster (dosage markers); order by chr,pos
mkref <- as.data.table(readRDS(here("data/zeal/zeal_snp50k_dosage.rds"))$markers)[
  , .(marker, chr = as.integer(chr), pos = as.integer(pos))
]
G <- G[marker %in% mkref$marker]
log_info("markers: %d VCF sites on the SNP50K panel (roster %d)", nrow(G), nrow(mkref))
mk <- mkref[match(G$marker, marker)]
ord <- order(mk$chr, mk$pos)

# GT -> 0/1/2 dosage (NA for ./., multiallelic, anything not biallelic REF/ALT)
lut <- c("0/0" = 0L, "0|0" = 0L, "0/1" = 1L, "1/0" = 1L, "0|1" = 1L, "1|0" = 1L, "1/1" = 2L, "1|1" = 2L)
Mgt <- as.matrix(G[, nil$skim_id, with = FALSE])
state <- matrix(lut[Mgt], nrow = nrow(Mgt), dimnames = list(G$marker, nil$pedigree))
state <- state[ord, , drop = FALSE]
mk <- mk[ord]
rownames(state) <- mk$marker

# QC: B73 checks ~ 0 teosinte; NIL panel low presence
comp <- prop.table(table(factor(state, levels = 0:2)))
b73 <- intersect(ss[is_B73 == TRUE, pedigree], colnames(state))
log_info(
  "hwe_post_gt: %d markers x %d lines | 0=%.1f%% 1=%.1f%% 2=%.1f%% NA=%.1f%% | presence=%.3f",
  nrow(state), ncol(state), 100 * comp["0"], 100 * comp["1"], 100 * comp["2"],
  100 * mean(is.na(state)), sum(comp[c("1", "2")])
)
if (length(b73)) {
  log_info("QC: B73 checks (n=%d) mean dosage = %.3f (expect ~0)", length(b73), mean(state[, b73], na.rm = TRUE))
}
saveRDS(
  list(markers = mk, state = state, lines = data.table(skim_id = nil$skim_id, pedigree = nil$pedigree)),
  here("data/zeal/zeal_hwe_post_gt.rds")
)
log_info("wrote data/zeal/zeal_hwe_post_gt.rds")
