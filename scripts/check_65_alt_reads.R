#!/usr/bin/env Rscript
# The 50K set (bzea_50K_cohort.vcf.gz) and the HQ_BZEA panel disagree on the ALT allele at a
# handful of sites. Before choosing which ALT to attach, ask: at those disagreement sites, do
# ANY reads land on ALT across our samples? If there are ~0 ALT reads, the choice makes no
# difference to any count.
#
#   step 1: sites where cohort-ALT (bzea_50K_cohort.vcf.gz) != HQ_BZEA-ALT
#   step 2: for the 11 skim samples with GATK counts, sum ALT reads at those sites, under BOTH
#           allele definitions (GATK infers alt from reads, so "alt reads" = reads matching the
#           attached alt base; we read the raw allelic_counts50K columns to get the actual base).
#
#   Rscript scripts/check_65_alt_reads.R    # sandbox DISABLED (rsstu automount)

suppressMessages(library(data.table))
here <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
BZ <- "/Volumes/rsstu/users/r/rrellan/BZea/bzeaseq"
COH <- file.path(BZ, "50K/results/joint/bzea_50K_cohort.vcf.gz") # the 50K set = cohort-called alleles
HQ <- file.path(BZ, "nilhmm/vcf/HQ_BZEA.vcf.gz") # HQ panel = teosinte-vs-B73 donor allele
RAW <- file.path(BZ, "50K/allelic_counts50K.tsv") # GATK, 51,991/sample, has REF/ALT base cols
BCF <- "/opt/homebrew/bin/bcftools"

query_alleles <- function(vcf, cols) {
  d <- fread(
    cmd = sprintf("%s query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\n' %s", BCF, shQuote(vcf)),
    header = FALSE, col.names = c("contig", "pos", cols[1], cols[2])
  )
  d[, chr := ifelse(grepl("^chr", contig), contig, paste0("chr", contig))]
  d[, pos := as.integer(pos)]
  d[, c("chr", "pos", cols[1], cols[2]), with = FALSE]
}

# --- 50K-set alleles (cohort VCF): the ALT the extractor now attaches ---
coh <- query_alleles(COH, c("c_ref", "c_alt"))
cat(sprintf("cohort VCF (50K set) records: %d\n", nrow(coh)))

# --- HQ_BZEA panel alleles: the designated donor ALT; the merge restricts it to the 50K positions.
#     Self-contained: queries HQ_BZEA directly (this comparison predates and does not depend on
#     data/zeal/snp50k_alleles.tsv, which is itself cohort-derived).
pan <- query_alleles(HQ, c("p_ref", "p_alt"))

m <- merge(coh, pan, by = c("chr", "pos"))
cat(sprintf("sites joined (cohort vs panel): %d\n", nrow(m)))
cat(sprintf("REF differs: %d ; ALT differs: %d\n", m[c_ref != p_ref, .N], m[c_alt != p_alt, .N]))
disc <- m[c_alt != p_alt]
cat(sprintf("\n%d ALT-disagreement sites. panel_alt -> cohort_alt:\n", nrow(disc)))
print(disc[, .N, by = .(p_alt, c_alt)][order(-N)])

# --- count reads: read the raw GATK table for the 11 skim samples, at the disagreement sites ---
disc[, pk := paste(sub("chr", "", chr), pos)]
skim_ids <- sub("\\.tsv$", "", list.files(file.path(here, "data/skimsweep/skim/counts_50k"), pattern = "\\.tsv$"))
cat(sprintf("\nscanning raw GATK table for %d skim samples at %d sites...\n", length(skim_ids), nrow(disc)))

# stream the 2.5G table once; keep only wanted samples & disc positions
disckey <- setNames(seq_len(nrow(disc)), disc$pk)
con <- file(RAW, "r")
invisible(readLines(con, 1)) # header
hits <- list()
nline <- 0L
repeat {
  ch <- readLines(con, n = 200000L)
  if (length(ch) == 0L) break
  nline <- nline + length(ch)
  parts <- tstrsplit(ch, "\t", fixed = TRUE)
  # cols: SAMPLE CONTIG POSITION REF_COUNT ALT_COUNT REF_NUCLEOTIDE ALT_NUCLEOTIDE
  keep <- parts[[1]] %in% skim_ids
  if (!any(keep)) next
  k <- paste(sub("chr", "", parts[[2]][keep]), parts[[3]][keep])
  sel <- k %in% names(disckey)
  if (!any(sel)) next
  hits[[length(hits) + 1L]] <- data.table(
    sample = parts[[1]][keep][sel], pk = k[sel],
    ref_count = as.integer(parts[[4]][keep][sel]),
    alt_count = as.integer(parts[[5]][keep][sel]),
    gatk_alt = parts[[7]][keep][sel]
  )
}
close(con)
H <- rbindlist(hits)
cat(sprintf("raw lines scanned: %d ; rows at disc sites (11 samples): %d\n", nline, nrow(H)))

H <- merge(H, disc[, .(pk, chr, pos, p_alt, c_alt)], by = "pk")
cat(sprintf(
  "\ntotal ALT reads (GATK, its inferred base) across all 11 samples at the %d sites: %d\n",
  nrow(disc), sum(H$alt_count)
))
cat(sprintf(
  "rows where GATK's inferred alt base == cohort_alt: %d (alt reads there: %d)\n",
  H[gatk_alt == c_alt, .N], H[gatk_alt == c_alt, sum(alt_count)]
))
cat(sprintf(
  "rows where GATK's inferred alt base == panel_alt : %d (alt reads there: %d)\n",
  H[gatk_alt == p_alt, .N], H[gatk_alt == p_alt, sum(alt_count)]
))
cat("\nVERDICT: if both 'alt reads there' totals are ~0, it makes no difference which ALT allele we attach.\n")
