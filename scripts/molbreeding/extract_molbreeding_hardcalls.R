#!/usr/bin/env Rscript
# Extract MolBreeding GBTS genotype HARD CALLS on the v5 wideseq-filtered grid
# (9,157 donor-informative sites) = the TRUTH genotype input for the ZEAL nnil
# calibration. The count pipeline (molbreeding_to_gatk_table.py) DROPPED the vendor
# genotype calls, keeping only read counts; this recovers them from the raw delivery
# and lands them on the SAME v5 wsfilt sites + REF=B73/ALT=donor polarity as the
# counts (verified by per-cell concordance at the end).
#
# Conventions matched to the count pipeline EXACTLY:
#   - REF/ALT taken from the wsfilt count table (authoritative site polarity), NOT recomputed
#   - v3->v5 remaps coordinates only, NO allele flip (molbreeding_gatk_table_to_v5.py)
#   - g: 0 hom-REF (recurrent/B73), 1 het, 2 hom-ALT (donor); NA = missing/other
#   - sample well_id -> PN#_SID# via the sample map (same names as the count files)
#
# In : data/zeal/molbreeding/source/All.Genotype.xls          (raw GBTS delivery, v3, genotype() cols)
#      data/zeal/molbreeding/gatk_table_SNP_wsfilt_v5.tsv      (v5 wsfilt sites + REF/ALT polarity)
#      data/zeal/molbreeding/sites_v5_SNP.tsv                  (v3 -> v5 liftover)
#      data/zeal/molbreeding/molbreeding_sample_map.tsv        (well_id -> PN#_SID# + pedigree)
#      data/zeal/molbreeding/counts_targetseq/<PN>.tsv         (110x counts, for validation)
# Out: data/zeal/molbreeding/molbreeding_hardcalls_wsfilt.tsv  (name, pedigree, marker, chr, pos, g)
suppressMessages({
  library(here)
  library(data.table)
})
MB <- here("data/zeal/molbreeding")
SRC <- file.path(MB, "source/All.Genotype.xls") # tab-separated despite .xls
WSF <- file.path(MB, "gatk_table_SNP_wsfilt_v5.tsv")
LIFT <- file.path(MB, "sites_v5_SNP.tsv")
SMAP <- file.path(MB, "molbreeding_sample_map.tsv")
COUNTS <- file.path(MB, "counts_targetseq")
OUT <- file.path(MB, "molbreeding_hardcalls_wsfilt.tsv")

# --- authoritative v5 wsfilt sites + REF=B73/ALT=donor polarity (from the count table)
wsf <- fread(WSF)
sites <- unique(wsf[, .(
  chr = as.integer(sub("chr", "", CONTIG)), pos = POSITION,
  ref = REF_NUCLEOTIDE, alt = ALT_NUCLEOTIDE
)])

# --- raw genotype delivery (v3), melt the per-well genotype() columns to long
src <- fread(SRC, sep = "\t", header = TRUE)
gcols <- grep("^genotype\\(", names(src), value = TRUE)
long <- melt(src,
  id.vars = c("chrom", "position"), measure.vars = gcols,
  variable.name = "gcol", value.name = "gt", variable.factor = FALSE
)
long[, well := sub("^genotype\\((.*)\\)$", "\\1", gcol)][, gcol := NULL]

# --- v3 -> v5 (coordinates only, no allele flip), then restrict to the wsfilt sites
lift <- fread(LIFT)[, .(chr_v3, pos_v3, chr_v5, pos_v5)]
long <- merge(long, lift, by.x = c("chrom", "position"), by.y = c("chr_v3", "pos_v3"))
long <- merge(long, sites, by.x = c("chr_v5", "pos_v5"), by.y = c("chr", "pos"))

# --- genotype string (e.g. "AA","AG","GG") vs the site REF/ALT -> g in {0,1,2}
gt_to_g <- function(gt, ref, alt) {
  gt <- toupper(trimws(as.character(gt)))
  a1 <- substr(gt, 1, 1)
  a2 <- substr(gt, 2, 2)
  miss <- is.na(gt) | gt %in% c("", "NN", "NA", "./.", ".", "--") | nchar(gt) != 2L
  nref <- (a1 == ref) + (a2 == ref)
  nalt <- (a1 == alt) + (a2 == alt)
  g <- rep(NA_integer_, length(gt))
  ok <- !miss & (nref + nalt) == 2L
  g[ok & nref == 2L] <- 0L
  g[ok & nref == 1L & nalt == 1L] <- 1L
  g[ok & nalt == 2L] <- 2L
  g
}
long[, g := gt_to_g(gt, ref, alt)]

# --- well_id -> PN#_SID# + canonical pedigree (drop registry ".B")
smap <- fread(SMAP)
long[, name := smap$sample_id[match(well, smap$well_id)]]
long[, pedigree := sub("\\.B$", "", smap$genotype[match(well, smap$well_id)])]

out <- long[!is.na(name), .(name, pedigree,
  marker = paste0("S", chr_v5, "_", pos_v5), chr = chr_v5, pos = pos_v5, g
)][order(name, chr, pos)]
fwrite(out, OUT, sep = "\t")
cat(sprintf(
  "wrote %s: %d samples x %d sites (%d rows); g missing/NA = %.2f%%\n",
  basename(OUT), uniqueN(out$name), uniqueN(out$marker), nrow(out), 100 * mean(is.na(out$g))
))
cat("per-sample g composition (0 REF / 1 HET / 2 ALT / NA):\n")
print(out[, .(
  REF = sum(g == 0, na.rm = TRUE), HET = sum(g == 1, na.rm = TRUE),
  ALT = sum(g == 2, na.rm = TRUE), NA_ = sum(is.na(g))
), by = name])

# --- VALIDATION: vendor g vs the 110x counts (proves polarity + site alignment)
count_call <- function(rc, ac, min_tot = 10L, band = 0.1) {
  tot <- rc + ac
  af <- ac / pmax(tot, 1)
  ifelse(tot < min_tot, NA_integer_, ifelse(af < band, 0L, ifelse(af > 1 - band, 2L, 1L)))
}
val <- rbindlist(lapply(unique(out$name), function(nm) {
  f <- file.path(COUNTS, paste0(nm, ".tsv"))
  if (!file.exists(f)) {
    return(NULL)
  }
  cf <- fread(f, header = FALSE, col.names = c("contig", "pos", "rb", "rc", "ab", "ac"))
  cf[, chr := as.integer(sub("chr", "", contig))]
  cf[, gc := count_call(rc, ac)]
  merge(out[name == nm, .(chr, pos, g)], cf[, .(chr, pos, gc)], by = c("chr", "pos"))[, name := nm]
}))
v <- val[!is.na(g) & !is.na(gc)]
cat(sprintf(
  "\nVALIDATION vs 110x counts: %d comparable cells | g-vs-count concordance = %.4f\n",
  nrow(v), mean(v$g == v$gc)
))
cat("confusion (rows = vendor g, cols = count-derived):\n")
print(table(vendor_g = v$g, count = v$gc))
