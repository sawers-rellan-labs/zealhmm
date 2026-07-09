#!/usr/bin/env Rscript
# ZEAL/BZea — binhmm/nilhmm ancestry mosaic from the REAL per-bin GENOTYPE calls in
# bzeaseq's all_samples_bin_genotypes.tsv (the pipeline's own per-1Mb-bin ancestry;
# presence 0.108 matches the established ~0.10 — the per-bin call recovers the sparse
# 0.4x teosinte signal that a per-SNP recall misses). We map the GENOTYPE column
# REF/HET/ALT -> 0/1/2 and expand bins onto the SNP50K markers. No re-calling.
# Output: data/zeal/zeal_binhmm_mosaic.rds  list(markers, state[marker x line], lines)
suppressMessages({
  library(here)
  library(data.table)
})
source(here("scripts/logging.R"))

bins <- fread(here("data/zeal/all_samples_bin_genotypes.tsv"),
  select = c("SAMPLE", "CONTIG", "BIN_START", "BIN_END", "GENOTYPE")
)
setnames(bins, c("name", "contig", "bin_start", "bin_end", "gt"))
bins[, chr := as.integer(sub("^chr", "", contig))]
bins[, state := fcase(gt == "REF", 0L, gt == "HET", 1L, gt == "ALT", 2L, default = NA_integer_)]
bins <- bins[!is.na(chr) & chr %in% 1:10 & !is.na(state)]

mk <- fread(here("data/zeal/markers_snp50k_v5.tsv"))[, .(marker, chr = as.integer(chr), pos = as.integer(pos))]
setorder(mk, chr, pos)
ss <- fread(here("data/zeal/samplesheet_3way.csv"))[gwas_nil == TRUE & !is.na(skim_id)]
panel <- ss[skim_id %in% unique(bins$name), skim_id]
ped <- ss[match(panel, skim_id), pedigree]
bins <- bins[name %in% panel]
log_info("binhmm (from bzeaseq per-bin GENOTYPE): %d panel lines x %d markers | %d bins", length(panel), nrow(mk), nrow(bins))

# expand per-bin GENOTYPE onto markers (findInterval on bin starts, per name x chr)
setkey(bins, name, chr, bin_start)
state <- matrix(NA_integer_, nrow = nrow(mk), ncol = length(panel), dimnames = list(mk$marker, ped))
mk_by_chr <- split(seq_len(nrow(mk)), mk$chr)
for (j in seq_along(panel)) {
  bj <- bins[.(panel[j])]
  for (ch in names(mk_by_chr)) {
    sc <- bj[chr == as.integer(ch)][order(bin_start)]
    if (!nrow(sc)) next
    ri <- mk_by_chr[[ch]]
    idx <- findInterval(mk$pos[ri], sc$bin_start)
    idx[idx < 1L] <- 1L
    state[ri, j] <- sc$state[idx]
  }
}
comp <- prop.table(table(factor(state, levels = 0:2)))
log_info(
  "binhmm mosaic: %d x %d | NA=%.2f%% | B73=%.1f%% het=%.1f%% teo=%.1f%% | PRESENCE=%.3f (established ~0.101)",
  nrow(state), ncol(state), 100 * mean(is.na(state)), 100 * comp["0"], 100 * comp["1"], 100 * comp["2"], sum(comp[c("1", "2")])
)
saveRDS(
  list(markers = mk, state = state, lines = data.table(skim_id = panel, pedigree = ped)),
  here("data/zeal/zeal_binhmm_mosaic.rds")
)
log_info("wrote data/zeal/zeal_binhmm_mosaic.rds")
