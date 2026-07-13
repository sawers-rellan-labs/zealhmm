#!/usr/bin/env Rscript
# =============================================================================
# Regenerate the per-trait candidate CONTRACT tables the notebooks read, from the
# normalized registry (reference/candidate_{genes,evidence}.csv — the tracked source
# of truth). Replaces the old per-trait build_*_candidate_tsv.R scripts.
#
# Emits, for each ZEAL trait present in the registry:
#   data/teonam/<trait>_candidate_genes.tsv   (contract: symbol gene_id chr start end
#       qtl_chen2019 v5_canonical_symbol pathway)  -- consumed by the GWAS + R/qtl notebooks
#   results/sim/zeal/<trait>_candidate_overlap.csv (symbol qtl gene_id chr start end)
#       -- the Manhattan/lollipop gene overlay
# Uncloned loci (locus:<symbol> keys, e.g. mhl1) emit a blank gene_id, as before.
# Run: Rscript scripts/build_candidate_tsvs.R
# =============================================================================
suppressMessages({
  library(data.table)
  library(here)
})

genes <- fread(here("reference/candidate_genes.csv"), encoding = "UTF-8")
ev <- fread(here("reference/candidate_evidence.csv"), encoding = "UTF-8")
reg <- merge(ev, genes, by = "gene_id", all.x = TRUE)
reg[, out_gene_id := fifelse(grepl("^locus:", gene_id), "", gene_id)]
reg[, qtl_out := fifelse(is.na(qtl) | qtl == "", '""', qtl)] # literal "" matches the legacy TSVs

dir.create(here("results/sim/zeal"), recursive = TRUE, showWarnings = FALSE)
for (tr in sort(unique(reg$trait))) {
  d <- reg[trait == tr][order(chr, start)]
  # Only placeable candidates reach the notebooks. Coordinate-less loci (e.g. rt1 —
  # a classical uncloned locus, Jenkins 1930, with no v5 coordinates) stay recorded in
  # the registry but are not emitted as plotted/tabled candidates.
  dp <- d[!is.na(chr) & !is.na(start)]
  dropped <- setdiff(d$symbol, dp$symbol)
  tsv <- dp[, .(
    symbol,
    gene_id = out_gene_id, chr, start, end,
    qtl_chen2019 = qtl_out, v5_canonical_symbol = symbol, pathway = functional_note
  )]
  fwrite(tsv, here(sprintf("data/teonam/%s_candidate_genes.tsv", tr)), sep = "\t", quote = FALSE)
  ov <- dp[, .(
    symbol,
    qtl = sprintf("%s(chr%d)", toupper(tr), chr), gene_id = out_gene_id, chr, start, end
  )]
  fwrite(ov, here(sprintf("results/sim/zeal/%s_candidate_overlap.csv", tr)))
  cat(sprintf(
    "%-10s %2d placed%s\n", tr, nrow(dp),
    if (length(dropped)) sprintf("  (registry-only, unmapped: %s)", paste(dropped, collapse = ", ")) else ""
  ))
}
