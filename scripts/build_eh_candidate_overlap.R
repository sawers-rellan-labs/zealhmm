#!/usr/bin/env Rscript
# =============================================================================
# Ear-height (EH) Manhattan/lollipop gene overlay, derived from the tracked,
# hand-maintained candidate TSV (EH has no canonical ref CSV / TSV builder).
#
# Source:      data/teonam/eh_candidate_genes.tsv   (tracked; edit this to add genes)
# Regenerated: results/sim/zeal/eh_candidate_overlap.csv
#   overlay schema: symbol  qtl  gene_id  chr  start  end
#
# Run: Rscript scripts/build_eh_candidate_overlap.R
# =============================================================================

suppressPackageStartupMessages(library(data.table))
library(here)

tsv <- fread(here("data/teonam/eh_candidate_genes.tsv"), encoding = "UTF-8")
tsv[, chr := as.integer(chr)]
setorder(tsv, chr, start)

overlap <- tsv[, .(
  symbol,
  qtl = sprintf("EH(chr%d)", chr), gene_id, chr, start, end
)]
ov_path <- here("results/sim/zeal/eh_candidate_overlap.csv")
fwrite(overlap, ov_path)
cat(sprintf("wrote %s (%d genes)\n", ov_path, nrow(overlap)))
