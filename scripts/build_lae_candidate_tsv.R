#!/usr/bin/env Rscript
# =============================================================================
# LAE (leaves above the uppermost ear = leaf/node-number count) candidate gene:
# derive the caller TSV input + the notebook overlap CSV from the canonical CSV.
#
# Source:
#   data/ref/leaf_number_candidate_genes_v5.csv   (hand-curated, TU1-centric leaf-number
#     panel: tu1 upstream + its plastochron/meristem/polarity targets + data-driven
#     positional hits, all with B73 v5 gene models + coordinates; self-contained, no bulk
#     ear-height overlay)
# Regenerated:
#   data/teonam/lae_candidate_genes.tsv        (caller/GWAS candidate contract)
#   results/sim/zeal/lae_candidate_overlap.csv (Manhattan/lollipop gene overlay)
#
# TSV schema (shared candidate-gene contract, cf. en/eh *_candidate_genes.tsv):
#   symbol  gene_id  chr  start  end  qtl_chen2019  v5_canonical_symbol  pathway
# Overlap schema: symbol  qtl  gene_id  chr  start  end
#
# `pathway` is derived from tier + protein + the first clause of trait_role.
# Run: Rscript scripts/build_lae_candidate_tsv.R
# =============================================================================

suppressPackageStartupMessages(library(data.table))
library(here)

csv <- here("data/ref/leaf_number_candidate_genes_v5.csv")
d <- fread(csv, encoding = "UTF-8")

# Canonical CSV (v6) schema:
#   label  gene  entrez_or_model_id  chr  start_v5  end_v5  strand  tier  mechanism
#   length_bp  v5_gene_model
tsv <- data.table(
  symbol              = d$gene,
  gene_id             = d$v5_gene_model,
  chr                 = as.integer(sub("^chr", "", d$chr)),
  start               = as.integer(d$start_v5),
  end                 = as.integer(d$end_v5),
  qtl_chen2019        = '""', # literal, matching the other candidate TSVs
  v5_canonical_symbol = d$gene,
  pathway             = paste0("Tier", d$tier, " — ", d$mechanism)
)
# Safety net: drop any candidate without a placeable v5 coordinate (the current panel
# places every row) so an unplaceable row can never reach the map.
n_all <- nrow(tsv)
tsv <- tsv[!is.na(chr) & !is.na(start) & !is.na(end)]
if (nrow(tsv) < n_all) {
  cat(sprintf("dropped %d coordinate-less candidate(s) from the placed panel\n", n_all - nrow(tsv)))
}
setorder(tsv, chr, start)

tsv_path <- here("data/teonam/lae_candidate_genes.tsv")
fwrite(tsv, tsv_path, sep = "\t", quote = FALSE)
cat(sprintf("wrote %s (%d gene(s))\n", tsv_path, nrow(tsv)))

overlap <- tsv[, .(symbol, qtl = paste0("LAE(chr", chr, ")"), gene_id, chr, start, end)]
ov_path <- here("results/sim/zeal/lae_candidate_overlap.csv")
fwrite(overlap, ov_path)
cat(sprintf("wrote %s (%d gene(s))\n", ov_path, nrow(overlap)))
