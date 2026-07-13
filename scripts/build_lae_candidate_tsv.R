#!/usr/bin/env Rscript
# =============================================================================
# LAE (leaves above the uppermost ear = leaf/node-number count) candidate gene:
# derive the caller TSV input + the notebook overlap CSV from the canonical CSV.
#
# Sources:
#   data/ref/leaf_number_candidate_genes_v5.csv   (hand-curated; tu1 / tunicate1)
#   data/teonam/eh_candidate_genes.tsv            (ear-height panel, overlaid: LAE = ear/
#     node position shares internode / GA / BR / auxin architecture with ear height)
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

first_clause <- function(x) {
  after <- trimws(sub("^[^—]*—\\s*", "", x)) # drop "<category> — " if present
  after <- sub("\\.\\s.*$", "", after) # keep first sentence
  sub("\\.$", "", after)
}

tsv <- data.table(
  symbol              = d$symbol,
  gene_id             = d$v5_gene_model,
  chr                 = as.integer(sub("^chr", "", d$chr)),
  start               = as.integer(d$start_bp),
  end                 = as.integer(d$end_bp),
  qtl_chen2019        = '""', # literal, matching the other candidate TSVs
  v5_canonical_symbol = d$symbol,
  pathway             = paste0("Tier", d$tier, " — ", d$protein, "; ", first_clause(d$trait_role))
)
# Overlay the ear-height candidate panel (LAE ear/node position shares internode / GA /
# BR / auxin architecture with ear height), deduping by gene id.
eh_path <- here("data/teonam/eh_candidate_genes.tsv")
if (file.exists(eh_path)) {
  eh <- fread(eh_path, encoding = "UTF-8")
  eh[, qtl_chen2019 := as.character(qtl_chen2019)]
  eh[is.na(qtl_chen2019), qtl_chen2019 := ""]
  tsv <- unique(rbind(tsv, eh, fill = TRUE), by = "gene_id")
  cat(sprintf("overlaid %d ear-height candidates\n", nrow(eh)))
} else {
  warning("eh_candidate_genes.tsv not found; LAE panel = leaf-number genes only")
}
setorder(tsv, chr, start)

tsv_path <- here("data/teonam/lae_candidate_genes.tsv")
fwrite(tsv, tsv_path, sep = "\t", quote = FALSE)
cat(sprintf("wrote %s (%d gene(s))\n", tsv_path, nrow(tsv)))

overlap <- tsv[, .(symbol, qtl = paste0("LAE(chr", chr, ")"), gene_id, chr, start, end)]
ov_path <- here("results/sim/zeal/lae_candidate_overlap.csv")
fwrite(overlap, ov_path)
cat(sprintf("wrote %s (%d gene(s))\n", ov_path, nrow(overlap)))
