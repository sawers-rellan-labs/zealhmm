#!/usr/bin/env Rscript
# =============================================================================
# Brace-root (NBR = number of brace roots) candidate genes: derive the caller
# TSV input + the notebook overlap CSV from the single canonical annotation CSV.
#
# Single source of truth:
#   data/ref/brace_root_genes_v5.csv   (hand-curated, rich; v5/v4/v3 ids + refs)
# Regenerated:
#   data/teonam/nbr_candidate_genes.tsv        (caller/GWAS candidate contract)
#   results/sim/zeal/nbr_candidate_overlap.csv (Manhattan/lollipop gene overlay)
#
# Genes with no v5 coordinates (e.g. rt1, a classical uncloned locus) cannot be
# placed on the genome and are dropped from both outputs (logged).
#
# TSV schema (shared candidate-gene contract, cf. dta/eh/en *_candidate_genes.tsv):
#   symbol  gene_id  chr  start  end  qtl_chen2019  v5_canonical_symbol  pathway
# Overlap schema (cf. *_candidate_overlap.csv):
#   symbol  qtl  gene_id  chr  start  end
#
# `pathway` is derived deterministically from tier + phenotype_effect + the first
# clause of the phenotype text, so editing the CSV is the only edit needed.
#
# Run: Rscript scripts/build_brace_root_candidate_tsv.R
# =============================================================================

suppressPackageStartupMessages(library(data.table))
library(here)

csv <- here("data/ref/brace_root_genes_v5.csv")
d <- fread(csv, encoding = "UTF-8", colClasses = list(character = c("v5_id", "chr")))

# drop uncloned / unplaced loci (no v5 id or coordinates)
placed <- d[v5_id != "" & !is.na(start) & !is.na(end)]
dropped <- setdiff(d$gene_symbol, placed$gene_symbol)
if (length(dropped)) {
  cat(sprintf("dropped (no v5 coordinates): %s\n", paste(dropped, collapse = ", ")))
}

tier_short <- function(tier) {
  fifelse(
    tier == 1, "Tier1 brace-root regulator",
    fifelse(
      tier == 2, "Tier2 shoot-borne root initiation/emergence",
      fifelse(tier == 3, "Tier3 lateral-root / modifier", as.character(tier))
    )
  )
}

# first clause of the phenotype text, up to the first sentence stop
first_clause <- function(x) {
  x <- sub("\\.\\s.*$", "", x) # keep first sentence
  sub("\\.$", "", x) # trim trailing period
}

tsv <- data.table(
  symbol = placed$gene_symbol,
  gene_id = placed$v5_id,
  chr = as.integer(placed$chr),
  start = as.integer(placed$start),
  end = as.integer(placed$end),
  qtl_chen2019 = '""', # literal, matching the other candidate TSVs
  v5_canonical_symbol = placed$gene_symbol,
  pathway = paste0(
    tier_short(placed$tier), " — ", placed$phenotype_effect,
    "; ", first_clause(placed$phenotype)
  )
)
setorder(tsv, chr, start)

tsv_path <- here("data/teonam/nbr_candidate_genes.tsv")
fwrite(tsv, tsv_path, sep = "\t", quote = FALSE)
cat(sprintf("wrote %s (%d genes)\n", tsv_path, nrow(tsv)))

overlap <- tsv[, .(symbol, qtl = paste0("NBR(chr", chr, ")"), gene_id, chr, start, end)]
ov_path <- here("results/sim/zeal/nbr_candidate_overlap.csv")
fwrite(overlap, ov_path)
cat(sprintf("wrote %s (%d genes)\n", ov_path, nrow(overlap)))
