#!/usr/bin/env Rscript
# =============================================================================
# Ear-number / prolificacy candidate genes: derive the caller TSV inputs from
# the single canonical annotation CSV.
#
# Single source of truth:
#   data/ref/maize_earnumber_prolificacy_cloned_genes.csv   (hand-curated, rich)
# Regenerated inputs (identical; two trait aliases used by the notebooks):
#   data/teonam/en_candidate_genes.tsv
#   data/teonam/prolif_candidate_genes.tsv
#
# The TSV schema is the shared candidate-gene contract used across the ZEAL
# notebooks (cf. dta/stpu/spad *_candidate_genes.tsv):
#   symbol  gene_id  chr  start  end  qtl_chen2019  v5_canonical_symbol  pathway
# `pathway` is derived deterministically from the CSV tier + protein + the first
# clause of trait_role, so editing the CSV is the only edit needed.
#
# Run: Rscript scripts/build_earnumber_candidate_tsv.R
# =============================================================================

suppressPackageStartupMessages(library(data.table))
library(here)

csv <- here("data/ref/maize_earnumber_prolificacy_cloned_genes.csv")
d <- fread(csv, encoding = "UTF-8")

# --- tier -> short label prefix used in the pathway string ------------------
tier_short <- function(tier) {
  fifelse(
    grepl("^1", tier), "Tier1 core",
    fifelse(grepl("^2", tier), "Tier2 infl. architecture", tier)
  )
}

# first clause of trait_role: text after the em-dash, up to the first sentence stop
first_clause <- function(x) {
  after <- trimws(sub("^[^—]*—\\s*", "", x)) # drop "<category> — "
  after <- sub("\\.\\s.*$", "", after) # keep first sentence
  sub("\\.$", "", after) # trim trailing period
}

out <- data.table(
  symbol = d$symbol,
  gene_id = d$v5_gene_model,
  chr = as.integer(sub("^chr", "", d$chr)),
  start = d$start_bp,
  end = d$end_bp,
  qtl_chen2019 = '""', # literal, matching the other candidate TSVs
  v5_canonical_symbol = d$symbol,
  pathway = paste0(
    tier_short(d$tier), " — ", d$protein,
    "; ", first_clause(d$trait_role)
  )
)
setorder(out, chr, start) # match the chr/position ordering of the inputs

for (trait in c("en", "prolif")) {
  path <- here(sprintf("data/teonam/%s_candidate_genes.tsv", trait))
  fwrite(out, path, sep = "\t", quote = FALSE)
  cat(sprintf("wrote %s (%d genes)\n", path, nrow(out)))
}
