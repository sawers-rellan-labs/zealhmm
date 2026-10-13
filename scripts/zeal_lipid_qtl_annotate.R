#!/usr/bin/env Rscript
# =============================================================================
# Annotate the lipid QTL summary with the v5 glycerolipid candidates (from
# scripts/zeal_lipid_candidates_v5.py) whose gene locus falls inside each QTL's
# LOD-support CI (same chromosome). Closes the loop: significant QTL -> candidate genes.
# In:  results/sim/zeal/rqtl/lipid/zeal_lipid_qtl_summary<sfx>.csv
#      data/zeal/lipid_candidates_v5.tsv
# Out: results/sim/zeal/rqtl/lipid/zeal_lipid_qtl_summary_annotated<sfx>.csv
# TAG env selects the scan: "" (default, log-scale main run) or e.g. "raw" -> "_raw" suffix.
# =============================================================================
suppressMessages({
  library(data.table)
  library(here)
})
TAG <- Sys.getenv("TAG", "")
sfx <- if (TAG %in% c("", "taxon")) "" else paste0("_", TAG)
q <- fread(here(sprintf("results/sim/zeal/rqtl/lipid/zeal_lipid_qtl_summary%s.csv", sfx)))
g <- fread(here("data/zeal/lipid_candidates_v5.tsv"))
g[, chr_i := as.integer(sub("^chr", "", chr))]
g[, mid_mb := (as.numeric(start) + as.numeric(end)) / 2e6]

ann <- function(ch, lo, hi) {
  h <- g[chr_i == ch & mid_mb >= lo & mid_mb <= hi]
  if (!nrow(h)) {
    return(list(0L, "", ""))
  }
  list(
    nrow(h), paste(h$gene_v5, collapse = ";"),
    paste(unique(unlist(strsplit(paste(h$pathway_names, collapse = ";"), ";"))), collapse = " | ")
  )
}
res <- q[,
  {
    a <- ann(chr, ci_left_mb, ci_right_mb)
    .(n_candidates = a[[1]], candidate_genes = a[[2]], candidate_pathways = a[[3]])
  },
  by = seq_len(nrow(q))
]
out <- cbind(q, res[, .(n_candidates, candidate_genes, candidate_pathways)])
fwrite(out, here(sprintf("results/sim/zeal/rqtl/lipid/zeal_lipid_qtl_summary_annotated%s.csv", sfx)))
cat(sprintf("QTL with >=1 candidate in CI: %d / %d\n", sum(out$n_candidates > 0), nrow(out)))
cat("=== QTL carrying candidate genes (trait, chr, peak Mb, LOD, n_cand, genes) ===\n")
print(
  out[n_candidates > 0][order(-lod), .(trait, class, chr,
    peak_mb = round(peak_mb, 1), lod = round(lod, 1),
    ci_mb = paste0(round(ci_left_mb, 1), "-", round(ci_right_mb, 1)), n_candidates, candidate_genes
  )],
  nrow = 60
)
