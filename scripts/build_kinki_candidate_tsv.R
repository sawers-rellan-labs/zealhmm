#!/usr/bin/env Rscript
# =============================================================================
# Kinki (zigzag-culm / kinked-stem severity) candidate loci: derive the caller
# TSV input + the notebook overlap CSV from the single canonical annotation CSV.
#
# Single source of truth:
#   data/ref/zigzag_gwas_loci_v5.csv   (hand-curated; Locus, GeneID_v5, coords, Relevance)
# Regenerated:
#   data/teonam/kinki_candidate_genes.tsv        (caller/GWAS candidate contract)
#   results/sim/zeal/kinki_candidate_overlap.csv (Manhattan/lollipop gene overlay)
#
# TSV schema (shared candidate-gene contract, cf. stpi/eh *_candidate_genes.tsv):
#   symbol  gene_id  chr  start  end  qtl_chen2019  v5_canonical_symbol  pathway
# Overlap schema (cf. *_candidate_overlap.csv):
#   symbol  qtl  gene_id  chr  start  end
#
# `symbol` is the short gene tag parsed from the free-text Locus column; `pathway`
# is the Relevance text (a working note "— YOUR comparison gene" is stripped for
# display). Editing the CSV is the only edit needed.
#
# Run: Rscript scripts/build_kinki_candidate_tsv.R
# =============================================================================

suppressPackageStartupMessages(library(data.table))
library(here)

csv <- here("data/ref/zigzag_gwas_loci_v5.csv")
d <- fread(csv, encoding = "UTF-8")

# short symbol from Locus: text before "(", then take the gene after any
# "... pathway: X" prefix, and the first of any "a / b" alias pair.
short_symbol <- function(locus) {
  lead <- trimws(sub("\\s*\\(.*$", "", locus)) # "ct2", "kn1 pathway: rs2", "gn1 / knox4"
  lead <- trimws(sub("^.*:\\s*", "", lead)) # "kn1 pathway: rs2" -> "rs2"
  trimws(sub("\\s*/.*$", "", lead)) # "gn1 / knox4" -> "gn1"
}

# Relevance -> pathway; drop the informal working annotation
clean_pathway <- function(x) trimws(sub("\\s*[—-]+\\s*YOUR comparison gene", "", x))

tsv <- data.table(
  symbol              = vapply(d$Locus, short_symbol, character(1)),
  gene_id             = d$GeneID_v5,
  chr                 = as.integer(d$Chr),
  start               = as.integer(d$Start_bp),
  end                 = as.integer(d$End_bp),
  qtl_chen2019        = '""', # literal, matching the other candidate TSVs
  v5_canonical_symbol = vapply(d$Locus, short_symbol, character(1)),
  pathway             = vapply(d$Relevance, clean_pathway, character(1))
)
setorder(tsv, chr, start)

tsv_path <- here("data/teonam/kinki_candidate_genes.tsv")
fwrite(tsv, tsv_path, sep = "\t", quote = FALSE)
cat(sprintf("wrote %s (%d loci)\n", tsv_path, nrow(tsv)))

overlap <- tsv[, .(symbol, qtl = paste0("KINKI(chr", chr, ")"), gene_id, chr, start, end)]
ov_path <- here("results/sim/zeal/kinki_candidate_overlap.csv")
fwrite(overlap, ov_path)
cat(sprintf("wrote %s (%d loci)\n", ov_path, nrow(overlap)))
