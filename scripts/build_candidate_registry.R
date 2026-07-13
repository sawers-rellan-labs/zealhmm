#!/usr/bin/env Rscript
# =============================================================================
# Normalize the fragmented candidate-gene tables into ONE registry (3NF).
#
# Sources (heterogeneous schemas), unioned + cross-filled by v5 gene_id:
#   data/ref/brace_root_genes_v5.csv            (NBR; richest: doi/phenotype/effect/v3v4)
#   data/ref/maize_earnumber_prolificacy_cloned_genes.csv (EN + Prolif; key_reference)
#   data/ref/leaf_number_candidate_genes_v5.csv (LAE; mechanism/genbank/entrez)
#   data/ref/zigzag_gwas_loci_v5.csv            (Kinki; Relevance)
#   data/teonam/*_candidate_genes.tsv           (all 14 traits; the derived contract)
#
# Output (tracked in reference/ — the single hand-curated source of truth going forward;
# one-time consolidator, run once from the legacy sources, then edit the registry directly):
#   reference/candidate_genes.csv     one row per distinct v5 gene (gene-intrinsic)
#   reference/candidate_evidence.csv  one row per gene x trait (the causal claim + evidence)
#
# The evidence table carries the DOI / evidence_type / causation_direction the causal
# claim comes from. Populated from the sources that already have it (brace_root fully;
# earnumber via key_reference; leaf_number via mechanism); everything else is left as
# an explicit TODO for curation. NO fabrication.
# Run: Rscript scripts/build_candidate_registry.R
# =============================================================================
suppressMessages({
  library(data.table)
  library(here)
})

norm_chr <- function(x) as.integer(sub("^chr", "", as.character(x)))
first_doi <- function(x) {
  m <- regmatches(x, regexpr("10\\.[0-9]{4,}/[^ ,;\"]+", x))
  if (length(m)) m[1] else NA_character_
}
blank <- function(n) rep(NA_character_, n)

# ---- 1. gene-intrinsic rows from each source (canonical columns) -------------
recs <- list()

# brace_root (NBR) — the richest schema; the causal model template
br <- fread(here("data/ref/brace_root_genes_v5.csv"), encoding = "UTF-8")
recs$nbr <- data.table(
  gene_id = br$v5_id, symbol = br$gene_symbol, full_name = br$full_name,
  entrez = NA_character_, genbank = NA_character_, v4_id = br$v4_id, v3_id = br$v3_id,
  chr = norm_chr(br$chr), start = as.integer(br$start), end = as.integer(br$end),
  strand = as.character(br$strand), length_bp = as.integer(br$length_bp), biotype = br$biotype,
  protein = NA_character_,
  trait = "nbr", tier = as.character(br$tier), qtl = NA_character_,
  functional_note = br$phenotype, causation_direction = br$phenotype_effect,
  evidence_type = NA_character_, reference = br$reference, doi = br$doi, pmid = NA_character_
)

# earnumber / prolificacy (feeds EN and Prolif)
en <- fread(here("data/ref/maize_earnumber_prolificacy_cloned_genes.csv"), encoding = "UTF-8")
en_base <- data.table(
  gene_id = en$v5_gene_model, symbol = en$symbol, full_name = en$gene_name,
  entrez = NA_character_, genbank = NA_character_, v4_id = NA_character_, v3_id = NA_character_,
  chr = norm_chr(en$chr), start = as.integer(en$start_bp), end = as.integer(en$end_bp),
  strand = as.character(en$strand), length_bp = NA_integer_, biotype = NA_character_,
  protein = en$protein, tier = as.character(en$tier), qtl = NA_character_,
  functional_note = en$trait_role, causation_direction = NA_character_,
  evidence_type = NA_character_, reference = en$key_reference,
  doi = vapply(en$key_reference, first_doi, ""), pmid = NA_character_
)
recs$en <- copy(en_base)[, trait := "en"]
recs$prolif <- copy(en_base)[, trait := "prolif"]

# leaf_number (LAE) — v7 schema
lae <- fread(here("data/ref/leaf_number_candidate_genes_v5.csv"), encoding = "UTF-8")
recs$lae <- data.table(
  gene_id = lae$v5_gene_model, symbol = lae$gene, full_name = lae$label,
  entrez = as.character(lae$entrez_or_model_id), genbank = as.character(lae$genbank_accession),
  v4_id = NA_character_, v3_id = NA_character_,
  chr = norm_chr(lae$chr), start = as.integer(lae$start_v5), end = as.integer(lae$end_v5),
  strand = as.character(lae$strand), length_bp = as.integer(lae$length_bp), biotype = NA_character_,
  protein = NA_character_, trait = "lae", tier = as.character(lae$tier), qtl = NA_character_,
  functional_note = lae$mechanism, causation_direction = NA_character_,
  evidence_type = NA_character_, reference = NA_character_, doi = NA_character_, pmid = NA_character_
)

# zigzag (Kinki)
zz <- fread(here("data/ref/zigzag_gwas_loci_v5.csv"), encoding = "UTF-8")
recs$kinki <- data.table(
  gene_id = zz$GeneID_v5, symbol = sub(" .*", "", zz$Locus), full_name = zz$Locus,
  entrez = NA_character_, genbank = NA_character_, v4_id = NA_character_, v3_id = NA_character_,
  chr = norm_chr(zz$Chr), start = as.integer(zz$Start_bp), end = as.integer(zz$End_bp),
  strand = NA_character_, length_bp = NA_integer_, biotype = NA_character_, protein = NA_character_,
  trait = "kinki", tier = NA_character_, qtl = NA_character_,
  functional_note = zz$Relevance, causation_direction = NA_character_,
  evidence_type = NA_character_, reference = NA_character_, doi = NA_character_, pmid = NA_character_
)

# ---- 2. the 14 derived contract TSVs (defines trait membership; thin) --------
tsvs <- list.files(here("data/teonam"), pattern = "_candidate_genes\\.tsv$", full.names = TRUE)
for (f in tsvs) {
  tr <- sub("_candidate_genes.tsv$", "", basename(f))
  d <- fread(f, encoding = "UTF-8")
  recs[[paste0("tsv_", tr)]] <- data.table(
    gene_id = d$gene_id, symbol = d$symbol, full_name = NA_character_,
    entrez = NA_character_, genbank = NA_character_, v4_id = NA_character_, v3_id = NA_character_,
    chr = norm_chr(d$chr), start = as.integer(d$start), end = as.integer(d$end),
    strand = NA_character_, length_bp = NA_integer_, biotype = NA_character_, protein = NA_character_,
    trait = tr, tier = NA_character_,
    qtl = if ("qtl_chen2019" %in% names(d)) as.character(d$qtl_chen2019) else NA_character_,
    functional_note = if ("pathway" %in% names(d)) d$pathway else NA_character_,
    causation_direction = NA_character_, evidence_type = NA_character_,
    reference = NA_character_, doi = NA_character_, pmid = NA_character_
  )
}

reg <- rbindlist(recs, use.names = TRUE, fill = TRUE)
reg[, qtl := ifelse(is.na(qtl) | qtl %in% c("", '""'), NA_character_, qtl)]
# Uncloned loci (e.g. mhl1/inv9f for StPu) have no v5 gene model — key them by
# locus:<symbol> so they survive the registry instead of being dropped.
reg[is.na(gene_id) | gene_id == "", gene_id := paste0("locus:", symbol)]
reg <- reg[!is.na(gene_id) & gene_id != ""]

# ---- 3. candidate_genes: one row per gene, cross-filled (first non-NA) -------
coalesce_by <- function(x) {
  ok <- !is.na(x)
  if (is.character(x)) ok <- ok & x != ""
  if (any(ok)) x[which(ok)[1]] else x[1] # x[1] keeps the column's type when all-NA
}
GCOLS <- c("symbol", "full_name", "entrez", "genbank", "v4_id", "v3_id", "chr", "start", "end", "strand", "length_bp", "biotype", "protein")
genes <- reg[, lapply(.SD, coalesce_by), by = gene_id, .SDcols = GCOLS]
setorder(genes, chr, start)
fwrite(genes, here("reference/candidate_genes.csv"))

# ---- 4. candidate_evidence: one row per gene x trait -------------------------
ECOLS <- c("gene_id", "trait", "tier", "qtl", "functional_note", "phenotype", "causation_direction", "evidence_type", "reference", "doi", "pmid", "curated")
ev <- unique(reg[, .(gene_id, trait, tier, qtl, functional_note, causation_direction, evidence_type, reference, doi, pmid)], by = c("gene_id", "trait"))
ev[, phenotype := NA_character_] # reserved: the explicit causal claim (Tier-1 curation)
ev[, curated := fifelse(!is.na(doi) & doi != "", "source", "TODO")]
setcolorder(ev, ECOLS)
setorder(ev, trait, gene_id)
fwrite(ev, here("reference/candidate_evidence.csv"))

cat(sprintf(
  "candidate_genes.csv: %d genes | candidate_evidence.csv: %d gene x trait rows\n  with DOI: %d/%d  | TODO (no evidence source): %d\n",
  nrow(genes), nrow(ev), sum(!is.na(ev$doi) & ev$doi != ""), nrow(ev), sum(ev$curated == "TODO")
))
