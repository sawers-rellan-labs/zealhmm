#!/usr/bin/env Rscript
# ZEAL/BZea Phase 1 — unified three-way sample sheet.
#
# Establishes, one row per line, the correspondence between:
#   pedigree string  <->  skim-seq library id  <->  BRB-seq library id
# plus the family factor (taxon), donor founder accession, external
# Sanchez/Holland collection code, field row, and BRB plate/well/passport.
#
# Anchor = the authored skim<->BRB crosswalk (zealtiger sample_correspondence.qmd);
# SNP50K-genotyped lines not present in that crosswalk are appended as skim-only rows
# so the sheet is a superset of the GWAS roster. See agent/zeal_dta_repro_plan.md (Phase 1).
#
# Output: data/zeal/samplesheet_3way.csv
#
# Cross-experiment key = the canonical pedigree (drop trailing ".B"); skim & MolBreeding
# share the skim id space, BRB-seq has its own (PN#_SID# is per-experiment, NOT a global key).

suppressMessages({
  library(here)
  library(data.table)
  library(readxl)
})
source(here("scripts/logging.R"))

# ---- source paths (override via env) ---------------------------------------
ZT <- Sys.getenv("ZEALTIGER", "/Users/fvrodriguez/Desktop/zealtiger")
P_CORR <- file.path(ZT, "results/sample_correspondence/skim_brbseq_correspondence.csv")
P_MASTER <- file.path(ZT, "data/sample_metadata_master.csv")
P_BRB <- file.path(ZT, "data/brbseq_metadata_master.csv")
P_SNP <- here("data/zeal/bzea_50K_cohort_ref_metadata.csv") # staged from Hazel
P_FB23 <- here("data/zeal/CLY23_D4_FieldBook.xlsx") # REF-all: accession -> external code
OUT <- here("data/zeal/samplesheet_3way.csv")

for (p in c(P_CORR, P_MASTER, P_BRB, P_SNP, P_FB23)) {
  if (!file.exists(p)) stop(sprintf("missing input: %s", p))
}

canon_ped <- function(x) sub("\\.B$", "", x) # drop the registry ".B" suffix
acc_of <- function(ped) sub("_.*$", "", ped) # founder accession = pedigree prefix
taxon_of <- function(acc) substr(acc, 1, 2) # taxon = accession prefix (Zx/Zv/Zd/Zl/Zh/Zm)

# ---- 1. anchor: skim<->BRB crosswalk ---------------------------------------
corr <- fread(P_CORR)
log_info("crosswalk: %d rows (%d skim, %d brbseq)", nrow(corr), sum(corr$in_skim), sum(corr$in_brbseq))
ss <- corr[, .(
  type, label,
  pedigree = canon_ped(pedigree),
  skim_id = skim_prefix, brbseq_id = brbseq_prefix,
  in_skim, in_brbseq
)]
ss[skim_id == "" | is.na(skim_id), skim_id := NA_character_]
ss[brbseq_id == "" | is.na(brbseq_id), brbseq_id := NA_character_]

# ---- 2. SNP50K genotyped roster (authoritative taxa/accession) --------------
snp <- fread(P_SNP)
snpc <- snp[is_cohort == TRUE, .(
  skim_id = sample, in_snp50k = TRUE,
  snp_taxon = maizegdb_prefix, snp_accession = donor_accession,
  snp_taxa_label = donor_taxa,
  is_B73 = is_B73, is_purple = is_purple, is_check = is_check,
  snp_pedigree = canon_ped(NIL_pedigree)
)]
log_info("SNP50K cohort: %d genotyped lines", nrow(snpc))

ss <- merge(ss, snpc, by = "skim_id", all.x = TRUE)
ss[is.na(in_snp50k), in_snp50k := FALSE]

# append SNP50K-genotyped lines absent from the crosswalk (skim-only) -----------
missing <- snpc[!skim_id %in% ss$skim_id]
if (nrow(missing)) {
  add <- missing[, .(
    type = "NIL", label = snp_pedigree, pedigree = snp_pedigree,
    skim_id, brbseq_id = NA_character_, in_skim = TRUE, in_brbseq = FALSE,
    in_snp50k = TRUE, snp_taxon, snp_accession, snp_taxa_label,
    is_B73, is_purple, is_check, snp_pedigree
  )]
  ss <- rbind(ss, add, fill = TRUE)
  log_warn("appended %d SNP50K-genotyped lines missing from the skim<->BRB crosswalk", nrow(add))
}

# ---- 3. field row + project (skim/MolB master registry) ---------------------
mst <- fread(P_MASTER)
mst <- mst[, .(
  skim_id = sample, field_row, project,
  mst_check = is_check, founder_group, donor_id, mst_pedigree = canon_ped(pedigree)
)]
ss <- merge(ss, mst, by = "skim_id", all.x = TRUE)

# ---- 4. BRB-seq plate/well + passport ---------------------------------------
brb <- fread(P_BRB)
setnames(brb, "teo-species", "teo_species", skip_absent = TRUE)
brbm <- brb[, .(
  brbseq_id = sample_id, brb_plate = plate, brb_well = plate_pos,
  brb_batch = batch, elevation, teo_species,
  brb_accession = accession, inv4m_genotype
)]
ss <- merge(ss, brbm, by = "brbseq_id", all.x = TRUE)

# ---- 5. external Sanchez/Holland accession code (REF-all bridge) ------------
ref <- as.data.table(read_excel(P_FB23, sheet = "REF-all"))
amap <- unique(ref[
  !is.na(accession_id) & !is.na(old_accession_id),
  .(accession_id, old_accession_id, taxa_code)
])
amap <- amap[, .SD[1], by = accession_id] # one external code per accession

# ---- 6. resolve family factor + accession (prefer SNP50K, else derive) ------
ss[, donor_accession := fifelse(!is.na(snp_accession), snp_accession, acc_of(pedigree))]
ss[, taxon := fifelse(!is.na(snp_taxon), snp_taxon, taxon_of(donor_accession))]
ss[, is_B73 := fcoalesce(is_B73, FALSE) | type %chin% "B73"]
ss[, is_purple := fcoalesce(is_purple, FALSE) | type %chin% "Purple"]
# a line is a check if any source says so, or it is a B73/Purple/non-NIL entry
ss[, is_check := fcoalesce(is_check, mst_check, FALSE) | is_B73 | is_purple |
  type %chin% c("B73", "Purple", "check", "Check")]
ss[, taxa_label := fcoalesce(snp_taxa_label, teo_species)]

# data-integrity check (before dropping the SNP50K columns): pedigree-derived
# taxon should agree with the authoritative SNP50K taxon on genotyped lines
n_disc <- ss[in_snp50k == TRUE & !is.na(snp_accession) &
  taxon_of(acc_of(pedigree)) != snp_taxon, .N]

ss <- merge(ss, amap, by.x = "donor_accession", by.y = "accession_id", all.x = TRUE)

# ---- 7. GWAS-panel definition -----------------------------------------------
# The 5 teosinte-donor taxa are the "Family" fixed factor (TeoNAM's 5 TIL analog).
# Zm = Zea mays (B73 / maize landrace), NOT a teosinte donor -> excluded.
# Zv (parviglumis) and Zx (mexicana) are Z. mays subspecies; Zd/Zl/Zh are distinct species.
TEO <- c("Zx", "Zv", "Zd", "Zl", "Zh")
ss[, donor_class := fifelse(taxon %chin% TEO, "teosinte", "maize")]
ss[, gwas_nil := in_snp50k & !is_check & project %chin% "bzea" &
  taxon %chin% TEO & !is.na(pedigree) & pedigree != ""]

# ---- 8. finalize ------------------------------------------------------------
keep <- c(
  "pedigree", "taxon", "donor_class", "donor_accession", "old_accession_id",
  "skim_id", "brbseq_id", "in_skim", "in_brbseq", "in_snp50k", "gwas_nil",
  "type", "label", "is_check", "is_B73", "is_purple",
  "field_row", "project", "brb_plate", "brb_well", "brb_batch",
  "elevation", "teo_species", "taxa_label", "inv4m_genotype", "taxa_code"
)
ss <- ss[, ..keep]
setorder(ss, -gwas_nil, taxon, donor_accession, pedigree)

fwrite(ss, OUT)
log_info("wrote %s: %d rows x %d cols", OUT, nrow(ss), ncol(ss))

# ---- QC report --------------------------------------------------------------
log_info("--- QC ---")
log_info(
  "in_skim=%d  in_brbseq=%d  in_snp50k=%d  all three=%d",
  sum(ss$in_skim, na.rm = TRUE), sum(ss$in_brbseq, na.rm = TRUE),
  sum(ss$in_snp50k, na.rm = TRUE),
  sum(ss$in_skim & ss$in_brbseq & ss$in_snp50k, na.rm = TRUE)
)
log_info("checks=%d (B73=%d purple=%d)", sum(ss$is_check), sum(ss$is_B73), sum(ss$is_purple))
log_info("SNP50K lines with field_row=%d / %d", sum(ss$in_snp50k & !is.na(ss$field_row)), sum(ss$in_snp50k))
log_info("SNP50K lines with external accession=%d / %d", sum(ss$in_snp50k & !is.na(ss$old_accession_id)), sum(ss$in_snp50k))
tx <- ss[gwas_nil == TRUE, .N, by = taxon][order(-N)]
log_info("GWAS NIL panel (gwas_nil==TRUE) = %d bzea teosinte-donor lines", sum(ss$gwas_nil))
log_info("family factor (taxon): %s", paste(sprintf("%s=%d", tx$taxon, tx$N), collapse = " "))
log_info(
  "excluded from panel: %d checks, %d maize (Zm), %d non-bzea project",
  sum(ss$is_check), sum(ss$in_snp50k & ss$donor_class == "maize" & !ss$is_check),
  sum(ss$in_snp50k & !ss$is_check & !(ss$project %chin% "bzea"))
)
if (n_disc) {
  log_warn("%d genotyped lines: pedigree-derived taxon != SNP50K taxon", n_disc)
} else {
  log_info("taxon integrity: pedigree-derived == SNP50K authoritative on all genotyped lines")
}
