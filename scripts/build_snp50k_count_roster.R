#!/usr/bin/env Rscript
# Build data/zeal/snp50k_count_roster.tsv: one row per unique physical SNP50K skim sample that any
# figure needs, with the canonical `sample` id and its role flags. This is the artifact the whole
# count-tree consolidation keys on (one canonical GATK-counts-with-real-alleles store per `sample`,
# each figure selecting its subset). See [[50k-set-terminology]].
#
# Columns:
#   sample       canonical id: PEDIGREE for NILs; synthetic B73_<PN> / Purple_<PN> for checks
#                (checks are the B73/Purple recurrent-parent controls, pedigree=NA)
#   skim_prefix  the PN#_SID# id GATK's allelic_counts50K.tsv is keyed on (the extractor's input key)
#   is_check     TRUE for B73/Purple controls
#   in_dnarna    DNA-skim + RNA-BRBseq base cohort   -> data/zeal/paired_cohort/skim/counts_50k
#   in_molb      14-NIL MolBreeding comparison cohort -> data/zeal/molb_cohort/skim/counts_50k
#   in_dnasweep  DNA coverage sweep (Fig 3 selection) -> data/skimsweep/skim/counts_50k
#
# Master PN<->pedigree map: data/zeal/correspondence/skim_brbseq_correspondence.csv (type, pedigree,
# skim_prefix, ...). Role flags come from the CURRENT tree memberships (what each figure actually
# reads today), so the roster reproduces exactly the current sample sets, just unified.
#
#   Rscript scripts/build_snp50k_count_roster.R

suppressMessages(library(data.table))
library(here)

corr <- fread(here("data/zeal/correspondence/skim_brbseq_correspondence.csv"))
stopifnot(all(c("type", "pedigree", "skim_prefix") %in% names(corr)))

# canonical id + check flag
corr[, is_check := type %in% c("B73", "Purple")]
corr[, sample := fifelse(is_check, paste(type, skim_prefix, sep = "_"), pedigree)]

# authoritative role memberships = which count directory each sample appears in
base_no_ext <- function(dir) sub("\\.tsv$", "", list.files(here(dir), pattern = "\\.tsv$"))
dnasweep_dir <- base_no_ext("data/skimsweep/skim/counts_50k") # PN-named  (coverage sweep)
dnarna_dir <- base_no_ext("data/zeal/paired_cohort/skim/counts_50k") # ped-named (skim+BrB base cohort)
molb_dir <- base_no_ext("data/zeal/molb_cohort/skim/counts_50k") # ped-named (MolBreeding cohort)

corr[, in_dnasweep := skim_prefix %in% dnasweep_dir]
corr[, in_dnarna := !is.na(pedigree) & pedigree %in% dnarna_dir]
corr[, in_molb := !is.na(pedigree) & pedigree %in% molb_dir]

roster <- corr[
  in_dnarna | in_molb | in_dnasweep,
  .(sample, skim_prefix, is_check, in_dnarna, in_molb, in_dnasweep)
]

# integrity checks: no dup PN, no dup canonical id, and every tree member accounted for
dup_pn <- roster[duplicated(skim_prefix), unique(skim_prefix)]
dup_id <- roster[duplicated(sample), unique(sample)]
if (length(dup_pn)) cat("WARN duplicate skim_prefix:", paste(dup_pn, collapse = ", "), "\n")
if (length(dup_id)) cat("WARN duplicate sample id:", paste(dup_id, collapse = ", "), "\n")
miss_dnasweep <- setdiff(dnasweep_dir, roster$skim_prefix)
miss_dnarna <- setdiff(dnarna_dir, roster$sample)
miss_molb <- setdiff(molb_dir, roster$sample)
if (length(miss_dnasweep)) cat("WARN dnasweep PN not mapped:", paste(miss_dnasweep, collapse = ", "), "\n")
if (length(miss_dnarna)) cat("WARN dnarna pedigree not in roster:", paste(miss_dnarna, collapse = ", "), "\n")
if (length(miss_molb)) cat("WARN molb pedigree not in roster:", paste(miss_molb, collapse = ", "), "\n")

setorder(roster, -in_dnarna, -in_molb, -in_dnasweep, sample)
fwrite(roster, here("data/zeal/snp50k_count_roster.tsv"), sep = "\t")

cat(sprintf("\nwrote data/zeal/snp50k_count_roster.tsv: %d samples\n", nrow(roster)))
cat(sprintf("  NILs: %d ; checks: %d\n", roster[!is_check == TRUE, .N], roster[is_check == TRUE, .N]))
cat(sprintf(
  "  in_dnarna: %d ; in_molb: %d ; in_dnasweep: %d\n",
  roster[in_dnarna == TRUE, .N], roster[in_molb == TRUE, .N], roster[in_dnasweep == TRUE, .N]
))
cat(sprintf("  molb-only (in_molb & !in_dnarna): %d\n", roster[in_molb == TRUE & in_dnarna == FALSE, .N]))
cat(sprintf("  dnasweep that are checks: %d\n", roster[in_dnasweep == TRUE & is_check == TRUE, .N]))
