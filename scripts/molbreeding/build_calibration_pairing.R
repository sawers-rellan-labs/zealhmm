#!/usr/bin/env Rscript
# Build the calibration-context truth<->test sample pairing for the MolBreeding
# (target-seq, ~110x, TRUTH) vs SNP50K-skim (~0.4x, TEST) nnil calibration.
#
# The CANONICAL cross-source table
# (data/zeal/correspondence/molbreeding_3way_correspondence.csv) is label-based and
# assumes NO mislabels -- it is left UNTOUCHED. This derives a calibration-ONLY pairing
# from it and then applies documented sample-swap overrides established by genotype
# identity (donor Jaccard), valid only in this small-n calibration context.
#
# Documented override (zealtiger pn4_sid330_mislabel.qmd): the MolBreeding well labelled
# PN4_SID330 (well 4B6) is a label swap -- its DNA is PN4_SID322 (skim donor Jaccard 0.84,
# concordance 0.91; vs its own label 0.10, rank ~506/1393). So the TRUTH counts file
# PN4_SID330.tsv pairs with the skim TEST sample PN4_SID322 (true pedigree
# Zx.0030_P2_P2_P1.2.1.1). Skim PN4_SID330 (Zx.0030_P2_P4_P2.1.1.1) is a distinct,
# correctly-labelled individual with NO MolBreeding truth -> excluded from the calibration.
#
# Output: data/zeal/correspondence/calibration_pairing.csv
#   truth_sample (molb PN), test_sample (skim PN), true_pedigree, is_check, note
suppressMessages({
  library(here)
  library(data.table)
})
corr <- fread(here("data/zeal/correspondence/molbreeding_3way_correspondence.csv"))

# documented, genotype-verified overrides (calibration context ONLY)
overrides <- data.table(
  truth_sample = "PN4_SID330",
  test_sample = "PN4_SID322",
  true_pedigree = "Zx.0030_P2_P2_P1.2.1.1",
  note = "MolBreeding well 4B6 label swap: DNA is PN4_SID322 (donor Jaccard 0.84 vs skim-322, 0.10 vs own label; pn4_sid330_mislabel.qmd)"
)

pair <- corr[, .(
  truth_sample = molb_prefix,
  test_sample = skim_prefix, # canonical label pairing (identity in the id space)
  true_pedigree = molb_ped,
  is_check = (type == "B73"),
  note = ""
)]
# B73-bulk checks are MolBreeding-only (no skim NIL partner)
pair[is_check == TRUE, `:=`(
  test_sample = NA_character_,
  note = "B73-bulk control (MolBreeding-only, no skim partner)"
)]

# apply the documented mislabel overrides (canonical file stays untouched)
for (i in seq_len(nrow(overrides))) {
  j <- which(pair$truth_sample == overrides$truth_sample[i])
  if (length(j) != 1L) stop(sprintf("override truth_sample %s matched %d rows", overrides$truth_sample[i], length(j)))
  pair[j, `:=`(
    test_sample = overrides$test_sample[i],
    true_pedigree = overrides$true_pedigree[i],
    note = overrides$note[i]
  )]
}

# USE ALL libraries: no DNA-quality (Hannah Pil `quality`) gate is applied. Every NIL
# with a skim test partner is in the calibration, including the Jaccard-identified pair
# (truth PN4_SID330 <-> skim PN4_SID322). B73-bulk checks have no skim partner -> truth-only controls.
pair[, in_calibration := !is_check & !is.na(test_sample)]

setorder(pair, is_check, truth_sample)
fwrite(pair, here("data/zeal/correspondence/calibration_pairing.csv"))
cat(sprintf(
  "wrote calibration_pairing.csv: %d rows | %d NIL calibration pairs (in_calibration) | %d B73 controls | overrides applied: %d\n",
  nrow(pair), sum(pair$in_calibration), sum(pair$is_check), nrow(overrides)
))
print(pair[note != ""])
