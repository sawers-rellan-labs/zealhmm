#!/usr/bin/env Rscript
# =============================================================================
# Regenerate the TeoNAM JLM marker set + hapmap DIRECTLY from the 51K GWAS pool
# (teonam_map_v5_gwas, 51,004 genotyped markers), replacing the old path that ran
# FastIndep on the deduped 47,750 grid (teonam_map_v5_gwas_nr). This is the
# traceable rebuild of the JLM-generating script that was lost.
#
# Verified faithful reconstruction (agent/teonam_jlm_verify_source.R):
#   markers  = select_independent(cM-distance, threshold=0.1, n_runs=1) per chr
#   genotype = per-family W22TIL0x files step-interpolated onto the selected cM
#             (mode="step") -> complete A/M/C  [reproduced old hapmap: 0 mismatch]
#
# Change here: FastIndep POOL = 51,004 (all genotyped markers, duplicate-cM kept)
# instead of the 47,750 dedup. cM-distance @0.1 subsumes the dedup (tied markers
# are at distance 0 < 0.1 -> one survivor per cluster), so no pre-dedup needed.
# Run: Rscript agent/teonam_jlm_build.R
# =============================================================================
suppressMessages({
  library(data.table)
  devtools::load_all("/Users/fvrodriguez/repos/nilhmm", quiet = TRUE)
})
setwd("/Users/fvrodriguez/repos/zealhmm")
THR <- 0.1
HMP <- "data/teonam/tassel/geno.hmp.txt"
HMP_OLD <- "data/teonam/tassel/geno_gwas_nr.hmp.txt.bak" # backup of the 47,750-pool JLM

mc <- fread("data/teonam/marker_info_v5_cm.tsv")
setnames(mc, "chr_v5", "chr")
cm_by <- setNames(mc$cm, mc$marker)
pos_by <- setNames(mc$pos_v5, mc$marker)

# --- 51,004 GWAS pool = annotation(with cM) ∩ genotype columns --------------
gcols <- names(fread("data/teonam/TeoNAM_genotype_clean.csv", nrows = 0))[-(1:3)]
pool <- mc[marker %in% gcols & is.finite(cm)][order(chr, cm)]
cat(sprintf("FastIndep pool (teonam_map_v5_gwas): %d markers, %d chr\n", nrow(pool), uniqueN(pool$chr)))

# --- FastIndep (deterministic greedy) per chromosome, cM-distance @0.1 -------
run_chr <- function(m) { # m: data.table(marker, cm), sorted
  m <- m[order(cm)]
  D <- abs(outer(m$cm, m$cm, `-`))
  dimnames(D) <- list(m$marker, m$marker)
  as.character(select_independent(D,
    threshold = THR, n_runs = 1L, sense = "distance",
    max_markers = nrow(D) + 1L
  ))
}
sel <- unlist(lapply(sort(unique(pool$chr)), function(ch) run_chr(pool[chr == ch, .(marker, cm)])))
cat(sprintf("selected JLM markers (51K pool @%.2f cM): %d\n", THR, length(sel)))

# --- compare to the old 47,750-pool JLM set ----------------------------------
if (file.exists(HMP) && !file.exists(HMP_OLD)) {
  old <- fread(HMP, select = "rs#")[["rs#"]]
  cat(sprintf(
    "old JLM (47,750-pool): %d | in both: %d | new-only: %d | old-only: %d | Jaccard: %.4f\n",
    length(old), length(intersect(sel, old)), length(setdiff(sel, old)),
    length(setdiff(old, sel)), length(intersect(sel, old)) / length(union(sel, old))
  ))
}

# --- genotypes: per-family step-interpolation onto the selected cM -----------
fams <- c(
  TIL01 = "W22TIL01_genotype.csv", TIL03 = "W22TIL03_genotype.csv",
  TIL11 = "W22TIL11_genotype.csv", TIL14 = "W22TIL14_genotype.csv", TIL25 = "W22TIL25_genotype.csv"
)
fam_data <- lapply(names(fams), function(fam) {
  g <- fread(file.path("data/teonam", fams[fam]))
  g <- g[!duplicated(g[[1]])]
  list(g = g, keys = paste0(fam, sub("^.*Line_", "", g[[1]])))
})
names(fam_data) <- names(fams)

seldt <- data.table(marker = sel, chr = mc[match(sel, marker), chr], cm = cm_by[sel])[order(chr, cm)]
G <- do.call(rbind, lapply(sort(unique(seldt$chr)), function(ch) {
  tgt <- seldt[chr == ch]
  tgt_df <- data.frame(chr = ch, cm = tgt$cm)
  blocks <- lapply(names(fam_data), function(fam) {
    g <- fam_data[[fam]]$g
    obsmk <- intersect(names(g)[-(1:3)], mc[chr == ch, marker])
    obs <- data.frame(marker = obsmk, cm = cm_by[obsmk])
    obs <- obs[order(obs$cm), ]
    obs <- obs[!duplicated(obs$cm), ]
    geno <- t(as.matrix(g[, obs$marker, with = FALSE]))
    storage.mode(geno) <- "double"
    dn <- interpolate_genotype(geno, data.frame(chr = ch, cm = obs$cm), tgt_df, mode = "step")
    colnames(dn) <- fam_data[[fam]]$keys
    dn
  })
  D <- do.call(cbind, blocks)
  rownames(D) <- tgt$marker
  D
}))
cat(sprintf("interpolated genotype matrix: %d markers x %d taxa; NA=%d\n", nrow(G), ncol(G), sum(is.na(G))))

# --- recode 0/1/2 -> A/M/C and write hapmap (sorted chrom,pos) ---------------
cn <- matrix("N", nrow(G), ncol(G))
cn[round(G) == 0] <- "A"
cn[round(G) == 1] <- "M"
cn[round(G) == 2] <- "C"
hd <- data.table(
  `rs#` = rownames(G), alleles = "A/C", chrom = as.integer(mc[match(rownames(G), marker), chr]),
  pos = as.integer(pos_by[rownames(G)]), strand = "+", `assembly#` = NA, center = NA,
  protLSID = NA, assayLSID = NA, panelLSID = NA, QCcode = NA
)
hmp <- cbind(hd, as.data.table(cn))
setnames(hmp, c(names(hd), colnames(G)))
setorder(hmp, chrom, pos)

if (file.exists(HMP) && !file.exists(HMP_OLD)) {
  file.copy(HMP, HMP_OLD)
  cat(sprintf("backed up old hapmap -> %s\n", basename(HMP_OLD)))
}
fwrite(hmp, HMP, sep = "\t", quote = FALSE, na = "N")
gaps <- seldt[, diff(cm), by = chr]$V1
cat(sprintf(
  "wrote %s: %d markers x %d taxa | min adjacent cM gap: %.4f (>= %.2f = %s)\n",
  basename(HMP), nrow(hmp), ncol(G), min(gaps), THR, min(gaps) >= THR - 1e-9
))
