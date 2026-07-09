#!/usr/bin/env Rscript
# ZEAL/BZea Phase 4 — build the JLM HapMap from the RTIGER mosaic.
# Analog of teonam_jlm_build.R. JLM (stepwise regression) is the ONE analysis that thins:
# markers >= 0.1 cM apart (exact 1-D greedy MDdIS per chr; Chen 2019 l.78, user-confirmed);
# GWAS/MLM stay full. Mosaic state 0/1/2 -> A/M/C HapMap, taxa = pedigree.
# Output: data/zeal/tassel/geno_jlm.hmp.txt
suppressMessages({
  library(here)
  library(data.table)
})
source(here("scripts/logging.R"))
THR <- 0.1

M0 <- readRDS(here("data/zeal/zeal_rtiger_mosaic.rds"))
state <- M0$state
mk <- M0$markers
inv <- tryCatch(readRDS(here("data/zeal/snp50k_invariant_markers.rds")), error = function(e) character(0))
cm <- fread(here("data/zeal/markers_snp50k_cm.tsv"))
mk[, cm := cm$cm[match(marker, cm$marker)]]
mk <- mk[!marker %in% inv & !is.na(cm)]

# 0.1-cM thin (JLM only): exact 1-D greedy per chr, sorted by cM
sel <- mk[order(chr, cm)][,
  {
    keep <- logical(.N)
    last <- -Inf
    for (i in seq_len(.N)) {
      if (cm[i] - last >= THR) {
        keep[i] <- TRUE
        last <- cm[i]
      }
    }
    .SD[keep]
  },
  by = chr
]
setorder(sel, chr, pos)
log_info("JLM thin: %d -> %d markers (>=%.1f cM apart)", nrow(mk), nrow(sel), THR)

G <- state[sel$marker, , drop = FALSE]
if (any(duplicated(colnames(G)))) G <- G[, !duplicated(colnames(G)), drop = FALSE]
cn <- matrix("N", nrow(G), ncol(G))
cn[round(G) == 0] <- "A"
cn[round(G) == 1] <- "M"
cn[round(G) == 2] <- "C" # B73/het/teo
hd <- data.table(
  `rs#` = sel$marker, alleles = "A/C", chrom = sel$chr, pos = sel$pos,
  strand = "+", `assembly#` = NA, center = NA, protLSID = NA, assayLSID = NA, panelLSID = NA, QCcode = NA
)
hmp <- cbind(hd, as.data.table(cn))
setnames(hmp, c(names(hd), colnames(G)))
dir.create(here("data/zeal/tassel"), showWarnings = FALSE, recursive = TRUE)
fwrite(hmp, here("data/zeal/tassel/geno_jlm.hmp.txt"), sep = "\t", quote = FALSE, na = "N")
log_info("wrote data/zeal/tassel/geno_jlm.hmp.txt: %d markers x %d taxa", nrow(hmp), ncol(G))
