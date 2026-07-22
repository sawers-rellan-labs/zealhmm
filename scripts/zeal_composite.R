#!/usr/bin/env Rscript
# ZEAL Phase 4 — DTA composite (A: JLM QTL lollipop · B: MLM Taxon+K Manhattan).
# Analog of the TeoNAM DTA composite (teonam-qtl-recovery-dta-mlm-118k.qmd). ZEAL has one
# ancestry mosaic (rtiger), so this is a 2-panel A/B (TeoNAM's Panel C was a 2nd caller).
# JLM (Panel A) is the model-independent DTA analysis (recovers zmcct10, chr10); MLM
# (Panel B) shows the ancestry-confound suppression. TRAIT via env (default DTA).
suppressMessages({
  library(here)
  library(data.table)
  library(ggplot2)
  library(fastman)
  library(ggrepel)
  library(ggtext)
  library(scales)
  library(cowplot)
})
source(here("scripts/teonam_notebook_plots.R"))
source(here("scripts/logging.R"))
TRAIT <- toupper(Sys.getenv("TRAIT", "DTA"))
TTAG <- tolower(TRAIT)
FAMCOL <- Sys.getenv("FAMILY_COL", "donor_accession") # best-calibrated MLM (lambda 1.19)
OUT <- here("results/sim/zeal")
overlap <- here("results/sim/zeal/dta_candidate_overlap.csv")

lam <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  round(median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1), 2)
}
gt <- fread(here(sprintf("data/zeal/%s_gwas_mlm_hwe_post_gt_%s_snp50k.csv", TTAG, FAMCOL))) # Panel B
rtiger <- fread(here(sprintf("data/zeal/%s_gwas_mlm_rtiger_mosaic_%s_snp50k.csv", TTAG, FAMCOL))) # Panel C

# A: JLM lollipop (model-independent QTL; context scan = the rtiger_mosaic MLM)
p_lolli <- plot_lollipop(
  here(sprintf("results/sim/zeal/tassel/%s_jlm_native1.txt", TTAG)), rtiger,
  sprintf("ZEAL %s — JLM QTL genomic distribution (nested in taxon)", TRAIT),
  overlap_csv = overlap, label_genes = c("zcn8", "zcn12", "zmcct9", "zmcct10", "dlf1", "tu1"),
  out_png = file.path(OUT, sprintf("%s_jlm_lollipop_snp50k.png", TTAG))
)
# B: per-SNP genotype (bcftools HWE-posterior) MLM (sharp SNP peaks; sparse -> lambda-deflated at 0.4x)
p_B <- plot_manhattan(gt,
  sprintf("ZEAL %s — MLM (Taxon+K), per-SNP genotype (bcftools HWE-posterior) (lambda_GC = %.2f)", TRAIT, lam(gt$P)),
  overlap_csv = overlap, out_png = file.path(OUT, sprintf("%s_gwas_mlm_hwe_post_gt_manhattan_snp50k.png", TTAG))
)
# C: rtiger ancestry mosaic MLM (block-smoothed; well-calibrated)
p_C <- plot_manhattan(rtiger,
  sprintf("ZEAL %s — MLM (Taxon+K), rtiger ancestry mosaic (lambda_GC = %.2f)", TRAIT, lam(rtiger$P)),
  overlap_csv = overlap, out_png = file.path(OUT, sprintf("%s_gwas_mlm_rtiger_mosaic_manhattan_snp50k.png", TTAG))
)

composite <- cowplot::plot_grid(p_lolli, p_B, p_C,
  ncol = 1, labels = c("A", "B", "C"),
  label_size = 20, align = "v", axis = "lr"
)
ggsave(file.path(OUT, sprintf("%s_composite_snp50k.png", TTAG)), composite,
  width = 9, height = 13, dpi = 200, bg = "white"
)
log_info("wrote %s/%s_composite_snp50k.png (A JLM · B per-SNP · C mosaic)", OUT, TTAG)
