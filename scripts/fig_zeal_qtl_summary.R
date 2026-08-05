#!/usr/bin/env Rscript
# scripts/fig_zeal_qtl_summary.R
#
# DRAFT. 4-panel linear whole-genome QTL summary figure for the ZEAL maize
# population paper (ZEAL = Zea Exotic Allele Library). Modeled on a linear
# per-chromosome QTL map (not a circos).
#
# Panels:
#   A  ten vertical chromosome ideograms (x = chr 1..10, y = position Mb v5);
#      one colored interval bar per QTL (ci_left..ci_right, Mb), bar linewidth
#      encodes LOD, color encodes trait; candidate-gene symbols labeled at the
#      overlapping peak. Source: ZEAL R/qtl taxon-covariate joint-scan QTL,
#      one confidence-interval peak per chromosome (ci_taxon table).
#   B  per-QTL percent variance explained (points), one column per trait,
#      colored by trait, annotated with the QTL count. Source: fitqtl drop-one.
#   C  two headline loci: phenotype BLUE (y) split by genotype class
#      A = B73 recurrent / H = het / B = teosinte donor (x), state palette.
#   D  per-trait total variance explained by the mapped QTL joint model.
#
# Conventions: figure text sans-serif; position in Mb (v5), never cM; panel
# tags bold capitals A,B,C,D; every series/axis labeled with dataset + caller +
# trait. No em/en dashes anywhere.

suppressMessages({
  library(qtl)
  library(ggplot2)
  library(patchwork)
  library(ggbeeswarm)
  library(ggrepel)
})

# ---- config -----------------------------------------------------------------
TRAITS <- c("dta", "dts", "ph", "prolif") # 4 TeoNAM-overlap traits
PHENO_COL <- c(dta = "DTA", dts = "DTS", ph = "PH", prolif = "Prolif")
TRAIT_LABEL <- c(
  dta = "DTA (days to anthesis)",
  dts = "DTS (days to silk)",
  ph = "PH (plant height)",
  prolif = "Prolif (ear prolificacy count)"
)

# colorblind-friendly trait palette (Okabe-Ito subset), distinct from the
# A/H/B genotype state palette used in Panel C.
TRAIT_PALETTE <- c(
  dta = "#0072B2", dts = "#56B4E9", ph = "#E69F00", prolif = "#009E73"
)

RQTL_DIR <- "results/sim/zeal/rqtl"
MARKERS_TSV <- "data/zeal/markers_snp50k_v5.tsv"
OUT_PNG <- "results/sim/zeal/fig_zeal_qtl_summary.png"

BASE_SIZE <- 8
BASE_FAMILY <- "sans"

# reuse the shared A/H/B genotype palette
source("R/plotting.R")
GENO_PAL <- state_palette() # REF=gold, HET=springgreen4, ALT=purple4

marker_bp <- function(m) as.numeric(sub("^S[0-9]+_", "", m))

# ---- load data --------------------------------------------------------------
# chromosome extents (max bp per chr) from the SNP50K v5 marker map
markers <- read.delim(MARKERS_TSV, stringsAsFactors = FALSE)
markers <- markers[markers$chr %in% 1:10, ]
chr_max <- tapply(markers$pos, markers$chr, max)
chr_ext <- data.frame(chr = as.integer(names(chr_max)), max_mb = as.numeric(chr_max) / 1e6)
chr_ext <- chr_ext[order(chr_ext$chr), ]

# peaks + effects + candidate overlaps, per trait
peaks_all <- list()
eff_all <- list()
cand_all <- list()
for (tr in TRAITS) {
  peaks_all[[tr]] <- read.csv(file.path(RQTL_DIR, sprintf("zeal_%s_peaks_ci_taxon.csv", tr)), stringsAsFactors = FALSE)
  eff_all[[tr]] <- read.csv(file.path(RQTL_DIR, sprintf("zeal_%s_qtl_effects.csv", tr)), stringsAsFactors = FALSE)
  cp <- file.path("results/sim/zeal", sprintf("%s_candidate_overlap.csv", tr))
  cand_all[[tr]] <- if (file.exists(cp)) read.csv(cp, stringsAsFactors = FALSE) else NULL
}

# tidy peaks with Mb coordinates and trait tag
peaks_df <- do.call(rbind, lapply(TRAITS, function(tr) {
  p <- peaks_all[[tr]]
  data.frame(
    trait = tr, chr = p$chr,
    peak_mb = marker_bp(p$marker) / 1e6,
    lo_mb = p$ci_left / 1e6, hi_mb = p$ci_right / 1e6,
    lod = p$lod, marker = p$marker, name = p$name,
    stringsAsFactors = FALSE
  )
}))

# effects tidy
eff_df <- do.call(rbind, eff_all)

# variance partition table (written by scripts/zeal_qtl_effects.R)
part_df <- read.csv(file.path(RQTL_DIR, "zeal_qtl_variance_partition.csv"), stringsAsFactors = FALSE)
part_df <- part_df[part_df$trait %in% TRAITS, ]

# ---- candidate-gene overlap: peak CI vs gene interval -----------------------
# for each peak, find candidate genes (same chr) whose [start,end] intersects
# the QTL support interval [ci_left, ci_right]; keep the closest symbol.
overlap_gene <- function(tr, chr, lo_bp, hi_bp) {
  cand <- cand_all[[tr]]
  if (is.null(cand)) {
    return(NA_character_)
  }
  cc <- cand[cand$chr == chr & cand$end >= lo_bp & cand$start <= hi_bp, ]
  if (nrow(cc) == 0) {
    return(NA_character_)
  }
  paste(cc$symbol, collapse = ", ")
}
peaks_df$symbol <- mapply(
  overlap_gene, peaks_df$trait, peaks_df$chr,
  peaks_df$lo_mb * 1e6, peaks_df$hi_mb * 1e6
)

# ---- shared theme -----------------------------------------------------------
theme_zeal <- function() {
  theme_bw(base_size = BASE_SIZE, base_family = BASE_FAMILY) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = BASE_SIZE + 1, face = "plain"),
      plot.tag = element_text(size = BASE_SIZE + 4, face = "bold"),
      legend.key.size = unit(0.35, "cm"),
      legend.title = element_text(size = BASE_SIZE),
      legend.text = element_text(size = BASE_SIZE - 1)
    )
}

trait_scale_col <- scale_color_manual(
  values = TRAIT_PALETTE, breaks = TRAITS,
  labels = TRAIT_LABEL[TRAITS], name = "Trait (ZEAL R/qtl joint scan)",
  guide = guide_legend(nrow = 2, override.aes = list(linewidth = 1.4, size = 1.6))
)
trait_scale_fill <- scale_fill_manual(
  values = TRAIT_PALETTE, breaks = TRAITS,
  labels = TRAIT_LABEL[TRAITS], name = "Trait (ZEAL R/qtl joint scan)"
)

# ==== Panel A: linear chromosome QTL map =====================================
# dodge multiple traits within a chromosome slot so overlapping intervals are
# visible: assign each trait a small horizontal offset around the chr integer.
n_tr <- length(TRAITS)
offsets <- setNames(seq(-0.28, 0.28, length.out = n_tr), TRAITS)
peaks_df$x <- peaks_df$chr + offsets[peaks_df$trait]

# candidate-gene label positions (only peaks that overlap a gene)
lab_df <- peaks_df[!is.na(peaks_df$symbol), ]

panelA <- ggplot() +
  # chromosome backbones
  geom_segment(
    data = chr_ext, aes(x = chr, xend = chr, y = 0, yend = max_mb),
    linewidth = 3.2, colour = "grey88", lineend = "round"
  ) +
  # QTL support-interval bars, linewidth ~ LOD, colour = trait
  geom_segment(
    data = peaks_df,
    aes(x = x, xend = x, y = lo_mb, yend = hi_mb, colour = trait, linewidth = lod),
    lineend = "round"
  ) +
  # peak points
  geom_point(
    data = peaks_df, aes(x = x, y = peak_mb, colour = trait),
    size = 0.5, show.legend = FALSE
  ) +
  # candidate-gene symbols at overlapping peaks
  geom_text_repel(
    data = lab_df, aes(x = x, y = peak_mb, label = symbol, colour = trait),
    size = 2, fontface = "italic", show.legend = FALSE,
    min.segment.length = 0, segment.size = 0.2, max.overlaps = 30,
    box.padding = 0.3
  ) +
  trait_scale_col +
  scale_linewidth_continuous(
    range = c(0.3, 2.2), breaks = c(5, 8, 11),
    name = "LOD", guide = guide_legend(nrow = 2)
  ) +
  scale_x_continuous(breaks = 1:10, expand = expansion(add = 0.5)) +
  scale_y_reverse(expand = expansion(mult = 0.02)) +
  labs(
    title = "ZEAL R/qtl joint-scan QTL, four TeoNAM-overlap traits",
    x = "Chromosome (maize v5)", y = "Position (Mb, v5)"
  ) +
  theme_zeal() +
  theme(panel.grid.major.x = element_blank())

# ==== Panel B: per-QTL percent variance explained ============================
eff_df$trait <- factor(eff_df$trait, levels = TRAITS)
nq <- table(factor(peaks_df$trait, levels = TRAITS))
xlab_b <- sprintf("%s\n(%d QTL)", TRAITS, as.integer(nq[TRAITS]))

panelB <- ggplot(eff_df, aes(x = trait, y = pct_var, colour = trait)) +
  geom_point(position = position_jitter(width = 0.12, seed = 1), size = 1.6, alpha = 0.9) +
  trait_scale_col +
  scale_x_discrete(labels = xlab_b) +
  guides(colour = "none") + # trait colour legend supplied by Panel A
  labs(
    title = "Per-QTL variance explained (ZEAL R/qtl fitqtl drop-one)",
    x = "Trait", y = "Per-QTL variance explained (%)"
  ) +
  theme_zeal() +
  theme(legend.position = "none")

# ==== Panel C: two headline loci, phenotype BLUE by genotype class ===========
# DTA headline: top-LOD DTA peak that overlaps a candidate gene.
# Prolif headline: prefer the ub3 locus (chr4) if it is among the ci peaks and
# overlaps its candidate; otherwise the top-LOD Prolif peak.
pick_peak <- function(tr, require_cand) {
  p <- peaks_df[peaks_df$trait == tr, ]
  if (require_cand) p <- p[!is.na(p$symbol), ]
  if (nrow(p) == 0) {
    return(NULL)
  }
  p[which.max(p$lod), ]
}
dta_head <- pick_peak("dta", TRUE)
if (is.null(dta_head)) dta_head <- pick_peak("dta", FALSE)

# Prolif: look for a chr4 peak whose overlapping candidate symbol includes ub3
prol <- peaks_df[peaks_df$trait == "prolif", ]
ub3_row <- prol[prol$chr == 4 & !is.na(prol$symbol) & grepl("\\bub3\\b", prol$symbol), ]
prol_head <- if (nrow(ub3_row) > 0) ub3_row[which.max(ub3_row$lod), ] else pick_peak("prolif", FALSE)

headline <- list(dta_head, prol_head)

geno_pheno_panel <- function(pk) {
  tr <- pk$trait
  cross <- readRDS(file.path(RQTL_DIR, sprintf("zeal_%s_cross.rds", tr)))
  g <- pull.geno(cross)[, pk$marker]
  y <- as.numeric(pull.pheno(cross, PHENO_COL[[tr]]))
  d <- data.frame(g = g, y = y)
  d <- d[!is.na(d$g) & !is.na(d$y), ]
  d$cls <- factor(d$g, levels = 1:3, labels = c("REF", "HET", "ALT"))
  ttl <- sprintf(
    "%s, chr%d @ %.1f Mb%s",
    toupper(tr), pk$chr, pk$peak_mb,
    if (!is.na(pk$symbol)) paste0(" (", pk$symbol, ")") else ""
  )
  ggplot(d, aes(x = cls, y = y, colour = cls)) +
    geom_quasirandom(size = 0.35, alpha = 0.5, width = 0.32) +
    geom_boxplot(outlier.shape = NA, fill = NA, colour = "grey25", linewidth = 0.3, width = 0.45) +
    scale_colour_manual(values = GENO_PAL, guide = "none") +
    scale_x_discrete(labels = c(REF = "A\nB73", HET = "H\nhet", ALT = "B\nteosinte")) +
    labs(title = ttl, x = "Genotype at peak marker", y = paste0(TRAIT_LABEL[[tr]], "\nBLUE")) +
    theme_zeal()
}
panelC1 <- geno_pheno_panel(headline[[1]])
panelC2 <- geno_pheno_panel(headline[[2]])
panelC <- panelC1 + panelC2 +
  plot_annotation(title = "Allele effect at headline loci (ZEAL, genotype A/H/B = B73/het/teosinte)") &
  theme(plot.title = element_text(size = BASE_SIZE, face = "plain"))

# ==== Panel D: total QTL-model variance explained per trait ==================
part_df$trait <- factor(part_df$trait, levels = TRAITS)
panelD <- ggplot(part_df, aes(x = trait, y = pct_var_qtl_model, fill = trait)) +
  geom_col(width = 0.68, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%", pct_var_qtl_model)),
    vjust = -0.4, size = 2.4
  ) +
  trait_scale_fill +
  scale_x_discrete(labels = TRAITS) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Variance explained by mapped QTL joint model (ZEAL R/qtl)",
    x = "Trait", y = "Full-model variance explained (%)"
  ) +
  theme_zeal() +
  theme(plot.title = element_text(size = BASE_SIZE))

# ==== compose ================================================================
# tag the four logical panels A,B,C,D. Panel C is a two-subpanel composite;
# wrap it so it carries a single "C" tag.
panelA <- panelA + labs(tag = "A")
panelB <- panelB + labs(tag = "B")
panelD <- panelD + labs(tag = "D")
panelC_wrapped <- wrap_elements(panelC) + labs(tag = "C")

layout <- (panelA) /
  (panelB | panelD) /
  (panelC_wrapped) +
  plot_layout(heights = c(1.5, 1, 1.1), guides = "collect") &
  theme(legend.position = "bottom")

ggsave(OUT_PNG, layout, width = 180, height = 235, units = "mm", dpi = 320)
cat("wrote", OUT_PNG, "\n")
cat("Panel C headline loci:\n")
for (h in headline) {
  cat(sprintf(
    "  %s chr%d @ %.1f Mb LOD %.2f symbol=%s\n",
    toupper(h$trait), h$chr, h$peak_mb, h$lod,
    if (is.na(h$symbol)) "none" else h$symbol
  ))
}
