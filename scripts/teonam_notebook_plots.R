# Shared plotting helpers for the TeoNAM QTL-recovery notebooks (OLS and MLM, STAM and
# DTA). Sourced by analysis/teonam-qtl-recovery*-118k.qmd so the function bodies live in
# one place instead of being pasted into every notebook's `helpers` chunk.
#   source(here::here("scripts/teonam_notebook_plots.R"))
# Trait-specific inputs (candidate-overlap path, lollipop gene labels) are ARGUMENTS,
# so the same functions serve STAM, DTA, etc. LOD = 5 (Chen 2019: P < 1e-5).
suppressMessages({
  library(data.table)
  library(ggplot2)
  library(fastman)
  library(ggrepel)
  library(ggtext)
  library(scales)
})

LOD <- 5 # fixed fallback only (Chen 2019: P < 1e-5). Prefer the per-trait perm threshold below.

#' Per-trait JLM permutation threshold (Chen-style empirical alpha), written by
#' scripts/zeal_jlm_run.R to data/zeal/gwas_perm_thresholds.csv. Trait is parsed from
#' the JLM output filename (`<trait>_jlm_native1.txt`). Falls back to the fixed global
#' LOD when the table/row is missing (e.g. a trait not yet permuted, or TeoNAM).
jlm_threshold <- function(jlm_txt, alpha = 0.05) {
  # locals named off the column names (trait, alpha) so the i-filter resolves them
  # from scope (data.table's `..` prefix only works in j, not in an i-filter).
  want_trait <- toupper(sub("_jlm.*$", "", basename(jlm_txt)))
  want_alpha <- alpha
  # TeoNAM (nilHMM) and ZEAL have separate JLM threshold tables — pick by the output
  # path so TeoNAM DTA does not collide with ZEAL DTA on the shared trait key.
  csv <- if (grepl("teonam", jlm_txt)) {
    here::here("data/teonam/jlm_perm_thresholds.csv")
  } else {
    here::here("data/zeal/gwas_perm_thresholds.csv")
  }
  if (file.exists(csv)) {
    thr <- fread(csv)
    hit <- thr[trait == want_trait & model == "jlm" & alpha == want_alpha, thr_neglog10p]
    if (length(hit) && is.finite(hit[1])) {
      return(hit[1])
    }
  }
  warning(sprintf("no JLM perm threshold for %s (alpha=%.2f); using fixed LOD %g", want_trait, alpha, LOD))
  LOD
}

#' Coordinate transformer replicating fastman_gg's internal x layout, so gene
#' labels / vertical lines can be placed at an arbitrary (chr, bp).
get_transformer <- function(m) {
  cmat <- do.call(rbind, lapply(
    split(m[order(m$CHR, m$BP), ], m$CHR[order(m$CHR, m$BP)]),
    function(d) {
      data.frame(
        CHR = d$CHR[1], min_bp = min(d$BP),
        width = max(d$BP) - min(d$BP),
        medgap = if (nrow(d) > 1) median(diff(sort(d$BP))) else 1
      )
    }
  ))
  cmat <- cmat[order(cmat$CHR), ]
  maxgap <- max(cmat$medgap, na.rm = TRUE)
  numc <- nrow(cmat)
  cmat$base <- 0
  cmat$midp <- 0
  cmat$midp[1] <- cmat$width[1] / 2
  for (i in 2:numc) {
    cmat$base[i] <- cmat$base[i - 1] + cmat$width[i - 1] + maxgap
    cmat$midp[i] <- cmat$base[i] + cmat$width[i] / 2
  }
  fac <- numc / cmat$midp[numc]
  cmat$basef <- fac * cmat$base
  function(chr, bp) {
    i <- match(chr, cmat$CHR)
    fac * (bp - cmat$min_bp[i]) + cmat$basef[i]
  }
}

#' Genomic inflation factor: median(chi^2_obs) / median(chi^2_null), 1-df scale.
lambda_gc <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  round(median(qchisq(p, df = 1, lower.tail = FALSE)) / qchisq(0.5, df = 1), 3)
}

#' Read a TASSEL MLM output into a scan (SNP, CHR, BP[v5], P) using the joint 2-df
#' additive+dominance p-value (column `p`), dropping the mean ("None") row.
read_mlm_tassel <- function(txt) {
  d <- fread(txt)
  d <- d[Marker != "None" & !is.na(Chr) & is.finite(p)]
  data.table(SNP = d$Marker, CHR = as.integer(d$Chr), BP = as.integer(d$Pos), P = as.numeric(d$p))
}

#' Per-(trait, model, geno) GWAS genome-wide FWER permutation threshold, written by
#' zeal_ols_taxon.R / zeal_mlm_taxon.R to data/zeal/gwas_perm_thresholds.csv. Replaces
#' the inline Benjamini-Hochberg fdr_lod() line in the recovery notebooks. model is
#' "ols" or "mlm"; geno is the caller ("rtiger_mosaic", "hwe_post_gt", ...). Returns NA
#' (no line drawn) if the row is missing.
gwas_thr <- function(trait, model, geno, alpha = 0.05) {
  want_trait <- toupper(trait)
  want_model <- model
  want_geno <- geno
  want_alpha <- alpha
  csv <- here::here("data/zeal/gwas_perm_thresholds.csv")
  if (file.exists(csv)) {
    d <- fread(csv)
    hit <- d[trait == want_trait & model == want_model & geno == want_geno & alpha == want_alpha, thr_neglog10p]
    if (length(hit) && is.finite(hit[1])) {
      return(hit[1])
    }
  }
  warning(sprintf("no GWAS perm threshold for %s / %s / %s (alpha=%.2f)", want_trait, model, geno, alpha))
  NA_real_
}

#' Per-(trait, model) TeoNAM GWAS genome-wide FWER permutation threshold for the raw
#' 118K point Manhattans, written by teonam_stam_gwas118k.R / teonam_mlm_family_118k.R
#' to data/teonam/gwas_perm_thresholds.csv (separate from ZEAL; TeoNAM DTA would collide).
#' model is "ols" or "mlm". Falls back to the fixed LOD if the row is missing.
teonam_gwas_thr <- function(trait, model, geno = "raw118k", alpha = 0.05) {
  want_trait <- toupper(trait)
  want_model <- model
  want_geno <- geno
  want_alpha <- alpha
  csv <- here::here("data/teonam/gwas_perm_thresholds.csv")
  if (file.exists(csv)) {
    d <- fread(csv)
    hit <- d[trait == want_trait & model == want_model & geno == want_geno & alpha == want_alpha, thr_neglog10p]
    if (length(hit) && is.finite(hit[1])) {
      return(hit[1])
    }
  }
  warning(sprintf("no TeoNAM GWAS perm threshold for %s / %s / %s; using fixed LOD %g", want_trait, model, geno, LOD))
  LOD
}

#' Manhattan with candidate genes starred at their true v5 positions.
#' `overlap_csv` = a candidate table with columns symbol, chr, start (v5).
plot_manhattan <- function(scan_csv, title, overlap_csv, out_png = NULL, mark = NULL, lod = LOD) {
  scan <- if (is.data.frame(scan_csv)) as.data.frame(scan_csv) else as.data.frame(fread(scan_csv))
  scan <- scan[is.finite(scan$P) & scan$P > 0, ]
  m <- scan[order(scan$CHR, scan$BP), c("SNP", "CHR", "BP", "P")]

  ov <- fread(overlap_csv)
  ov[, y := vapply(seq_len(.N), function(i) {
    w <- m$P[m$CHR == chr[i] & abs(m$BP - start[i]) <= 5e5]
    if (length(w)) -log10(min(w)) else NA_real_
  }, numeric(1))]
  tr <- get_transformer(m)
  ov[, BPn := tr(chr, start)]
  ov <- ov[is.finite(y)]

  mp <- -log10(min(m$P))
  ticks <- pretty(c(0, mp))
  unit <- diff(ticks)[1]
  ytop <- max(ticks) + unit

  p <- fastman_gg(
    m = m, snp = "SNP", col = c("black", "gray75"), maxP = mp,
    genomewideline = NULL, suggestiveline = NULL,
    ylab = expression(-log[10](italic(P))), xlab = ""
  ) +
    geom_hline(yintercept = if (is.finite(lod)) lod else Inf, linetype = "dotted", linewidth = 0.7, color = "black") +
    geom_point(data = ov, aes(x = BPn, y = y), inherit.aes = FALSE, shape = 8, size = 3, stroke = 1, color = "red") +
    ggrepel::geom_text_repel(
      data = ov, aes(x = BPn, y = y, label = symbol), inherit.aes = FALSE,
      fontface = "italic", size = 5, segment.color = "grey50", min.segment.length = 0,
      box.padding = 0.5, max.overlaps = Inf, ylim = c(mp + 0.3 * unit, ytop)
    ) +
    ggtitle(title) +
    theme_classic(base_size = 18) +
    scale_y_continuous(breaks = ticks, limits = c(0, ytop), expand = expansion(mult = c(0.01, 0))) +
    theme(
      plot.title = element_markdown(hjust = 1, size = 15, face = "bold"),
      legend.position = "none",
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2)
    )
  if (!is.null(mark)) {
    xm <- tr(mark$chr, mark$bp)
    p <- p + geom_vline(xintercept = xm, linetype = "dashed", color = "red", linewidth = 0.5)
  }
  if (!is.null(out_png)) ggsave(out_png, p, width = 9, height = 4.6, dpi = 200, bg = "white")
  p
}

#' Chen 2019 Fig 4A: JLM QTL genomic distribution as a lollipop. `overlap_csv` =
#' candidate table (symbol, chr, start); `label_genes` = symbols to label in bold-italic.
plot_lollipop <- function(jlm_txt, scan_csv, title, overlap_csv, label_genes, out_png = NULL, label_dist = 5e5, lod = NULL) {
  if (is.null(lod)) lod <- jlm_threshold(jlm_txt) # per-trait perm threshold (Chen-style)
  scan <- if (is.data.frame(scan_csv)) as.data.frame(scan_csv) else as.data.frame(fread(scan_csv))
  scan <- scan[is.finite(scan$P) & scan$P > 0, ]
  maxis <- scan[order(scan$CHR, scan$BP), c("SNP", "CHR", "BP", "P")]

  j <- fread(jlm_txt)
  j <- j[!(Name %in% c("mean", "Family", "Error"))]
  q <- data.table(SNP = j$Name, CHR = as.integer(j$Locus), BP = as.integer(j$Position), logP = -log10(as.numeric(j$`pr>F`)))
  q <- q[logP >= lod]
  ov <- fread(overlap_csv)
  ov <- ov[symbol %in% label_genes]
  # Label a QTL with a candidate gene only when the gene actually co-locates with it
  # (nearest QTL within label_dist bp = the same +/-500 kb window the candidate/OLS tables
  # use). Without this, a distant gene is snapped onto an unrelated QTL on the same
  # chromosome, so its lollipop position disagrees with its true position in the tables
  # (e.g. chi1 at chr1 300 Mb pulled onto an 86 Mb QTL).
  q[, gene := ""]
  for (i in seq_len(nrow(ov))) {
    idx <- which(q$CHR == ov$chr[i])
    if (length(idx)) {
      d <- abs(q$BP[idx] - ov$start[i])
      if (min(d) <= label_dist) q$gene[idx[which.min(d)]] <- ov$symbol[i]
    }
  }
  tr <- get_transformer(maxis)
  q[, BPn := tr(CHR, BP)]
  maxlogP <- if (nrow(q)) max(q$logP) else lod # no QTL >= threshold: keep axis finite, draw an empty panel
  brks <- seq(0, floor(maxlogP / 10) * 10, by = 10)
  ytop <- maxlogP * 1.18

  p <- fastman_gg(
    m = maxis, snp = "SNP", col = c("#00000000", "#00000000"), maxP = maxlogP,
    genomewideline = NULL, suggestiveline = NULL,
    ylab = expression(-log[10](italic(P))), xlab = ""
  ) +
    geom_hline(yintercept = lod, linetype = "dotted", linewidth = 0.6, color = "black") +
    geom_segment(data = q, aes(x = BPn, xend = BPn, y = 0, yend = logP), inherit.aes = FALSE, color = "blue", linewidth = 0.5) +
    geom_point(data = q, aes(x = BPn, y = logP), inherit.aes = FALSE, color = "blue", size = 2.5) +
    ggrepel::geom_text_repel(
      data = q[gene != ""], aes(x = BPn, y = logP, label = gene), inherit.aes = FALSE,
      fontface = "bold.italic", size = 5, segment.color = NA, box.padding = 0.4,
      point.padding = 0.3, nudge_y = 0.05 * maxlogP, seed = 1, max.overlaps = Inf
    ) +
    ggtitle(title) + theme_classic(base_size = 18) +
    scale_y_continuous(breaks = brks, limits = c(0, ytop), expand = expansion(mult = c(0.01, 0))) +
    theme(
      plot.title = element_markdown(hjust = 1, size = 15, face = "bold"),
      legend.position = "none",
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2)
    )
  if (!is.null(out_png)) ggsave(out_png, p, width = 9, height = 4, dpi = 200, bg = "white")
  p
}

#' Coverage-sweep line-Manhattan (one geom_line per coverage; viridis, infinity black).
#' Adaptive y-limit fits the tallest peak; infinity/0.5x LOD-5 peak loci marked as
#' below-axis / top triangles. Coverage legend top-left.
plot_sweep_line <- function(sweep_csv, title, legend = TRUE, base_size = 14,
                            legend_corner = c("topleft", "topright")) {
  legend_corner <- match.arg(legend_corner)
  sw <- as.data.table(fread(sweep_csv))
  sw <- sw[is.finite(P) & P > 0]
  sw[, logP := -log10(P)]
  tr <- get_transformer(as.data.frame(sw[coverage == sw$coverage[1], .(CHR, BP)]))
  sw[, BPn := tr(CHR, BP)]
  cov_levels <- c("∞", "20", "10", "5", "1", "0.5", "0.2", "0.1")
  sw[, cov_lab := factor(ifelse(is.finite(coverage), formatC(coverage), "∞"), levels = cov_levels)]
  cov_asc <- c("0.1", "0.2", "0.5", "1", "5", "10", "20")
  pal <- setNames(grDevices::hcl.colors(7, "viridis"), cov_asc)
  pal["∞"] <- "black"
  axis_df <- sw[, .(mid = (min(BPn) + max(BPn)) / 2), by = CHR][order(CHR)]
  ytop <- ceiling(max(sw$logP, na.rm = TRUE) * 1.02) # fit the data (peaks vary by trait; a fixed limit clips)
  sw[, dord := ifelse(is.finite(coverage), -coverage, Inf)]

  # LOD-5 peak loci (above-threshold markers clumped within 1 Mb -> top marker per
  # clump): infinity = black up-triangles just below the axis (caller ceiling), 0.5x =
  # down-triangles in the 0.5x legend colour near the top.
  peak_loci <- function(cvsub, gap = 1e6) {
    s <- cvsub[logP > LOD]
    if (!nrow(s)) {
      return(s[0])
    }
    setorder(s, CHR, BP)
    s[, pk := cumsum(CHR != shift(CHR, fill = -1L) | BP - shift(BP, fill = -1L) > gap)]
    s[, .SD[which.max(logP)], by = pk]
  }
  pk_inf <- peak_loci(sw[!is.finite(coverage)])
  pk_05 <- peak_loci(sw[coverage == 0.5])

  ggplot(sw[order(dord, CHR, BP)], aes(x = BPn, y = logP, color = cov_lab, group = interaction(CHR, coverage))) +
    geom_line(linewidth = 0.4) +
    geom_hline(yintercept = LOD, linetype = "dotted", linewidth = 0.7, color = "grey30") +
    geom_point(
      data = pk_inf, aes(x = BPn, y = -0.045 * ytop), inherit.aes = FALSE,
      shape = 24, fill = "black", colour = "black", size = 2.4
    ) + # infinity LOD-5 loci: black up-triangles BELOW the x-axis (shape 24 mirrors 25)
    geom_point(
      data = pk_05, aes(x = BPn, y = 0.98 * ytop), inherit.aes = FALSE,
      shape = 25, fill = pal[["0.5"]], colour = pal[["0.5"]], size = 2.4
    ) + # 0.5x LOD-5 loci: down-triangles in the 0.5x legend colour, at the top
    scale_color_manual(values = pal, name = "coverage (×)", drop = FALSE) +
    scale_x_continuous(breaks = axis_df$mid, labels = axis_df$CHR, expand = c(0.01, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0.01, 0))) +
    coord_cartesian(ylim = c(0, ytop), clip = "off") + # allow the below-axis infinity triangles to render in the margin
    labs(x = "chromosome", y = expression(-log[10](italic(P))), title = title) +
    theme_classic(base_size = base_size) +
    theme(
      plot.margin = margin(t = 5, r = 6, b = 16, l = 6), # room below the axis for the infinity triangles
      plot.title = element_markdown(hjust = 0, size = base_size, face = "bold"),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2),
      legend.position = if (legend) c(if (legend_corner == "topright") 0.995 else 0.005, 0.97) else "none",
      legend.justification = if (legend_corner == "topright") c(1, 1) else c(0, 1),
      legend.background = element_rect(fill = scales::alpha("white", 0.7), colour = NA),
      legend.key.size = unit(0.85, "lines"),
      legend.title = element_text(size = base_size - 3),
      legend.text = element_text(size = base_size - 4),
      legend.margin = margin(1, 3, 1, 3)
    ) +
    guides(color = guide_legend(override.aes = list(linewidth = 1.2)))
}
