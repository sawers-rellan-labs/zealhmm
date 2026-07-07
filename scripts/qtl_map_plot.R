#!/usr/bin/env Rscript
# PNG of the composite bcsft(1,4) genetic map (cm_qtl) built by teonam_qtl_map.R.
# Builds an R/qtl `map` object (per-chr cM vectors, origin shifted to 0) and plots
# it two ways: qtl::plot.map (marker-density bars) + a ggplot cM-length panel.
suppressMessages({
  library(data.table)
  library(qtl)
  library(ggplot2)
})
setwd("/Users/fvrodriguez/repos/zealhmm")

d <- fread("data/teonam/marker_info_v5_cm_qtl.tsv")
d <- d[!is.na(cm_qtl) & !is.na(chr_v5)][order(chr_v5, cm_qtl)]

# R/qtl map object: named list of per-chromosome cM vectors, each origin at 0
chrs <- sort(unique(d$chr_v5))
map <- lapply(chrs, function(ch) {
  z <- d[chr_v5 == ch]
  v <- z$cm_qtl - min(z$cm_qtl) # shift origin to 0 (your plotting fix)
  names(v) <- z$marker
  v
})
names(map) <- as.character(chrs)
class(map) <- "map"

tot <- sum(vapply(map, max, numeric(1)))
cat(sprintf("composite map: %d markers, %d chr, %.0f cM total\n", nrow(d), length(map), tot))

# gap check (your d[d>18] diagnostic)
gaps <- unlist(lapply(map, function(x) diff(x)))
cat(sprintf(
  "inter-marker gap: mean %.3f, median %.3f cM; gaps >18 cM: %d\n",
  mean(gaps), median(gaps), sum(gaps > 18)
))

# ---- 1. R/qtl plot.map (the plot(consensus_map) look) -----------------------
png("results/sim/teonam/qtl_map_plot.png", width = 1400, height = 1000, res = 150)
plot.map(map, main = sprintf("TeoNAM composite bcsft(1,4) map — %d markers, %.0f cM", nrow(d), tot))
dev.off()

# ---- 2. ggplot: chromosome cM bars with marker ticks ------------------------
mk <- rbindlist(lapply(chrs, function(ch) data.table(chr = ch, cm = map[[as.character(ch)]])))
ends <- mk[, .(cm = max(cm)), by = chr]
p <- ggplot(mk, aes(x = factor(chr), y = cm)) +
  geom_segment(aes(xend = factor(chr), y = 0, yend = cm), linewidth = 6, colour = "grey85") +
  geom_segment(aes(xend = factor(chr), y = cm - 0.15, yend = cm + 0.15),
    linewidth = 6,
    colour = "grey30", alpha = 0.20
  ) + # marker ticks
  geom_text(data = ends, aes(label = sprintf("%.0f", cm)), vjust = -0.6, size = 3.5) +
  scale_y_reverse(expand = expansion(mult = c(0.02, 0.06))) +
  labs(
    x = "chromosome (v5)", y = "position (cM)",
    title = sprintf("TeoNAM composite bcsft(1,4) genetic map — %d markers, %.0f cM", nrow(d), tot)
  ) +
  theme_classic(base_size = 14)
ggsave("results/sim/teonam/qtl_map_cm_bars.png", p, width = 9, height = 6, dpi = 200, bg = "white")
cat("wrote results/sim/teonam/qtl_map_plot.png and qtl_map_cm_bars.png\n")
