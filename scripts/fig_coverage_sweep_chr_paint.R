#!/usr/bin/env Rscript
# Figure 3 (fig:paint): chromosome painting of the coverage-sweep NILs, six dispatcher-correct
# caller lanes per row at their recorded operating points (reference/zeal_skim_calibration_
# operating_points.csv). Rows = B73 control (must paint all-REF) + 10 vary_skim NILs, ordered by
# skim coverage lambda. Lanes (top->bottom): Skim {nnil, bbnil, rtiger, binhmm} + BrB {googa, atlas}.
# Uses the NATIVE nilHMM::paint_calls() with a `method` track; this script only prepares the
# multi-track common-schema calls (relabel skim/brb samples to the NIL id) from cached segments.
#   Rscript scripts/coverage_sweep_chr_paint.R

suppressMessages({
  library(devtools)
  library(data.table)
  library(here)
  library(ggplot2)
})
devtools::load_all(path.expand("~/repos/nilhmm"), quiet = TRUE) # paint_calls
OUT <- here::here("results/sim/zeal_nil")
KEEP <- c("name", "chr", "start_bp", "end_bp", "state")

# ---- cached caller segments --------------------------------------------------
SEG <- file.path(OUT, "paint_seg_cache")
seg <- list(
  "Skim-nnil"   = readRDS(file.path(SEG, "nnil.rds")),
  "Skim-bbnil"  = readRDS(file.path(SEG, "bbnil.rds")),
  "Skim-rtiger" = readRDS(file.path(SEG, "rtiger.rds")),
  "Skim-binhmm" = readRDS(file.path(SEG, "binhmm.rds")),
  "BrB-googa"   = readRDS(file.path(OUT, "googa_nir_seg_cache", "nir_0.200.rds")),
  "BrB-atlas"   = readRDS(file.path(OUT, "atlas_nir_seg_cache", "nir_0.100.rds"))
)
methods <- names(seg)
skim_methods <- c("Skim-nnil", "Skim-bbnil", "Skim-rtiger", "Skim-binhmm")
LANE_ORDER <- c("Skim-binhmm", "Skim-nnil", "Skim-bbnil", "Skim-rtiger", "BrB-googa", "BrB-atlas") # top->bottom

# ---- sample -> NIL correspondence (B73 control + 10 vary_skim NILs) ----------
mem <- fread(here::here("data/skimsweep/coverage_sweep_members.csv"))[sweep == "vary_skim"]
corr <- rbind(
  data.table(nil = "B73", skim = "PN10_SID893", brb = "PN3_SID213", skim_cov = NA_real_),
  mem[, .(nil, skim = skim_name, brb = brb_name, skim_cov)]
)
skim2nil <- setNames(corr$nil, corr$skim)
brb2nil <- setNames(corr$nil, corr$brb)

relabel <- function(s, map, method) {
  s <- as.data.table(s)[, ..KEEP][name %in% names(map)]
  s[, `:=`(name = map[name], method = method)]
  s
}
calls <- rbindlist(lapply(methods, function(m) {
  relabel(seg[[m]], if (m %in% skim_methods) skim2nil else brb2nil, m)
}))

# ---- ordered factors: rows (B73 top, then coverage) + track lanes -----------
ord <- c("B73", corr[nil != "B73"][order(skim_cov), nil])
row_lab <- setNames(ifelse(ord == "B73", "B73\ncontrol",
  sprintf("%.2fx\n%s", corr$skim_cov[match(ord, corr$nil)], gsub("_", "\n", ord))
), ord) # wrap at underscores
calls[, name := factor(row_lab[name], levels = row_lab[ord])]
calls[, method := factor(method, levels = LANE_ORDER)] # lane order: top = first level

# ---- native painter + styling: no x-axis, bold larger left/right labels ------
p <- paint_calls(as.data.frame(calls), track = "method") +
  labs(
    x = NULL,
    title = "Chromosome painting: coverage-sweep NILs x six calibrated callers",
    subtitle = "rows = B73 control + 10 vary_skim NILs (skim coverage / pedigree); lanes = 4 Skim (DNA) + 2 BrB (RNA-seq) callers"
  ) +
  theme(
    axis.title.x = element_blank(), axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), axis.line.x = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 12, face = "bold", lineheight = 0.85),
    axis.text.y.right = element_text(size = 12, face = "bold")
  )
ggsave(file.path(OUT, "coverage_sweep_chr_paint.png"), p, width = 14, height = 13, dpi = 150)
cat(sprintf(
  "rows=%d (incl B73) | lanes=%d | native paint_calls | wrote coverage_sweep_chr_paint.png\n",
  length(ord), length(methods)
))
