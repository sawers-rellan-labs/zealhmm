#!/usr/bin/env Rscript
# Two-panel donor-introgression fragment-size validation in cM, ONE ROW x TWO COLUMNS.
# Replaces the 11-line vary_skim `fig_zeal_fragment_size_cm_validation.png`. Panel order is
# calibration first, then application at scale (Fausto):
#
#   A  CALIBRATION AGAINST MolBreeding TRUTH. Reference = molb target-seq (~110x) donor fragments
#      (nnil on the hard calls at the fixed truth config), NOT the simulation. Callers = the four
#      skim callers on the same 14 molb NILs, plus googa/atlas for the 4 molb NILs that have BrB.
#   B  ALL COHORT vs SIMULATION. The full skim-x-BrB paired cohort (330 NILs after the coverage-QC
#      exclusions) vs the BC2S3 sim latent donor-tract cM law. NO GAMMA FIT (Fausto).
#
# It is a fragment-SIZE DISTRIBUTION, so fragments are POOLED per series and NO per-line pairing is
# required: series legitimately differ in n, and every label states its own dataset, caller,
# parameters and n ([[label-every-series-explicitly]]).
#
# cM via nilHMM::bp_to_cm (Hyman monotone Marey spline) on the NATIVE TeoNAM v5 map
# (markers_snp50k_cm.tsv, 1559 cM), identical in both panels. Linear cM x so the ~exponential
# BC2S3 law reads straight.
#
# Inputs (run these first):
#   scripts/stage_zeal_paired_cohort.R + scripts/zeal_paired_cohort_decode.R   -> panel B
#   scripts/zeal_paired_cohort_coverage_qc.R                                  -> QC exclusions
#   scripts/stage_zeal_molb_cohort.R  + scripts/zeal_molb_cohort_decode.R     -> panel A
#
#   Rscript scripts/fig_zeal_fragment_size_cm_two_panel.R

suppressMessages({
  library(devtools)
  library(data.table)
  library(here)
  library(ggplot2)
  library(patchwork)
})
devtools::load_all(path.expand("~/repos/nilhmm"), quiet = TRUE)
source(here::here("R/metrics.R"))
source(here::here("scripts/logging.R"))

OUT <- here::here("results/sim/zeal_nil")
PB <- file.path(OUT, "paired_cohort")
PA <- file.path(OUT, "molb_panel")
KEEP <- c("name", "chr", "start_bp", "end_bp", "state")
BASE <- 16

# ---- shared cM spline (same map, both panels) --------------------------------
nmap <- fread(here::here("data/zeal/markers_snp50k_cm.tsv"))[, .(chr = as.integer(chr), bp = as.integer(pos), cm = as.numeric(cm))]
to_cm <- bp_to_cm(as.data.frame(nmap))

# ---- QC exclusions apply to ALL SIX callers, never per caller ----------------
excl <- fread(here::here("reference/zeal_paired_cohort_exclusions.csv"))
excl[, pedigree := sub("\\.B$", "", pedigree)]
qc_out <- excl[reason_code == "failed_coverage_qc", pedigree]
log_info("[fig2] coverage-QC exclusions applied to every caller: %s", paste(qc_out, collapse = ", "))

read_cm <- function(dir, caller, drop = character(0)) {
  f <- file.path(dir, sprintf("%s_donor_cm.csv", caller))
  if (!file.exists(f)) {
    log_warn("[fig2] missing %s", f)
    return(data.table(name = character(0), cm = numeric(0)))
  }
  d <- fread(f)[, .(name, cm)]
  n0 <- uniqueN(d$name)
  d <- d[!name %in% drop]
  if (n0 != uniqueN(d$name)) log_info("[fig2] %s: dropped %d QC line(s)", caller, n0 - uniqueN(d$name))
  d
}

# ---- palette: one colour per caller, shared across panels --------------------
# LEGEND CARRIES THE METHOD NAME ONLY (Fausto). Dataset, emission/duration parameters and the
# per-series n belong in the figure caption in main.tex, not crammed into the key. The series
# levels are identical in both panels (reference + six callers) so patchwork collects ONE legend;
# each panel's subtitle says what its reference is. Caller names stay lowercase = our implementation
# ([[caller-casing-lowercase-vs-namesake]]).
CAL <- c("nnil", "bbnil", "rtiger", "binhmm", "googa", "atlas")
pal_caller <- c(
  nnil = "#009E73", bbnil = "#0072B2", rtiger = "#56B4E9",
  binhmm = "#999999", googa = "#D55E00", atlas = "#E69F00"
)
L_REF <- "reference"
# Operating points are NOT plotted and NOT put in the figure caption (Fausto: tabular data belongs
# in a table). They are emitted as a LaTeX table fragment that main.tex \input's, so no number is
# ever hand-copied into the paper.
SPEC <- data.table(
  caller = c("nnil", "bbnil", "rtiger", "binhmm", "googa", "atlas"),
  macro = c("\\nnil", "\\bbnil", "\\rtiger", "\\binhmm", "\\googa", "\\atlas"),
  # kept terse: the table must fit \textwidth in a 9pt twocolumn layout, so symbols are
  # defined in the table caption rather than spelled out in every cell.
  data_source = c("skim", "skim", "skim", "skim bins", "BRB-seq", "BRB-seq"),
  emission = c("categorical", "BetaBin.", "BetaBin.", "Gaussian", "categorical", "categorical"),
  duration = c(
    "$r=6.4$e$-4$", "$r=6.4$e$-4$", "$\\ell=5$", "$\\sigma=0.995$", "$r=8.9$e$-4$", "$\\ell=50$"
  ),
  extra = c("$0.70$", "---", "---", "---", "$0.20$", "$0.10$")
)
DESC <- setNames(sprintf("%s, %s / %s", SPEC$data_source, SPEC$emission, SPEC$duration), SPEC$caller)

# ================================ PANEL A ===================================
# reference = molb target-seq truth; callers = 4 skim on 14 NILs + googa/atlas on the 4 with BrB
molb_ped <- fread(here::here("data/zeal/molb_cohort/molb_cohort_manifest.csv"))$pedigree
A <- list()
tr <- fread(file.path(PA, "molb_truth_donor_cm.csv"))[, .(name, cm)]
A[[L_REF]] <- tr
for (cl in c("nnil", "bbnil", "rtiger", "binhmm")) A[[cl]] <- read_cm(PA, cl)
# BrB callers: only the molb NILs that have BrB, taken from the paired-cohort decode
for (cl in c("googa", "atlas")) A[[cl]] <- read_cm(PB, cl, drop = qc_out)[name %in% molb_ped]
dtA <- rbindlist(lapply(names(A), function(k) data.table(cm = A[[k]]$cm, series = k)))
lvlA <- c(L_REF, CAL)
dtA[, series := factor(series, levels = lvlA)]
ksA <- dtA[series != L_REF, .(D = fragment_size_ks(cm, tr$cm), n_frag = .N, med_cm = median(cm)), by = series]
ksA[, `:=`(n_nil = vapply(as.character(series), function(k) uniqueN(A[[k]]$name), integer(1)), descriptor = DESC[as.character(series)])]
log_info(
  "[fig2] PANEL A reference = MolBreeding target-seq ~110x, nnil (rrate=map 1.7e-3, nir=0.70): %d NILs, %d fragments, median %.1f cM",
  uniqueN(tr$name), nrow(tr), median(tr$cm)
)
log_info("[fig2] PANEL A: cM KS vs molb truth (caption material)")
print(ksA[, .(series, D = round(D, 3), n_nil, n_frag, med_cm = round(med_cm, 1))])

# ================================ PANEL B ===================================
# reference = BC2S3 sim latent donor tracts; callers = all six on the 330-NIL cohort. NO gamma fit.
sim_truth <- as.data.table(readRDS(file.path(OUT, "zeal_nil_bc2s3_full.rds"))$truth)
db <- as.data.table(.donor_blocks(sim_truth[, ..KEEP]))
db[, cm := to_cm(chr, end_bp) - to_cm(chr, start_bp)]
ref_cm <- db[is.finite(cm) & cm > 0, cm]
B <- list()
B[[L_REF]] <- data.table(name = NA_character_, cm = ref_cm)
for (cl in CAL) B[[cl]] <- read_cm(PB, cl, drop = qc_out)
dtB <- rbindlist(lapply(names(B), function(k) data.table(cm = B[[k]]$cm, series = k)))
lvlB <- c(L_REF, CAL)
dtB[, series := factor(series, levels = lvlB)]
ksB <- dtB[series != L_REF, .(D = fragment_size_ks(cm, ref_cm), n_frag = .N, med_cm = median(cm)), by = series]
ksB[, `:=`(n_nil = vapply(as.character(series), function(k) uniqueN(B[[k]]$name), integer(1)), descriptor = DESC[as.character(series)])]
log_info(
  "[fig2] PANEL B reference = BC2S3 simulated latent ancestry: %d NILs, %d fragments, median %.1f cM",
  uniqueN(sim_truth$name), length(ref_cm), median(ref_cm)
)
log_info("[fig2] PANEL B: cM KS vs BC2S3 sim latent truth (caption material)")
print(ksB[, .(series, D = round(D, 3), n_nil, n_frag, med_cm = round(med_cm, 1))])

# shared scales: identical series levels in both panels -> patchwork collects ONE legend
PAL <- setNames(c("black", unname(pal_caller[CAL])), c(L_REF, CAL))
LTY <- setNames(c("dotted", rep("solid", length(CAL))), c(L_REF, CAL))

# ================================= plot =====================================
xmax <- as.numeric(quantile(c(dtA$cm, dtB$cm), 0.99))
# Axis and legend carry the MINIMUM (Fausto): no map provenance on the axis, no reference
# restatement in a subtitle. The map (native TeoNAM v5), the dotted-line convention, the operating
# points and every per-series n live in the figure caption in main.tex.
mk <- function(d, ttl) {
  ggplot(d, aes(cm, colour = series, linetype = series)) +
    stat_ecdf(linewidth = 0.9) +
    scale_colour_manual(values = PAL, name = NULL, drop = FALSE, guide = guide_legend(nrow = 1)) +
    scale_linetype_manual(values = LTY, name = NULL, drop = FALSE, guide = guide_legend(nrow = 1)) +
    coord_cartesian(xlim = c(0, xmax)) +
    labs(x = "donor introgression size (cM)", y = "ECDF", title = ttl) +
    theme_bw(base_size = BASE) +
    theme(plot.title = element_text(size = BASE * 0.95))
}
pA <- mk(dtA, "Calibration against MolBreeding truth")
pB <- mk(dtB, "Full paired cohort against simulation")

fig <- (pA | pB) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 20, face = "bold"), plot.tag.location = "plot", plot.tag.position = "topleft",
    legend.position = "bottom", legend.direction = "horizontal",
    legend.text = element_text(size = BASE * 0.85), legend.key.width = grid::unit(26, "pt")
  )
ggsave(file.path(OUT, "zeal_fragment_size_cm_two_panel.png"), fig, width = 15, height = 7, dpi = 150) # name matches main.tex \includegraphics (NO fig_ prefix)
fwrite(
  rbind(
    ksA[, .(panel = "A", reference = "MolBreeding target-seq ~110x truth", series, descriptor, D, n_nil, n_frag, med_cm)],
    ksB[, .(panel = "B", reference = "BC2S3 simulated latent ancestry", series, descriptor, D, n_nil, n_frag, med_cm)]
  ),
  file.path(OUT, "fig_zeal_fragment_size_cm_two_panel_ks.csv")
)

# ---- LaTeX table: operating points + per-panel agreement --------------------
# The paper keeps this table INLINE in main.tex (Fausto's call, matching tab:ranking), so this file
# is not \input by the paper. It is the regenerable reference: re-run this script and diff it
# against main.tex to catch the inline copy drifting from the data.
TABDIR <- OUT
# The all-recurrent baseline, computed not asserted: what fraction of MolBreeding truth calls is REF,
# i.e. the per-marker accuracy a caller gets for free by painting the whole genome recurrent.
hc_g <- fread(here::here("data/zeal/molbreeding/molbreeding_hardcalls_wsfilt.tsv"))[!is.na(g), .N, by = g]
ref_frac <- hc_g[g == 0L, N] / sum(hc_g$N)
log_info(
  "[fig2] all-recurrent baseline: %.1f%% of molb truth calls are REF (HET %.1f%%, ALT %.1f%%) -> free per-marker accuracy %.3f, zero donor fragments",
  100 * ref_frac, 100 * hc_g[g == 1L, N] / sum(hc_g$N), 100 * hc_g[g == 2L, N] / sum(hc_g$N), ref_frac
)
tab <- merge(SPEC, ksA[, .(caller = as.character(series), nA = n_nil, fA = n_frag, mA = med_cm, dA = D)], by = "caller")
tab <- merge(tab, ksB[, .(caller = as.character(series), nB = n_nil, fB = n_frag, mB = med_cm, dB = D)], by = "caller")
tab <- tab[match(SPEC$caller, caller)] # keep the dispatcher-ordered caller sequence
# LaTeX thousands separator, so 12624 prints as 12{,}624
tex_num <- function(x) gsub(",", "{,}", formatC(as.integer(x), format = "d", big.mark = ","), fixed = TRUE)
rows <- tab[, sprintf(
  "%s & %s & %.1f & %.3f & %s & %.1f & %.3f \\\\",
  macro, tex_num(fA), mA, dA, tex_num(fB), mB, dB
)]
tex <- c(
  "% GENERATED by scripts/fig_zeal_fragment_size_cm_two_panel.R -- do not edit by hand.",
  # ONE job: how well each caller's DONOR (HET or ALT) fragments match each reference. The
  # operating points are Methods material and live in the text, not here. Column names say
  # "donor" explicitly: ~90% of these NIL genomes is REF, so an all-recurrent painter scores ~0.90
  # per-marker accuracy while recovering ZERO donor fragments, and any metric not restricted to
  # donor fragments is close to meaningless (Fausto; cf. [[zeal-mismatch-is-degenerate-use-fragdsc]]).
  "\\begin{table}[tb]",
  "\\centering",
  "\\small",
  "\\begin{tabular}{@{}l rrr rrr@{}}",
  "\\toprule",
  " & \\multicolumn{3}{c}{A: vs MolBreeding} & \\multicolumn{3}{c}{B: vs simulation} \\\\",
  "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}",
  "caller & \\multicolumn{1}{c}{$n$} & median & KS & \\multicolumn{1}{c}{$n$} & median & KS \\\\",
  "\\midrule",
  rows,
  "\\bottomrule",
  "\\end{tabular}",
  sprintf(
    paste0(
      "\\caption{\\textbf{Donor-fragment recovery against each reference.} ",
      "Every column is restricted to \\emph{donor} (HET or ALT) fragments, because in these NILs ",
      "%.1f\\%% of truth calls are REF: a caller that painted the whole genome recurrent would score ",
      "$%.3f$ per-marker accuracy while recovering zero donor fragments, so unrestricted accuracy ",
      "cannot separate callers. $n$ is the number of donor fragments the caller calls, \\emph{median} ",
      "their median length in cM, and KS the Kolmogorov--Smirnov distance between the caller's ",
      "donor-fragment length distribution and the reference's. Lower KS is better; $n$ far above the ",
      "reference's indicates over-fragmentation. The panel A reference is the MolBreeding target-seq ",
      "(${\\sim}110\\times$) mosaic on $9{,}157$ sites: %s donor fragments from %d NILs, median ",
      "%.1f\\,cM. The panel B reference is the BC2S3 simulated latent ancestry: %s donor fragments ",
      "from %s simulated NILs, median %.1f\\,cM. Skim callers use all %d MolBreeding NILs in panel A, ",
      "while \\googa and \\atlas act on the BRB-seq pangene grid and so contribute only the %d that ",
      "also carry BRB-seq data; in panel B the cohort is %d NILs, of which %d and %d respectively ",
      "yield at least one donor fragment. Operating points are given in the text.}"
    ),
    100 * ref_frac, ref_frac,
    tex_num(nrow(tr)), uniqueN(tr$name), median(tr$cm),
    tex_num(length(ref_cm)), tex_num(uniqueN(sim_truth$name)), median(ref_cm),
    tab[caller == "nnil", nA], tab[caller == "googa", nA],
    tab[caller == "nnil", nB], tab[caller == "googa", nB], tab[caller == "atlas", nB]
  ),
  "\\label{tab:fragsize}",
  "\\end{table}"
)
writeLines(tex, file.path(TABDIR, "zeal_fragment_size_operating_points.tex"))
log_info("[fig2] wrote zeal_fragment_size_cm_two_panel.png + _ks.csv + tables/zeal_fragment_size_operating_points.tex")
