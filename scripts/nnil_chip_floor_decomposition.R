# Decompose the ~1.9% chip-calibration marker-mismatch floor (nnil states decoded
# from real GBS vs chip calls, 24 both-platform NILs) by the underlying GBS genotype
# call. Answers: is the floor the heterozygous calls expected from low-coverage ML
# (TASSEL-GBS) calling? Result (2026-07-22): no. HET is 0.2% of scored genotypes and
# 5.5% of mismatches; the floor is dominated by REF-hom / missing calls where the
# caller reads B73 but the chip reads donor (donor-allele dropout), not hets.
# Mirrors the data prep in scripts/fig_nnil_sim_calibration.R. Run: Rscript this file.
suppressMessages({
  library(data.table)
  library(here)
  library(nilHMM)
})
source(here::here("R/metrics.R"))
SIMDIR <- here::here("results/sim/nested_nnil")
DESIGN <- "BC5S4"
NIR_CHIP <- 0.80 # chip-calibration optimum

mv5 <- fread(here::here("data/nnil_foil/markers_v5.tsv"))
v4_to_v5 <- setNames(mv5$marker, mv5$marker_v4)
v5_to_v4 <- setNames(mv5$marker_v4, mv5$marker)

sim <- readRDS(file.path(SIMDIR, sprintf("nested_nnil_%s_nnil_gbs.rds", tolower(DESIGN))))
sim_mk <- paste0("S", sim$grid$chr, "_", sim$grid$pos)
common <- intersect(sim_mk, mv5$marker)
total_cM <- sum(mv5[marker %in% common, max(cm, na.rm = TRUE), by = chr]$V1)

chip <- fread(here::here("data/nnil_foil/chip_truth_projected.csv"))
geno <- fread(here::here("data/nnil_equiv/geno_recoded.csv"))
idcol <- names(geno)[1]
lines24 <- intersect(chip$Line, geno[[idcol]])
mk_common <- Reduce(intersect, list(
  sim_mk, setdiff(names(chip), "Line"),
  unname(v4_to_v5[setdiff(names(geno), idcol)])
))
RRATE <- 2 * total_cM / (100 * length(mk_common))
cat(sprintf("lines24=%d  mk_common=%d  rrate=%.2e\n", length(lines24), length(mk_common), RRATE))

# chip states (per-marker segments)
chip_long <- melt(chip[Line %in% lines24, c("Line", mk_common), with = FALSE],
  id.vars = "Line", variable.name = "marker", value.name = "state"
)
chip_long[, `:=`(
  name = Line, chr = as.integer(sub("S(\\d+)_.*", "\\1", marker)),
  pos = as.integer(sub("S\\d+_", "", marker))
)]
chip_seg <- chip_long[order(name, chr, pos),
  .(start_bp = pos, end_bp = pos, state = as.integer(state)),
  by = .(name, chr)
]

# real GBS genotypes g (0=REF-hom, 1=HET, 2=ALT/donor-hom, 3=missing)
cols_v4 <- v5_to_v4[mk_common]
gsub24 <- geno[get(idcol) %in% lines24, c(idcol, cols_v4), with = FALSE]
setnames(gsub24, c(idcol, cols_v4), c("name", mk_common))
real_long <- melt(gsub24, id.vars = "name", variable.name = "marker", value.name = "g")
real_long[, `:=`(
  chr = as.integer(sub("S(\\d+)_.*", "\\1", marker)),
  pos = as.integer(sub("S\\d+_", "", marker)), g = as.integer(g)
)]
real_long <- real_long[order(name, chr, pos)]

cat("\n=== raw GBS genotype composition on 24 lines x mk_common ===\n")
gtab <- real_long[, .N, by = g][order(g)]
gtab[, frac := N / sum(N)]
gtab[, label := c("0 REF-hom", "1 HET", "2 ALT/donor-hom", "3 missing")[g + 1L]]
print(gtab)
cat(sprintf(
  "HET rate (of all): %.3f%%   HET rate (of scored, g<3): %.3f%%\n",
  100 * gtab[g == 1, N] / sum(gtab$N),
  100 * gtab[g == 1, N] / gtab[g < 3, sum(N)]
))

cat("\n=== chip state composition ===\n")
print(chip_long[, .N, by = .(state = as.integer(state))][order(state)])

# nnil decode at chip optimum
nnil_seg <- as.data.table(call_ancestry(
  real_long[, .(name, chr, pos, g)],
  caller = "nnil", design = DESIGN, rrate = RRATE, nir = NIR_CHIP
))

# per-marker states on the chip grid, joined with raw g
grid_chip <- unique(chip_seg[, .(chr, pos = start_bp)])
per <- rbindlist(lapply(lines24, function(nm) {
  data.table(
    name = nm, chr = grid_chip$chr, pos = grid_chip$pos,
    nnil = rasterize_states(nnil_seg[name == nm], grid_chip)$state,
    chip = rasterize_states(chip_seg[name == nm], grid_chip)$state
  )
}))
per <- merge(per, real_long[, .(name, chr, pos, g)], by = c("name", "chr", "pos"), all.x = TRUE)
ev <- per[!is.na(nnil) & !is.na(chip)] # bins scored by marker_dice
ev[, mismatch := nnil != chip]

cat(sprintf(
  "\n=== overall chip-calibration mismatch (should ~1.9%%): %.3f%% (n=%d bins) ===\n",
  100 * mean(ev$mismatch), nrow(ev)
))

cat("\n=== genotype composition of MISMATCHED bins ===\n")
mm <- ev[mismatch == TRUE]
mmtab <- mm[, .N, by = g][order(g)]
mmtab[, frac_of_mismatches := N / sum(N)]
mmtab[, label := c("0 REF-hom", "1 HET", "2 ALT/donor-hom", "3 missing")[g + 1L]]
print(mmtab)

cat("\n=== mismatch rate WITHIN each genotype class ===\n")
byg <- ev[, .(n = .N, mism = sum(mismatch), rate = mean(mismatch)), by = g][order(g)]
byg[, label := c("0 REF-hom", "1 HET", "2 ALT/donor-hom", "3 missing")[g + 1L]]
print(byg)

cat("\n=== how much of the total mismatch do HET calls explain? ===\n")
cat(sprintf(
  "HET bins among mismatches: %d / %d = %.2f%% of all mismatches\n",
  mm[g == 1, .N], nrow(mm), 100 * mm[g == 1, .N] / nrow(mm)
))

cat("\n=== nnil vs chip confusion among mismatches (rows=nnil, cols=chip) ===\n")
print(table(nnil = mm$nnil, chip = mm$chip))
