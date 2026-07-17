#!/usr/bin/env Rscript
# =============================================================================
# ZEAL Phase 3 — LB-Impute ancestry mosaic on the REAL SNP50K counts.
#
# The gateway genotype artifact: LB-Impute (map-aware, unit="cm") decodes teosinte
# ancestry state (0=B73, 1=het, 2=teo) at EVERY SNP50K marker, per chromosome, for the
# full NIL panel, from the observed allele counts. recombdist (cM) is taken from the
# BC2S3 simulation calibration at the coverage nearest the real per-site depth
# (scripts/zeal_lbimpute_calib_bycov.R). No thinning (all markers). The lbimpute HMM
# decodes each line independently given the BC2S3 start seed (f_1,f_2), so the panel
# runs per chromosome in one call (threads-parallel) — no per-family grouping needed.
#
# Feeds: JLM genotypes (then 0.1-cM thinned), Panel C MLM, and VanRaden K.
# Input : data/zeal/zeal_snp50k_dosage.rds  (n_ref, n_alt, markers, samples)
#         data/zeal/markers_snp50k_cm.tsv    (marker, chr, pos, cm)
#         data/zeal/samplesheet_3way.csv     (gwas_nil skim_id -> pedigree, donor_accession)
#         results/sim/zeal/lbimpute_calib_bycov.csv  (recombdist* per coverage)
# Output: data/zeal/zeal_lbimpute_mosaic.rds list(markers, state[marker x line], lines, recombdist, lambda)
# =============================================================================
suppressMessages({
  library(data.table)
  devtools::load_all("/Users/fvrodriguez/repos/nilhmm", quiet = TRUE)
})
ROOT <- "/Users/fvrodriguez/repos/zealhmm"
setwd(ROOT)
for (f in list.files(file.path(ROOT, "R"), "\\.R$", full.names = TRUE)) source(f)
source(file.path(ROOT, "scripts/logging.R"))

GENOTYPEERR <- 0.05
DRP <- TRUE
ERR <- 0.01
THREADS <- as.integer(Sys.getenv("LBIMPUTE_THREADS", as.character(max(1L, parallel::detectCores() - 2L))))
EXP <- breeding_prior("BC2S3") # BC2S3 start seed
F1 <- as.numeric(EXP["HET"])
F2 <- as.numeric(EXP["ALT"])

# --- inputs ------------------------------------------------------------------
D <- readRDS("data/zeal/zeal_snp50k_dosage.rds")
mk <- copy(D$markers)[, .(marker, chr = as.integer(chr), pos = as.integer(pos))]
cm <- fread("data/zeal/markers_snp50k_cm.tsv")
mk[, cm := cm$cm[match(marker, cm$marker)]]
mk <- mk[!is.na(cm)]
ss <- fread("data/zeal/samplesheet_3way.csv")[gwas_nil == TRUE & !is.na(skim_id)]
panel <- ss[skim_id %in% colnames(D$n_ref), unique(skim_id)]
ped <- ss[match(panel, skim_id), pedigree]
log_info("panel: %d NIL lines x %d markers (chr1-10)", length(panel), nrow(mk))

# recombdist at the coverage nearest the real per-site depth ------------------
obs_depth <- mean(D$cov[mk$marker, panel])
calib <- fread("results/sim/zeal/lbimpute_calib_bycov.csv")[coverage != "Inf"]
calib[, lam := as.numeric(coverage)]
sel <- calib[which.min(abs(lam - obs_depth))]
RD <- sel$recombdist
log_info(
  "real per-site depth = %.3f -> nearest calib coverage lambda=%s, recombdist* = %.2f cM (FDR %.3f)",
  obs_depth, sel$coverage, RD, sel$fdr
)

# --- decode per chromosome (whole panel; lbimpute decodes each line) ----------
nref <- D$n_ref[mk$marker, panel]
nalt <- D$n_alt[mk$marker, panel]
t0 <- Sys.time()
blocks <- list()
for (ch in 1:10) {
  idx <- which(mk$chr == ch)
  long <- data.table(
    name = rep(panel, each = length(idx)),
    chr = ch, pos = rep(mk$pos[idx], times = length(panel)),
    cm = rep(mk$cm[idx], times = length(panel)),
    n_ref = as.vector(nref[idx, ]), n_alt = as.vector(nalt[idx, ])
  )
  st <- as.data.table(call_states(long,
    caller = "lbimpute", unit = "cm", recombdist = RD,
    err = ERR, genotypeerr = GENOTYPEERR, drp = DRP, f_1 = F1, f_2 = F2,
    min_reads = 1L, threads = THREADS
  ))
  blocks[[ch]] <- st[, .(name, chr, pos, state)]
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  log_info(">>> chr%d done (%d markers) | elapsed %.1f min | ETA ~%.1f min", ch, length(idx), el, el / ch * (10 - ch))
}
S <- rbindlist(blocks)

# --- assemble marker x line state matrix (columns keyed by pedigree) ----------
mat <- dcast(S, chr + pos ~ name, value.var = "state")
setkey(mat, chr, pos)
mk2 <- mk[order(chr, pos)]
mat <- mat[mk2[, .(chr, pos)], on = c("chr", "pos")]
state <- as.matrix(mat[, panel, with = FALSE])
rownames(state) <- mk2$marker
colnames(state) <- ped # GWAS line id = pedigree

teo_frac <- mean(state, na.rm = TRUE) / 2
# breakpoints/line QC: within-chr state transitions, averaged over lines (BC2S3 sim ~14.7)
chr_of <- mk2$chr
bp_per_line <- apply(state, 2, function(s) {
  sum(vapply(1:10, function(ch) {
    v <- s[chr_of == ch]
    sum(v[-1] != v[-length(v)], na.rm = TRUE)
  }, integer(1)))
})
log_info(
  "mosaic: %d markers x %d lines | mean teosinte fraction = %.3f (BC2S3 expect ~0.11) | NA=%d",
  nrow(state), ncol(state), teo_frac, sum(is.na(state))
)
log_info(
  "QC breakpoints/line: mean %.1f (BC2S3 sim ~14.7; far higher => REF-noise over-fragmentation)",
  mean(bp_per_line)
)
# 102 invariant-REF markers flagged for downstream K / per-SNP GWAS exclusion (kept here as backbone)
inv <- rownames(state)[rowSums(D$n_alt[mk2$marker, panel]) == 0]
log_info("invariant-REF markers (drop from K/per-SNP GWAS, keep in mosaic backbone): %d", length(inv))
saveRDS(inv, "data/zeal/snp50k_invariant_markers.rds")
saveRDS(
  list(
    markers = mk2, state = state, lines = data.table(skim_id = panel, pedigree = ped),
    recombdist = RD, lambda = sel$coverage, obs_depth = obs_depth
  ),
  "data/zeal/zeal_lbimpute_mosaic.rds"
)
log_info("wrote data/zeal/zeal_lbimpute_mosaic.rds")
