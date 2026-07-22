#!/usr/bin/env Rscript
# =============================================================================
# ZEAL Phase 3 — rtiger ancestry mosaic (the usable one; supersedes the
# LB-Impute mosaic, which under-called teosinte ~8x — see agent/zeal_dta_repro_plan.md).
#
# Reproduce-don't-approximate: use the EXISTING, validated per-accession rtiger-SNP50K
# calls from zealtiger (fit_rtiger_by_donor.R, rigidity=8 autotune pick on the SNP50K
# sim, benchmarked on real NILs) rather than a fresh sim calibration. Those calls are
# ancestry SEGMENTS (name, chr, start_bp, end_bp, state 0/1/2); here we back-project them
# onto the full SNP50K marker grid to get the marker x line state matrix used by the GWAS.
# Presence (het+teo) ~10.7% matches the established BzeaSeq block estimate (~10.1%).
#
# Input : data/zeal/rtiger_50K_calls.csv     (rtiger segments; staged from zealtiger)
#         data/zeal/markers_snp50k_v5.tsv    (marker, chr, pos)
#         data/zeal/samplesheet_3way.csv     (gwas_nil skim_id -> pedigree)
# Output: data/zeal/zeal_rtiger_mosaic.rds   list(markers, state[marker x line], lines)
# =============================================================================
suppressMessages({
  library(here)
  library(data.table)
})
source(here("scripts/logging.R"))

seg <- fread(here("data/zeal/rtiger_50K_calls.csv"))[, .(name,
  chr = as.integer(chr),
  start_bp = as.integer(start_bp), end_bp = as.integer(end_bp), state = as.integer(state)
)]
mk <- fread(here("data/zeal/markers_snp50k_v5.tsv"))[, .(marker, chr = as.integer(chr), pos = as.integer(pos))]
setorder(mk, chr, pos)
ss <- fread(here("data/zeal/samplesheet_3way.csv"))[gwas_nil == TRUE & !is.na(skim_id)]

panel <- ss[skim_id %in% unique(seg$name), skim_id]
ped <- ss[match(panel, skim_id), pedigree]
miss <- ss[!skim_id %in% seg$name, .N]
log_info("panel lines with rtiger calls: %d (%d gwas_nil lines lack calls, dropped)", length(panel), miss)
log_info("markers=%d | segments=%d over %d samples", nrow(mk), nrow(seg), uniqueN(seg$name))

# back-project segments -> per-marker state, per line (segments tile each chromosome;
# findInterval on segment starts gives the covering segment for each marker position)
setkey(seg, name, chr, start_bp)
state <- matrix(NA_integer_, nrow = nrow(mk), ncol = length(panel), dimnames = list(mk$marker, ped))
mk_by_chr <- split(seq_len(nrow(mk)), mk$chr)
t0 <- Sys.time()
for (j in seq_along(panel)) {
  sg <- seg[.(panel[j])]
  for (ch in names(mk_by_chr)) {
    sc <- sg[chr == as.integer(ch)]
    if (!nrow(sc)) next
    ri <- mk_by_chr[[ch]]
    idx <- findInterval(mk$pos[ri], sc$start_bp) # 0 if before first start -> clamp to 1
    idx[idx < 1L] <- 1L
    state[ri, j] <- sc$state[idx]
  }
  if (j %% 300 == 0) {
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    log_info(">>> %d/%d lines | elapsed %.1f min | ETA ~%.1f min", j, length(panel), el, el / j * (length(panel) - j))
  }
}

na_frac <- mean(is.na(state))
comp <- prop.table(table(factor(state, levels = 0:2)))
presence <- sum(comp[c("1", "2")])
dosage <- (comp["1"] * 0.5 + comp["2"]) # teosinte dosage fraction (het=0.5)
log_info(
  "mosaic: %d markers x %d lines | NA=%.2f%% | B73=%.1f%% het=%.1f%% teo=%.1f%%",
  nrow(state), ncol(state), 100 * na_frac, 100 * comp["0"], 100 * comp["1"], 100 * comp["2"]
)
log_info("teosinte PRESENCE (het+teo) = %.3f (established ~0.101) | DOSAGE (het=0.5) = %.3f", presence, dosage)
# breakpoints/line QC (transitions along sorted markers within chr)
chr_of <- mk$chr
bpl <- apply(state, 2, function(s) {
  sum(vapply(unique(chr_of), function(c) {
    v <- s[chr_of == c]
    v <- v[!is.na(v)]
    sum(v[-1] != v[-length(v)])
  }, integer(1)))
})
log_info("breakpoints/line: mean %.1f", mean(bpl))

saveRDS(
  list(markers = mk, state = state, lines = data.table(skim_id = panel, pedigree = ped)),
  here("data/zeal/zeal_rtiger_mosaic.rds")
)
log_info("wrote data/zeal/zeal_rtiger_mosaic.rds")
