#!/usr/bin/env Rscript
# ZEAL — ancestry mosaic from any nilhmm caller on the real SNP50K counts.
# CALLER env: nnil | binhmm | lbimpute (rtiger uses zealtiger's per-donor calls via
# zeal_rtiger_mosaic.R). Runs call_ancestry per chromosome over the whole NIL panel
# (each line decoded independently given the BC2S3 priors), assembles the marker x line
# state matrix, and VALIDATES teosinte presence vs the established BzeaSeq ~0.101.
# This is the caller-comparison: which callers recover the right ancestry at 0.4x.
#
# Output: data/zeal/zeal_<caller>_mosaic.rds  list(markers, state[marker x line], lines)
suppressMessages({
  library(here)
  library(data.table)
  library(SpATS)
})
devtools::load_all("/Users/fvrodriguez/repos/nilhmm", quiet = TRUE)
source(here("scripts/logging.R"))
CALLER <- Sys.getenv("CALLER", "nnil")
DESIGN <- "BC2S3"
ERR <- 0.01
THREADS <- as.integer(Sys.getenv("LBIMPUTE_THREADS", as.character(max(1L, parallel::detectCores() - 2L))))

D <- readRDS(here("data/zeal/zeal_snp50k_dosage.rds"))
mk <- copy(D$markers)[, .(marker, chr = as.integer(chr), pos = as.integer(pos))]
cm <- fread(here("data/zeal/markers_snp50k_cm.tsv"))
mk[, cm := cm$cm[match(marker, cm$marker)]]
mk <- mk[!is.na(cm)]
ss <- fread(here("data/zeal/samplesheet_3way.csv"))[gwas_nil == TRUE & !is.na(skim_id)]
panel <- ss[skim_id %in% colnames(D$n_ref), unique(skim_id)]
ped <- ss[match(panel, skim_id), pedigree]
nref <- D$n_ref[mk$marker, panel]
nalt <- D$n_alt[mk$marker, panel]
log_info("caller=%s | panel %d lines x %d markers", CALLER, length(panel), nrow(mk))

# per-caller call_states args
args_for <- function(sub) {
  base <- list(caller = CALLER, design = DESIGN, err = ERR, min_reads = 1L, threads = 1L)
  if (CALLER == "nnil") base$rrate <- 3.3e-5 # settled nNIL rate (memory nnil-rrate-flat-plateau)
  if (CALLER == "binhmm") {
    base$bin_size <- 1e6
    base$cluster_method <- "gauss"
  }
  if (CALLER == "lbimpute") {
    base$unit <- "cm"
    base$recombdist <- 15.69
    base$genotypeerr <- 0.05
    base$drp <- TRUE
  }
  c(list(data = sub), base)
}

t0 <- Sys.time()
blocks <- list()
for (ch in 1:10) {
  idx <- which(mk$chr == ch)
  long <- data.table(
    name = rep(panel, each = length(idx)), chr = ch,
    pos = rep(mk$pos[idx], times = length(panel)), cm = rep(mk$cm[idx], times = length(panel)),
    n_ref = as.vector(nref[idx, ]), n_alt = as.vector(nalt[idx, ])
  )
  st <- as.data.table(do.call(call_states, args_for(long)))
  blocks[[ch]] <- st[, .(name, chr, pos, state)]
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  log_info(">>> chr%d (%d mk) | %.1f min | ETA ~%.1f min", ch, length(idx), el, el / ch * (10 - ch))
}
S <- rbindlist(blocks)
mat <- dcast(S, chr + pos ~ name, value.var = "state")[mk[order(chr, pos), .(chr, pos)], on = c("chr", "pos")]
state <- as.matrix(mat[, panel, with = FALSE])
rownames(state) <- mk[order(chr, pos)]$marker
colnames(state) <- ped

comp <- prop.table(table(factor(state, levels = 0:2)))
log_info(
  "%s mosaic: %d x %d | B73=%.1f%% het=%.1f%% teo=%.1f%% | PRESENCE=%.3f (established ~0.101)",
  CALLER, nrow(state), ncol(state), 100 * comp["0"], 100 * comp["1"], 100 * comp["2"], sum(comp[c("1", "2")])
)
saveRDS(
  list(markers = mk[order(chr, pos)], state = state, lines = data.table(skim_id = panel, pedigree = ped)),
  here(sprintf("data/zeal/zeal_%s_mosaic.rds", CALLER))
)
log_info("wrote data/zeal/zeal_%s_mosaic.rds", CALLER)
