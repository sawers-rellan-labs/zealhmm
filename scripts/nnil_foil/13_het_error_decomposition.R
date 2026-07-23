#!/usr/bin/env Rscript
# Calibration foil, step 13: TEST the het-error-floor hypothesis.
#
# The real 24-line GBS-vs-chip marker mismatch bottoms at ~1.8% (nir ~= 0.95),
# far above Holland's File_S04 0.07%. Hypothesis (Fausto): that ~2% floor is the
# HET-calling-error minimum -- disagreements that involve a HET call on one side
# (spurious/dropped hets), not genuine REF<->ALT ancestry disagreement. Test it by
# decomposing the 3x3 REF/HET/ALT confusion (chip truth x nnil call) at a ladder of
# nir into HET-involving vs hom<->hom (REF<->ALT) mismatch. If HET-involving
# dominates and the hom<->hom residual is ~Holland's 0.07%, the floor is het error.
#
#   Rscript scripts/nnil_foil/13_het_error_decomposition.R
# Output: data/nnil_foil/het_error_decomposition.csv

suppressMessages({
  library(nilHMM)
  library(BEDMatrix)
  library(data.table)
  library(jsonlite)
})
root <- here::here()
for (f in list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)) source(f) # rasterize_states/state_accuracy
source(file.path(root, "scripts/logging.R"))
FOIL <- file.path(root, "data/nnil_foil")
EQUIV <- file.path(root, "data/nnil_equiv")

# ---- data: identical load to scripts/nnil_foil/08_nir_sweep.R -----------------
xw <- fread(file.path(FOIL, "markers_v5.tsv"))
setkey(xw, marker)
chip <- fread(file.path(FOIL, "chip_truth_projected.csv"))
gbs_lines <- readLines(file.path(EQUIV, "lines.csv"))
both <- intersect(chip$Line, gbs_lines)
stopifnot(length(both) == 24)
chip <- chip[Line %in% both]
setorder(chip, Line)
chip_markers <- setdiff(names(chip), "Line")
chip_markers <- chip_markers[chip_markers %in% xw$marker]
cm_info <- xw[chip_markers, .(marker, chr, pos = pos_v5)]
setorder(cm_info, chr, pos)
grid_eval <- cm_info[, .(chr, pos)]

states_to_segments <- function(mat, lines, marker_pos) {
  mp <- marker_pos[match(colnames(mat), marker)]
  out <- vector("list", nrow(mat))
  for (i in seq_len(nrow(mat))) {
    dt <- data.table(chr = mp$chr, pos = mp$pos, state = as.integer(mat[i, ]))[!is.na(state) & state != 3L]
    if (!nrow(dt)) next
    setorder(dt, chr, pos)
    dt[, run := rleid(chr, state)]
    seg <- dt[, .(start_bp = min(pos), end_bp = max(pos), state = state[1]), by = .(chr, run)]
    seg[, name := lines[i]]
    out[[i]] <- seg[, .(name, chr, start_bp, end_bp, state)]
  }
  rbindlist(out)
}
chip_mat <- as.matrix(chip[, ..chip_markers])
rownames(chip_mat) <- chip$Line
tr <- states_to_segments(chip_mat, chip$Line, cm_info[, .(marker, chr, pos)])

geno <- BEDMatrix(file.path(EQUIV, "geno.bed"))
md <- fread(file.path(EQUIV, "markers.csv"))
row_idx <- match(both, gbs_lines)
xw_by_v4 <- xw[match(md$marker, marker_v4)]
keep_col <- which(!is.na(xw_by_v4$marker))
g_raw <- geno[row_idx, keep_col, drop = FALSE]
v5_chr <- xw_by_v4$chr[keep_col]
v5_pos <- xw_by_v4$pos_v5[keep_col]
ord <- order(v5_chr, v5_pos)
g_raw <- g_raw[, ord, drop = FALSE]
v5_chr <- v5_chr[ord]
v5_pos <- v5_pos[ord]
data <- rbindlist(lapply(seq_along(both), function(i) {
  g <- as.integer(g_raw[i, ])
  g[is.na(g)] <- 3L
  data.table(name = both[i], chr = v5_chr, pos = v5_pos, g = g)
}))
setorder(data, name, chr, pos)
hp <- fromJSON(file.path(EQUIV, "params.json"))
map_r <- fromJSON(file.path(FOIL, "chip_calib.json"))$map_r
log_info("loaded %d NILs vs chip | %d shared markers", length(both), length(chip_markers))

# ---- accumulate the chip x call confusion over all lines, per nir ------------
lab <- c("REF", "HET", "ALT")
confusion_at <- function(nir) {
  called <- as.data.table(call_ancestry(
    data = data, caller = "nnil", rrate = map_r,
    germ = hp$germ, gert = hp$gert, p = hp$p, nir = nir, mr = hp$mr, f_1 = hp$f_1, f_2 = hp$f_2
  ))
  C <- matrix(0, 3, 3, dimnames = list(lab, lab))
  for (L in both) {
    a <- state_accuracy(called[name == L], tr[name == L], grid_eval)$confusion
    C <- C + as.matrix(a)
  }
  C
}

nir_ladder <- c(0.594, 0.9, 0.95)
res <- rbindlist(lapply(nir_ladder, function(nir) {
  C <- confusion_at(nir)
  N <- sum(C)
  diagsum <- sum(diag(C))
  mismatch <- 1 - diagsum / N
  # off-diagonal decomposition
  het_involving <- (C["HET", "REF"] + C["HET", "ALT"] + C["REF", "HET"] + C["ALT", "HET"]) / N
  hom_hom <- (C["REF", "ALT"] + C["ALT", "REF"]) / N
  log_info(
    "nir=%.3f | mismatch=%.4f | het-involving=%.4f (%.0f%% of mismatch) | hom<->hom=%.4f (%.0f%%)",
    nir, mismatch, het_involving, 100 * het_involving / mismatch, hom_hom, 100 * hom_hom / mismatch
  )
  log_info(
    "  directions: chipREF->callHET=%.4f  chipHET->callREF=%.4f  chipALT->callHET=%.4f  chipHET->callALT=%.4f  chipREF<->ALT=%.4f",
    C["REF", "HET"] / N, C["HET", "REF"] / N, C["ALT", "HET"] / N, C["HET", "ALT"] / N, (C["REF", "ALT"] + C["ALT", "REF"]) / N
  )
  data.table(
    nir = nir, mismatch = mismatch, het_involving = het_involving, hom_hom = hom_hom,
    het_pct_of_mismatch = het_involving / mismatch,
    chipREF_callHET = C["REF", "HET"] / N, chipHET_callREF = C["HET", "REF"] / N,
    chipALT_callHET = C["ALT", "HET"] / N, chipHET_callALT = C["HET", "ALT"] / N,
    chipREF_ALT = (C["REF", "ALT"] + C["ALT", "REF"]) / N
  )
}))
fwrite(res, file.path(FOIL, "het_error_decomposition.csv"))
# full 3x3 confusion at the optimum (genome fraction), for the notebook
C95 <- confusion_at(0.95)
C95f <- round(C95 / sum(C95), 6)
fwrite(as.data.table(C95f, keep.rownames = "chip_truth"), file.path(FOIL, "het_error_confusion.csv"))
log_info(
  "VERDICT (nir 0.95): mismatch %.4f = %.0f%% hom<->hom (REF<->ALT, donor recall) + %.0f%% HET-involving -> het-error-floor hypothesis %s.",
  res[nir == 0.95]$mismatch, 100 * (1 - res[nir == 0.95]$het_pct_of_mismatch), 100 * res[nir == 0.95]$het_pct_of_mismatch,
  ifelse(res[nir == 0.95]$het_pct_of_mismatch < 0.5, "REFUTED", "supported")
)
cat("\n=== confusion at the optimum (nir=0.95), row=chip truth, col=nnil call ===\n")
print(C95f)
