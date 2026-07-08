#!/usr/bin/env Rscript
# Assemble the FSFHap-block-smoothed 118K truth. For each family, read FSFHap's
# PARENTAL (A/C) output (the clean ancestry mosaic, ~9 breakpoints/chr), convert
# to dosage (A=0 W22, M=1 het, C=2 teo), and step-interpolate that family's clean
# blocks onto the FULL 118K union cM grid -> a DENSE, clean-block 118K truth for
# every line. Combining the 5 families gives the drop-in replacement for
# teonam_gwas118k_dosage_polar.rds: same shape (markers x 1257), but block-clean
# instead of the choppy per-SNP imputed genotypes.
#
# Output: data/teonam/teonam_gwas118k_dosage_fsfhap.rds (list: dos, markers, lines)
# Run: Rscript scripts/teonam_gwas118k_truth_assemble.R
suppressMessages({
  library(data.table)
  devtools::load_all("/Users/fvrodriguez/repos/nilhmm", quiet = TRUE)
})
ROOT <- "/Users/fvrodriguez/repos/zealhmm"
setwd(ROOT)
TDIR <- file.path(ROOT, "data/teonam/tassel/fsfhap118k")

# --- 118K union cM grid (same as the sweeps) ---------------------------------
mc <- fread("data/teonam/markers_v5_gwas118k_cm.tsv") # marker, chr, pos_v5, cm
setnames(mc, "pos_v5", "pos")
setorder(mc, chr, cm)
target_df <- data.frame(chr = as.integer(mc$chr), cm = as.numeric(mc$cm))
union_markers <- mc$marker
cm_by <- setNames(mc$cm, mc$marker)
chr_by <- setNames(mc$chr, mc$marker)
FAMS <- c("TIL01", "TIL03", "TIL11", "TIL14", "TIL25")

hetcodes <- c("M", "R", "Y", "S", "K", "W")
blocks <- list()
for (fam in FAMS) {
  f <- file.path(TDIR, sprintf("imputed_%s1.hmp.txt", fam)) # ...1 = imputed_parents (A/C)
  d <- fread(f, colClasses = "character")
  keys <- names(d)[-(1:11)]
  mk <- d[["rs#"]]
  keep <- mk %in% union_markers & !is.na(cm_by[mk]) # parental markers on the cM grid
  d <- d[keep]
  mk <- mk[keep]
  G <- as.matrix(d[, ..keys]) # markers x lines, A/C/M/N
  dos <- matrix(NA_real_, nrow(G), ncol(G))
  dos[G == "A"] <- 0
  dos[G %in% hetcodes] <- 1
  dos[G == "C"] <- 2
  # order by cM (dedup ties per chr -> strictly increasing for interpolation)
  mt <- data.table(marker = mk, chr = as.integer(chr_by[mk]), cm = as.numeric(cm_by[mk]))
  o <- order(mt$chr, mt$cm)
  mt <- mt[o]
  dos <- dos[o, , drop = FALSE]
  dup <- mt[, .(dup = duplicated(cm)), by = chr]$dup
  mt <- mt[!dup]
  dos <- dos[!dup, , drop = FALSE]
  # PER-LINE step-interpolation onto the union grid (handles the few scattered N in
  # the parental output by using only each line's observed markers, like the control sweep)
  b <- matrix(NA_real_, nrow(mc), ncol(dos), dimnames = list(union_markers, keys))
  for (j in seq_len(ncol(dos))) {
    obs <- which(!is.na(dos[, j]))
    if (length(obs) < 2) next
    okc <- unique(mt$chr[obs])
    tsel <- target_df$chr %in% okc
    b[tsel, j] <- interpolate_genotype(
      matrix(dos[obs, j], ncol = 1L),
      data.frame(chr = mt$chr[obs], cm = mt$cm[obs]),
      target_df[tsel, , drop = FALSE],
      mode = "step"
    )[, 1L]
  }
  blocks[[fam]] <- b
  cat(sprintf("  %s: %d parental markers -> %d union rows x %d lines\n", fam, nrow(mt), nrow(b), ncol(b)))
}
dos <- do.call(cbind, blocks)
rownames(dos) <- union_markers
storage.mode(dos) <- "integer"
lines <- colnames(dos)
saveRDS(
  list(dos = dos, markers = union_markers, lines = lines),
  "data/teonam/teonam_gwas118k_dosage_fsfhap.rds"
)
cat(sprintf(
  "wrote clean-block truth: %d markers x %d lines (%.2f%% NA)\n",
  nrow(dos), ncol(dos), 100 * mean(is.na(dos))
))

# --- verify choppiness collapsed (chr1, TIL01A001) ---------------------------
c1 <- which(chr_by[union_markers] == 1)
v <- dos[c1, "TIL01A001"][order(mc$pos[c1])]
v <- v[!is.na(v)]
cat(sprintf(
  "choppiness check chr1 TIL01A001: %d transitions (was 473 pre-FSFHap)\n",
  sum(v[-1] != v[-length(v)])
))
