#!/usr/bin/env Rscript
# Pin down the genotype pipeline behind the committed JLM hapmap (geno.hmp.txt).
# Candidate source = per-family W22TIL0x_genotype.csv step-interpolated onto the
# selected markers' cM (same densify as agent/teonam_mlm_interp.R), NOT raw clean.
# Read-only. Run: Rscript agent/teonam_jlm_verify_source.R
suppressMessages({
  library(data.table)
  devtools::load_all("/Users/fvrodriguez/repos/nilhmm", quiet = TRUE)
})
setwd("/Users/fvrodriguez/repos/zealhmm")

hmp <- fread("data/teonam/tassel/geno.hmp.txt")
jlm_mk <- hmp[["rs#"]]
taxa <- names(hmp)[-(1:11)]

mc <- fread("data/teonam/map_v5_coe2008.tsv")
setnames(mc, "chr_v5", "chr")
cm_by <- setNames(mc$cm, mc$marker)

fams <- c(
  TIL01 = "W22TIL01_genotype.csv", TIL03 = "W22TIL03_genotype.csv",
  TIL11 = "W22TIL11_genotype.csv", TIL14 = "W22TIL14_genotype.csv", TIL25 = "W22TIL25_genotype.csv"
)
fam_data <- lapply(names(fams), function(fam) {
  g <- fread(file.path("data/teonam", fams[fam]))
  g <- g[!duplicated(g[[1]])]
  list(g = g, keys = paste0(fam, sub("^.*Line_", "", g[[1]])))
})
names(fam_data) <- names(fams)
fam_markers <- unique(unlist(lapply(fam_data, function(x) names(x$g)[-(1:3)])))

cat(sprintf("JLM markers: %d\n", length(jlm_mk)))
cat(sprintf("  in annotation cM table  : %d\n", sum(jlm_mk %in% mc$marker)))
cat(sprintf("  in per-family geno files: %d\n", sum(jlm_mk %in% fam_markers)))
cat(sprintf(
  "  all taxa reconstructible from fam keys: %s\n",
  all(taxa %in% unlist(lapply(fam_data, `[[`, "keys")))
))

# Reconstruct genotypes for the committed JLM markers by step-interpolation onto
# THEIR cM (grid = jlm markers as target; obs = per-family observed), then compare.
recode <- c("A", "M", "C") # 0/1/2
build <- function(mk_set) {
  mk <- data.table(marker = mk_set, chr = mc[match(mk_set, marker), chr], cm = cm_by[mk_set])
  mk <- mk[!is.na(cm)][order(chr, cm)]
  blocks_by_chr <- lapply(sort(unique(mk$chr)), function(ch) {
    tgt <- mk[chr == ch]
    tgt_df <- data.frame(chr = ch, cm = tgt$cm)
    fam_blocks <- lapply(names(fam_data), function(fam) {
      g <- fam_data[[fam]]$g
      obsmk <- intersect(names(g)[-(1:3)], mc[chr == ch, marker])
      obs <- data.frame(marker = obsmk, cm = cm_by[obsmk])
      obs <- obs[order(obs$cm), ]
      obs <- obs[!duplicated(obs$cm), ]
      geno <- t(as.matrix(g[, obs$marker, with = FALSE]))
      storage.mode(geno) <- "double"
      dn <- interpolate_genotype(geno, data.frame(chr = ch, cm = obs$cm), tgt_df, mode = "step")
      colnames(dn) <- fam_data[[fam]]$keys
      dn
    })
    D <- do.call(cbind, fam_blocks)
    rownames(D) <- tgt$marker
    D
  })
  do.call(rbind, blocks_by_chr)
}
G <- build(jlm_mk)
cat(sprintf("\nreconstructed matrix: %d markers x %d taxa\n", nrow(G), ncol(G)))
common_mk <- intersect(rownames(G), jlm_mk)
common_tx <- intersect(colnames(G), taxa)
rec <- matrix(recode[round(G[common_mk, common_tx]) + 1], length(common_mk)) # markers x taxa
obs <- as.matrix(hmp[match(common_mk, jlm_mk), ..common_tx]) # markers x taxa
cat(sprintf(
  "compare on %d markers x %d taxa: match=%s  mismatches=%d (%.4f%%)\n",
  length(common_mk), length(common_tx), all(rec == obs), sum(rec != obs),
  100 * mean(rec != obs)
))
