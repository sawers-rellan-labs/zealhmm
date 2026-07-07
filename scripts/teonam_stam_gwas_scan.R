#!/usr/bin/env Rscript
# TeoNAM STAM GWAS scan (for Chen 2019 Fig 4C reproduction).
# Densify the FULL union of markers (no LD prune -- GWAS uses all markers), then
# a per-marker across-family GLM scan: STAM ~ Family + marker(additive), 1 df.
# Output: data/teonam/stam_gwas_scan_interpolated.csv (SNP, CHR[v5], BP[v5], P).
# Run: Rscript agent/teonam_stam_gwas_scan.R
suppressMessages({
  library(data.table)
  library(parallel)
  library(readxl)
  devtools::load_all("/Users/fvrodriguez/repos/nilhmm", quiet = TRUE)
})

mc <- fread("data/teonam/marker_info_v5_cm.tsv")
setnames(mc, "chr_v5", "chr")
cm_by <- setNames(mc$cm, mc$marker)
pos_by <- setNames(mc$pos_v5, mc$marker)
# FULL 51,004 GWAS target (teonam_map_v5_gwas). interpolate_genotype accepts a
# duplicate-cM target, so target all genotyped markers directly (no unique-cM dedup).
gcols <- names(fread("data/teonam/TeoNAM_genotype_clean.csv", nrows = 0))[-(1:3)]
GWAS_MK <- intersect(gcols, mc$marker)

fams <- c(
  TIL01 = "W22TIL01_genotype.csv", TIL03 = "W22TIL03_genotype.csv",
  TIL11 = "W22TIL11_genotype.csv", TIL14 = "W22TIL14_genotype.csv",
  TIL25 = "W22TIL25_genotype.csv"
)
fam_data <- lapply(names(fams), function(fam) {
  g <- fread(file.path("data/teonam", fams[fam]))
  g <- g[!duplicated(g[[1]])]
  list(g = g, keys = paste0(fam, sub("^.*Line_", "", g[[1]]))) # robust key: TIL01 + A001
})
names(fam_data) <- names(fams)

# --- densify the FULL unique-cM union per chromosome (no prune) --------------
densify_chr <- function(ch) {
  mch <- mc[chr == ch]
  setorder(mch, cm)
  tgt <- mch[marker %in% GWAS_MK]
  setorder(tgt, cm)
  tgt_df <- data.frame(chr = ch, cm = tgt$cm)
  blocks <- lapply(names(fam_data), function(fam) {
    g <- fam_data[[fam]]$g
    mk <- intersect(names(g)[-(1:3)], mch$marker)
    obs <- data.frame(marker = mk, cm = cm_by[mk])
    obs <- obs[order(obs$cm), ]
    obs <- obs[!duplicated(obs$cm), ]
    cols <- obs$marker
    geno <- t(as.matrix(g[, ..cols]))
    storage.mode(geno) <- "double"
    dn <- interpolate_genotype(geno, data.frame(chr = ch, cm = obs$cm), tgt_df, mode = "step")
    colnames(dn) <- fam_data[[fam]]$keys
    dn
  })
  D <- do.call(cbind, blocks)
  rownames(D) <- tgt$marker
  list(geno = D, markers = tgt[, .(marker, chr, cm)])
}
res <- mclapply(1:10, densify_chr, mc.cores = 5)
G <- do.call(rbind, lapply(res, `[[`, "geno"))
markers <- rbindlist(lapply(res, `[[`, "markers"))
cat("full densified union:", nrow(G), "markers x", ncol(G), "lines\n")

# --- GWAS scan: STAM ~ Family + marker (single additive effect across families) ---
ph <- as.data.frame(read_excel("data/teonam/9250682/TeoNAM_1257RILs_22traits_phenotype_data.xlsx"))
names(ph)[1] <- "line"
y <- suppressWarnings(as.numeric(setNames(ph$STAM, ph$line)[colnames(G)]))
fam <- factor(substr(colnames(G), 1, 5))
ok <- !is.na(y)
y <- y[ok]
fam <- fam[ok]
Gm <- G[, ok, drop = FALSE]
Xr <- model.matrix(~fam)
n <- length(y)
RSS0 <- sum(lm.fit(Xr, y)$residuals^2)
scan1 <- function(i) {
  g <- Gm[i, ]
  if (sd(g) == 0) {
    return(NA_real_)
  }
  fit <- lm.fit(cbind(Xr, g), y)
  RSS1 <- sum(fit$residuals^2)
  df2 <- n - fit$rank
  pf(((RSS0 - RSS1) / 1) / (RSS1 / df2), 1, df2, lower.tail = FALSE)
}
P <- unlist(mclapply(seq_len(nrow(Gm)), scan1, mc.cores = 6))

scan <- data.table(
  SNP = rownames(G), CHR = as.integer(markers$chr),
  BP = as.integer(pos_by[rownames(G)]), P = P
)[order(CHR, BP)]
fwrite(scan, "data/teonam/stam_gwas_scan_interpolated.csv")
cat("scan markers:", nrow(scan), " max -log10P:", round(-log10(min(P, na.rm = TRUE)), 1), "\n")
print(head(scan[order(P)], 5))
