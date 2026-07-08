#!/usr/bin/env Rscript
# TeoNAM STAM GWAS scan on the AUTHENTIC 118,838-SNP Chen 2019 GWAS panel
# (Qiuyue Chen's Drive release), lifted to v5 -- the faithful Fig 4C reproduction.
# Replaces the earlier 51K-map-panel approximation (stam_gwas_scan_*.csv).
#   per-marker model: STAM ~ Family + marker(additive), F-test on marker (1 df),
#   Family levels dropped to those present for that marker -- identical model to
#   scripts/teonam_stam_gwas_family_imputed.R.
# Output: data/teonam/stam_gwas_scan_118k.csv (SNP, CHR[v5], BP[v5], P, n).
# Run: Rscript scripts/teonam_stam_gwas118k.R
suppressMessages({
  library(data.table)
  library(parallel)
  library(readxl)
})
source("/Users/fvrodriguez/repos/zealhmm/scripts/logging.R")

g118 <- readRDS("data/teonam/teonam_gwas118k_dosage.rds")
dos <- g118$dos # integer matrix [markers x lines], NA where untyped
lines <- colnames(dos)
fam <- factor(substr(lines, 1, 5)) # TIL01A001 -> TIL01

ph <- as.data.frame(read_excel("data/teonam/9250682/TeoNAM_1257RILs_22traits_phenotype_data.xlsx"))
names(ph)[1] <- "line"
y <- suppressWarnings(as.numeric(setNames(ph$STAM, ph$line)[lines]))

mc <- fread("data/teonam/markers_v5_gwas118k.tsv") # 118K v2->v5 liftover roster
chr_by <- setNames(mc$chr_v5, mc$marker)
pos_by <- setNames(mc$pos_v5, mc$marker)
mk <- intersect(rownames(dos), mc$marker) # lifted markers only
mk_idx <- match(mk, rownames(dos))
log_info("%s", paste("lines:", length(lines), " markers (lifted):", length(mk), " STAM non-NA:", sum(!is.na(y))))

scan1 <- function(i) {
  g <- dos[i, ] # dosage across lines, NA where untyped
  ok <- !is.na(g) & !is.na(y)
  n <- sum(ok)
  if (n < 20) {
    return(c(NA_real_, n))
  }
  gg <- as.numeric(g[ok])
  if (stats::sd(gg) == 0) {
    return(c(NA_real_, n))
  }
  yy <- y[ok]
  ff <- droplevels(fam[ok])
  Xr <- if (nlevels(ff) > 1) model.matrix(~ff) else matrix(1, n, 1)
  RSS0 <- sum(lm.fit(Xr, yy)$residuals^2)
  fit <- lm.fit(cbind(Xr, gg), yy)
  RSS1 <- sum(fit$residuals^2)
  df2 <- n - fit$rank
  if (df2 <= 0 || RSS1 <= 0) {
    return(c(NA_real_, n))
  }
  c(pf(((RSS0 - RSS1) / 1) / (RSS1 / df2), 1, df2, lower.tail = FALSE), n)
}

out <- mclapply(mk_idx, scan1, mc.cores = max(1L, detectCores() - 2L))
P <- vapply(out, `[`, numeric(1), 1)
N <- vapply(out, `[`, numeric(1), 2)

scan <- data.table(
  SNP = mk, CHR = as.integer(chr_by[mk]), BP = as.integer(pos_by[mk]),
  P = P, n = as.integer(N)
)[order(CHR, BP)]
fwrite(scan, "data/teonam/stam_gwas_scan_118k.csv")
log_info("%s", paste(
  "scan markers:", nrow(scan), " tested:", sum(!is.na(scan$P)),
  " max -log10P:", round(-log10(min(scan$P, na.rm = TRUE)), 1),
  " median n/marker:", median(scan$n)
))
print(head(scan[order(P)], 10))
