#!/usr/bin/env Rscript
# STAM GWAS scan on the UNIMPUTED (raw, block-sparse) genotypes -- no
# interpolation. Each marker is tested only on the lines that actually carry a
# call at it (the families that typed it); NA lines are dropped. This avoids the
# step-interpolation plateau artifact seen with the densified matrix.
#   per-marker model: STAM ~ Family + marker(additive), F-test on marker (1 df),
#   Family levels dropped to those present for that marker.
# Output: data/teonam/stam_gwas_scan_family_imputed.csv (SNP, CHR[v5], BP[v5], P, n).
# Run: Rscript agent/teonam_stam_gwas_family_imputed.R
suppressMessages({
  library(data.table)
  library(parallel)
  library(readxl)
})

geno <- fread("data/teonam/TeoNAM_genotype_clean.csv") # 1237 lines x 51482 markers, NA = untyped
key <- paste0(sub("^W22", "", geno[[2]]), sub("^.*Line_", "", geno[[1]])) # TIL01 + A001
fam <- factor(sub("^W22", "", geno[[2]]))

ph <- as.data.frame(read_excel("data/teonam/9250682/TeoNAM_1257RILs_22traits_phenotype_data.xlsx"))
names(ph)[1] <- "line"
y <- suppressWarnings(as.numeric(setNames(ph$STAM, ph$line)[key]))

mc <- fread("data/teonam/markers_v5.tsv") # map-neutral v2->v5 liftover: roster + v5 chr/bp
chr_by <- setNames(mc$chr_v5, mc$marker)
pos_by <- setNames(mc$pos_v5, mc$marker)
mk <- intersect(names(geno)[-(1:3)], mc$marker) # lifted markers only
cat("lines:", nrow(geno), " markers (lifted):", length(mk), " STAM non-NA:", sum(!is.na(y)), "\n")

scan1 <- function(m) {
  g <- geno[[m]] # numeric, NA where untyped
  ok <- !is.na(g) & !is.na(y)
  if (sum(ok) < 20) {
    return(c(NA_real_, sum(ok)))
  }
  gg <- g[ok]
  if (stats::sd(gg) == 0) {
    return(c(NA_real_, sum(ok)))
  }
  yy <- y[ok]
  ff <- droplevels(fam[ok])
  Xr <- if (nlevels(ff) > 1) model.matrix(~ff) else matrix(1, sum(ok), 1)
  n <- sum(ok)
  RSS0 <- sum(lm.fit(Xr, yy)$residuals^2)
  fit <- lm.fit(cbind(Xr, gg), yy)
  RSS1 <- sum(fit$residuals^2)
  df2 <- n - fit$rank
  if (df2 <= 0 || RSS1 <= 0) {
    return(c(NA_real_, n))
  }
  c(pf(((RSS0 - RSS1) / 1) / (RSS1 / df2), 1, df2, lower.tail = FALSE), n)
}
out <- mclapply(mk, scan1, mc.cores = 6)
P <- vapply(out, `[`, numeric(1), 1)
N <- vapply(out, `[`, numeric(1), 2)

scan <- data.table(
  SNP = mk, CHR = as.integer(chr_by[mk]), BP = as.integer(pos_by[mk]),
  P = P, n = as.integer(N)
)[order(CHR, BP)]
fwrite(scan, "data/teonam/stam_gwas_scan_family_imputed.csv")
cat(
  "scan markers:", nrow(scan), " tested:", sum(!is.na(scan$P)),
  " max -log10P:", round(-log10(min(scan$P, na.rm = TRUE)), 1),
  " median n/marker:", median(scan$n), "\n"
)
print(head(scan[order(P)], 6))
