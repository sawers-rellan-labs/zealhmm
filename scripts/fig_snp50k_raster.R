#!/usr/bin/env Rscript
# Raster of the SNP50K genotype-call matrix (the raw HWE/pileup calls, pre-HMM).
# One pixel per cell: markers on x (genomic order), lines on y (grouped by taxon).
# Palette: REF gold, HET green, ALT purple, missing grey.
# Source: data/zeal/zeal_hwe_post_gt.rds  ($state = 49,002 markers x 1,403 lines,
# values 0/1/2/NA). Same object exported to release/bzea_genotypes/snp50k/.

suppressPackageStartupMessages(library(png))

GOLD <- "#E0A81C"
GREEN <- "#2E9B57"
PURPLE <- "#6B3FA0"
GREY <- "#9AA0A6"
OUT <- "nilhmm-paper/figures/snp50k_genotype_raster.png"

x <- readRDS("data/zeal/zeal_hwe_post_gt.rds")
st <- x$state # markers x lines, 0/1/2/NA
mk <- x$markers

mord <- order(mk$chr, mk$pos) # markers by genome position (x)
lord <- order(colnames(st)) # lines grouped by taxon (Zd/Zh/Zl/Zv/Zx) (y)
st <- st[mord, lord, drop = FALSE]

idx <- t(st) # rows = lines (y), cols = markers (x)
h <- nrow(idx)
w <- ncol(idx)
code <- idx
code[is.na(code)] <- 3L # 0 REF, 1 HET, 2 ALT, 3 missing
code <- code + 1L # 1..4 for palette lookup
rm(idx, st)
gc()

rgbcol <- function(hex) as.vector(col2rgb(hex)) / 255
palR <- sapply(c(GOLD, GREEN, PURPLE, GREY), function(h) rgbcol(h)[1])
palG <- sapply(c(GOLD, GREEN, PURPLE, GREY), function(h) rgbcol(h)[2])
palB <- sapply(c(GOLD, GREEN, PURPLE, GREY), function(h) rgbcol(h)[3])

write_raster <- function(code_mat, path) {
  hh <- nrow(code_mat)
  ww <- ncol(code_mat)
  a <- array(0, dim = c(hh, ww, 3))
  a[, , 1] <- palR[code_mat]
  a[, , 2] <- palG[code_mat]
  a[, , 3] <- palB[code_mat]
  writePNG(a, path)
  cat(sprintf("  %-52s %d x %d px\n", path, ww, hh))
}

cat("wrote:\n")
write_raster(code, OUT) # full resolution

# display version: thin markers so ~1 marker == 1 screen pixel (no blending on a slide)
step <- ceiling(w / 2400)
OUT2 <- sub("\\.png$", "_display.png", OUT)
write_raster(code[, seq(1, w, by = step), drop = FALSE], OUT2)
cat(sprintf("  (display = every %dth marker, %d markers)\n", step, length(seq(1, w, by = step))))

n <- length(code)
cat(sprintf(
  "REF %.1f%%  HET %.2f%%  ALT %.4f%%  missing %.1f%%\n",
  100 * sum(code == 1) / n, 100 * sum(code == 2) / n,
  100 * sum(code == 3) / n, 100 * sum(code == 4) / n
))
