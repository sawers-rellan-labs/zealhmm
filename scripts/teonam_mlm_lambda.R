#!/usr/bin/env Rscript
# OLS vs MLM(Q+K) on the SAME interpolated 47,750-marker matrix: genomic inflation
# factor lambda_GC, QQ plots, and the ts2 peak under each model. Answers the
# original ts2-inflation question at the MODEL level (holding genotypes fixed).
suppressMessages({
  library(data.table)
  library(ggplot2)
})
setwd("/Users/fvrodriguez/repos/zealhmm")
OUT <- "results/sim/teonam"

lam <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  round(median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(.5, 1), 3)
}

ols <- fread("data/teonam/stam_gwas_scan_interpolated.csv") # interpolated OLS (47,750): SNP,CHR,BP,P
mlm <- fread("data/teonam/tassel/mlm_interp2.txt") # MLM: Trait,Marker,Chr,Pos,df,F,p,...
mlm <- mlm[Marker != "None"]
mlm[, p := as.numeric(p)]
mlm <- mlm[is.finite(p)]

cat(sprintf(
  "OLS-interp : n=%d  lambda_GC=%.3f  max -log10P=%.2f\n",
  nrow(ols[is.finite(P) & P > 0]), lam(ols$P), max(-log10(ols$P[is.finite(ols$P) & ols$P > 0]))
))
cat(sprintf(
  "MLM  Q+K   : n=%d  lambda_GC=%.3f  max -log10P=%.2f\n",
  nrow(mlm), lam(mlm$p), max(-log10(mlm$p[mlm$p > 0]))
))

# ts2 peak (candidate overlap has symbol,chr,start[v5])
ov <- fread("results/sim/teonam/stam_candidate_overlap.csv")
ts2 <- ov[symbol == "ts2"]
if (nrow(ts2)) {
  peak <- function(chr, bp, p, g) {
    w <- which(chr == g$chr & abs(bp - g$start) <= 5e5 & is.finite(p) & p > 0)
    if (length(w)) round(max(-log10(p[w])), 2) else NA
  }
  cat(sprintf(
    "\nts2 (chr%d ~%.1f Mb) peak -log10P:  OLS=%s   MLM=%s\n", ts2$chr[1], ts2$start[1] / 1e6,
    peak(ols$CHR, ols$BP, ols$P, ts2), peak(mlm$Chr, mlm$Pos, mlm$p, ts2)
  ))
}

# QQ (both models)
qq <- function(p, nm) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  o <- -log10(sort(p))
  e <- -log10(ppoints(length(p)))
  k <- c(which(o > 2), sample(which(o <= 2), min(sum(o <= 2), 6000)))
  data.table(model = nm, exp = e[k], obs = o[k])
}
lab_ols <- sprintf("OLS (no K), λ=%.2f", lam(ols$P))
lab_mlm <- sprintf("MLM Q+K, λ=%.2f", lam(mlm$p))
qd <- rbind(qq(ols$P, lab_ols), qq(mlm$p, lab_mlm))
qd[, model := factor(model, levels = c(lab_ols, lab_mlm))]
p <- ggplot(qd[order(-as.integer(model))], aes(exp, obs, colour = model)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey40", linetype = "dashed") +
  geom_point(size = .6, alpha = .5) +
  scale_colour_manual(values = setNames(c("firebrick", "steelblue4"), c(lab_ols, lab_mlm)), name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  labs(
    x = expression(Expected ~ -log[10](italic(P))), y = expression(Observed ~ -log[10](italic(P))),
    title = "STAM GWAS QQ (interpolated 47,750 markers): OLS vs MLM (Q=5PC + K=cIBS)",
    subtitle = "One panel; dashed = null. MLM (blue) hugs the null through the bulk; OLS (red) is inflated."
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
    legend.background = element_rect(fill = scales::alpha("white", 0.7), colour = NA)
  )
ggsave(file.path(OUT, "stam_ols_vs_mlm_qq.png"), p, width = 10, height = 5, dpi = 170, bg = "white")
cat(sprintf("\nwrote %s\n", file.path(OUT, "stam_ols_vs_mlm_qq.png")))
