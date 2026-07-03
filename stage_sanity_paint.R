#!/usr/bin/env Rscript
# Stage the skim-sweep subset (11 samples: B73 control + 10 vary-skim NILs) into
# data/, organized by (source x input type). Consumed by
# analysis/nilhmm_sanity_paint.qmd. Idempotent. ~ tens of MB.
#
# Layout:
#   data/ref/                                   genome-level, dataset-independent
#     gene_to_pangene.tsv, b73_gene_coords.tsv  (ATLAS references)
#   data/skimsweep/                             the vary-skim cohort
#     coverage_sweep_members.csv                (full sweep table; note filters it)
#     skim/counts_50k/<skim>.tsv                -> Skim-nNIL, Skim-rtiger
#     skim/bins/<skim>_bin_genotypes.tsv        -> Skim-binhmm (PER-SAMPLE only)
#     brb/counts_wideseq/<brb>.tsv              -> BrB-nNIL, BrB-rtiger
#     brb/pangene/<species>/<brb>.pangene_counts.tsv  -> BrB-atlas
#
# NOTE: touches the rsstu network mount -> run with the sandbox DISABLED so the
# mount is visible. Provenance recorded in DATA.md.

suppressMessages(library(data.table))

# real sources
ZT_SRC <- "/Users/fvrodriguez/Desktop/zealtiger"
CAS_SRC <- "/Volumes/rsstu/users/r/rrellan/tlaloc/cassini"
MNT_SRC <- "/Volumes/rsstu/users/r/rrellan/BZea/bzeaseq/ancestry"

# local roots
SS <- "data/skimsweep"
SKIM <- file.path(SS, "skim")
BRB <- file.path(SS, "brb")
REF <- "data/ref"

copy1 <- function(src, dest) {
  if (is.na(src) || !file.exists(src)) {
    return(FALSE)
  }
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  file.copy(src, dest, overwrite = TRUE)
}
report <- function(tag, ok, n) cat(sprintf("%-14s %d / %d\n", tag, sum(ok), n))

# ---- sample list (species-only dir; note maps tax -> species the same way) ---
mem_src <- file.path(ZT_SRC, "results/sim_calibration/coverage_sweep_members.csv")
mem <- fread(mem_src)[sweep == "vary_skim"][order(vary_cov)]
mem[, tax := sub("\\..*", "", nil)]
tax2sp <- c(
  Zx = "mexicana", Zl = "luxurians", Zd = "diploperennis",
  Zv = "parviglumis", Zh = "huehuetenangensis"
)
mem[, sp := tax2sp[tax]]
rows <- rbind(
  data.table(skim = "PN10_SID893", brb = "PN3_SID213", sp = "mexicana"),
  mem[, .(skim = skim_name, brb = brb_name, sp)]
)
# the note reads + filters the full sweep table, so stage it verbatim
copy1(mem_src, file.path(SS, "coverage_sweep_members.csv"))
cat("samples:", nrow(rows), "\n")

# ---- skim 50K counts ---------------------------------------------------------
ok <- vapply(unique(rows$skim), function(s) {
  f <- list.files(file.path(ZT_SRC, "data/rtiger_50K"), paste0("^", s, "\\.tsv$"),
    full.names = TRUE, recursive = TRUE
  )
  length(f) > 0 && copy1(f[1], file.path(SKIM, "counts_50k", paste0(s, ".tsv")))
}, logical(1))
report("skim counts", ok, length(unique(rows$skim)))

# ---- skim PER-SAMPLE bins (never the aggregate) ------------------------------
ok <- vapply(unique(rows$skim), function(s) {
  copy1(
    file.path(MNT_SRC, paste0(s, "_bin_genotypes.tsv")),
    file.path(SKIM, "bins", paste0(s, "_bin_genotypes.tsv"))
  )
}, logical(1))
report("skim bins", ok, length(unique(rows$skim)))
if (!all(ok)) cat("  MISSING bins:", paste(unique(rows$skim)[!ok], collapse = ", "), "\n")

# ---- BRB wideseq counts ------------------------------------------------------
ok <- vapply(unique(rows$brb), function(b) {
  copy1(
    file.path(ZT_SRC, "results/sim_calibration/brbseq_ks_wideseq/counts", paste0(b, ".tsv")),
    file.path(BRB, "counts_wideseq", paste0(b, ".tsv"))
  )
}, logical(1))
report("brb counts", ok, length(unique(rows$brb)))

# ---- ATLAS pangene counts (brb/pangene/<species>/) ---------------------------
ok <- vapply(seq_len(nrow(rows)), function(i) {
  copy1(
    file.path(CAS_SRC, "results", rows$sp[i], "pangene", paste0(rows$brb[i], ".pangene_counts.tsv")),
    file.path(BRB, "pangene", rows$sp[i], paste0(rows$brb[i], ".pangene_counts.tsv"))
  )
}, logical(1))
report("atlas pangene", ok, nrow(rows))

# ---- genome-level ATLAS references (dataset-independent) ---------------------
copy1(file.path(CAS_SRC, "data/pangene/gene_to_pangene.tsv"), file.path(REF, "gene_to_pangene.tsv"))
copy1(file.path(CAS_SRC, "data/meta/b73_gene_coords.tsv"), file.path(REF, "b73_gene_coords.tsv"))
cat("ref gene maps  ", sum(file.exists(file.path(REF, c("gene_to_pangene.tsv", "b73_gene_coords.tsv")))), "/ 2\n")

sz <- sum(file.info(list.files("data", recursive = TRUE, full.names = TRUE))$size, na.rm = TRUE)
cat(sprintf("\nstaged under data/ - %.1f MB\n", sz / 1e6))
