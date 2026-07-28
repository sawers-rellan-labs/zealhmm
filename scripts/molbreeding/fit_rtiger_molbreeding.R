#!/usr/bin/env Rscript
# Fit RTIGER on the MolBreeding target sequencing data using REAL allele counts (replaces the
# old molbreeding_to_rtiger.R, which faked 0/1/2 pseudo-counts at fixed depth).
#
# target sequencing is target-capture SEQUENCING (~110x), so ref_depth/alt_depth are real read
# counts — fed to RTIGER like the SNP50K skim (Source 2). RTIGER's BetaBinomial EM
# cost scales with the number of distinct (k,n) count pairs, i.e. with coverage.
#
# DEFAULTS (the validated run config; see docs/rtiger-coverage-benchmark.md):
#   --estimator=mle  --cap=20    -> MLE E-step, binomial-thin to ~20x.
# At <=20x, MLE is fast and robust and validated against simulated truth (marker acc
# 0.9999, het seg F1 0.95 at 20x), so capping makes the high-coverage cost — and the
# MLE-vs-MoM choice — moot. Both knobs are overridable:
#   --estimator : mle (default) | mom (closed-form MoM, faster; only trust on capped data —
#                 it diverges on raw high coverage)
#   --cap=N     : binomial-thin counts to ~Nx (default 20; --cap=0 disables). No-ops when
#                 the set's mean depth is already below N (e.g. the 0.4x skim).
#   --r=N       : rigidity (default 8, matching the skim)
#
# Input : data/molbreeding_45k/gatk_table_<set>_v5.tsv  (GATK table, v5)
#         data/molbreeding_45k_sample_map.tsv, data/rtiger_50K/seqlengths.csv
# Output: data/molbreeding_targetseq/<set>/{counts/,expDesign.csv,fit/,calls_common_schema.csv}
#           source = MolBreeding_<set>_r<r>_<METHOD>[_cap<N>]
#
# Usage:  Rscript fit_rtiger_molbreeding.R [set ...] [--estimator=mle|mom] [--cap=N] [--r=N]

suppressMessages({
  library(data.table)
  library(RTIGER)
  library(readr)
  library(tidyverse)
  library(GenomicRanges)
  library(S4Vectors)
})
source("R/06_score_rtiger.R") # segments_from_rtiger_object()
set.seed(42) # reproducible binomial thinning

## ---- args ------------------------------------------------------------------
argv <- commandArgs(trailingOnly = TRUE)
getopt <- function(flag, default) {
  hit <- grepl(paste0("^", flag, "="), argv)
  if (any(hit)) sub(paste0("^", flag, "="), "", argv[hit][1]) else default
}
r <- as.integer(getopt("--r", "8"))
estimator <- tolower(getopt("--estimator", "mle")) # mle (default) | mom
cap_arg <- as.integer(getopt("--cap", "20")) # default 20x; --cap=0 disables
user_cap <- if (is.na(cap_arg) || cap_arg <= 0) NA_integer_ else cap_arg
stopifnot(estimator %in% c("mle", "mom"))
sets <- argv[!grepl("^--", argv) & nzchar(argv)]
if (length(sets) == 0) sets <- "SNP"

## ---- Julia -----------------------------------------------------------------
JULIA_HOME <- Sys.getenv("JULIA_HOME", "")
if (!nzchar(JULIA_HOME) || !file.exists(file.path(JULIA_HOME, "julia"))) {
  launcher <- path.expand("~/.juliaup/bin/julia")
  if (file.exists(launcher)) {
    JULIA_HOME <- system(paste(
      shQuote(launcher), "--startup-file=no -e",
      shQuote("print(Sys.BINDIR)")
    ), intern = TRUE)
  }
}
setupJulia(JULIA_HOME = JULIA_HOME)
sourceJulia() # loads stock (MLE) emissionUpdateState

# Swap the BetaBinomial M-step in Main between MLE (stock) and MoM (override).
# Re-sourcing is cheap and lets one run mix regimes across sets.
.mstep <- "mle"
set_mstep <- function(method) {
  if (method == "mom" && .mstep != "mom") {
    if (!file.exists("mom_mstep.jl")) stop("mom_mstep.jl not found (required for --estimator=mom)")
    JuliaCall::julia_source(normalizePath("mom_mstep.jl"))
    .mstep <<- "mom"
  } else if (method == "mle" && .mstep != "mle") {
    sourceJulia()
    .mstep <<- "mle" # restores stock MLE emissionUpdateState
  }
}

## ---- cap helper ------------------------------------------------------------
cap_counts <- function(ref, alt, cap) { # binomial down-thinning to ~cap x
  tot <- ref + alt
  p <- pmin(1, cap / pmax(tot, 1L))
  list(ref = rbinom(length(ref), ref, p), alt = rbinom(length(alt), alt, p))
}

## ---- shared inputs ---------------------------------------------------------
smap <- fread("data/molbreeding_45k_sample_map.tsv")
sl <- read_csv("data/rtiger_50K/seqlengths.csv", show_col_types = FALSE)
seqv <- setNames(sl$len, sl$chr_label)

fit_set <- function(setname) {
  cat(sprintf("\n=============== %s (rigidity=%d) ===============\n", setname, r))
  gt_path <- sprintf("data/molbreeding_45k/gatk_table_%s_v5.tsv", setname)
  out_root <- file.path(getopt("--outdir", "data/molbreeding_targetseq"), setname)
  cnt_dir <- file.path(out_root, "counts")
  fit_dir <- file.path(out_root, "fit")
  dir.create(cnt_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)

  gt <- fread(gt_path) # SAMPLE CONTIG POSITION REF_NUCLEOTIDE ALT_NUCLEOTIDE REF_COUNT ALT_COUNT TOTAL_COUNT
  mean_depth <- mean(gt$REF_COUNT + gt$ALT_COUNT)

  method <- estimator
  capv <- user_cap
  capped <- !is.na(capv) && mean_depth > capv # cap only bites above target depth
  if (capped) {
    th <- cap_counts(gt$REF_COUNT, gt$ALT_COUNT, capv)
    gt[, `:=`(REF_COUNT = th$ref, ALT_COUNT = th$alt)]
  }
  cat(sprintf(
    "mean depth %.1fx | estimator %s | cap %s%s\n",
    mean_depth, toupper(method),
    if (is.na(capv)) "off" else paste0(capv, "x"),
    if (capped) {
      sprintf(" -> realized %.1fx", mean(gt$REF_COUNT + gt$ALT_COUNT))
    } else if (!is.na(capv)) " (below cap, no-op)" else ""
  ))
  set_mstep(method)

  # well_id -> PN#_SID# and accession (donor)
  si <- match(gt$SAMPLE, smap$well_id)
  gt[, name := ifelse(is.na(si), SAMPLE, smap$sample_id[si])]
  donor_of <- function(w) {
    g <- smap$genotype[match(w, smap$well_id)]
    ifelse(is.na(g), w, sub("_.*$", "", g))
  }

  ed <- data.table(files = character(), name = character(), donor = character())
  for (w in unique(gt$SAMPLE)) {
    sub <- gt[SAMPLE == w]
    six <- data.table(
      chr = sub$CONTIG, pos = sub$POSITION,
      RefBase = sub$REF_NUCLEOTIDE, RefCount = sub$REF_COUNT,
      AltBase = sub$ALT_NUCLEOTIDE, AltCount = sub$ALT_COUNT
    )
    f <- normalizePath(file.path(cnt_dir, paste0(sub$name[1], ".tsv")), mustWork = FALSE)
    fwrite(six, f, sep = "\t", col.names = FALSE)
    ed <- rbind(ed, data.table(files = f, name = sub$name[1], donor = donor_of(w)))
  }
  fwrite(ed, file.path(out_root, "expDesign.csv"))

  ## ---- fit (one joint EM over the samples) ---------------------------------
  prog <- file.path(fit_dir, "fit_progress.log")
  if (file.exists(prog)) file.remove(prog)
  t0 <- Sys.time()
  res <- RTIGER(
    expDesign = as.data.frame(ed[, .(files, name)]),
    outputdir = fit_dir, seqlengths = seqv, rigidity = r,
    autotune = FALSE, progress_log = prog, save.results = FALSE,
    threads = 1L, verbose = TRUE
  )
  saveRDS(res, file.path(fit_dir, "rtiger_result.rds"))
  cat(sprintf(
    "Fit done in %.2f min (%d samples, M-step=%s).\n",
    as.numeric(difftime(Sys.time(), t0, units = "mins")), length(res@Viterbi), toupper(method)
  ))

  ## ---- common schema -------------------------------------------------------
  tag <- sprintf(
    "MolBreeding_%s_r%d_%s%s", setname, r, toupper(method),
    if (capped) paste0("_cap", capv) else ""
  )
  seg <- segments_from_rtiger_object(res)
  seg$donor <- ed$donor[match(seg$name, ed$name)]
  seg$source <- tag
  seg <- seg |>
    select(source, donor, name, chr, start_bp, end_bp, state) |>
    arrange(name, chr, start_bp)
  out_csv <- file.path(out_root, "calls_common_schema.csv")
  write_csv(seg, out_csv)

  l <- seg$end_bp - seg$start_bp + 1
  dose <- seg |>
    mutate(len = l) |>
    group_by(name) |>
    summarise(d = sum(len * state / 2) / sum(len), .groups = "drop")
  cat(sprintf(
    "Wrote %s [%s]: %d samples, %d segments, mean dosage %.3f\n",
    out_csv, tag, n_distinct(seg$name), nrow(seg), mean(dose$d)
  ))
}

for (s in sets) fit_set(s)
cat("\nALL DONE.\n")
