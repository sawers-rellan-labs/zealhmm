#!/usr/bin/env Rscript
# Decode panel A of the two-panel fragment-size cM figure: "calibration against molb truth".
#
# REFERENCE  = molb truth donor fragments: nnil run on the MolBreeding target-seq (~110x) hard calls
#              at the fixed high-quality truth config (this is Jim's chip -> HMM step, NOT a re-call).
# CALLERS    = the four skim callers on the SAME 14 molb NILs, at the panel-B operating points, from
#              the mount-frame counts staged by scripts/stage_zeal_molb_cohort.R.
#              googa/atlas are NOT run here: only 4 of the 14 molb NILs have BrB, and all 4 sit inside
#              the staged 332 paired cohort, so panel A reuses their segments from that decode.
#
# It is a fragment-SIZE DISTRIBUTION, so fragments are POOLED and no per-line pairing is required
# (Fausto, explicit). Series therefore differ in n, and every label states its own n.
#
#   Rscript scripts/zeal_molb_cohort_decode.R
#   ZEAL_MOLB_RECOMPUTE=1 Rscript scripts/zeal_molb_cohort_decode.R

suppressMessages({
  library(devtools)
  library(data.table)
  library(here)
})
devtools::load_all(path.expand("~/repos/nilhmm"), quiet = TRUE)
source(here::here("R/metrics.R"))
source(here::here("scripts/logging.R"))

COH <- here::here("data/zeal/molb_cohort")
OUT <- here::here("results/sim/zeal_nil/molb_panel")
KEEP <- c("name", "chr", "start_bp", "end_bp", "state")
threads <- min(parallel::detectCores() - 2L, 8L)
RECOMPUTE <- Sys.getenv("ZEAL_MOLB_RECOMPUTE") != ""
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

man <- fread(file.path(COH, "molb_cohort_manifest.csv"))
setorder(man, pedigree)
peds <- man$pedigree
log_info("[molb-decode] %d molb calibration NILs", length(peds))

# ---- native TeoNAM v5 Marey spline (bp -> cM), same as panel B ---------------
nmap <- fread(here::here("data/zeal/markers_snp50k_cm.tsv"))[, .(chr = as.integer(chr), bp = as.integer(pos), cm = as.numeric(cm))]
to_cm <- bp_to_cm(as.data.frame(nmap))
db_cm <- function(seg) {
  db <- .donor_blocks(as.data.table(seg)[, ..KEEP])
  if (!nrow(db)) {
    return(data.table(name = character(0), chr = integer(0), start_bp = integer(0), end_bp = integer(0), cm = numeric(0)))
  }
  db <- as.data.table(db)
  db[, cm := to_cm(chr, end_bp) - to_cm(chr, start_bp)]
  db[is.finite(cm) & cm > 0, .(name, chr, start_bp, end_bp, cm)]
}

# ---- operating points -------------------------------------------------------
OPS <- fread(here::here("reference/zeal_skim_calibration_operating_points.csv"))
RRATE_SKIM <- 6.36e-4
RIG_RTIGER <- as.integer(OPS[caller == "rtiger", value])
# molb TRUTH config: the fixed high-quality nnil configuration on the 9,157-site molbreeding panel.
# rrate is panel-specific (2L/(100M) on 9,157 sites, not the SNP50K value).
MOLB_RRATE <- 1.67e-3
MOLB_NIR <- 0.70
MOLB_GERM <- 1e-3
MOLB_GERT <- 1e-4
MOLB_P <- 0.9
log_info(
  "[molb-decode] truth config: nnil on molb hardcalls, rrate=%.2e nir=%.2f germ=%.0e gert=%.0e p=%.2f | skim callers at rrate=%.2e, rtiger rigidity=%d",
  MOLB_RRATE, MOLB_NIR, MOLB_GERM, MOLB_GERT, MOLB_P, RRATE_SKIM, RIG_RTIGER
)

# ---- loaders ----------------------------------------------------------------
load_counts <- function(pp) {
  rbindlist(lapply(pp, function(p) {
    cf <- fread(file.path(COH, "skim/counts_50k", sprintf("%s.tsv", p)),
      header = FALSE, col.names = c("contig", "pos", "rb", "rc", "ab", "ac")
    )
    data.table(
      name = p, chr = as.integer(sub("chr", "", cf$contig)), pos = as.integer(cf$pos),
      n_ref = as.integer(cf$rc), n_alt = as.integer(cf$ac)
    )
  }))
}
load_bins <- function(pp) {
  rbindlist(lapply(pp, function(p) {
    b <- fread(file.path(COH, "skim/bins", sprintf("%s.tsv", p)))
    data.table(
      name = b$SAMPLE, chr = as.integer(sub("chr", "", b$CONTIG)), pos = as.integer(b$BIN_POS),
      alt_freq = as.numeric(b$ALT_FREQ), start_bp = as.integer(b$BIN_START),
      end_bp = as.integer(b$BIN_END), ninf = as.integer(b$INFORMATIVE_VARIANT_COUNT)
    )
  }))
}

# ---- molb TRUTH: nnil on the target-seq hard calls, keyed on pedigree --------
truth_csv <- file.path(OUT, "molb_truth_donor_cm.csv")
if (RECOMPUTE || !file.exists(truth_csv)) {
  hc <- fread(here::here("data/zeal/molbreeding/molbreeding_hardcalls_wsfilt.tsv"))
  log_info("[molb-decode] hardcalls: %d rows, %d samples, %d sites", nrow(hc), uniqueN(hc$name), uniqueN(hc$marker))
  # Key truth on the MolBreeding WELL (`name` = truth_sample) mapped through calibration_pairing to
  # the CORRECTED pedigree. The hardcall table's own `pedigree` column carries the MISLABEL: well
  # PN4_SID330 is labelled Zx.0030_P2_P4_P2.1.1.1 but its DNA is PN4_SID322 = Zx.0030_P2_P2_P1.2.1.1.
  # Joining on that column silently drops the corrected line (13 of 14 instead of 14).
  key <- fread(here::here("data/zeal/correspondence/calibration_pairing.csv"))[
    in_calibration == TRUE, .(name = trimws(truth_sample), ped_corrected = sub("\\.B$", "", true_pedigree))
  ]
  hc <- merge(hc, key, by = "name") # inner join also drops the 2 B73-bulk controls (not NILs)
  mism <- unique(hc[sub("\\.B$", "", pedigree) != ped_corrected, .(name, ped_labelled = pedigree, ped_corrected)])
  if (nrow(mism)) {
    log_info("[molb-decode] %d well(s) relabelled via calibration_pairing:", nrow(mism))
    print(mism)
  }
  log_info("[molb-decode] truth after restricting to the %d calibration NILs: %d samples", length(peds), uniqueN(hc$ped_corrected))
  stopifnot(uniqueN(hc$ped_corrected) == length(peds))
  hc[is.na(g), g := 3L]
  tr <- call_ancestry(
    data.frame(name = hc$ped_corrected, chr = as.integer(hc$chr), pos = as.integer(hc$pos), g = as.integer(hc$g)),
    caller = "nnil", rrate = MOLB_RRATE, nir = MOLB_NIR, germ = MOLB_GERM, gert = MOLB_GERT, p = MOLB_P,
    design = "BC2S3", threads = threads
  )
  d <- db_cm(tr)
  d[, series := "molb_truth"]
  fwrite(d, truth_csv)
  log_info("[molb-decode] molb truth: %d NILs, %d donor blocks, median %.1f cM", uniqueN(d$name), nrow(d), median(d$cm))
} else {
  log_info("[molb-decode] molb truth cached")
}

# ---- skim callers on the same 14 NILs ---------------------------------------
CALLERS <- list(
  nnil = function(pp) {
    sk <- load_counts(pp)
    g <- call_gt(sk$n_ref, sk$n_alt, prior = breeding_prior("BC2S3"), error = 0.01, return = "call")
    g[is.na(g)] <- 3L
    call_ancestry(data.frame(name = sk$name, chr = sk$chr, pos = sk$pos, g = as.integer(g)),
      caller = "nnil", rrate = RRATE_SKIM, nir = 0.70, germ = 1e-4, gert = 1e-2, p = 0.1,
      design = "BC2S3", threads = threads
    )
  },
  bbnil = function(pp) {
    call_ancestry(as.data.frame(load_counts(pp)),
      caller = "bbnil", rrate = RRATE_SKIM, fit_means = TRUE, conc = 20, err = 0.01,
      design = "BC2S3", parallel = TRUE, threads = threads
    )
  },
  rtiger = function(pp) {
    call_ancestry(as.data.frame(load_counts(pp)),
      caller = "rtiger", rigidity = RIG_RTIGER, design = "BC2S3", threads = threads
    )
  },
  binhmm = function(pp) {
    call_ancestry(as.data.frame(load_bins(pp)), caller = "binhmm", design = "BC2S3")
  }
)

for (cl in names(CALLERS)) {
  f <- file.path(OUT, sprintf("%s_donor_cm.csv", cl))
  if (!RECOMPUTE && file.exists(f)) {
    log_info("[molb-decode] %s cached", cl)
    next
  }
  t0 <- Sys.time()
  s <- tryCatch(as.data.table(CALLERS[[cl]](peds))[, ..KEEP],
    error = function(e) {
      # same QC discipline as panel B: retry per NIL so one thin line does not drop the cohort,
      # and name the lines the caller's internal QC refuses (that refusal is a real finding).
      log_warn("[molb-decode] %s aborted (%s); retrying per NIL", cl, sub("\n.*$", "", conditionMessage(e)))
      ok <- Filter(Negate(is.null), lapply(peds, function(p) {
        tryCatch(as.data.table(CALLERS[[cl]](p))[, ..KEEP],
          error = function(e2) {
            log_warn("[molb-decode] %s DROPPED %s: %s", cl, p, sub("\n.*$", "", conditionMessage(e2)))
            NULL
          }
        )
      }))
      if (!length(ok)) NULL else rbindlist(ok)
    }
  )
  if (is.null(s)) {
    log_warn("[molb-decode] %s produced nothing", cl)
    next
  }
  d <- db_cm(s)
  d[, series := cl]
  fwrite(d, f)
  log_info(
    "[molb-decode] %s: %d NILs, %d donor blocks, median %.1f cM (%.0fs)",
    cl, uniqueN(d$name), nrow(d), median(d$cm), as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )
}
log_info("[molb-decode] done; outputs in %s", OUT)
