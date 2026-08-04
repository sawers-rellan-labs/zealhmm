#!/usr/bin/env Rscript
# Decode the six paint callers across the staged 332-NIL ZEAL skim-x-BrB paired cohort, at the
# operating points recorded in reference/zeal_skim_calibration_operating_points.csv, and write each
# caller's donor-introgression sizes in cM (native TeoNAM v5 Marey spline) for panel 1 of the
# two-panel fragment-size validation figure.
#
# Panel 1 replaces the old 11-line vary_skim version: the FULL paired cohort compared against the
# BC2S3 sim latent donor-tract cM law (no gamma fit). See [[zeal-paired-cohort-332]].
#
# INCREMENTAL + RESUMABLE, by (caller, chunk): each chunk's segments are cached to RDS and its
# donor-block cM rows appended to the caller's CSV, so an interrupt costs at most one chunk.
# Inputs are the pedigree-named staged files from scripts/stage_zeal_paired_cohort.R, so the
# series `name` taken from a basename IS the NIL ID (never a platform-local PN#_SID# prefix).
#
#   Rscript scripts/zeal_paired_cohort_decode.R
#   ZEAL_DECODE_CALLERS=nnil,bbnil Rscript scripts/zeal_paired_cohort_decode.R  # subset
#   ZEAL_DECODE_RECOMPUTE=1 Rscript scripts/zeal_paired_cohort_decode.R         # clear caches

suppressMessages({
  library(devtools)
  library(data.table)
  library(here)
})
devtools::load_all(path.expand("~/repos/nilhmm"), quiet = TRUE)
source(here::here("R/metrics.R"))
source(here::here("scripts/logging.R"))

COH <- here::here("data/zeal/paired_cohort")
OUT <- here::here("results/sim/zeal_nil/paired_cohort")
KEEP <- c("name", "chr", "start_bp", "end_bp", "state")
CHUNK <- 20L
threads <- min(parallel::detectCores() - 2L, 8L)
RECOMPUTE <- Sys.getenv("ZEAL_DECODE_RECOMPUTE") != ""
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

man <- fread(file.path(COH, "cohort_manifest.csv"))
setorder(man, pedigree)
peds <- man$pedigree
chunks <- split(peds, ceiling(seq_along(peds) / CHUNK))
log_info("[decode] %d NILs in %d chunks of <=%d | %d threads", length(peds), length(chunks), CHUNK, threads)

# ---- native TeoNAM v5 Marey spline (bp -> cM) --------------------------------
nmap <- fread(here::here("data/zeal/markers_snp50k_cm.tsv"))[, .(chr = as.integer(chr), bp = as.integer(pos), cm = as.numeric(cm))]
to_cm <- bp_to_cm(as.data.frame(nmap))
db_cm <- function(seg) {
  db <- .donor_blocks(as.data.table(seg)[, ..KEEP])
  if (!nrow(db)) {
    return(data.table(name = character(0), cm = numeric(0)))
  }
  db <- as.data.table(db)
  db[, cm := to_cm(chr, end_bp) - to_cm(chr, start_bp)]
  db[is.finite(cm) & cm > 0, .(name, chr, start_bp, end_bp, cm)]
}

# ---- per-chunk input loaders (pedigree-named staged files) -------------------
STORE <- here::here("data/zeal/snp50k_counts") # unified GATK store (canonical id = pedigree for NILs)
load_counts <- function(pp) {
  rbindlist(lapply(pp, function(p) {
    cf <- fread(file.path(STORE, sprintf("%s.tsv", p)), # was data/zeal/paired_cohort/skim/counts_50k (bcftools)
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

# BrB pangene: positions are a fixed property of the pangene <-> B73 gene map, not of the samples
gtp <- fread(here::here("data/ref/gene_to_pangene.tsv"))
gc <- fread(here::here("data/ref/b73_gene_coords.tsv"), header = FALSE, col.names = c("gene", "chr", "start", "end"))
gc[, chrn := as.integer(sub("chr", "", chr))]
pg <- merge(gtp[species == "B73", .(gene, pan_gene)], gc[!is.na(chrn), .(gene, chrn, start)], by = "gene")
pg <- pg[pan_gene %in% pg[, .(nchr = uniqueN(chrn)), by = pan_gene][nchr == 1L, pan_gene]]
pgpos <- pg[, .(chr = as.integer(chrn[1]), pos = as.integer(stats::median(start))), by = pan_gene]
load_pangene <- function(pp) {
  brb <- rbindlist(lapply(pp, function(p) {
    fread(file.path(COH, "brb/pangene", sprintf("%s.pangene_counts.tsv", p)))[, .(name = sample, pan_gene, n_recur, n_donor)]
  }))
  b <- merge(brb, pgpos, by = "pan_gene")
  setorder(b, name, chr, pos)
  b[, .(name, chr, pos, n_ref = as.integer(n_recur), n_alt = as.integer(n_donor))]
}

# ---- operating points -------------------------------------------------------
# skim rrate = 2L/(100M) on the SNP50K grid; googa rrate = the pangene-grid analog. Both are the
# PINNED map-derived values recorded in the operating-points CSV, not swept here.
OPS <- fread(here::here("reference/zeal_skim_calibration_operating_points.csv"))
RRATE_SKIM <- 6.36e-4
RIG_ATLAS <- as.integer(OPS[caller == "atlas", value])
RIG_RTIGER <- as.integer(OPS[caller == "rtiger", value])
STAY_BINHMM <- as.numeric(OPS[caller == "binhmm", value])
ATLAS_THRESH <- 0.95
ATLAS_HET <- 0.25
ATLAS_MIN_READS <- 5L
GERM_BRB <- 0.05
GERT_BRB <- 0.10
P_BRB <- 0.5
MR_BRB <- 0.10
# googa rrate is PINNED TO THE MAP, so it must be RECOMPUTED for the grid actually decoded rather
# than inherited: 2L/(100M) depends on M. The recorded 1.32e-3 operating point was pinned to the
# ~23K-pan_gene grid covered by the 4 BrB samples staged at calibration time; this cohort traverses
# a larger grid, so that value is pinned to a grid we no longer decode. Fausto's call: use the
# recomputed value. Interpolate the native TeoNAM v5 cM onto the pangene positions, Holland 2L/(100M).
gr <- copy(pgpos)[, cm := as.numeric(NA)]
for (cc in unique(gr$chr)) {
  m <- nmap[chr == cc][order(bp)]
  i <- gr$chr == cc
  gr$cm[i] <- stats::approx(m$bp, m$cm, xout = gr$pos[i], rule = 2)$y
}
L_PANGENE <- sum(vapply(split(gr$cm, gr$chr), function(x) max(x) - min(x), numeric(1)))
RRATE_GOOGA <- 2 * L_PANGENE / (100 * nrow(gr))
log_info(
  "[decode] pangene grid %d markers, %.0f cM -> map-pinned googa rrate = %.3e (recorded operating point %.3e was pinned to the smaller ~23K calibration grid; NOT used)",
  nrow(gr), L_PANGENE, RRATE_GOOGA, as.numeric(OPS[caller == "googa", value])
)
log_info(
  "[decode] operating points: skim rrate=%.3e | googa rrate=%.3e nir=0.20 | atlas rigidity=%d nir=0.10 | rtiger rigidity=%d | binhmm stay=%.3f",
  RRATE_SKIM, RRATE_GOOGA, RIG_ATLAS, RIG_RTIGER, STAY_BINHMM
)

# ---- caller definitions: chunk of pedigrees -> segments ---------------------
CALLERS <- list(
  nnil = function(pp) { # categorical: hard-call first via the BC2S3 DESIGN-prior MAP
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
  },
  googa = function(pp) {
    call_ancestry(as.data.frame(load_pangene(pp)),
      caller = "googa", rrate = RRATE_GOOGA, nir = 0.20,
      germ = GERM_BRB, gert = GERT_BRB, p = P_BRB, mr = MR_BRB,
      atlas_thresh = ATLAS_THRESH, atlas_het = ATLAS_HET, atlas_min_reads = ATLAS_MIN_READS,
      design = "BC2S3", threads = threads
    )
  },
  atlas = function(pp) {
    call_ancestry(as.data.frame(load_pangene(pp)),
      caller = "atlas", rigidity = RIG_ATLAS, nir = 0.10,
      germ = GERM_BRB, gert = GERT_BRB, p = P_BRB, mr = MR_BRB,
      atlas_thresh = ATLAS_THRESH, atlas_het = ATLAS_HET, atlas_min_reads = ATLAS_MIN_READS,
      design = "BC2S3", threads = threads
    )
  }
)

want <- Sys.getenv("ZEAL_DECODE_CALLERS")
if (want != "") CALLERS <- CALLERS[trimws(strsplit(want, ",")[[1]])]
log_info("[decode] callers: %s", paste(names(CALLERS), collapse = ", "))

# ---- decode, chunk by chunk, appending as we go -----------------------------
for (cl in names(CALLERS)) {
  seg_dir <- file.path(OUT, sprintf("%s_seg_cache", cl))
  dir.create(seg_dir, showWarnings = FALSE, recursive = TRUE)
  cm_csv <- file.path(OUT, sprintf("%s_donor_cm.csv", cl))
  if (RECOMPUTE) {
    unlink(list.files(seg_dir, full.names = TRUE))
    unlink(cm_csv)
  }
  t0 <- Sys.time()
  n_this_run <- 0L
  for (k in seq_along(chunks)) {
    f <- file.path(seg_dir, sprintf("chunk_%03d.rds", k))
    if (file.exists(f)) next
    pp <- chunks[[k]]
    s <- tryCatch(as.data.table(CALLERS[[cl]](pp))[, ..KEEP],
      error = function(e) {
        # A whole-chunk abort must not cost every NIL in the chunk: rtiger refuses a chunk outright
        # when ANY (sample, chr) chain has < 2*rigidity covered markers. Retry sample by sample so
        # only the genuinely undecodable lines are lost, and say which ones.
        log_warn("[decode] %s chunk %d aborted (%s); retrying per NIL", cl, k, sub("\n.*$", "", conditionMessage(e)))
        one <- lapply(pp, function(p) {
          tryCatch(as.data.table(CALLERS[[cl]](p))[, ..KEEP],
            error = function(e2) {
              log_warn("[decode] %s DROPPED %s: %s", cl, p, sub("\n.*$", "", conditionMessage(e2)))
              NULL
            }
          )
        })
        ok <- Filter(Negate(is.null), one)
        log_info("[decode] %s chunk %d recovered %d/%d NILs per-NIL", cl, k, length(ok), length(pp))
        if (!length(ok)) NULL else rbindlist(ok)
      }
    )
    if (is.null(s)) next
    saveRDS(s, f)
    d <- db_cm(s)
    d[, `:=`(caller = cl, chunk = k)]
    fwrite(d, cm_csv, append = file.exists(cm_csv))
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    n_this_run <- n_this_run + 1L # chunks decoded in THIS run (cached ones cost no time)
    log_info(
      ">>> %s chunk %d/%d (%d NILs, %d donor blocks) | elapsed %.1f min | avg %.2f min/chunk | ETA ~%.1f min",
      cl, k, length(chunks), length(pp), nrow(d), el, el / n_this_run,
      (el / n_this_run) * (length(chunks) - k)
    )
  }
  nseg <- length(list.files(seg_dir, pattern = "\\.rds$"))
  if (file.exists(cm_csv)) {
    tot <- fread(cm_csv)
    log_info(
      "[decode] %s DONE: %d/%d chunks | %d NILs | %d donor blocks | median %.1f cM",
      cl, nseg, length(chunks), uniqueN(tot$name), nrow(tot), median(tot$cm)
    )
  } else {
    log_warn("[decode] %s produced no output", cl)
  }
}
log_info("[decode] all requested callers finished; outputs in %s", OUT)
