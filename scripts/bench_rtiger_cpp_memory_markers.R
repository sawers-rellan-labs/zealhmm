#!/usr/bin/env Rscript
# C++ RTIGER peak-RSS across the Arabidopsis shared-3 marker sweep (same panel and
# odd-index thinning as the time figure), so memory is reported on the SAME axis as
# rtiger_marker_scaling_cpp.png. Each size fits in its own process under
# /usr/bin/time -l.
#   Rscript scripts/bench_rtiger_cpp_memory_markers.R            # driver
#   Rscript scripts/bench_rtiger_cpp_memory_markers.R --worker L # one size (internal)
# Writes results/bench/rtiger_memory_markers.csv.
ROOT <- "/Users/fvrodriguez/repos/zealhmm"
PANEL <- file.path(ROOT, "data/rtiger_shared3_input")
OUTDIR <- file.path(ROOT, "results/bench")
SCRIPT <- file.path(ROOT, "scripts/bench_rtiger_cpp_memory_markers.R")
files <- c(S1 = "sampleBN.txt", S2 = "sampleZ.txt", S3 = "sampleAU.txt")
chrs <- paste0("Chr", 1:5)
LEVELS <- 0:4 # 109,703 -> 6,857 markers/sample
RIG <- as.integer(Sys.getenv("RTIGER_RIG", "2")) # operating rigidity (env-set)
OUTCSV <- if (RIG == 2L) "rtiger_memory_markers.csv" else sprintf("rtiger_memory_markers_r%d.csv", RIG)
args <- commandArgs(trailingOnly = TRUE)

# Build obs from a PRE-THINNED panel directory (thin_L<level>/, written by
# scripts/materialize_thinned_panel.R). The measured worker reads ONLY the thinned
# set, so its peak RSS is not floored by a full-panel read.
build_obs_from <- function(dir) {
  raw <- lapply(files, function(f) {
    read.table(file.path(dir, f),
      sep = "\t",
      col.names = c("chr", "pos", "refA", "refC", "altA", "altF"), stringsAsFactors = FALSE
    )
  })
  obs <- list()
  for (s in names(files)) {
    d <- raw[[s]]
    obs[[s]] <- list()
    for (c in chrs) {
      rows <- which(d$chr == c)
      if (!length(rows)) next
      k <- as.integer(round(d$refC[rows]))
      n <- k + as.integer(round(d$altF[rows]))
      obs[[s]][[c]] <- list(k = k, n = n)
    }
  }
  obs
}

if ("--worker" %in% args) { # fit at one marker level, quit
  suppressMessages(library(nilHMM))
  level <- as.integer(args[which(args == "--worker") + 1L])
  obs <- build_obs_from(file.path(PANEL, sprintf("thin_L%d", level)))
  mps <- sum(vapply(obs$S1, function(x) length(x$k), 0L))
  invisible(nilHMM:::.rtiger_fit(obs,
    r = RIG, nstates = 3L, eps = 0.01, max_iter = 6L,
    init_alpha = c(20, 20, 1), init_beta = c(1, 20, 20)
  ))
  cat(sprintf("MEMMARK_DONE level=%d markers=%d\n", level, mps))
  quit(save = "no")
}

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
rows <- list()
for (lv in LEVELS) {
  out <- suppressWarnings(system2("/usr/bin/time",
    c("-l", "Rscript", SCRIPT, "--worker", lv),
    stdout = TRUE, stderr = TRUE
  ))
  mk <- as.numeric(sub(".*markers=([0-9]+).*", "\\1", grep("MEMMARK_DONE", out, value = TRUE)))
  by <- as.numeric(sub("^\\s*([0-9]+).*", "\\1", grep("maximum resident set size", out, value = TRUE)))
  rows[[length(rows) + 1L]] <- data.frame(markers = mk, peak_rss_mib = round(by / 1024^2, 1))
  cat(sprintf("[mem-markers] %6.0f markers  %.0f MiB\n", mk, by / 1024^2))
}
d <- do.call(rbind, rows)
write.csv(d, file.path(OUTDIR, OUTCSV), row.names = FALSE)
cat(sprintf("wrote %s (rigidity=%d)\n", OUTCSV, RIG))
