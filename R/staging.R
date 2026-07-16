# Shared staging helpers: mount detection, file lookup, common call schema.
#
# Functions only — analysis notes source this file and own the paths/sample
# lists. Mirrors the mount-drift handling in the zealtiger sanity paint
# (`rsstu-mount-drift`): mount points differ across machines, so resolve each
# root to the first path that exists.

#' First existing path among candidates (mount-drift guard)
#'
#' @param ... Candidate absolute paths, most-preferred first.
#' @return The first path that exists, or the first candidate if none exist
#'   (so downstream `file.exists()` guards fail loudly rather than erroring here).
#' @export
pick <- function(...) {
  cands <- c(...)
  for (p in cands) {
    if (file.exists(p)) {
      return(p)
    }
  }
  cands[[1]]
}

#' Resolve the standard ZEAL mount roots
#'
#' Returns the cassini (`CAS`) and ancestry (`MNT`) roots, each resolved through
#' [pick()] across the known rsstu / BZea-share / local-volume layouts. Notes
#' that touch mount data should guard their chunks on `file.exists()` so the
#' site still renders off-mount.
#'
#' @return A named list with `cassini` and `ancestry` paths.
#' @export
bzea_mounts <- function() {
  list(
    cassini = pick(
      "/Volumes/rsstu/users/r/rrellan/tlaloc/cassini",
      "/Volumes/tlaloc/cassini", "/Volumes/cassini"
    ),
    ancestry = pick(
      "/Volumes/rsstu/users/r/rrellan/BZea/bzeaseq/ancestry",
      "/Volumes/BZea/bzeaseq/ancestry"
    )
  )
}

#' Find a single `<sample>.tsv` count file under a directory (recursive)
#'
#' @param dir Directory to search.
#' @param sample Sample id (matched as `^<sample>\\.tsv$`).
#' @return The first matching path, or `NA_character_` if none.
#' @export
find1 <- function(dir, sample) {
  f <- list.files(dir, paste0("^", sample, "\\.tsv$"),
    full.names = TRUE, recursive = TRUE
  )
  if (length(f)) f[1] else NA_character_
}

#' The common ancestry-call schema used across all analyses
#'
#' Every caller/source is staged to one tidy segment table so they can be
#' compared directly. State is integer-coded REF = 0, HET = 1, ALT = 2.
#'
#' @param dt A `data.table`/`data.frame` of segments.
#' @param source,name,donor Track/sample/donor labels (recycled).
#' @return A `data.table` with columns
#'   `source, donor, name, chr, start_bp, end_bp, state`.
#' @export
as_common_schema <- function(dt, source, name, donor = NA_character_) {
  dt <- data.table::as.data.table(dt)
  data.table::data.table(
    source = source, donor = donor, name = name,
    chr = as.integer(dt$chr),
    start_bp = as.numeric(dt$start_bp),
    end_bp = as.numeric(dt$end_bp),
    state = as.integer(dt$state)
  )[!is.na(chr) & !is.na(state)][order(chr, start_bp)]
}
