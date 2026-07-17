# TeoNAM QTL-recovery helpers: ID reconciliation, (later) cM attach, densification,
# QTL mapping, recovery scoring. Functions only — scripts/notes own the paths.
# See agent/teonam-qtl-recovery-plan.md.

#' Normalise a TeoNAM genotype ID to (family, line)
#'
#' Genotype IDs look like `"W22 x TIL01 BC1S4# Line_A001"`.
#' @param id Character vector of genotype IDs.
#' @return data.frame(family, line, key) where key = paste0(family, line).
#' @export
teonam_geno_key <- function(id) {
  family <- sub(".*\\bTIL([0-9]+)\\b.*", "TIL\\1", id)
  line <- sub(".*Line_([A-Za-z]?[0-9]+).*", "\\1", id)
  data.frame(
    family = family, line = line,
    key = paste0(family, line), stringsAsFactors = FALSE
  )
}

#' Normalise a TeoNAM phenotype Line code to (family, line)
#'
#' Phenotype Line codes look like `"TIL01A001"`.
#' @param line Character vector of phenotype Line codes.
#' @return data.frame(family, line, key).
#' @export
teonam_pheno_key <- function(line) {
  family <- sub("^(TIL[0-9]+)([A-Za-z][0-9]+)$", "\\1", line)
  ln <- sub("^(TIL[0-9]+)([A-Za-z][0-9]+)$", "\\2", line)
  data.frame(
    family = family, line = ln,
    key = paste0(family, ln), stringsAsFactors = FALSE
  )
}

#' Reconcile TeoNAM genotypes and phenotypes to the shared canonical RIL set
#'
#' Inner-joins the (deduped, cleaned) genotype IDs to the 1,257-RIL phenotype
#' table on the normalised `(family, line)` key. Reports the join and the
#' phenotyped RILs that lack genotypes in this (EasiGP) copy.
#'
#' @param geno_path Path to the cleaned genotype CSV (ID, population, factor,
#'   markers...). Only the ID + population columns are read.
#' @param pheno_path Path to the 1,257-RIL phenotype xlsx.
#' @return list(matched = data.frame(key, family, line, geno_id, pheno_line,
#'   population, <traits>), n_geno, n_pheno, n_matched, missing_geno = data.frame).
#' @export
reconcile_teonam <- function(geno_path, pheno_path) {
  stopifnot(file.exists(geno_path), file.exists(pheno_path))
  gid <- data.table::fread(geno_path,
    select = 1:2,
    header = TRUE, data.table = FALSE
  )
  names(gid)[1:2] <- c("geno_id", "population")
  gk <- teonam_geno_key(gid$geno_id)
  geno <- cbind(gk, gid) # key, family, line, geno_id, population

  ph <- as.data.frame(readxl::read_excel(pheno_path))
  names(ph)[1:2] <- c("pheno_line", "pheno_pop")
  pk <- teonam_pheno_key(ph$pheno_line)
  pheno <- cbind(pk, ph) # key, family, line, pheno_line, ...

  matched <- merge(geno[, c("key", "family", "line", "geno_id", "population")],
    pheno[, setdiff(names(pheno), c("family", "line"))],
    by = "key"
  )
  matched <- matched[order(matched$family, matched$line), ]

  missing_geno <- pheno[
    !(pheno$key %in% geno$key),
    c("key", "family", "line", "pheno_line")
  ]
  list(
    matched = matched,
    n_geno = nrow(geno), n_pheno = nrow(pheno), n_matched = nrow(matched),
    missing_geno = missing_geno[order(missing_geno$family, missing_geno$line), ]
  )
}

#' Attach consensus-map cM to lifted v5 markers (interpolation, not a join)
#'
#' Evaluates the consensus-map Marey spline at each marker's `pos_v5`
#' (nilHMM::bp_to_cm, Hyman monotone, clamped). The consensus map
#' (Ed Coe 2008 composite "Genetic", v5-anchored) defines the spline; markers are
#' query points — NOT a coordinate join. Requires `R/simulate.R` sourced.
#'
#' @param v5 data.frame with `chr_v5`, `pos_v5` (e.g. from `liftover_teonam()`).
#' @param map Consensus map data.table (chr, cm, bp); e.g. `load_consensus_map()`.
#' @return `v5` with an added `cm` column.
#' @export
attach_cm_v5 <- function(v5, map) {
  map <- data.table::as.data.table(map)
  to_cm <- nilHMM::bp_to_cm(map) # monotone bp -> cM Marey spline (nilHMM::bp_to_cm, R/map.R)
  v5$cm <- NA_real_
  for (ch in unique(v5$chr_v5)) {
    idx <- v5$chr_v5 == ch
    v5$cm[idx] <- to_cm(as.integer(ch), v5$pos_v5[idx])
  }
  v5
}
