# TeoNAM marker liftover AGPv2 -> v4 -> v5 (two hops; no direct v2->v5 chain).
# Adapts Desktop/zealtiger/molbreeding_liftover.R. Chains live in
# data/ref/chain_files/ and MUST be tab-delimited for rtracklayer::import.chain
# (verified). Functions only; a script owns the paths. See
# agent/teonam-qtl-recovery-plan.md §4.2.

#' One liftOver hop keeping only unique 1:1 mappings
#'
#' @param gr GRanges (mcol `marker` carried through).
#' @param chain_path Path to a tab-delimited UCSC chain file.
#' @return GRanges of inputs that mapped to exactly one range (marker preserved).
#' @export
lift_unique <- function(gr, chain_path) {
  stopifnot(file.exists(chain_path))
  ch <- rtracklayer::import.chain(chain_path)
  lifted <- rtracklayer::liftOver(gr, ch) # GRangesList parallel to gr
  n1 <- S4Vectors::elementNROWS(lifted) == 1L # exactly one target range
  out <- unlist(lifted[n1])
  out$marker <- gr$marker[n1]
  out
}

#' Liftover TeoNAM AGPv2 markers to B73 NAM v5 (v2 -> v4 -> v5)
#'
#' @param marker_info_path `markers_v2.csv` (chromosome, name = `S<chr>_<pos>`
#'   in AGPv2, start, end).
#' @param chain_dir Directory holding the two tab-delimited chains.
#' @return data.frame(marker, chr_v2, pos_v2, chr_v5, pos_v5) for markers that
#'   lift 1:1 through both hops, stay on chr 1-10, and do not change chromosome.
#' @export
liftover_teonam <- function(marker_info_path,
                            chain_dir = "data/ref/chain_files") {
  mi <- data.table::fread(marker_info_path, data.table = FALSE)
  # AGPv2 position lives in the marker name S<chr>_<pos>; chromosome col is 1..10
  pos_v2 <- as.integer(sub("^S[0-9]+_", "", mi$name))
  chr_v2 <- as.character(mi$chromosome)
  ok <- !is.na(pos_v2) & chr_v2 %in% as.character(1:10)
  gr <- GenomicRanges::GRanges(
    seqnames = chr_v2[ok],
    ranges   = IRanges::IRanges(start = pos_v2[ok], width = 1L),
    marker   = mi$name[ok]
  )
  n_in <- length(gr)

  gr_v4 <- lift_unique(gr, file.path(chain_dir, "AGPv2_to_B73_RefGen_v4.chain"))
  gr_v5 <- lift_unique(
    gr_v4,
    file.path(chain_dir, "B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain")
  )

  res <- data.frame(
    marker = gr_v5$marker,
    chr_v5 = as.character(GenomicRanges::seqnames(gr_v5)),
    pos_v5 = GenomicRanges::start(gr_v5),
    stringsAsFactors = FALSE
  )
  # attach original v2 coords; drop chromosome-changers; keep chr 1-10
  src <- data.frame(
    marker = mi$name[ok], chr_v2 = chr_v2[ok], pos_v2 = pos_v2[ok],
    stringsAsFactors = FALSE
  )
  res <- merge(src, res, by = "marker")
  res <- res[res$chr_v5 %in% as.character(1:10) & res$chr_v5 == res$chr_v2, ]
  res <- res[
    order(as.integer(res$chr_v5), res$pos_v5),
    c("marker", "chr_v2", "pos_v2", "chr_v5", "pos_v5")
  ]

  attr(res, "n_in") <- n_in
  res
}
