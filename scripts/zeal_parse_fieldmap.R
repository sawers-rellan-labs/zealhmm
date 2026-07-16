#!/usr/bin/env Rscript
# ZEAL Phase 2a — reconstruct plot Row/Range from the field-layout maps.
#
# Two field books, two map formats (both are 2-D field grids; body cells = plot ids
# laid serpentine, borders = "B", gaps = "x", serpentine markers = arrows/</>):
#
#   CLY25-B5  (data/zeal/CLY25.xlsx)      "label" mode: range labels in column 1,
#             field column index in the header row; interior range-gutter columns.
#   CLY23-D4  (data/zeal/RR-23-Fields.xlsx, sheet D4-Map_eval_3_blocks_Grid)
#             "grid" mode: a clean serpentine grid, no labels -> position IS the
#             coordinate (physical row = range, physical col = row/column axis).
#
# cf. the canonical target Manifest_map_template.xlsx (Plot-ID, Genotype, row, column).
# Output per field: data/zeal/<field>_fieldmap.csv (plot_id, block, range, col),
# validated against that field's phenotype sheet (100% of phenotyped plots must map).

suppressMessages({
  library(here)
  library(data.table)
  library(readxl)
})
source(here("scripts/logging.R"))

FIELDS <- list(
  cly25_b5 = list(
    xlsx = here("data/zeal/CLY25.xlsx"),
    sheets = c("CLY25-B5-blocks-1-4-map", "CLY25-B5-block-5-map"),
    mode = "label", plot_min = 1L, plot_max = 5588L,
    out = here("data/zeal/cly25_b5_fieldmap.csv"),
    pheno_xlsx = here("data/zeal/CLY25-Fieldbook.xlsx"), pheno_sheet = "B5_BZea_eval"
  ),
  cly23_d4 = list(
    xlsx = here("data/zeal/RR-23-Fields.xlsx"),
    sheets = "D4-Map_eval_3_blocks_Grid",
    mode = "grid", plot_min = 3001L, plot_max = 7432L, # full 3-block D4 eval (3001-7432)
    out = here("data/zeal/cly23_d4_fieldmap.csv"),
    pheno_xlsx = here("data/zeal/CLY23_D4_FieldBook.xlsx"), pheno_sheet = "UPDATED_CLY23_D4_FieldBook"
  )
)

read_grid <- function(xlsx, sheet) {
  m <- suppressMessages(as.matrix(read_excel(xlsx, sheet = sheet, col_names = FALSE)))
  num <- suppressWarnings(matrix(as.numeric(m), nrow = nrow(m)))
  list(m = m, num = num)
}

# "label" mode: range axis in column 1, plots in the body, drop gutters/headers ----
melt_label <- function(xlsx, sheet, pmin, pmax) {
  g <- read_grid(xlsx, sheet)
  is_plot <- !is.na(g$num) & g$num == floor(g$num) & g$num >= pmin & g$num <= pmax
  is_plot[, 1:2] <- FALSE # range-axis + "B" border columns are never plots
  rng <- suppressWarnings(as.integer(g$m[, 1]))
  is_plot[is.na(rng), ] <- FALSE # header rows (col index, arrows) have empty column 1
  for (j in which(colSums(is_plot) > 0)) { # drop range-gutter cols (value == range)
    cells <- which(is_plot[, j])
    if (mean(g$num[cells, j] == rng[cells]) >= 0.5) is_plot[, j] <- FALSE
  }
  plot_cols <- which(colSums(is_plot) > 0)
  cidx <- setNames(seq_along(plot_cols), plot_cols)
  idx <- which(is_plot, arr.ind = TRUE)
  data.table(
    plot_id = as.integer(g$num[is_plot]), block = sheet,
    range = rng[idx[, "row"]], col = unname(cidx[as.character(idx[, "col"])])
  )[!is.na(range)]
}

# "grid" mode: clean serpentine grid, physical position = field coordinate ----------
melt_grid <- function(xlsx, sheet, pmin, pmax) {
  g <- read_grid(xlsx, sheet)
  is_plot <- !is.na(g$num) & g$num == floor(g$num) & g$num >= pmin & g$num <= pmax
  idx <- which(is_plot, arr.ind = TRUE)
  # range = physical row (top->bottom), col = physical column (left->right)
  data.table(
    plot_id = as.integer(g$num[is_plot]), block = sheet,
    range = idx[, "row"], col = idx[, "col"]
  )
}

build_field <- function(name, cfg) {
  log_info("=== %s (%s mode) ===", name, cfg$mode)
  fm <- rbindlist(lapply(cfg$sheets, function(s) {
    d <- if (cfg$mode == "label") melt_label(cfg$xlsx, s, cfg$plot_min, cfg$plot_max) else melt_grid(cfg$xlsx, s, cfg$plot_min, cfg$plot_max)
    log_info("  %-26s -> %4d plots | range %d..%d | col 1..%d", s, nrow(d), min(d$range), max(d$range), max(d$col))
    d
  }))
  dup <- fm[, .N, by = plot_id][N > 1]
  if (nrow(dup)) log_warn("  %d plot ids in >1 cell", nrow(dup))
  npos <- uniqueN(fm[, .(block, range, col)])

  ph <- as.data.table(read_excel(cfg$pheno_xlsx, sheet = cfg$pheno_sheet))
  setnames(ph, 1, "plot_id")
  ph[, plot_id := suppressWarnings(as.integer(plot_id))]
  ph <- ph[!is.na(plot_id) & plot_id >= cfg$plot_min & plot_id <= cfg$plot_max]
  covered <- length(intersect(ph$plot_id, fm$plot_id))
  log_info(
    "  parsed=%d unique=%d positions=%d | pheno(in map id-range)=%d covered=%d (%.1f%%)",
    nrow(fm), uniqueN(fm$plot_id), npos, nrow(ph), covered, 100 * covered / max(1, nrow(ph))
  )
  miss <- setdiff(ph$plot_id, fm$plot_id)
  if (length(miss)) log_warn("  %d phenotyped plots lack a coord (e.g. %s)", length(miss), paste(head(sort(miss)), collapse = ","))

  setorder(fm, block, range, col)
  fwrite(fm, cfg$out)
  log_info("  wrote %s (%d rows)", cfg$out, nrow(fm))
}

for (nm in names(FIELDS)) build_field(nm, FIELDS[[nm]])
