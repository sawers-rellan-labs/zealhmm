#!/usr/bin/env Rscript
# ZEAL — Kinki (zigzag-culm / kinked-stem severity) phenotype, empirical logit.
# Kinki is a CLY23-ONLY ordinal 0/1/2 severity score for the classic maize "zigzag culm"
# phenotype (dwarf, shortened/thickened arching ear-region internodes; Eyster 1920,
# J Hered 11:349). Like StPi/StPu it is not continuous, so no SpATS: the score is binarized
# to "kinked" (>=1) and modeled as the Haldane-Anscombe empirical logit of the per-genotype
# kinked-plot proportion  log((k+0.5)/(n-k+0.5)). CLY23 only (absent in CLY25), single field.
# Candidate gene(s) TBD (under investigation from the zigzag-culm literature).
# Outputs: data/zeal/pheno_kinki_elogit.csv (Genotype, k, n, prop, Kinki_mean)
#          data/zeal/tassel/pheno_kinki_all.txt (TASSEL: Taxa | Kinki | Family=taxon)
suppressMessages({
  library(here)
  library(data.table)
  library(readxl)
})
source(here("scripts/logging.R"))
canon_ped <- function(x) sub("\\.B$", "", x)

manifest_cly23 <- function() {
  fm <- fread(here("data/zeal/cly23_d4_fieldmap.csv"))
  ph <- as.data.table(read_excel(here("data/zeal/CLY23_D4_FieldBook.xlsx"), sheet = "UPDATED_CLY23_D4_FieldBook"))
  setnames(ph, "CLY23_D4", "plot_id")
  ph[, plot_id := suppressWarnings(as.integer(plot_id))]
  gc <- as.data.table(read_excel(here("data/zeal/CLY23_D4_FieldBook.xlsx"), sheet = "GENOTYPE-CONVERSION"))
  o2n <- unique(rbind(
    gc[!is.na(oldold_genotype), .(old = oldold_genotype, new = new_genotype)],
    gc[!is.na(old_genotype), .(old = old_genotype, new = new_genotype)]
  ))[, .SD[1], by = old]
  ph <- merge(ph, o2n, by.x = "old_genotype", by.y = "old", all.x = TRUE)
  ph[, Genotype := fcase(
    Species == "B73" | old_genotype == "B73", "B73",
    Species == "Check", "Purple",
    !is.na(new), canon_ped(new),
    default = old_genotype
  )]
  # binarize the 0/1/2 severity score to kinked presence (>=1)
  ph[, Kinki := as.integer(suppressWarnings(as.numeric(Kinki)) >= 1)]
  merge(fm, ph[, .(plot_id, Genotype, Kinki)], by = "plot_id")
}

plots <- manifest_cly23()[is.finite(Kinki) & !is.na(Genotype)]
el <- plots[, .(k = sum(Kinki), n = .N), by = Genotype]
el[, prop := k / n][, Kinki_mean := log((k + 0.5) / (n - k + 0.5))]
fwrite(el[, .(Genotype, k, n, prop, Kinki_mean)], here("data/zeal/pheno_kinki_elogit.csv"))
log_info(
  "elogit Kinki (CLY23 only): %d genotypes | %d ever-kinked (k>0) | elogit range [%.2f, %.2f]",
  nrow(el), sum(el$k > 0), min(el$Kinki_mean), max(el$Kinki_mean)
)

ss <- fread(here("data/zeal/samplesheet_3way.csv"))
tass <- merge(ss[gwas_nil == TRUE, .(pedigree, taxon)], el[, .(pedigree = Genotype, y = Kinki_mean)], by = "pedigree")[is.finite(y)]
ph_out <- here("data/zeal/tassel/pheno_kinki_all.txt")
writeLines(c("<Phenotype>", "taxa\tdata\tfactor", "Taxa\tKinki\tFamily"), ph_out)
fwrite(tass[, .(pedigree, round(y, 4), taxon)], ph_out, sep = "\t", append = TRUE, col.names = FALSE)
log_info("wrote %s (elogit, %d gwas_nil lines)", ph_out, nrow(tass))
