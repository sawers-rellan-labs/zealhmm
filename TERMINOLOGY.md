# Terminology (ZEAL/BZea — paper & data)

Paper-facing and general/biological vocabulary for the ZEAL/BZea work: phenotype
definitions, germplasm/taxa, datasets, the population-structure covariate, and the
ZEAL analysis object names. The analysis notebooks, `DATA.md`, the scripts, and the
genotype browser follow this.

> **Scope.** Package **coding/naming conventions** — the genotype↔ancestry-mosaic
> wall, `call_gt` / `call_ancestry` / `call_states`, the caller list, the
> `_gt`/`_mosaic` suffix convention, and GL/GP/MAP — live in
> **`nilhmm/design/TERMINOLOGY.md`** and are **not repeated here**. This doc is the
> ZEAL side: what the biological and data terms mean for the paper.

## Datasets

| name | what | states | provenance |
|------|------|--------|-----------|
| **SNP50K** | current teosinte-vs-B73 panel (~52 K sites, B73 v5), read counts → callers | 3-state 0/1/2 | Fausto; `bzea-genotype-call-sets` |
| **200K** (inv4m / introfinder) | *previous* RTIGER introgression calls | **2-state** (B73 / introgression) | Nirwan; always label "previous" |
| **wideseq** (~27.6 M) | dense teosinte panel (Schnable 2023, MAF ≥ 0.05) — the expansion base | counts | for the Poisson-QC union with SNP50K |

Spelling: **SNP50K** (not `50K`/`snp50k`/`SNP 50K`), **wideseq** (not
`27M`/`wideseq_ref`), **200K** (always tagged "previous"; it is 2-state).

## Germplasm / taxa

Always **taxon / taxa** — **never** species, subspecies, population, or race. The BZea
donors are of *mixed taxonomic rank* (some species, e.g. *Z. diploperennis* /
*Z. luxurians*; some subspecies of *Z. mays*, e.g. *mexicana* / *parviglumis* /
*huehuetenangensis*), so no single rank word is correct — "taxon" is the rank-agnostic
umbrella. Taxa: **Zv / Zx / Zl / Zd / Zh**. Dotted codes `Zx.####` are individual
**founder accessions** (82 across the 5 taxa; see the `bzea-taxa-naming` memory).

## Population-structure covariate (Family vs. taxon)

The MLM structure covariate is **not** named the same across datasets:

- **TeoNAM (Chen 2019):** Q (first 5 PCs of the kinship K) ≈ a 5-level **Family** factor
  = the 5 donor parents (1 parent = 1 family = 1 taxon draw). My TeoNAM *validation*
  substitutes an explicit 5-level "Family" factor for Q.
- **ZEAL/BZea:** the covariate is **taxon (5)**, *not* per-founder family — the 82
  founders collapse into 5 taxa. "Family" carried over from TeoNAM is misleading here.

Correction pending (relabel ZEAL "Family" → taxon; explain the distinction in the
notebooks/paper) — see the `zeal-q-analog-is-taxon-not-family` memory.

## Phenotypes / trait panel

The **fieldbooks are the authoritative phenotype source**, and **definitions can vary by
field/year** — always cite the field when the protocol or units matter. Two BZea evaluation
fields feed the panel (both Clayton, NC):

- **CLY23-D4** (2023) — `data/zeal/CLY23_D4_FieldBook.xlsx`, sheet `UPDATED_CLY23_D4_FieldBook`.
  Trait columns present, **no glossary sheet**. `DTA`/`DTS` are precomputed columns; the binary
  stem traits (`StPi`, `StPu`) and `Kinki` are here.
- **CLY25-B5** (2025) — `data/zeal/CLY25-Fieldbook.xlsx`, sheet `B5_BZea_eval`. Its
  **`Column_definitions` sheet is the canonical glossary** (definitions below are quoted from it).
  Silking/anthesis are the raw date columns **DOA/DOS**; the codes **DTA/DTS** are *derived*
  (Excel-serial date − planting 2025-04-03).

Traits enter the GWAS/JLM as **per-genotype values**: **continuous** traits as a **per-field
SpATS spatial BLUE** (unioned across the fields a trait was scored in); **binary/ordinal** traits
as an **empirical logit** of the per-genotype proportion of positive plots, with **no spatial
correction**.

| code | trait | definition | type | per-genotype value | field |
|------|-------|-----------|------|--------------------|-------|
| **DTA** | days to anthesis | from **DOA**: date 50% of plants shed pollen on central + lateral tassel spikes | continuous (d) | SpATS BLUE | CLY23 (precomputed), CLY25 (DOA − planting) |
| **DTS** | days to silk | from **DOS**: date 50% of plants had visible silks | continuous (d) | SpATS BLUE | CLY23, CLY25 |
| **PH** | plant height | base of plant → tip of highest tassel | continuous (cm) | SpATS BLUE | CLY23, CLY25 |
| **EH** | ear height | base of plant → primary ear-bearing node | continuous (cm) | SpATS BLUE | CLY23, CLY25 |
| **EN** | ear number | number of **nodes** with ears per plant | continuous (count) | SpATS BLUE | CLY23, CLY25 |
| **Prolif** | prolificacy | **total** number of ears per plant* | continuous (count) | SpATS BLUE | CLY23, CLY25 |
| **NBR** | brace-root number | number of nodes with brace roots per plant | continuous (count) | SpATS BLUE | CLY23, CLY25 |
| **LAE** | leaves above ear | number of leaves above the primary ear | continuous (count) | SpATS BLUE | CLY25 |
| **SPAD** | leaf greenness | SPAD-meter greenness index on one leaf above the primary ear-bearing node | continuous | SpATS BLUE (see variants) | CLY25 |
| **StPi** | stem pigment | binary stem-**anthocyanin** (pigmented) score | binary 0/1 | empirical logit, no spatial corr. | CLY23 |
| **StPu** | stem pubescence | binary stem-**macrohair** (pubescence-present) score | binary 0/1 | empirical logit, no spatial corr. | CLY23 |
| **Kinki** | zigzag culm | kinked/zigzag stem (Eyster 1920); CLY23 ordinal severity, binarized | binary (from ordinal) | empirical logit, no spatial corr. | **CLY23 only** |

\*Prolif protocol (CLY25 glossary): count all ears where a cob can be felt; if none, count ears
≥ 50% the size of the primary ear.

**Naming / variants (important):**

- **DOA/DOS** are the *raw date* fieldbook columns; **DTA/DTS** are the *derived* days-to counts.
  Growing-degree-day variants (`GDDTA`/`GDDTS`) also exist in the CLY23 fieldbook.
- **SPAD is three notebooks**: `spad` = CLY25 `SPAD` SpATS BLUE (spatially corrected);
  `spad20das` / `spad36das` = per-line means from the zealbrowser `Bzea_merged`
  `SPAD_leaf_greenness_{20,36}DAS` columns, **as-is, no spatial correction**. (CLY23 `SPAD1` is a
  placeholder of all 1.0; the real CLY23 reading is `SPAD2`.)
- **StPi ≠ StPu**: StPi = stem **pigment** (anthocyanin); StPu = stem **pubescence** (macrohairs).
  Do not conflate. **Kinki is CLY23-only.**
- Definitions above are the **CLY25** `Column_definitions` glossary; a trait scored in both fields
  may differ in protocol by year.

## ZEAL analysis object naming

ZEAL's genotype and mosaic *objects* (used by `scripts/`, `DATA.md`, the browser, and the
release). For the genotype↔mosaic distinction, the `_gt`/`_mosaic` convention, and GL/GP/MAP,
see **`nilhmm/design/TERMINOLOGY.md`**.

- **Genotype object — `hwe_post_gt`** (`_gt` layer): the authoritative cohort VCF
  `bzea_50K_cohort.vcf.gz` (`bcftools mpileup … | call -mv`, HWE-prior MAP = "HWE-posterior");
  REF = B73, so `0/1/2 = B73 / het / teosinte` allele dosage. Built into `zeal_hwe_post_gt.rds`
  (`scripts/zeal_hwe_post_gt.R`); GWAS Panel B; shipped in the release. Single-sample
  GL/argmax-GL genotypes and the mistaken `ml_gt` / `persnp` are **retired**, not used.
- **Mosaic objects — `<caller>_mosaic`** (`_mosaic` layer): our HMM ancestry inference via
  `call_ancestry()` (e.g. `zeal_rtiger_mosaic.rds`); callers `rtiger` / `nnil` / `binhmm` /
  `lbimpute`. Schema `list(markers, state[marker × line, 0/1/2], lines)`. `mosaic` is a noun,
  never a caller name.
- **`GENO` / predictor selection** takes the **object name directly** — the suffix carries the
  type, no prefix or alias: `GENO=hwe_post_gt` (genotype); `GENO=rtiger_mosaic` / `nnil_mosaic`
  / `binhmm_mosaic` / `lbimpute_mosaic` (mosaic). Dead/removed: the `mosaic:` prefix, the bare
  `mosaic`=rtiger alias, `persnp`, `ml_gt`.
