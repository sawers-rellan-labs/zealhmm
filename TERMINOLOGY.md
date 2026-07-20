# Terminology (ZEAL — paper & data)

Paper-facing and general/biological vocabulary for the ZEAL work: phenotype
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

Always **taxon / taxa** — **never** species, subspecies, population, or race. The ZEAL
donors are of *mixed taxonomic rank* (some species, e.g. *Z. diploperennis* /
*Z. luxurians*; some subspecies of *Z. mays*, e.g. *mexicana* / *parviglumis* /
*huehuetenangensis*), so no single rank word is correct — "taxon" is the rank-agnostic
umbrella. Taxa: **Zv / Zx / Zl / Zd / Zh**. Dotted codes `Zx.####` are individual
**founder accessions** (82 across the 5 taxa; see the `bzea-taxa-naming` memory).

## Population-structure covariate (Family vs. taxon)

The structure covariate differs by dataset **and by analysis** — and, importantly, the ZEAL
scans as *run* are not what the older "Taxon" labels implied:

- **TeoNAM (Chen 2019):** Q (first 5 PCs of the kinship K) ≈ a 5-level factor = the 5 donor
  parents (1 parent = 1 family). The TeoNAM validation uses an explicit 5-level **Family**
  factor in place of Q.
- **ZEAL — what the code actually did:** the **MLM and per-marker OLS** scans were
  computed with the **82-donor-accession** factor (`FAMILY_COL=donor_accession`,
  parent-plant-level correction) — this is what "**Family**" denotes in the notebooks. The
  **JLM** instead nests marker effects within the **5 donor taxa** (the TASSEL pheno-file
  `Family` column holds taxon codes `Zd`/`Zx`/…). So MLM/OLS (82) and JLM (5) genuinely use
  **different** structure factors; the notebooks say "**Family**" for the 82-accession scans
  and "**taxon**" for the JLM nesting, each matching what it computed.

The 82-accession correction is a strong (parent-plant-level) structure control that likely
over-corrects and costs power. **Pending comparison** (`zeal-q-analog-is-taxon-not-family`
memory): rerun the MLM/OLS at the **taxon (5)** level and decide — if over-inflated, discard;
if it recovers power without inflation, it *replaces* the donor-accession scans in the paper.

## Phenotypes / trait panel

The **fieldbooks are the authoritative phenotype source**, and **definitions can vary by
field/year** — always cite the field when the protocol or units matter. Two ZEAL evaluation
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

## nNIL calibration foil (nilHMM methods)

These are **nilHMM methods** terms, not ZEAL biology: the calibration *procedure*, the
*reference* it is scored against, the *objective metric*, and the two HMM *knobs* used by
the nNIL calibration foil (`scripts/nnil_foil/`, `analysis/nnil-calibration-foil.qmd`),
which diagnoses how Zhong et al. (2025) calibrated the `nnil` caller. Keeping them distinct
matters, because conflating them caused real confusion. (Package coding/naming conventions
still live in **`nilhmm/design/TERMINOLOGY.md`**.)

| term | meaning | notes / do not confuse with |
|------|---------|-----------------------------|
| **chip-supervised GBS calibration** (of the HMM) | Tuning the GBS caller's parameters to minimize disagreement between its GBS calls and the chip calls on the 24 both-platform lines. **GBS is the sole HMM input; the chip supplies reference labels.** This is Holland's File_S16. | Deprecated labels: "joint calibration", "GBS ∪ chip", "GBS-vs-chip calibration". It is **not** a union co-fit (chip genotypes are never stacked as input) and **not** Kennedy-O'Hagan (no simulator). |
| **in-sample calibration** | The calibration lines are a subset of the prediction set (GBS_calib ⊂ GBS_predict); the tuned parameters are then applied to the same and larger GBS set. | Not held out. The clean nilHMM analysis keeps calibration (on simulation) separate from evaluation (on held-out GBS). |
| **chip calls** | The chip-only HMM introgression calls (Holland File_S14), used as the reference labels for chip-supervised calibration. | Always "calls", **never "chip truth"**: it is a chip-only fit, not ground truth. |
| **simcross truth** | The latent generating ancestry in the `simcross` simulation. The one genuine ground truth in the foil. | Only the simulation has it; the real nNIL has only the sparse chip calls. |
| **GBS-vs-chip calls mismatch** | Per-marker disagreement between the GBS caller's calls and the chip calls on shared markers. Holland's original objective. | The metric everywhere the target is the chip calls. Distinct from **marker mismatch**. |
| **marker mismatch** | The analogous per-marker disagreement against the **simcross truth** in the simulation. | Kept as a distinct name on purpose: the target differs (latent truth, not chip calls). |
| **nir** (emission non-informative rate) | Probability that a donor haplotype still shows the REF (B73) allele at a marker. An *emission* parameter, **not** a segment-length prior. | Carries the calibration signal (97.7% of File_S04's mismatch variance). |
| **GBS nir** | The *caller's* nir as applied to the GBS data; grid-tuned to ~0.9. Labels the GBS-caller sweep (e.g. panel D of the fragment-size figure). | Distinguish from the founder value below. |
| **founder non-informative rate (f0)** | The *biological/structural* non-informative rate from the nested B73 × NAM-founder design (~0.59, measured on the NAM-founder chip genotypes). | The gap **GBS nir (0.9) minus founder nir (0.59)** is the GBS data-quality effect, not biology. |
| **r** / **map r** (intermarker recombination fraction) | Per-adjacent-marker transition probability (off-diagonal HMM mass), unitless and per-marker-interval (not per bp or per cM, not realized recombination). **map r** = 2L/(100N) from the genetic map (native TeoNAM v5). | A benign hyperparameter here (0.03% of File_S04's variance); a genetic map sets it, no chip or simulation needed. |

## Model calibration: Kennedy-O'Hagan (KOH)

**Kennedy-O'Hagan (KOH)** is the orthodox Bayesian framework for calibrating a simulator
against real data:

> z(x) = rho * eta(x, theta) + delta(x) + epsilon

where **eta** is the simulator, **theta** the calibration parameters, **delta(x)** a
model-discrepancy term absorbing what the simulator systematically gets wrong, and
**epsilon** observation error. **theta** and **delta** are confounded: with no explicit
discrepancy term a calibration parameter silently absorbs model inadequacy and stops
meaning what you think (Brynjarsdottir and O'Hagan 2014).

Mapping to the nNIL foil:

- **eta** = `simcross` ancestry + Holland's emission (a genuine generative model).
- **theta** = {`r` (identified: the map value, insensitive), `nir`}.
- **delta** = the GBS data-quality gap, made identifiable by the paired 24-line design
  (chip and GBS on the same individuals).
- The **nir 0.594 to 0.9 gap is the discrepancy** leaking into `theta`.

What is and is not KOH here:

- **Sim-only nir recovery** (`09_sim_nir_sweep.R`): generate at a known nir = 0.594, score
  against `simcross` truth, recover 0.594. Legitimate simulation-based parameter recovery,
  and the **eta-only, delta = 0 corner** of KOH, not a full KOH calibration.
- **chip-supervised calibration** (Holland File_S16): fits `theta` to the chip calls with
  **no discrepancy term** and in-sample. **Not KOH** (no simulator; the target is a second
  inference on a second platform).
- **KOH proper** (not yet done): infer `theta` and an explicit `delta` jointly from the
  simulator plus the paired chip data. Only reach for `delta` on genuinely noisy sets
  (nNIL GBS, ZEAL skim/BrB), never on the clean TeoNAM GBS.

References: Kennedy and O'Hagan (2001), *J. R. Stat. Soc. B* 63(3):425-464,
doi:10.1111/1467-9868.00294; Brynjarsdottir and O'Hagan (2014), *Inverse Problems*
30(11):114007, doi:10.1088/0266-5611/30/11/114007.
