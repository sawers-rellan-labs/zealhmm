# Terminology & naming conventions

The canonical vocabulary for ZEAL/BZea genotype and ancestry work. Scripts, `DATA.md`, the
analysis notebooks, and the genotype browser all follow this. The governing idea is a hard wall
between **observed genotypes** (inputs) and **inferred ancestry** (our output).

## The two layers

| | **Genotype** | **Ancestry mosaic** |
|-|--------------|---------------------|
| Role | **input / evidence** | **our inference / output** |
| Source | actual observations, or calls from a prior caller | inferred by us from the evidence |
| Method | per-site, **no linkage / no HMM** | HMM **across** sites (recombination + design prior) |
| Verb / fn | `call_gt()` | `call_ancestry()` |
| Output suffix | **`_gt`** | **`_mosaic`** |
| States 0/1/2 | genotype dosage: 0=REF(B73), 1=HET, 2=ALT(donor/teosinte) | ancestry: 0=B73, 1=het, 2=teosinte |

**Why the wall matters:** a mosaic *overrides* the genotype — a teosinte-ancestry block reports
state `2` even at an invariant site (allele identical to B73), and it collapses observed hets into
an ancestry state. So a mosaic is **not** a genotype and must never be filed, labeled, or served as
one. "Imputed genotype" is also banned: it's ambiguous (the 200K set was *already* Beagle-imputed)
and it blurs the wall.

Never call our output "genotypes." Genotypes are inputs; the mosaic is what we infer.

## Genotype layer — `_gt` (`hwe_post`)

A genotype is a per-site call from the reads (the evidence), **not** an HMM inference. The
genotype ZEAL actually uses and shares is **`hwe_post`** — the authoritative cohort VCF
`bzea_50K_cohort.vcf.gz`, produced by

```
bcftools mpileup -f <B73 v5> -R <SNP50K sites> | bcftools call -mv     # HWE-prior MAP = "HWE-posterior"
```

REF = B73, so its `0/1/2 = B73 / het / teosinte` allele dosage. Built into `zeal_hwe_post_gt.rds`
by `scripts/zeal_hwe_post_gt.R` (the real bcftools calls, extracted for the panel); used as GWAS
Panel B and shipped in the release on the common callable panel — values are the bcftools calls,
not a reconstruction (no ambiguity); same call set as the VCF sent to Jim Holland.

**Not used / not shared:** single-sample **GL / argmax-GL** genotypes (`bcftools`'s underlying GL,
or nilhmm `call_gt(prior="flat")`). They were never part of the pipeline. An earlier `ml_gt`
(single-sample argmax-GL) was created by mistake during the persnp retirement and has been removed
from the repo and the release.

`GL/GP/MAP` terms (reference): **GL/PL** = likelihood `P(reads|G)`; **GP** = posterior `P(G|reads)`
(VCF field); **MAP** = argmax-GP estimator. `call -m` applies an HWE genotype prior → its `GT` is
the MAP (hence *HWE-posterior*). See the `gl-gp-map-caller-terminology` memory. `call_gt()` exists in nilhmm
(a per-site caller with a swappable prior) but the ZEAL genotype comes from bcftools, not `call_gt`.

- **`persnp` is retired** and so is the mistaken `ml_gt`; the genotype layer is `hwe_post` (above).
- `design_priors("BC2S3")` → `c(.859, .031, .109)` (55:2:7 / 64) remains available as a prior
  vector for anyone who wants a design-informed per-site caller, but ZEAL does not use it.

## Ancestry-mosaic layer — `<caller>_mosaic`

Our HMM ancestry inference, via `call_ancestry()` / `call_states()`. The **caller** is the method
(HMM engine + parameter set):

| caller | engine |
|--------|--------|
| `rtiger` | RTIGER |
| `nnil` | nilHMM nNIL |
| `binhmm` | bin-HMM |
| `lbimpute` | LB-Impute |

- Output object / file: `<caller>_mosaic` (e.g. `zeal_rtiger_mosaic.rds`), schema
  `list(markers, state[marker × line, 0/1/2], lines)`.
- **`mosaic` is a noun (the ancestry-state matrix), never a caller name.** The retired
  `GENO=mosaic` alias silently meant "rtiger" — that is gone.

## GENO / predictor selection

`GENO` (and any predictor selector) takes the **object name directly** — the suffix carries the
type. No prefix, no alias:

- `GENO=hwe_post_gt`  (genotype layer — the bcftools cohort genotypes)
- `GENO=rtiger_mosaic`, `GENO=nnil_mosaic`, `GENO=binhmm_mosaic`, `GENO=lbimpute_mosaic`  (mosaic layer)

Dead and removed: the `mosaic:` prefix, the bare `mosaic`=rtiger alias, `persnp`.

## GL / GP / MAP (reference)

- **GL / PL** — genotype *likelihood* `P(reads|G)` (phred-scaled = PL). GATK's `GT` = argmax GL (ML).
- **GP** — genotype *posterior* `P(G|reads)` (VCF ≥4.3 FORMAT field). Beagle/GLIMPSE/bcftools.
- **MAP** — reserved for the *estimator* only: argmax over GP. Not a token, not a filename.
  (Kept distinct from `map` the genetic/recombination map and `Map()` the function — hence GP, not
  MAP, in the `_gt` tokens.)

## Datasets

| name | what | states | provenance |
|------|------|--------|-----------|
| **SNP50K** | current teosinte-vs-B73 panel (~52 K sites, B73 v5), read counts → `call_gt`/callers | 3-state 0/1/2 | Fausto; `bzea-genotype-call-sets` |
| **200K** (inv4m / introfinder) | *previous* RTIGER introgression calls | **2-state** (B73 / introgression) | Nirwan; always label "previous" |
| **wideseq** (~27.6 M) | dense teosinte panel (Schnable 2023, MAF ≥ 0.05) — the expansion base | counts | for the Poisson-QC union with SNP50K |

Spelling: **SNP50K** (not `50K`/`snp50k`/`SNP 50K`), **wideseq** (not `27M`/`wideseq_ref`),
**200K** (always tagged "previous"; it is 2-state).

## Germplasm

**taxa**, never "species": Zv / Zx / Zl / Zd / Zh. Dotted codes `Zx.####` are individual founder
accessions. See the `bzea-taxa-naming` memory.

## Verbs vs nouns (quick reference)

- **genotype** (`hwe_post`, `_gt`) → `bcftools mpileup | call -mv` (HWE-prior MAP); built into
  `zeal_hwe_post_gt.rds` by `scripts/zeal_hwe_post_gt.R`. (`call_gt()` exists in nilhmm ≥ 0.3.0 as a
  per-site caller, but ZEAL's genotype comes from bcftools, not `call_gt`.)
- `call_ancestry()` / `call_states()` — verb, mosaic layer → an **ancestry mosaic** (`_mosaic`).
- **caller** = a method (engine + params). **genotype** = input. **mosaic** = output. **GENO** =
  the object name.