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
one. "Imputed genotype" is also banned: it's ambiguous (the 250K set was *already* Beagle-imputed)
and it blurs the wall.

Never call our output "genotypes." Genotypes are inputs; the mosaic is what we infer.

## Genotype layer — `_gt`

A genotype is a per-site call from read counts (the evidence), via `call_gt()` (nilhmm). The
estimator is always `argmax`; the **method = the prior**:

| token | estimator | prior | `call_gt(...)` |
|-------|-----------|-------|----------------|
| `gl` | argmax **GL** (ML) | flat / none | `prior = "flat"` |
| `gphwe` | argmax **GP** (MAP) | HWE at panel AF | `prior = "hwe", af = <AF>` |
| `gpdesign` | argmax **GP** (MAP) | Mendelian single-locus expectation of the cross | `prior = design_prior("BC2S3")` |
| `gp` | argmax **GP** (MAP) | custom | `prior = <length-3 vector>` |

Notes:
- Tokens mirror the **VCF FORMAT fields**: `GL` (likelihood) vs `GP` (posterior). See the
  `gl-gp-map-caller-terminology` memory + `genotype_likelihoods_and_hmm.qmd §2`.
- Bare **`gp`** = "MAP with a custom prior you pass in"; `gphwe`/`gpdesign` are the two named
  priors. (`gphwe` and `gpdesign` *are* also MAP — `gp` alone just means custom.)
- **`gpdesign` is per-site and marginal** — it applies the design's single-locus expectation
  independently at each site, *no linkage*. That's the deliberate contrast to `call_ancestry`,
  which adds the recombination/linkage prior. The design prior is **derived** from
  `design_priors("BC2S3")` → `c(.859, .031, .109)` (55:2:7 / 64), never hand-typed.
- **Representation:** the canonical `_gt` is the **hardcall** (`return="call"`, 0/1/2 → the VCF
  `GT`). A continuous **dosage** `E[G|reads]` for GWAS power is derived from the GP array
  (`return="post"`): `gp[,,2]*1 + gp[,,3]*2`. Do **not** use `return="dosage"` for this — it
  returns the hardcall as a double, not the posterior mean.
- **`persnp` is retired** (replaced by the `_gt` scheme). The old ad-hoc `round(2·alt/cov)` hardcall
  is now built by `scripts/zeal_gt.R` (`METHOD=gl|gphwe|gpdesign`) via `call_gt` → `zeal_<method>_gt.rds`;
  `gl_gt` is the canonical per-SNP genotype (GWAS Panel B / the per-SNP control).

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

- `GENO=gl_gt`, `GENO=gphwe_gt`, `GENO=gpdesign_gt`, `GENO=gp_gt`  (genotype layer)
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
| **250K** (inv4m / introfinder) | *previous* RTIGER introgression calls | **2-state** (B73 / introgression) | Nirwan; always label "previous" |
| **wideseq** (~27.6 M) | dense teosinte panel (Schnable 2023, MAF ≥ 0.05) — the expansion base | counts | for the Poisson-QC union with SNP50K |

Spelling: **SNP50K** (not `50K`/`snp50k`/`SNP 50K`), **wideseq** (not `27M`/`wideseq_ref`),
**250K** (always tagged "previous"; it is 2-state).

## Germplasm

**taxa**, never "species": Zv / Zx / Zl / Zd / Zh. Dotted codes `Zx.####` are individual founder
accessions. See the `bzea-taxa-naming` memory.

## Verbs vs nouns (quick reference)

- `call_gt()` — verb, genotype layer → a **genotype** (`_gt`); nilhmm ≥ 0.3.0 (the old `call_gl`
  name was removed, no alias). Builder: `scripts/zeal_gt.R`.
- `call_ancestry()` / `call_states()` — verb, mosaic layer → an **ancestry mosaic** (`_mosaic`).
- **caller** = a method (engine + params). **genotype** = input. **mosaic** = output. **GENO** =
  the object name.