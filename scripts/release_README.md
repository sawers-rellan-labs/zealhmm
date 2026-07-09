# BZea NIL genotype release

Ancestry and genotype calls for the BZea BC2S3 near-isogenic lines (~1,400 lines,
~80 teosinte donors × recurrent B73), on the **SNP50K** teosinte-vs-B73 panel
(~49,002 sites, **B73 AGPv5**). Built by `scripts/zeal_export_release.R` in the
[`zealhmm`](https://github.com/sawers-rellan-labs/zealhmm) repo. Terminology and
conventions: see `TERMINOLOGY.md` there.

## The hard distinction: ancestry mosaic vs genotype

- **`*_mosaic`** = an **ancestry** inference (which parental genome each segment came
  from), produced by an HMM. States `0/1/2 = B73 / het / teosinte` **ancestry**.
- **`ml_gt`** = an actual **genotype** called per-site from the read counts (no HMM).
  States `0/1/2 = REF / het / ALT` allele dosage.

> ⚠️ A mosaic is **not** a genotype. Its `0/1/2` is *ancestry* dosage — at an invariant
> site inside a teosinte block it reports `2` even though the allele equals B73. The
> PLINK/`012` files encode that ancestry dosage on the SNP's ref/alt alleles **only for
> tooling compatibility**; do not read a `*_mosaic` file as allele genotypes. Only
> `ml_gt` reports true alleles.

## Contents

```
markers/snp50k_markers.tsv        marker, chr, pos, ref, alt, cM  (B73 v5)
lines/snp50k_lines.tsv            pedigree, taxon, donor_accession, taxa_code, skim_id
snp50k/
  bzea_snp50k_<name>.{bed,bim,fam}   PLINK 1 binary (canonical; "binary 012")
  bzea_snp50k_<name>_012.tsv.gz      tidy marker × line 0/1/2 matrix
  bzea_snp50k_<name>.rds             native R list(markers, state, lines)
250k/
  bzea_250k_rtiger_introgression_segments.tsv.gz   previous 2-state segments
  bzea_250k_rtiger_introgression.rds
MANIFEST.tsv                       file, bytes, sha256
```

### 50K variations (`<name>`)

| name | layer | caller / method | notes |
|------|-------|-----------------|-------|
| `rtiger_mosaic` | ancestry | RTIGER | reference mosaic; recovers teosinte presence ~0.10 |
| `nnil_mosaic` | ancestry | nilHMM nNIL | per-SNP recall; under-calls at 0.4× |
| `binhmm_mosaic` | ancestry | bin-HMM (1 Mb) | block-aggregated |
| `lbimpute_mosaic` | ancestry | LB-Impute | |
| `ml_gt` | **genotype** | `call_gt(prior="flat")` = argmax-GL (ML) | per-site; the observed-evidence layer |

### 250K (previous)

The genotypes behind the **inv4m paper** — an older RTIGER call set on a larger
GATK/QC variant panel (Nirwan Tandukar's introgression finder). **2-state only**
(`B73` / `Introgression`), shipped as per-line segments. Kept for cross-comparison;
prefer the SNP50K 3-state calls above for new work.

## Encoding

`0/1/2` with `NA`/`./.` for no-call. PLINK: `A1`=ALT (teosinte-informative / donor),
`A2`=REF (B73). Coordinates are **B73 AGPv5** bp; `cM` is the native genetic map.

## Provenance & citation

SNP50K counts and callers: `zealhmm` (`DATA.md`). Package: `nilHMM` (≥0.3.0). 250K:
`nirwan1265/BZea_Introgression_Finder`. _Cite: <DOI pending CyVerse deposit>._
