# ZEAL NIL genotype release

Ancestry and genotype calls for the ZEAL BC2S3 near-isogenic lines (~1,400 lines,
~80 teosinte donors × recurrent B73) on the **SNP50K** teosinte-vs-B73 panel
(~49,002 sites, **B73 AGPv5**). Built by `scripts/zeal_export_release.R` in the
[`zealhmm`](https://github.com/sawers-rellan-labs/zealhmm) repo. Terminology and
conventions: `TERMINOLOGY.md` there.

## The hard distinction: ancestry mosaic vs genotype

- **`<caller>_mosaic`** = an **ancestry** inference (which parental genome each segment came
  from), produced by an HMM. States `0/1/2 = B73 / het / teosinte` **ancestry**.
- **`hwe_post_gt`** = the actual **genotype** — the real `bcftools mpileup | call -mv`
  (HWE-prior MAP) calls from the cohort VCF. REF = B73, so `0/1/2 = B73 / het / teosinte`
  allele dosage.

> ⚠️ A mosaic is **not** a genotype. Its `0/1/2` is *ancestry* dosage — at an invariant
> site inside a teosinte block it reports `2` even though the allele equals B73. The
> PLINK/`012` files encode that ancestry dosage on the SNP's ref/alt alleles **only for
> tooling compatibility**; do not read a `*_mosaic` file as allele genotypes. Only
> `hwe_post_gt` reports true alleles.
>
> No single-sample GL / argmax-GL genotypes are distributed — only the cohort-called
> `hwe_post_gt`, which is what the project actually produced and used.

## One shared panel

Every object uses the **same** line panel — the intersection of all objects' callable lines.
rtiger's per-chromosome "≥ 2×rigidity informative markers" QC is the binding constraint: 8
near-empty libraries (< 0.3% of SNP50K sites covered) that rtiger cannot call are dropped, so
every object carries the identical set of ~1,395 lines in the same column order.

## Contents

```
markers/snp50k_markers.tsv        marker, chr, pos, ref, alt, cM  (B73 v5)
lines/snp50k_lines.tsv            pedigree, taxon, donor_accession, taxa_code, skim_id
snp50k/
  bzea_snp50k_<name>.{bed,bim,fam}   PLINK 1 binary (canonical; "binary 012")
  bzea_snp50k_<name>_012.tsv.gz      tidy marker × line 0/1/2 matrix
  bzea_snp50k_<name>.rds             native R list(markers, state, lines)
MANIFEST.tsv                       file, bytes, sha256
```

### The objects (`<name>`)

| name | layer | caller / method |
|------|-------|-----------------|
| `rtiger_mosaic` | ancestry | rtiger (reference mosaic; recovers teosinte presence ~0.10) |
| `nnil_mosaic` | ancestry | nilHMM nnil |
| `binhmm_mosaic` | ancestry | bin-HMM (1 Mb) |
| `lbimpute_mosaic` | ancestry | LB-Impute |
| `hwe_post_gt` | **genotype** | `bcftools mpileup \| call -mv` (HWE-prior MAP) — real cohort calls, panel-subset |

The legacy **200K** rtiger introgression set (inv4m paper; a separate GATK → Beagle-imputed
lineage, not from this pipeline) is intentionally **not** in this release — regenerate with
recalibrated rtiger if it is ever needed.

## Encoding

`0/1/2` with `NA`/`./.` for no-call. PLINK: `A1`=ALT (teosinte-informative / donor),
`A2`=REF (B73). Coordinates are **B73 AGPv5** bp (`chr1`..`chr10`); `cM` is the native map.

## Provenance & citation

SNP50K counts and callers: `zealhmm` (`DATA.md`). Genotype (`hwe_post_gt`) from
`bzea_50K_cohort.vcf.gz` (`bcftools call -mv`). Package: `nilHMM` (≥0.3.0).
_Cite: <DOI pending CyVerse deposit>._
