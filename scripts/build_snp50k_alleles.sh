#!/usr/bin/env bash
# Build the single site+allele authority for the SNP50K set:
#   data/zeal/snp50k_alleles.tsv    chr  pos  ref  alt     (49,002 rows, headerless)
#
# SOURCE = the 50K set itself: bzea_50K_cohort.vcf.gz, the biallelic-SNP VCF produced from the
# joint call cohort.vcf.gz by results/joint/get_cohort_and_reference_vcf.sh Step 1
#   (bcftools view -v snps -m2 -M2). One file defines BOTH which sites are in the 50K set AND
# what their REF/ALT alleles are, so there is a single authority for site membership and alleles.
#
# WHY THIS SOURCE (not HQ_BZEA + a cM marker list): the 50K set is DEFINED by this biallelic call,
# so its own REF/ALT is the self-consistent choice. It differs from the HQ_BZEA panel ALT at only
# 65 of 49,002 sites, and a read count at those sites (scripts/check_65_alt_reads.R) found just 6 ALT
# reads total across 11 skim samples, so it makes no difference which allele we pick. See DATA.md and
# [[50k-set-terminology]].
#
# The per-sample GATK count files carry a MADE-UP alt base (CollectAllelicCounts infers alt from
# reads and writes N at 0.4x), so scripts/extract_gatk_counts_50k.sh takes REF/ALT from THIS file,
# never from the count tables.
#
#   bash scripts/build_snp50k_alleles.sh
#   env: BZEASEQ_DIR (default /Volumes/rsstu/users/r/rrellan/BZea/bzeaseq)
set -uo pipefail
BZ=${BZEASEQ_DIR:-/Volumes/rsstu/users/r/rrellan/BZea/bzeaseq}
V="$BZ/50K/results/joint/bzea_50K_cohort.vcf.gz"
OUT=data/zeal/snp50k_alleles.tsv
BCFTOOLS=${BCFTOOLS:-$(command -v bcftools || echo /opt/homebrew/bin/bcftools)}

[ -r "$V" ] || { echo "ERROR: cannot read 50K-set VCF $V (mount down?)" >&2; exit 1; }
[ -x "$BCFTOOLS" ] || { echo "ERROR: bcftools not found (set BCFTOOLS)" >&2; exit 1; }

# CHROM is already chr-prefixed and every record is biallelic (0 multiallelic), so a straight query
# is the whole job. Normalize CHROM to chrN defensively so the output always matches the count
# tables' CONTIG column.
"$BCFTOOLS" query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$V" \
  | awk -F'\t' 'BEGIN{OFS="\t"}{ c=$1; if (c !~ /^chr/) c="chr"c; if ($4 ~ /,/) next; print c,$2,$3,$4 }' \
  > "$OUT"

n=$(wc -l < "$OUT" | tr -d ' ')
echo "wrote $OUT: $n rows (expect 49002)"
echo "REF composition: $(cut -f3 "$OUT" | sort | uniq -c | sort -rn | awk '{printf "%s=%s ", $2, $1}')"
echo "ALT composition: $(cut -f4 "$OUT" | sort | uniq -c | sort -rn | awk '{printf "%s=%s ", $2, $1}')"
[ "$n" = "49002" ] || { echo "ERROR: expected 49002 rows" >&2; exit 1; }