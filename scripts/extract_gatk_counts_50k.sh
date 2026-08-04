#!/usr/bin/env bash
# Extract per-sample SNP50K read counts from the aggregated GATK table into the RTIGER 6-column
# layout, filtered to the 49,002-site biallelic panel.
#
# WHY THIS SOURCE: `CollectAllelicCounts` tallies reads straight off each BAM at the panel
# positions (`-I bam`, `-L HQ_BZEA.vcf.gz` used only as an interval list), so it never consults a
# genotype call and reports a row for every panel site. The alternative
# `50K/results/allelic_counts/*_allele_counts.tsv` is `AD` read out of the joint-called
# `cohort.vcf.gz`, i.e. mediated by `bcftools call`, and the two disagree on 19.5% of covered sites
# (measured, scripts/measure_lineage_disagreement.R; GATK counts ~10% more reads because bcftools
# mpileup drops marginal reads via BAQ/quality filters). Everything in this repo standardizes on the
# BAM tally. See DATA.md "The SNP50K panel".
#
# Logic follows zealtiger `extract_check_counts_50K.sh` (same awk shape, same ALT filler rule),
# extended with the 49,002-site filter and an explicit sample-grouping assertion.
#
# Verify output with scripts/verify_snp50k_extractor.sh (asserts REF counts are byte-identical to the
# GATK-counts-with-filler-alleles tree; only allele labels + a few zeroed non-panel ALT reads change).
#
# ONE streaming pass over a 2.5 GB table, so call it ONCE with the union of every sample you need.
# Idempotent: samples whose output already exists are skipped, and if all are present the 2.5 GB
# read is skipped entirely.
#
#   bash scripts/extract_gatk_counts_50k.sh <sample_list.txt> <out_dir>
#     sample_list.txt : one PN#_SID# skim prefix per line
#     out_dir         : receives <prefix>.tsv (6 cols, headerless, 49,002 rows)
#   env: BZEASEQ_DIR (default /Volumes/rsstu/users/r/rrellan/BZea/bzeaseq)
#        PANEL       (default data/zeal/snp50k_alleles.tsv; the 50K-set site+allele authority
#                     built from bzea_50K_cohort.vcf.gz by scripts/build_snp50k_alleles.sh)
set -uo pipefail

LIST=${1:?usage: extract_gatk_counts_50k.sh <sample_list.txt> <out_dir>}
DEST=${2:?usage: extract_gatk_counts_50k.sh <sample_list.txt> <out_dir>}
BZ=${BZEASEQ_DIR:-/Volumes/rsstu/users/r/rrellan/BZea/bzeaseq}
IN="$BZ/50K/allelic_counts50K.tsv"
PANEL=${PANEL:-data/zeal/snp50k_alleles.tsv}   # chr pos ref alt, 49,002 rows; scripts/build_snp50k_alleles.sh

[ -r "$LIST" ]  || { echo "ERROR: cannot read sample list $LIST" >&2; exit 1; }
[ -r "$PANEL" ] || { echo "ERROR: cannot read $PANEL -- run scripts/build_snp50k_alleles.sh" >&2; exit 1; }
mkdir -p "$DEST"

# which samples are still missing?
missing=$(while read -r s; do [ -n "$s" ] && [ ! -s "$DEST/$s.tsv" ] && echo "$s"; done < "$LIST")
n_want=$(grep -cve '^[[:space:]]*$' "$LIST")
if [ -z "$missing" ]; then
  echo "[gatk50k] all $n_want samples already extracted in $DEST; skipping the 2.5 GB read"
  exit 0
fi
n_miss=$(printf '%s\n' "$missing" | grep -cve '^[[:space:]]*$')
echo "[gatk50k] $n_miss of $n_want samples missing; streaming $IN"
[ -r "$IN" ] || { echo "ERROR: cannot read $IN (mount down? set BZEASEQ_DIR)" >&2; exit 1; }

WANT=$(mktemp); printf '%s\n' "$missing" > "$WANT"
trap 'rm -f "$WANT"' EXIT

# pass 1: wanted samples; pass 2: the 49,002 panel sites; pass 3: the big table.
# The table is grouped by sample; we assert that rather than trust it, because a regrouped file
# would silently truncate every per-sample output.
# The alleles come from the PANEL (= the 50K-set authority, bzea_50K_cohort.vcf.gz), not from GATK.
# CollectAllelicCounts infers the alt base from the reads, so at 0.4x it writes N almost everywhere
# (and those rows always have ALT_COUNT 0), and where it does report a base that base is sometimes
# not the set's alt allele at all -- roughly 16 sites per sample, 1 to 2 reads each, i.e. sequencing
# error counted as donor evidence. So:
#   REF, ALT  <- panel  (no filler letter anywhere)
#   REF_COUNT <- as measured
#   ALT_COUNT <- as measured only when GATK's base IS the panel's alt; otherwise 0, because the
#                reads it counted were of some other base and the true count is not in this table.
awk -F'\t' -v DEST="$DEST" '
  FILENAME == ARGV[1] { if ($1 != "") want[$1] = 1; next }
  FILENAME == ARGV[2] { pref[$1 SUBSEP $2] = $3; palt[$1 SUBSEP $2] = $4; next }
  FNR == 1 { next }                                    # header of the big table
  {
    s = $1
    if (!(s in want)) next
    if (s != cur) {
      if (cur != "") { close(curfile); done[cur] = 1 }
      if (s in done) { print "ERROR: table is not grouped by sample (" s " seen again)" > "/dev/stderr"; exit 3 }
      cur = s; curfile = DEST "/" s ".tsv"; printf "" > curfile; nsamp++
    }
    k = $2 SUBSEP $3
    if (!(k in palt)) next                             # off-panel site
    if ($6 != pref[k]) { nrefbad++ }                   # should never happen; reported below
    ac = $5
    if ($7 != "N" && $7 != palt[k]) { if (ac > 0) { zeroed++; zreads += ac }; ac = 0 }
    print $2 "\t" $3 "\t" pref[k] "\t" $4 "\t" palt[k] "\t" ac >> curfile
  }
  END {
    printf("[gatk50k] wrote %d sample files\n", nsamp + 0) > "/dev/stderr"
    printf("[gatk50k] alleles taken from the panel; zeroed %d site(s) totalling %d read(s) where GATK counted a non-panel base\n", zeroed + 0, zreads + 0) > "/dev/stderr"
    if (nrefbad + 0 > 0) printf("[gatk50k] WARNING: %d row(s) where GATK REF != panel REF\n", nrefbad) > "/dev/stderr"
  }
' "$WANT" "$PANEL" "$IN" || exit $?

bad=0
while read -r s; do
  [ -z "$s" ] && continue
  r=$(wc -l < "$DEST/$s.tsv" 2>/dev/null | tr -d ' ')
  if [ "${r:-0}" != "49002" ]; then echo "  WARN $s: ${r:-0} rows (expected 49002)"; bad=$((bad + 1)); fi
done < "$WANT"
echo "[gatk50k] verified; $bad file(s) off the expected 49002 rows"