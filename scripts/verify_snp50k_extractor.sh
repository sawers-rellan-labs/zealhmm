#!/usr/bin/env bash
# Traceability check for scripts/extract_gatk_counts_50k.sh.
#
# CLAIM the extractor must satisfy: attaching real alleles from the 50K-set authority
# (data/snp50k_alleles.tsv, built from bzea_50K_cohort.vcf.gz) changes the ALLELE LABELS and zeroes
# a few non-panel ALT reads, but NEVER changes the read COUNTS that GATK measured. So versus the
# GATK-counts-with-filler-alleles tree (zealtiger's reformat, which carried filler ALT letters), on
# the shared samples, joined on (chr,pos):
#
#     n_ref differs      == 0        <- INVARIANT. GATK REF depth is untouched. Non-zero => a bug.
#     n_alt differs      == small    <- only the non-panel-base reads the extractor deliberately
#                                       zeroed (residual sequencing error; ~10-30 sites/sample).
#     alt-label differs  == large    <- EXPECTED & desired: filler letters (A, or C when REF=A)
#                                       replaced by the true allele from the 50K set (~34k/sample).
#
# This is the join that validated the extractor on 2026-08-03 (n_ref 0, n_alt 17, alt-label 33638 on
# PN3_SID225). Kept as a script so the check is reproducible, not a lost terminal one-liner.
#
#   bash scripts/verify_snp50k_extractor.sh <new_extractor_dir> [filler_allele_dir]
#     new_extractor_dir : output of extract_gatk_counts_50k.sh (<sample>.tsv, 6-col headerless)
#     filler_allele_dir : default data/skimsweep/skim/counts_50k (GATK counts with filler alleles)
#
# Both inputs are 6-col headerless:  chr  pos  ref  n_ref  alt  n_alt
# Exit 0 only if n_ref differs == 0 for every shared sample.
set -uo pipefail

NEW=${1:?usage: verify_snp50k_extractor.sh <new_extractor_dir> [filler_allele_dir]}
OLD=${2:-data/skimsweep/skim/counts_50k}
[ -d "$NEW" ] || { echo "ERROR: no such dir $NEW" >&2; exit 1; }
[ -d "$OLD" ] || { echo "ERROR: no such dir $OLD" >&2; exit 1; }

printf '%-16s %10s %10s %10s %12s\n' sample sites n_ref_diff n_alt_diff alt_label_diff
fail=0; nshared=0
for f in "$NEW"/*.tsv; do
  [ -e "$f" ] || continue
  s=$(basename "$f" .tsv)
  o="$OLD/$s.tsv"
  [ -f "$o" ] || continue                      # only samples present in BOTH trees
  nshared=$((nshared + 1))
  read -r sites nref nalt altlab < <(
    join -t$'\t' \
      <(awk -F'\t' '{print $1"_"$2"\t"$4"\t"$5"\t"$6}' "$f" | sort) \
      <(awk -F'\t' '{print $1"_"$2"\t"$4"\t"$5"\t"$6}' "$o" | sort) \
    | awk -F'\t' '
        { n++
          if ($2 != $5) nref++       # n_ref  (col 4 in the 6-col files)
          if ($4 != $7) nalt++       # n_alt  (col 6)
          if ($3 != $6) altlab++ }   # alt base label (col 5)
        END { printf "%d %d %d %d\n", n+0, nref+0, nalt+0, altlab+0 }'
  )
  printf '%-16s %10s %10s %10s %12s' "$s" "$sites" "$nref" "$nalt" "$altlab"
  if [ "$nref" != 0 ]; then printf '   <-- FAIL: n_ref must be 0\n'; fail=$((fail + 1)); else printf '\n'; fi
done

echo
if [ "$nshared" = 0 ]; then
  echo "ERROR: no samples shared between $NEW and $OLD" >&2; exit 1
fi
if [ "$fail" = 0 ]; then
  echo "PASS: n_ref identical on all $nshared shared sample(s); counts preserved, only alleles corrected."
else
  echo "FAIL: $fail of $nshared sample(s) changed GATK REF depth -- the extractor altered counts." >&2
  exit 2
fi
