#!/usr/bin/env bash
# Step 3 of the count-tree consolidation: relabel the PN-keyed GATK-counts-with-real-alleles cache
# (from scripts/extract_gatk_counts_50k.sh) into the canonical per-`sample` store.
#
#   data/zeal/gatk_counts_50k_cache/<skim_prefix>.tsv   ->   data/zeal/snp50k_counts/<sample>.tsv
#
# where <sample> is the canonical id from data/zeal/snp50k_count_roster.tsv (pedigree for NILs,
# B73_<PN>/Purple_<PN> for checks). One file per unique physical sample; every figure selects its
# subset from here by the roster's role flags. See [[50k-set-terminology]].
#
#   bash scripts/build_snp50k_count_store.sh
set -uo pipefail
ROSTER=${ROSTER:-data/zeal/snp50k_count_roster.tsv}
CACHE=${CACHE:-data/zeal/gatk_counts_50k_cache}
STORE=${STORE:-data/zeal/snp50k_counts}

[ -r "$ROSTER" ] || { echo "ERROR: no roster $ROSTER (run build_snp50k_count_roster.R)" >&2; exit 1; }
[ -d "$CACHE" ]  || { echo "ERROR: no cache $CACHE (run extract_gatk_counts_50k.sh)" >&2; exit 1; }
mkdir -p "$STORE"

n=0; miss=0; badrows=0
# columns: sample skim_prefix is_check in_dnarna in_molb in_dnasweep
while IFS=$'\t' read -r sample skim_prefix _rest; do
  [ "$sample" = "sample" ] && continue          # header
  src="$CACHE/$skim_prefix.tsv"
  dst="$STORE/$sample.tsv"
  if [ ! -s "$src" ]; then echo "  MISS $sample ($skim_prefix): no cache file"; miss=$((miss+1)); continue; fi
  cp -f "$src" "$dst"
  r=$(wc -l < "$dst" | tr -d ' ')
  [ "$r" = "49002" ] || { echo "  WARN $sample: $r rows (expected 49002)"; badrows=$((badrows+1)); }
  n=$((n+1))
done < "$ROSTER"

echo "[store] wrote $n files to $STORE ; missing cache: $miss ; wrong-row-count: $badrows"
echo "[store] store now holds $(ls "$STORE"/*.tsv 2>/dev/null | wc -l | tr -d ' ') files"
[ "$miss" = 0 ] && [ "$badrows" = 0 ] || { echo "[store] FAILED integrity" >&2; exit 2; }
echo "[store] OK: every roster sample has a 49,002-row canonical count file"
