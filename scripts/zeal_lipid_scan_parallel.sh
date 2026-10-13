#!/usr/bin/env bash
# Process-level parallel driver for the ZEAL lipid QTL scan (avoids R/qtl's
# fork-unsafe n.cluster on macOS): prep -> N independent chunk processes -> merge.
# Each chunk is a separate Rscript (n.cluster=1) scanning a slice of the 133 traits
# with a fixed seed (shared permutation shuffles). Usage: bash scripts/zeal_lipid_scan_parallel.sh [N]
set -euo pipefail
cd "$(dirname "$0")/.."
SCAN=scripts/zeal_lipid_rqtl_scan.R
TAG="${TAG:-taxon}"
if [ "$TAG" = "taxon" ]; then CH=results/sim/zeal/rqtl/lipid/chunks; else CH=results/sim/zeal/rqtl/lipid/chunks_$TAG; fi
export TAG
NCPU=$(sysctl -n hw.ncpu 2>/dev/null || nproc)
N="${1:-$(( NCPU > 3 ? NCPU - 2 : 1 ))}"; [ "$N" -gt 12 ] && N=12
echo "=== parallel lipid scan: N=$N chunks on $NCPU cpus ==="
rm -f "$CH"/scanone_*.csv "$CH"/perms_*.rds 2>/dev/null || true

echo "--- prep ---"
MODE=prep Rscript "$SCAN"

echo "--- $N chunks (parallel) ---"
pids=()
for i in $(seq 1 "$N"); do
  CHUNK="$i/$N" MODE=chunk Rscript "$SCAN" > "agent/lipid_chunk_${i}.log" 2>&1 &
  pids+=("$!")
done
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=1; done
if [ "$fail" -ne 0 ]; then echo "ERROR: a chunk failed; see agent/lipid_chunk_*.log"; exit 1; fi

echo "--- merge ---"
MODE=merge Rscript "$SCAN"
echo "=== parallel lipid scan DONE ==="
