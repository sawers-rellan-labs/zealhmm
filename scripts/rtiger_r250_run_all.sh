#!/bin/bash
# Reproduce the ENTIRE RTIGER r=250 benchmark dataset and its two figures
# (nilhmm-paper/figures/rtiger_{equivalence,marker_scaling}_r250.png, supplement Sec. 2)
# end-to-end from the staged shared-3 Arabidopsis panel. Every step is a tracked
# script; run order matters, because the C++ bench reads the original-Julia dumps
# this driver generates first. This REPLACES the ephemeral agent/ drivers so the
# whole pipeline is traceable from tracked code.
#
#   bash scripts/rtiger_r250_run_all.sh
#
# Inputs (gitignored, staged from the fork; see DATA.md):
#   data/rtiger_shared3_input/          shared-3 Col x Ler panel
#   ~/repos/rtiger-fork-assets/.../orig_base   the upstream-original Julia core
# Runtime: dominated by the original Julia at r=250 (the full 109,703 panel is ~43 min
# to convergence); the C++ steps are seconds. Expect ~1.5 h total.
set -euo pipefail
cd "$(dirname "$0")/.."                       # repo root
ROOT="$(pwd)"
ORIG="${RTIGER_ORIG_BASE:-$HOME/repos/rtiger-fork-assets/agent/scale_check/orig_base}"
RIG=250
BENCH="$ROOT/results/bench"
mkdir -p "$BENCH/orig_conv_r$RIG"

echo "== 1/6 SPLIT: materialize thinned panels (every benchmark reads only its split) =="
Rscript scripts/materialize_thinned_panel.R

echo "== 2/6 original Julia to convergence at r=$RIG, all five sizes =="
for lv in 0 1 2 3 4; do
  echo ">>> conv orig L$lv (r=$RIG)"
  julia scripts/rtiger_julia_conv_worker.jl "$ORIG" "$lv" "$RIG"
done

echo "== 3/6 original Julia peak RSS at r=$RIG, all five sizes =="
JMEM="$BENCH/rtiger_julia_memory_markers_r$RIG.csv"
echo "markers,peak_rss_mib" > "$JMEM"
for lv in 4 3 2 1 0; do
  echo ">>> mem orig L$lv (r=$RIG)"
  out="$(/usr/bin/time -l julia scripts/rtiger_julia_mem_worker.jl "$lv" "$RIG" 2>&1)"
  mk="$(printf '%s\n' "$out" | sed -n 's/.*markers=\([0-9]*\).*/\1/p' | head -1)"
  by="$(printf '%s\n' "$out" | grep 'maximum resident set size' | sed -E 's/^[[:space:]]*([0-9]+).*/\1/')"
  awk -v m="$mk" -v b="$by" 'BEGIN{printf "%s,%.1f\n", m, b/1048576}' >> "$JMEM"
done
echo "wrote $JMEM"

echo "== 4/6 C++ core: equivalence + throughput at r=$RIG (reads step-1 orig dumps) =="
Rscript scripts/bench_rtiger_r250.R

echo "== 5/6 C++ peak RSS at r=$RIG (thinned inputs) =="
RTIGER_RIG="$RIG" Rscript scripts/bench_rtiger_cpp_memory_markers.R

echo "== 6/6 figures =="
Rscript scripts/rtiger_r250_plot.R

echo "== DONE: RTIGER r=$RIG dataset + figures reproduced from tracked scripts =="
