#!/bin/bash
# Build the minimal R + data.table conda env that the missing-data summarization
# job (make_missing_data_summaries.lsf) activates.
#
# Run on the Hazel LOGIN node — `conda create` needs internet and compute nodes
# are isolated. The env lives in fast /share scratch (ephemeral, ~1 month) — that
# is fine, it only needs to exist while the summaries are (re)built. The full
# bzeaseq.yml spec is no longer dependency-solvable, so we build just what the
# reduction needs (data.table).
#
# Parameters (env vars, with defaults):
#   CONDA_BASE  conda install to source (default: Hazel system miniconda)
#   RDT_ENV     env prefix to create   (default: /share/<group>/<user>/conda/env/rdt)
#
# Usage:  bash scripts/build_missing_data_env.sh
set -uo pipefail

CONDA_BASE="${CONDA_BASE:-/usr/local/apps/conda/miniconda3/26.3.2}"
RDT_ENV="${RDT_ENV:-/share/$(id -gn)/$(whoami)/conda/env/rdt}"

source "$CONDA_BASE/etc/profile.d/conda.sh"
echo "start: $(date)  host: $(hostname)"
echo "conda: $(conda --version)   env prefix: $RDT_ENV"

# clean any partial prefix from a previous failed attempt
[ -e "$RDT_ENV" ] && conda env remove -y -p "$RDT_ENV" 2>/dev/null

t0=$(date +%s)
conda create -y -p "$RDT_ENV" -c conda-forge r-base r-data.table
rc=$?
echo "=== rc=$rc  elapsed=$(( $(date +%s) - t0 ))s ==="
if [ "$rc" -eq 0 ]; then
  "$RDT_ENV"/bin/Rscript -e 'cat(R.version.string, "| data.table", as.character(packageVersion("data.table")), "\n")'
fi
echo "FINISHED: $(date)"