#!/usr/bin/env bash
set -euo pipefail
# Copy to config/analysis.local.sh, edit paths, run from the repository root.
: "${SCHRODINGER:?Load the licensed Schrodinger environment first}"
GLIDE_INPUT="/private/docking/truncated_hits..maegz"
GLIDE_OUTPUT="results/glide100_run1"
GLIDE_EF_INPUT="/private/docking/full_poses_pv.maegz"
GLIDE_WITH_EF=0
args=(--input "$GLIDE_INPUT" --out "$GLIDE_OUTPUT"
      --n-clusters 100 --score-top 2000 --le-top 500)
if [[ "$GLIDE_WITH_EF" == "1" ]]; then
  args+=(--with-ef --ef-input "$GLIDE_EF_INPUT")
  # Add --chain A or --receptor /private/receptor.mae only after checking inputs.
fi
"$SCHRODINGER/run" python3 scripts/run_analysis.py "${args[@]}"
