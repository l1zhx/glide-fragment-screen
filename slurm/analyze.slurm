#!/usr/bin/env bash
#SBATCH --job-name=glide-analysis
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=glide-analysis-%j.log
set -euo pipefail
# Supply your site's partition/account/module requirements at submission time.
# The 16 GB request is an initial budget, not a measured requirement.
: "${SCHRODINGER:?Load your site's Schrodinger environment first}"
: "${GLIDE_PROJECT:?Absolute repository directory is required}"
: "${GLIDE_INPUT:?Absolute input MAEGZ path is required}"
: "${GLIDE_OUTPUT:?New absolute output directory is required}"
cd "$GLIDE_PROJECT"
export MPLBACKEND=Agg
"$SCHRODINGER/run" python3 scripts/run_analysis.py \
  --input "$GLIDE_INPUT" --out "$GLIDE_OUTPUT" \
  --n-clusters "${GLIDE_N_CLUSTERS:-100}" \
  --score-top "${GLIDE_SCORE_TOP:-2000}" --le-top "${GLIDE_LE_TOP:-500}"
