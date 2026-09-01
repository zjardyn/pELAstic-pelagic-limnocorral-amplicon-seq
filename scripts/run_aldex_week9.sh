#!/usr/bin/env bash
# ALDEx2 week-9 WS prev2of9 — plastic_level + MP retention.
# Default: Genus. For ASV (RF-matched prev2of9): ALDEX_TAX_LEVEL=ASV bash scripts/run_aldex_week9.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
IMAGE="${ALDEX_DOCKER_IMAGE:-zjardyn/plastic-amplicon:latest}"
MC_SAMPLES="${ALDEX_MC_SAMPLES:-128}"
DATASET="${ALDEX_DATASET:-both}"
TAX_LEVEL="${ALDEX_TAX_LEVEL:-Genus}"
MIN_READS="${ALDEX_MIN_READS:-3}"
MIN_SAMPLES="${ALDEX_MIN_SAMPLES:-2}"
OUT_DIR="${ALDEX_OUT_DIR:-output}"

if [[ "${TAX_LEVEL^^}" == "ASV" ]]; then
  FIG_DIR="${ALDEX_FIG_DIR:-figures/aldex_asv}"
  GLM_FIG_DIR="${ALDEX_GLM_FIG_DIR:-figures/aldex_glm_asv}"
else
  FIG_DIR="${ALDEX_FIG_DIR:-figures/aldex_genus}"
  GLM_FIG_DIR="${ALDEX_GLM_FIG_DIR:-figures/aldex_glm}"
fi

docker_aldex() {
  docker run --rm \
    -v "$ROOT:$ROOT" \
    -w "$ROOT" \
    -e "ALDEX_MC_SAMPLES=$MC_SAMPLES" \
    -e "ALDEX_DATASET=$DATASET" \
    -e "ALDEX_TAX_LEVEL=$TAX_LEVEL" \
    -e "ALDEX_MIN_READS=$MIN_READS" \
    -e "ALDEX_MIN_SAMPLES=$MIN_SAMPLES" \
    -e "ALDEX_OUT_DIR=$OUT_DIR" \
    -e "ALDEX_FIG_DIR=$FIG_DIR" \
    -e "ALDEX_GLM_FIG_DIR=$GLM_FIG_DIR" \
    "$IMAGE" \
    "$@"
}

echo "=== ALDEx2 week-9 WS | tax=$TAX_LEVEL | mc.samples=$MC_SAMPLES | filter=>=${MIN_READS} in >=${MIN_SAMPLES}/9 | dataset=$DATASET ==="
docker_aldex Rscript R/15_aldex_week9_ws_prev2of9.R

echo "=== ALDEx2 summary figures ($FIG_DIR) ==="
docker_aldex Rscript R/17_aldex_figures_week9_ws_prev2of9.R

echo "=== ALDEx2 GLM / sample-level figures ($GLM_FIG_DIR) ==="
docker_aldex Rscript R/18_aldex_glm_sample_figures_week9_ws_prev2of9.R

echo "Done."
