#!/usr/bin/env bash
# ALDEx2 week-9 WS prev2of9 — genus aggregate, plastic_level + MP retention.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
IMAGE="${ALDEX_DOCKER_IMAGE:-zjardyn/plastic-amplicon:latest}"
MC_SAMPLES="${ALDEX_MC_SAMPLES:-128}"
DATASET="${ALDEX_DATASET:-both}"
TAX_LEVEL="${ALDEX_TAX_LEVEL:-Genus}"
OUT_DIR="${ALDEX_OUT_DIR:-output}"
N_CL="${ALDEX_N_CL:-8}"

docker_aldex() {
  docker run --rm \
    -v "$ROOT:$ROOT" \
    -w "$ROOT" \
    -e "ALDEX_MC_SAMPLES=$MC_SAMPLES" \
    -e "ALDEX_DATASET=$DATASET" \
    -e "ALDEX_TAX_LEVEL=$TAX_LEVEL" \
    -e "ALDEX_OUT_DIR=$OUT_DIR" \
    "$IMAGE" \
    "$@"
}

echo "=== ALDEx2 week-9 WS | tax=$TAX_LEVEL | mc.samples=$MC_SAMPLES | dataset=$DATASET ==="
docker_aldex Rscript R/15_aldex_week9_ws_prev2of9.R

echo "=== ALDEx2 figures ==="
docker_aldex Rscript R/17_aldex_figures_week9_ws_prev2of9.R

echo "Done."
