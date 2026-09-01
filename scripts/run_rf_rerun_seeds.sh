#!/usr/bin/env bash
# Rerun RF (10k trees, log10 response) + regenerate figures for both suites.
# Same seeds: LOO fold i = 424242+i, full-n = 524242.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
IMAGE="${RF_DOCKER_IMAGE:-zjardyn/plastic-amplicon:latest}"
LOO_SEED="${RF_LOO_SEED_BASE:-424242}"
FULL_SEED="${RF_FULL_SEED:-524242}"
N_CL="${RF_N_CL:-8}"
TREES="${RF_NUM_TREES:-10000}"

docker_rf() {
  docker run --rm \
    -v "$ROOT:$ROOT" \
    -w "$ROOT" \
    -e "RF_N_CL=$N_CL" \
    -e "RF_NUM_TREES=$TREES" \
    -e "RF_LOO_SEED_BASE=$LOO_SEED" \
    -e "RF_FULL_SEED=$FULL_SEED" \
    "$IMAGE" \
    "$@"
}

ensure_rf_fig_packages() {
  docker_rf R -q -e '
    pkgs <- c("ComplexUpset", "ggVennDiagram")
    miss <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
    if (length(miss)) {
      install.packages(miss, repos = "https://cloud.r-project.org", quiet = TRUE)
    }
    stopifnot(vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1)))
  '
}

fig_suite() {
  ensure_rf_fig_packages
  local out_dir="$1"
  local fig_dir="$2"
  local suite="$3"
  docker_rf env \
    RF_OUT_DIR="$out_dir" \
    RF_10K_FIG_SUITE="$suite" \
    RF_10K_FIG_OUT_DIR="$fig_dir" \
    RF_10K_FIG_MIRROR_DIR= \
    Rscript R/23_rf_loocv_importance_10k_figures.R
  docker_rf env \
    RF_OUT_DIR="$out_dir" \
    RF_RANK_FIG_DIR="$fig_dir" \
    Rscript R/24_rf_clr1_rank_lineplots.R
  docker_rf env \
    RF_OUT_DIR="$out_dir" \
    RF_IMP_FIG_DIR="$fig_dir" \
    Rscript R/25_rf_clr_pc_importance_tracks.R
  docker_rf env \
    RF_OUT_DIR="$out_dir" \
    RF_IMP_DIFF_FIG_DIR="$fig_dir" \
    Rscript R/26_rf_clr_pc_importance_diff_histograms.R
  docker_rf env \
    RF_OUT_DIR="$out_dir" \
    RF_IMP_DIFF_RANK_FIG_DIR="$fig_dir" \
    Rscript R/27_rf_clr_pc_importance_diff_by_rank.R
}

echo "=== RF suite: all (clr 1/0.5/0.1, rclr x2, pa) ==="
docker_rf env \
  RF_OUT_DIR=output/rf_loocv_importance_all \
  RF_TRANSFORMS=clr_1,clr_0.5,clr_0.1,rclr,rclr_optspace,pa \
  Rscript R/14_random_forests_week9_ws_prev2of9.R

echo "=== RF suite: clr_pc (1, 0.5, 0.1) ==="
docker_rf env \
  RF_OUT_DIR=output/rf_loocv_importance_clr_pc \
  RF_TRANSFORMS=clr_1,clr_0.5,clr_0.1 \
  Rscript R/14_random_forests_week9_ws_prev2of9.R

echo "=== Figures: rf_loocv_importance_all ==="
fig_suite output/rf_loocv_importance_all figures/rf_loocv_importance_all all

echo "=== Figures: rf_loocv_importance_clr_pc ==="
fig_suite output/rf_loocv_importance_clr_pc figures/rf_loocv_importance_clr_pc clr_pc

echo "Done."
