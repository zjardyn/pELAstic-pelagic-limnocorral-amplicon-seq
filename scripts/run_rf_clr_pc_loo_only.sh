#!/usr/bin/env bash
# LOOCV-only RF (median fold importance) + clr_pc figures (clr_1, clr_0.5, clr_0.1).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
IMAGE="${RF_DOCKER_IMAGE:-zjardyn/plastic-amplicon:latest}"
LOO_SEED="${RF_LOO_SEED_BASE:-424242}"
N_CL="${RF_N_CL:-8}"
TREES="${RF_NUM_TREES:-10000}"

docker_rf() {
  docker run --rm \
    -v "$ROOT:$ROOT" \
    -w "$ROOT" \
    -e "RF_N_CL=$N_CL" \
    -e "RF_NUM_TREES=$TREES" \
    -e "RF_LOO_SEED_BASE=$LOO_SEED" \
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

echo "=== RF clr_pc (LOOCV only, median importance) ==="
docker_rf env \
  RF_OUT_DIR=output/rf_loocv_importance_clr_pc \
  RF_TRANSFORMS=clr_1,clr_0.5,clr_0.1 \
  Rscript R/14_random_forests_week9_ws_prev2of9.R

fig_suite() {
  ensure_rf_fig_packages
  docker_rf env \
    RF_OUT_DIR=output/rf_loocv_importance_clr_pc \
    RF_10K_FIG_SUITE=clr_pc \
    RF_10K_FIG_OUT_DIR=figures/rf_loocv_importance_clr_pc \
    RF_10K_FIG_MIRROR_DIR= \
    Rscript R/23_rf_loocv_importance_10k_figures.R
  docker_rf env \
    RF_OUT_DIR=output/rf_loocv_importance_clr_pc \
    RF_RANK_FIG_DIR=figures/rf_loocv_importance_clr_pc \
    Rscript R/24_rf_clr1_rank_lineplots.R
  docker_rf env \
    RF_OUT_DIR=output/rf_loocv_importance_clr_pc \
    RF_IMP_FIG_DIR=figures/rf_loocv_importance_clr_pc \
    Rscript R/25_rf_clr_pc_importance_tracks.R
  docker_rf env \
    RF_OUT_DIR=output/rf_loocv_importance_clr_pc \
    RF_IMP_DIFF_FIG_DIR=figures/rf_loocv_importance_clr_pc \
    Rscript R/26_rf_clr_pc_importance_diff_histograms.R
  docker_rf env \
    RF_OUT_DIR=output/rf_loocv_importance_clr_pc \
    RF_IMP_DIFF_RANK_FIG_DIR=figures/rf_loocv_importance_clr_pc \
    Rscript R/27_rf_clr_pc_importance_diff_by_rank.R
}

echo "=== Figures: rf_loocv_importance_clr_pc ==="
fig_suite

echo "Done."
