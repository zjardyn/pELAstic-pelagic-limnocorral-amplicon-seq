# Rebuild full RF analysis suite for 10k-tree runs.
#
# Default (8 transforms): figures/rf_loocv_importance_10k_8transforms/{16S,18S}/
#   Does not write ranks/. Does not mirror into figures/rf_loocv_importance/.
#
# 5-transform suite (clr_1, clr_0.5, clr_0.1, rclr, rclr_optspace):
#   RF_10K_FIG_SUITE=5 Rscript R/23_rf_loocv_importance_10k_figures.R
#   Writes figures/rf_loocv_importance_10k/{16S,18S}/ and mirrors marker
#   figures into figures/rf_loocv_importance/{16S,18S}/ (not ranks/).
#
# CLR 1/0.5/0.1 suite:
#   RF_10K_FIG_SUITE=4 Rscript R/23_rf_loocv_importance_10k_figures.R
#   Writes figures/rf_loocv_importance_10k_clr1_0.5_0.1/{16S,18S}/.
#
# Named suites (set RF_OUT_DIR to matching RDS folder):
#   RF_10K_FIG_SUITE=all RF_OUT_DIR=output/rf_loocv_importance_all \
#     RF_10K_FIG_OUT_DIR=figures/rf_loocv_importance_all
#   RF_10K_FIG_SUITE=clr_pc RF_OUT_DIR=output/rf_loocv_importance_clr_pc \
#     RF_10K_FIG_OUT_DIR=figures/rf_loocv_importance_clr_pc
#
# Optional env:
#   RF_10K_FIG_OUT_DIR, RF_10K_FIG_MIRROR_DIR (empty = no mirror)
#   RF_OUT_DIR (RDS directory, default output)
#   RF_10K_FIG_FLAT=1|0  (default 0; flat = marker in filename, no subdirs)
#
# Run: Rscript R/23_rf_loocv_importance_10k_figures.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(glue)
})

fig_suite <- Sys.getenv("RF_10K_FIG_SUITE", unset = "8")
if (!fig_suite %in% c("4", "5", "8", "all", "clr_pc")) {
  stop("RF_10K_FIG_SUITE must be 4, 5, 8, all, or clr_pc, got: ", fig_suite)
}

PANEL_TITLES <- c(
  "clr_1" = "CLR pseudocount = 1",
  "clr_0.5" = "CLR pseudocount = 0.5",
  "clr_0.1" = "CLR pseudocount = 0.1",
  "clr_0.01" = "CLR pseudocount = 0.01",
  "rclr" = "rCLR — no imputation",
  "rclr_optspace" = "rCLR — optspace",
  "gbm" = "GBM",
  "czm" = "CZM",
  "pa" = "Presence/absence"
)
UPSET_LABELS <- c(
  "clr_1" = "CLR pc=1",
  "clr_0.5" = "CLR pc=0.5",
  "clr_0.1" = "CLR pc=0.1",
  "clr_0.01" = "CLR pc=0.01",
  "rclr" = "rCLR no imp.",
  "rclr_optspace" = "rCLR optspace",
  "gbm" = "GBM",
  "czm" = "CZM",
  "pa" = "PA"
)
# Short x-axis labels for performance boxplots
BOX_LABELS <- c(
  "clr_1" = "CLR 1",
  "clr_0.5" = "CLR 0.5",
  "clr_0.1" = "CLR 0.1",
  "clr_0.01" = "CLR 0.01",
  "rclr" = "rCLR",
  "rclr_optspace" = "rCLR optspace",
  "gbm" = "GBM",
  "czm" = "CZM",
  "pa" = "PA"
)

if (identical(fig_suite, "8")) {
  TRANSFORM_IDS <- c(
    "clr_1", "clr_0.5", "clr_0.1", "clr_0.01",
    "rclr", "rclr_optspace",
    "gbm", "czm"
  )
  default_out <- "figures/rf_loocv_importance_10k_8transforms"
  default_mirror <- ""
  default_flat <- FALSE
  inter_all_tag <- "all8"
  inter_exact_all_tag <- "exact8of8"
  inter_exact_nm1_tag <- "exact7of8"
  inter_ge_tag <- "ge6of8"
  inter_ge_n <- 6L
  imp_width <- 28
} else if (identical(fig_suite, "4")) {
  TRANSFORM_IDS <- c("clr_1", "clr_0.5", "clr_0.1")
  default_out <- "figures/rf_loocv_importance_10k_clr1_0.5_0.1"
  default_mirror <- ""
  default_flat <- FALSE
  inter_all_tag <- "all3"
  inter_exact_all_tag <- "exact3of3"
  inter_exact_nm1_tag <- "exact2of3"
  inter_ge_tag <- "ge2of3"
  inter_ge_n <- 2L
  imp_width <- 14
} else if (identical(fig_suite, "all")) {
  TRANSFORM_IDS <- c(
    "clr_1", "clr_0.5", "clr_0.1",
    "rclr", "rclr_optspace",
    "pa"
  )
  default_out <- "figures/rf_loocv_importance_all"
  default_mirror <- ""
  default_flat <- FALSE
  inter_all_tag <- "all6"
  inter_exact_all_tag <- "exact6of6"
  inter_exact_nm1_tag <- "exact5of6"
  inter_ge_tag <- "ge4of6"
  inter_ge_n <- 4L
  imp_width <- 22
} else if (identical(fig_suite, "clr_pc")) {
  TRANSFORM_IDS <- c("clr_1", "clr_0.5", "clr_0.1")
  default_out <- "figures/rf_loocv_importance_clr_pc"
  default_mirror <- ""
  default_flat <- FALSE
  inter_all_tag <- "all3"
  inter_exact_all_tag <- "exact3of3"
  inter_exact_nm1_tag <- "exact2of3"
  inter_ge_tag <- "ge2of3"
  inter_ge_n <- 2L
  imp_width <- 14
} else {
  TRANSFORM_IDS <- c("clr_1", "clr_0.5", "clr_0.1", "rclr", "rclr_optspace")
  default_out <- "figures/rf_loocv_importance_10k"
  default_mirror <- "figures/rf_loocv_importance"
  default_flat <- FALSE
  inter_all_tag <- "all5"
  inter_exact_all_tag <- "exact5of5"
  inter_exact_nm1_tag <- "exact4of5"
  inter_ge_tag <- "ge3of5"
  inter_ge_n <- 3L
  imp_width <- 18
}
N_TRANSFORMS <- length(TRANSFORM_IDS)

out_dir <- Sys.getenv("RF_10K_FIG_OUT_DIR", unset = default_out)
mirror_dir <- Sys.getenv("RF_10K_FIG_MIRROR_DIR", unset = default_mirror)
flat_env <- Sys.getenv("RF_10K_FIG_FLAT", unset = "")
fig_flat <- if (nzchar(flat_env)) {
  flat_env %in% c("1", "true", "TRUE", "yes", "YES")
} else {
  default_flat
}
rf_dir <- Sys.getenv("RF_OUT_DIR", unset = "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
if (!fig_flat) {
  for (mk in c("16S", "18S")) {
    dir.create(file.path(out_dir, mk), showWarnings = FALSE, recursive = TRUE)
  }
}
if (nzchar(mirror_dir)) {
  dir.create(mirror_dir, showWarnings = FALSE, recursive = TRUE)
  if (!fig_flat) {
    for (mk in c("16S", "18S")) {
      dir.create(file.path(mirror_dir, mk), showWarnings = FALSE, recursive = TRUE)
    }
  }
}
marker_dir <- function(marker) {
  if (fig_flat) out_dir else file.path(out_dir, marker)
}

rf_path <- function(marker, transform) {
  file.path(rf_dir, sprintf("rf_%s_9_ws_n9_prev2of9_%s.rds", marker, transform))
}

needed_rds <- unlist(lapply(c("16S", "18S"), function(mk) {
  vapply(TRANSFORM_IDS, function(tr) rf_path(mk, tr), character(1))
}))
missing_rds <- needed_rds[!file.exists(needed_rds)]
if (length(missing_rds)) {
  stop(
    "Missing RF RDS (pull from hydra if needed):\n  ",
    paste(missing_rds, collapse = "\n  ")
  )
}

load_imp <- function(marker, transform) {
  x <- readRDS(rf_path(marker, transform))
  stopifnot(identical(as.integer(x$design$num_trees), 10000L) ||
              identical(x$design$num_trees, 10000))
  imp <- x$result$importance
  stopifnot("median_loo_importance" %in% names(imp) ||
              "mean_importance" %in% names(imp))
  if (!"mean_importance" %in% names(imp)) imp$mean_importance <- NA_real_
  if (!"full_n_importance" %in% names(imp)) imp$full_n_importance <- NA_real_
  imp %>%
    mutate(
      transform = transform,
      median_loo_importance = dplyr::coalesce(
        median_loo_importance, mean_importance, full_n_importance
      )
    ) %>%
    arrange(desc(median_loo_importance))
}

pos_set <- function(marker, transform) {
  imp <- load_imp(marker, transform)
  imp$asv[is.finite(imp$median_loo_importance) & imp$median_loo_importance > 0]
}

jaccard <- function(a, b) {
  a <- unique(a)
  b <- unique(b)
  u <- length(union(a, b))
  if (u == 0) return(NA_real_)
  length(intersect(a, b)) / u
}

# ---- 1xn importance panels ----
n_keep_rows <- list()
make_panel <- function(imp, panel_title, tag) {
  n_keep <- sum(is.finite(imp$median_loo_importance) & imp$median_loo_importance > 0)
  half_iqr <- dplyr::coalesce(imp$iqr_importance, 0) / 2
  df <- imp %>%
    mutate(
      rank = row_number(),
      ymin = median_loo_importance - half_iqr,
      ymax = median_loo_importance + half_iqr
    )
  df_gt0 <- dplyr::filter(df, is.finite(median_loo_importance), median_loo_importance > 0)
  ggplot(df, aes(x = rank, y = median_loo_importance)) +
    geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "grey80", alpha = 0.85) +
    geom_line(linewidth = 0.35, colour = "grey25") +
    geom_point(
      data = df_gt0,
      aes(x = rank, y = median_loo_importance),
      inherit.aes = FALSE,
      size = 0.45,
      colour = "grey20",
      alpha = 0.55
    ) +
    geom_hline(yintercept = 0, colour = "red", linetype = "dotted", linewidth = 0.5) +
    geom_vline(xintercept = n_keep + 0.5, colour = "red", linetype = "dotted", linewidth = 0.5) +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = sprintf("n = %d", n_keep),
      hjust = 1.1, vjust = 1.5, size = 2.6, colour = "red3"
    ) +
    labs(
      title = panel_title, tag = tag, x = "Rank",
      y = "Median LOO permutation importance ± IQR/2"
    ) +
    theme_bw(base_size = 8) +
    theme(
      plot.tag = element_text(face = "bold", size = 11),
      plot.title = element_text(size = 8)
    )
}

for (marker in c("16S", "18S")) {
  plots <- vector("list", length(TRANSFORM_IDS))
  n_asv <- NULL
  for (i in seq_along(TRANSFORM_IDS)) {
    tr <- TRANSFORM_IDS[[i]]
    imp <- load_imp(marker, tr)
    if (is.null(n_asv)) n_asv <- nrow(imp)
    n_keep <- sum(is.finite(imp$median_loo_importance) & imp$median_loo_importance > 0)
    n_keep_rows[[length(n_keep_rows) + 1L]] <- tibble(
      marker = marker, transform = tr, n_asv = n_asv, n_keep_gt0 = n_keep
    )
    plots[[i]] <- make_panel(imp, PANEL_TITLES[[tr]], LETTERS[[i]])
  }
  p <- wrap_plots(plots, nrow = 1) +
    plot_annotation(
      title = sprintf("%s rRNA | n = %d ASVs | %d transforms", marker, n_asv, N_TRANSFORMS),
      theme = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5))
    )
  ggsave(
    file.path(marker_dir(marker), sprintf("rf_imp_%s_all_transforms_all_taxa.pdf", marker)),
    p, width = imp_width, height = 4.5, useDingbats = FALSE
  )
  message("Wrote 1x", N_TRANSFORMS, " importance for ", marker)
}
n_keep_tbl <- bind_rows(n_keep_rows)
write.csv(
  n_keep_tbl,
  file.path(out_dir, "n_keep_median_loo_imp_gt0.csv"),
  row.names = FALSE
)

# ---- performance metrics ----
perf_rows <- list()
for (marker in c("16S", "18S")) {
  for (tr in TRANSFORM_IDS) {
    x <- readRDS(rf_path(marker, tr))
    lm <- x$result$loo_metrics
    perf_rows[[length(perf_rows) + 1L]] <- tibble(
      marker = marker,
      transform = tr,
      n_trees = x$design$num_trees,
      loo_rmse = lm$rmse,
      loo_mae = lm$mae,
      loo_r2 = lm$r2,
      baseline_rmse = lm$baseline_rmse,
      baseline_mae = lm$baseline_mae,
      baseline_r2 = lm$baseline_r2
    )
  }
}
perf <- bind_rows(perf_rows)
perf_mean <- perf %>%
  group_by(marker) %>%
  summarise(
    transform = "mean_across_transforms",
    n_trees = first(n_trees),
    loo_rmse = mean(loo_rmse),
    loo_mae = mean(loo_mae),
    loo_r2 = mean(loo_r2),
    baseline_rmse = mean(baseline_rmse, na.rm = TRUE),
    baseline_mae = mean(baseline_mae, na.rm = TRUE),
    baseline_r2 = mean(baseline_r2, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(
  bind_rows(perf, perf_mean),
  file.path(out_dir, "rf_loo_performance_metrics.csv"),
  row.names = FALSE
)

# ---- LOO fold performance boxplots ----
# Each LOO fold holds out one sample, so fold RMSE = fold MAE = |y_obs - y_pred|.
# Aggregate R² (all 9 predictions) cannot be split per fold; shown as bars.
fold_rows <- list()
for (marker in c("16S", "18S")) {
  for (tr in TRANSFORM_IDS) {
    x <- readRDS(rf_path(marker, tr))
    pred <- x$result$loo_predictions
    stopifnot(
      is.data.frame(pred),
      all(c("fold", "sample_id", "y_obs", "y_pred") %in% names(pred))
    )
    err <- as.numeric(pred$y_obs) - as.numeric(pred$y_pred)
    fold_rows[[length(fold_rows) + 1L]] <- tibble(
      marker = marker,
      transform = tr,
      fold = pred$fold,
      sample_id = pred$sample_id,
      y_obs = as.numeric(pred$y_obs),
      y_pred = as.numeric(pred$y_pred),
      abs_err = abs(err),
      fold_rmse = abs(err),
      fold_mae = abs(err)
    )
  }
}
fold_perf <- bind_rows(fold_rows)
write.csv(
  fold_perf,
  file.path(out_dir, "rf_loo_fold_performance.csv"),
  row.names = FALSE
)

transform_lab_levels <- unname(BOX_LABELS[TRANSFORM_IDS])
for (marker in c("16S", "18S")) {
  mk <- marker
  fd <- fold_perf %>%
    filter(marker == mk) %>%
    mutate(
      transform_lab = factor(
        transform,
        levels = TRANSFORM_IDS,
        labels = transform_lab_levels
      )
    )
  long_err <- bind_rows(
    fd %>% transmute(transform_lab, fold, sample_id, metric = "RMSE", value = fold_rmse),
    fd %>% transmute(transform_lab, fold, sample_id, metric = "MAE", value = fold_mae)
  ) %>%
    mutate(metric = factor(metric, levels = c("RMSE", "MAE")))

  p_err <- ggplot(long_err, aes(x = transform_lab, y = value)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, fill = "grey90", colour = "grey30") +
    geom_jitter(width = 0.12, height = 0, size = 1.6, alpha = 0.8, colour = "grey20") +
    facet_wrap(~metric, nrow = 1, scales = "free_y") +
    labs(
      title = sprintf("%s: LOOCV fold error by transform (n = %d folds)", mk, dplyr::n_distinct(fd$fold)),
      subtitle = "LOO holds out 1 sample/fold, so fold RMSE = fold MAE = |residual|",
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8))

  r2_df <- perf %>%
    filter(marker == mk) %>%
    mutate(
      transform_lab = factor(
        transform,
        levels = TRANSFORM_IDS,
        labels = transform_lab_levels
      )
    )
  p_r2 <- ggplot(r2_df, aes(x = transform_lab, y = loo_r2)) +
    geom_col(width = 0.65, fill = "grey80", colour = "grey30") +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
    labs(
      title = sprintf("%s: overall LOOCV R² by transform", mk),
      subtitle = "Single aggregate R² over all 9 LOO predictions (not fold-wise)",
      x = NULL,
      y = expression(R^2)
    ) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8))

  p_perf <- p_err / p_r2 + plot_layout(heights = c(1.15, 0.85))
  out_pdf <- file.path(
    marker_dir(mk),
    sprintf("rf_loo_performance_boxplot_%s.pdf", mk)
  )
  ggsave(out_pdf, p_perf, width = max(10, 1.15 * N_TRANSFORMS + 2), height = 8, useDingbats = FALSE)
  message("Wrote performance boxplots: ", out_pdf)
}

# ---- Jaccard triangle (>0) ----
for (marker in c("16S", "18S")) {
  sets <- lapply(TRANSFORM_IDS, function(tr) pos_set(marker, tr))
  names(sets) <- TRANSFORM_IDS
  mat <- outer(
    TRANSFORM_IDS, TRANSFORM_IDS,
    Vectorize(function(i, j) jaccard(sets[[i]], sets[[j]]))
  )
  dimnames(mat) <- list(TRANSFORM_IDS, TRANSFORM_IDS)
  write.csv(mat, file.path(marker_dir(marker), sprintf("jaccard_%s_fulln_imp_gt0.csv", marker)))

  lab <- unname(PANEL_TITLES[TRANSFORM_IDS])
  df <- as.data.frame(as.table(mat)) %>%
    rename(a = Var1, b = Var2, jaccard = Freq) %>%
    mutate(
      a = factor(a, levels = TRANSFORM_IDS, labels = lab),
      b = factor(b, levels = TRANSFORM_IDS, labels = lab),
      keep = as.integer(a) >= as.integer(b)
    ) %>%
    filter(keep)

  p_jac <- ggplot(df, aes(a, b, fill = jaccard)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", jaccard)), size = 2.6) +
    scale_fill_viridis_c(limits = c(0, 1), name = "Jaccard") +
    labs(
      title = sprintf("%s: Jaccard of ASVs with full-n importance > 0", marker),
      x = NULL, y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8))
  ggsave(
    file.path(marker_dir(marker), sprintf("jaccard_%s_gt0_triangle.pdf", marker)),
    p_jac, width = 11, height = 9, useDingbats = FALSE
  )
  message(
    marker, " set sizes >0: ",
    paste(sprintf("%s=%d", names(sets), lengths(sets)), collapse = ", ")
  )
}

# ---- intersections / UpSet / Venn ----
has_complexupset <- requireNamespace("ComplexUpset", quietly = TRUE)
has_ggvenn <- requireNamespace("ggVennDiagram", quietly = TRUE)

# Classical Venn: distinct set borders/labels (saturated Okabe–Ito so circles
# are not alike) + pale translucent region fills (not the default dark navy).
VENN_SET_PALETTE <- c(
  "#E69F00", "#56B4E9", "#009E73", "#CC79A7",
  "#0072B2", "#D55E00", "#F0E442", "#999999"
)
VENN_FILL_LOW <- "#FFFFFF"
VENN_FILL_HIGH <- "#B8D4EA" # light blue-grey; lighter than default #2C5D86
VENN_FILL_ALPHA <- 0.35
VENN_EDGE_SIZE <- 1.1
VENN_UPSET_SET_PALETTE <- VENN_SET_PALETTE
VENN_UPSET_TOP_BAR <- "#A8C5DE"
VENN_UPSET_MATRIX <- "#7EB8A0"

summary_rows <- list()
for (marker in c("16S", "18S")) {
  sets <- lapply(TRANSFORM_IDS, function(tr) pos_set(marker, tr))
  names(sets) <- TRANSFORM_IDS
  tax <- load_imp(marker, "clr_1") %>% select(asv, Genus)

  all_asvs <- unique(unlist(sets))
  memb <- sapply(sets, function(s) all_asvs %in% s)
  rownames(memb) <- all_asvs
  n_tr <- rowSums(memb)

  write_inter <- function(keep, tag) {
    keep <- unname(as.character(keep))
    if (!length(keep)) {
      wide <- tax %>% filter(FALSE)
      write.csv(
        wide,
        file.path(marker_dir(marker), sprintf("intersection_fulln_gt0_%s_%s.csv", tag, marker)),
        row.names = FALSE
      )
      summary_rows[[length(summary_rows) + 1L]] <<- tibble(
        marker = marker, set = tag, n_asv = 0L, n_genus = 0L
      )
      return(wide)
    }
    wide <- tax %>% filter(asv %in% keep)
    for (tr in TRANSFORM_IDS) {
      imp <- load_imp(marker, tr) %>% select(asv, !!tr := median_loo_importance)
      wide <- wide %>% left_join(imp, by = "asv")
    }
    which_tr <- vapply(keep, function(a) {
      paste(TRANSFORM_IDS[memb[a, ]], collapse = ";")
    }, character(1))
    wide <- wide %>%
      mutate(
        n_transforms = n_tr[match(asv, all_asvs)],
        transforms = which_tr[match(asv, keep)]
      ) %>%
      arrange(desc(n_transforms), Genus, asv)
    write.csv(
      wide,
      file.path(marker_dir(marker), sprintf("intersection_fulln_gt0_%s_%s.csv", tag, marker)),
      row.names = FALSE
    )
    summary_rows[[length(summary_rows) + 1L]] <<- tibble(
      marker = marker,
      set = tag,
      n_asv = nrow(wide),
      n_genus = n_distinct(wide$Genus)
    )
    wide
  }

  write_inter(all_asvs[n_tr == N_TRANSFORMS], inter_all_tag)
  write_inter(all_asvs[n_tr == N_TRANSFORMS], inter_exact_all_tag)
  write_inter(all_asvs[n_tr == (N_TRANSFORMS - 1L)], inter_exact_nm1_tag)
  write_inter(all_asvs[n_tr >= inter_ge_n], inter_ge_tag)

  # UpSet (ComplexUpset): full min_size>=1 intersection matrix
  if (has_complexupset) {
    upset_df <- as.data.frame(memb)
    colnames(upset_df) <- unname(UPSET_LABELS[TRANSFORM_IDS])
    for (nm in colnames(upset_df)) upset_df[[nm]] <- as.integer(upset_df[[nm]])
    upset_w <- max(14, 1.8 * N_TRANSFORMS + 4)
    pdf(
      file.path(marker_dir(marker), sprintf("upset_fulln_gt0_%s.pdf", marker)),
      width = upset_w, height = 7, useDingbats = FALSE
    )
    print(
      ComplexUpset::upset(
        upset_df,
        intersect = colnames(upset_df),
        name = "transform",
        width_ratio = 0.12,
        min_size = 1,
        stripes = ComplexUpset::upset_stripes(
          colors = c("white", "#EEF3F7")
        ),
        base_annotations = list(
          `Intersection size` = ComplexUpset::intersection_size(
            counts = TRUE,
            fill = VENN_UPSET_TOP_BAR
          )
        )
      )
    )
    dev.off()
    message("Wrote ComplexUpset for ", marker)
  } else {
    message("ComplexUpset not installed; skipped UpSet for ", marker)
  }

  # Venn / UpSet-Venn (ggVennDiagram):
  # - n<=7: classical region Venn -> venn_fulln_gt0_*.pdf
  # - any n: force_upset UpSet-style plot -> venn_upset_fulln_gt0_*.pdf
  # - n>7: classical Venn unsupported; also write upset into venn_fulln_gt0_*.pdf
  # aplot/upset_plot: ggplot2::ggsave on the raw object can write empty;
  # print() to pdf (or as.patchwork then ggsave) is required.
  if (has_ggvenn) {
    vlist <- sets
    names(vlist) <- unname(UPSET_LABELS[TRANSFORM_IDS])
    # Distinct non-empty membership patterns (avoids trailing size-0 bars)
    n_intersects_show <- {
      memb_pos <- memb[rowSums(memb) > 0, , drop = FALSE]
      if (!nrow(memb_pos)) {
        1L
      } else {
        length(unique(apply(memb_pos, 1L, function(r) paste(as.integer(r), collapse = ""))))
      }
    }

    save_aplot_pdf <- function(plot_obj, path, width, height) {
      pdf(path, width = width, height = height, useDingbats = FALSE)
      on.exit(dev.off(), add = TRUE)
      # Convert aplot -> patchwork when possible for more reliable drawing
      to_draw <- tryCatch(aplot::as.patchwork(plot_obj), error = function(e) plot_obj)
      print(to_draw)
      invisible(path)
    }

    set_pal <- VENN_SET_PALETTE[seq_along(vlist)]

    # Classical Venn when supported — distinct set borders, light region fills
    if (N_TRANSFORMS <= 7L) {
      ok_v <- tryCatch({
        p_v <- ggVennDiagram::ggVennDiagram(
          vlist,
          set_color = set_pal,
          label_alpha = 0,
          edge_size = VENN_EDGE_SIZE
        ) +
          ggplot2::scale_fill_gradient(
            low = VENN_FILL_LOW,
            high = VENN_FILL_HIGH,
            na.value = "white",
            guide = "none"
          )
        if (length(p_v$layers) >= 1L) {
          p_v$layers[[1]]$aes_params$alpha <- VENN_FILL_ALPHA
        }
        ggplot2::ggsave(
          file.path(marker_dir(marker), sprintf("venn_fulln_gt0_%s.pdf", marker)),
          p_v,
          width = 12,
          height = 10,
          useDingbats = FALSE
        )
        TRUE
      }, error = function(e) {
        message("ggVennDiagram classical Venn failed for ", marker, ": ", conditionMessage(e))
        FALSE
      })
      if (!ok_v) message("Skipped classical Venn for ", marker)
    }

    # UpSet-style Venn (explicit force_upset) — the "upset venn" figure
    ok_u <- tryCatch({
      p_u <- ggVennDiagram::ggVennDiagram(
        vlist,
        force_upset = TRUE,
        nintersects = n_intersects_show,
        order.intersect.by = "size",
        label_alpha = 0,
        relative_height = 3,
        relative_width = 0.28,
        sets.bar.color = VENN_UPSET_SET_PALETTE[seq_along(vlist)],
        top.bar.color = VENN_UPSET_TOP_BAR,
        intersection.matrix.color = VENN_UPSET_MATRIX
      )
      upset_venn_path <- file.path(
        marker_dir(marker), sprintf("venn_upset_fulln_gt0_%s.pdf", marker)
      )
      save_aplot_pdf(p_u, upset_venn_path, width = 14, height = 10)
      # For n>7, classical path is unavailable — keep legacy venn_* as a copy
      if (N_TRANSFORMS > 7L) {
        legacy_venn <- file.path(
          marker_dir(marker), sprintf("venn_fulln_gt0_%s.pdf", marker)
        )
        if (!file.copy(upset_venn_path, legacy_venn, overwrite = TRUE)) {
          stop("Failed to copy UpSet-Venn to ", legacy_venn)
        }
      }
      message("Wrote ggVennDiagram UpSet-Venn for ", marker)
      TRUE
    }, error = function(e) {
      message("ggVennDiagram UpSet-Venn failed for ", marker, ": ", conditionMessage(e))
      FALSE
    })
    if (!ok_u) message("Skipped UpSet-Venn for ", marker)
  }
}

write.csv(
  bind_rows(summary_rows),
  file.path(out_dir, "intersection_fulln_gt0_summary.csv"),
  row.names = FALSE
)

# Optional mirror of marker-level figures (not ranks/). Default: off for 8-transform.
if (nzchar(mirror_dir)) {
  mirror_pat <- paste(
    c(
      "^rf_imp_", "^jaccard_", "^upset_fulln_gt0_", "^venn_fulln_gt0_",
      "^venn_upset_fulln_gt0_",
      "^rf_loo_performance_boxplot_",
      sprintf(
        "intersection_fulln_gt0_(%s|%s|%s|%s)_",
        inter_all_tag, inter_exact_all_tag, inter_exact_nm1_tag, inter_ge_tag
      )
    ),
    collapse = "|"
  )
  for (mk in c("16S", "18S")) {
    src <- marker_dir(mk)
    dst <- file.path(mirror_dir, mk)
    files <- list.files(src, full.names = TRUE)
    files <- files[file.info(files)$isdir %in% FALSE]
    files <- files[grepl(mirror_pat, basename(files))]
    if (length(files)) {
      ok <- file.copy(files, file.path(dst, basename(files)), overwrite = TRUE)
      if (!all(ok)) warning("Some mirror copies failed for ", mk)
    }
  }
  shared <- c(
    "rf_loo_performance_metrics.csv",
    "rf_loo_fold_performance.csv",
    "intersection_fulln_gt0_summary.csv",
    "n_keep_fulln_imp_gt0.csv"
  )
  for (fn in shared) {
    src <- file.path(out_dir, fn)
    if (file.exists(src)) {
      file.copy(src, file.path(mirror_dir, fn), overwrite = TRUE)
    }
  }
  message("Mirrored marker figures to ", normalizePath(mirror_dir, mustWork = TRUE))
} else {
  message("No mirror (RF_10K_FIG_MIRROR_DIR empty).")
}

message("Done. Outputs in ", normalizePath(out_dir, mustWork = TRUE))
print(n_keep_tbl)
print(bind_rows(summary_rows))
print(perf_mean)
