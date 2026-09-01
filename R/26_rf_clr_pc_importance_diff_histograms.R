# Pairwise full-n importance differences across CLR pseudocounts (pc = 1, 0.5, 0.1).
#
# Writes:
#   figures/rf_loocv_importance_clr_pc/{16S,18S}/rf_imp_{marker}_clr_pc_diff_histograms.pdf
#   figures/rf_loocv_importance_clr_pc/{16S,18S}/rf_imp_{marker}_clr_pc_diff_pairs.csv
#
# Optional env:
#   RF_OUT_DIR           — RDS root (default output/rf_loocv_importance_clr_pc)
#   RF_IMP_DIFF_FIG_DIR  — figure root (default figures/rf_loocv_importance_clr_pc)
#   RF_IMP_DIFF_ASVS     — any_gt0 (default) | exact3of3 | all_finite
#
# Run: Rscript R/26_rf_clr_pc_importance_diff_histograms.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(glue)
})

TRANSFORMS <- c("clr_1", "clr_0.5", "clr_0.1")
TRANSFORM_LABELS <- c(
  "clr_1" = "pc = 1",
  "clr_0.5" = "pc = 0.5",
  "clr_0.1" = "pc = 0.1"
)
PAIRS <- list(
  c("clr_1", "clr_0.5"),
  c("clr_1", "clr_0.1"),
  c("clr_0.5", "clr_0.1")
)
PAIR_LABELS <- vapply(
  PAIRS,
  function(p) glue("{TRANSFORM_LABELS[[p[1]]]} − {TRANSFORM_LABELS[[p[2]]]}"),
  character(1)
)

rf_dir <- Sys.getenv("RF_OUT_DIR", unset = "output/rf_loocv_importance_clr_pc")
fig_dir <- Sys.getenv("RF_IMP_DIFF_FIG_DIR", unset = "figures/rf_loocv_importance_clr_pc")
asv_mode <- Sys.getenv("RF_IMP_DIFF_ASVS", unset = "any_gt0")

for (mk in c("16S", "18S")) {
  dir.create(file.path(fig_dir, mk), showWarnings = FALSE, recursive = TRUE)
}

rf_path <- function(marker, transform) {
  file.path(rf_dir, sprintf("rf_%s_9_ws_n9_prev2of9_%s.rds", marker, transform))
}

load_imp <- function(marker, transform) {
  path <- rf_path(marker, transform)
  if (!file.exists(path)) stop("Missing RF RDS: ", path)
  imp <- readRDS(path)$result$importance
  if (!"median_loo_importance" %in% names(imp)) {
    imp$median_loo_importance <- dplyr::coalesce(imp$mean_importance, imp$full_n_importance)
  }
  if (!"mean_importance" %in% names(imp)) imp$mean_importance <- NA_real_
  if (!"full_n_importance" %in% names(imp)) imp$full_n_importance <- NA_real_
  imp %>%
    dplyr::mutate(
      transform = transform,
      importance = dplyr::coalesce(
        median_loo_importance, full_n_importance, mean_importance
      )
    ) %>%
    dplyr::select(asv, Genus, transform, importance)
}

select_asvs <- function(wide) {
  if (asv_mode == "all_finite") {
    return(wide$asv)
  }
  mat <- as.matrix(wide[, TRANSFORMS, drop = FALSE])
  gt0_n <- rowSums(is.finite(mat) & mat > 0, na.rm = TRUE)
  if (asv_mode == "any_gt0") {
    return(wide$asv[gt0_n > 0])
  }
  if (asv_mode == "exact3of3") {
    return(wide$asv[gt0_n == length(TRANSFORMS)])
  }
  stop("RF_IMP_DIFF_ASVS must be any_gt0, exact3of3, or all_finite, got: ", asv_mode)
}

build_wide <- function(marker) {
  parts <- lapply(TRANSFORMS, function(tr) load_imp(marker, tr))
  wide <- parts[[1]] %>%
    dplyr::select(asv, Genus) %>%
    dplyr::left_join(
      dplyr::select(parts[[1]], asv, clr_1 = importance),
      by = "asv"
    ) %>%
    dplyr::left_join(
      dplyr::select(parts[[2]], asv, clr_0.5 = importance),
      by = "asv"
    ) %>%
    dplyr::left_join(
      dplyr::select(parts[[3]], asv, clr_0.1 = importance),
      by = "asv"
    )
  keep <- select_asvs(wide)
  wide %>% dplyr::filter(asv %in% keep)
}

pairwise_diffs <- function(wide) {
  rows <- vector("list", length(PAIRS))
  for (i in seq_along(PAIRS)) {
    a <- PAIRS[[i]][1]
    b <- PAIRS[[i]][2]
    d <- wide %>%
      dplyr::filter(is.finite(.data[[a]]), is.finite(.data[[b]])) %>%
      dplyr::transmute(
        asv,
        Genus,
        transform_a = a,
        transform_b = b,
        importance_a = .data[[a]],
        importance_b = .data[[b]],
        diff = importance_a - importance_b,
        abs_diff = abs(importance_a - importance_b),
        comparison = PAIR_LABELS[[i]]
      )
    rows[[i]] <- d
  }
  dplyr::bind_rows(rows)
}

HIST_FILL <- "#7EC8E3"
HIST_COLOUR <- "grey35"

hist_panel <- function(
  df,
  x_var,
  title,
  xlab,
  facet = FALSE,
  center_zero = FALSE,
  show_zero_line = FALSE
) {
  x_vals <- df[[x_var]]
  if (center_zero) {
    xlim <- max(abs(x_vals), na.rm = TRUE)
    if (!is.finite(xlim) || xlim <= 0) xlim <- 1
    xlim <- xlim * 1.05
    x_scale <- ggplot2::scale_x_continuous(
      limits = c(-xlim, xlim),
      breaks = pretty(c(-xlim, xlim), n = 7),
      expand = ggplot2::expansion(mult = c(0, 0))
    )
  } else {
    xlim <- max(x_vals, na.rm = TRUE)
    if (!is.finite(xlim) || xlim <= 0) xlim <- 1
    xlim <- xlim * 1.05
    x_scale <- ggplot2::scale_x_continuous(
      limits = c(0, xlim),
      breaks = pretty(c(0, xlim), n = 7),
      expand = ggplot2::expansion(mult = c(0, 0))
    )
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_var]])) +
    ggplot2::geom_histogram(
      fill = HIST_FILL,
      colour = HIST_COLOUR,
      linewidth = 0.15,
      bins = 40
    ) +
    x_scale +
    ggplot2::labs(title = title, x = xlab, y = "Count") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (show_zero_line) {
    p <- p + ggplot2::geom_vline(xintercept = 0, colour = "grey25", linewidth = 0.55)
  }

  if (facet) {
    p <- p +
      ggplot2::facet_wrap(~comparison, ncol = 1, scales = "free_y") +
      ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", size = 9))
  }

  p
}

build_hist_plot <- function(marker, diff_df) {
  diff_df <- diff_df %>%
    dplyr::mutate(
      comparison = factor(comparison, levels = PAIR_LABELS)
    )

  p_signed_all <- hist_panel(
    diff_df,
    x_var = "diff",
    title = "Signed difference — all pairs combined",
    xlab = "Importance difference",
    center_zero = TRUE,
    show_zero_line = TRUE
  )

  p_signed_split <- hist_panel(
    diff_df,
    x_var = "diff",
    title = "Signed difference — by pair",
    xlab = "Importance difference",
    facet = TRUE,
    center_zero = TRUE,
    show_zero_line = TRUE
  )

  p_abs_all <- hist_panel(
    diff_df,
    x_var = "abs_diff",
    title = "Absolute difference — all pairs combined",
    xlab = "|Importance difference|",
    center_zero = FALSE,
    show_zero_line = FALSE
  )

  p_abs_split <- hist_panel(
    diff_df,
    x_var = "abs_diff",
    title = "Absolute difference — by pair",
    xlab = "|Importance difference|",
    facet = TRUE,
    center_zero = FALSE,
    show_zero_line = FALSE
  )

  patchwork::wrap_plots(
    p_signed_all,
    p_signed_split,
    p_abs_all,
    p_abs_split,
    ncol = 1,
    heights = c(1, 1.4, 1, 1.4)
  ) +
    patchwork::plot_annotation(
      title = glue(
        "{marker} | CLR pseudocount importance differences | ASVs: {asv_mode} | n = {dplyr::n_distinct(diff_df$asv)}"
      ),
      subtitle = paste(
        "Top: all pairwise differences pooled; bottom: same data split by pc pair;",
        "signed panels centred at 0",
        sep = " "
      ),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 12),
        plot.subtitle = ggplot2::element_text(size = 8.5)
      )
    )
}

for (marker in c("16S", "18S")) {
  wide <- build_wide(marker)
  diff_df <- pairwise_diffs(wide)

  csv_out <- file.path(
    fig_dir, marker,
    sprintf("rf_imp_%s_clr_pc_diff_pairs.csv", marker)
  )
  utils::write.csv(diff_df, csv_out, row.names = FALSE)

  p <- build_hist_plot(marker, diff_df)
  pdf_out <- file.path(
    fig_dir, marker,
    sprintf("rf_imp_%s_clr_pc_diff_histograms.pdf", marker)
  )
  ggplot2::ggsave(pdf_out, p, width = 9, height = 14, useDingbats = FALSE)
  message("Wrote ", pdf_out, " (", nrow(diff_df), " pair-rows)")
  message("Wrote ", csv_out)
}

message("Done.")
