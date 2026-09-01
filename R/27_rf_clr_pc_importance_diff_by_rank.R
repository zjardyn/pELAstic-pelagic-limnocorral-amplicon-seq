# Importance difference by rank for each CLR pc pair (gt0 in both transforms of pair).
# Importance is normalized per transform: value / sum(importance > 0) before diffs.
#
# Writes figures/rf_loocv_importance_clr_pc/diff_by_rank/{16S,18S}/
#   rf_imp_{marker}_{pair}_diff_by_rank.pdf
#   rf_imp_{marker}_abs_diff_3panel.pdf
#   rf_imp_{marker}_{pair}_diff_by_importance.pdf
#   rf_imp_{marker}_abs_diff_3panel_by_importance.pdf
#
# Rank is by importance in the left transform of each pair (descending).
#
# Optional env:
#   RF_OUT_DIR              — RDS root (default output/rf_loocv_importance_clr_pc)
#   RF_IMP_DIFF_NORM        — sum_gt0 (default) | max_gt0 | none
#
# Run: Rscript R/27_rf_clr_pc_importance_diff_by_rank.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
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
  clr_1_clr_0.5 = c("clr_1", "clr_0.5"),
  clr_1_clr_0.1 = c("clr_1", "clr_0.1"),
  clr_0.5_clr_0.1 = c("clr_0.5", "clr_0.1")
)
PAIR_LABELS <- stats::setNames(
  vapply(
    PAIRS,
    function(p) glue("{TRANSFORM_LABELS[[p[1]]]} − {TRANSFORM_LABELS[[p[2]]]}"),
    character(1)
  ),
  names(PAIRS)
)

POINT_COLOUR <- "#7EC8E3"
BOX_FILL <- "#D9EEF7"
ABS_LABEL_THRESH <- as.numeric(Sys.getenv("RF_ABS_DIFF_LABEL_THRESH", unset = "0.005"))
imp_norm_mode <- Sys.getenv("RF_IMP_DIFF_NORM", unset = "sum_gt0")

rf_dir <- Sys.getenv("RF_OUT_DIR", unset = "output/rf_loocv_importance_clr_pc")
fig_root <- Sys.getenv("RF_IMP_DIFF_RANK_FIG_DIR", unset = "figures/rf_loocv_importance_clr_pc")
out_dir <- file.path(fig_root, "diff_by_rank")

marker_out_dir <- function(marker) {
  dir.create(file.path(out_dir, marker), showWarnings = FALSE, recursive = TRUE)
  file.path(out_dir, marker)
}

for (mk in c("16S", "18S")) {
  marker_out_dir(mk)
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
    dplyr::mutate(importance = dplyr::coalesce(
      median_loo_importance, full_n_importance, mean_importance
    )) %>%
    dplyr::select(asv, Genus, importance)
}

build_wide <- function(marker) {
  parts <- lapply(TRANSFORMS, function(tr) {
    load_imp(marker, tr) %>% dplyr::rename(!!tr := importance)
  })
  wide <- parts[[1]] %>%
    dplyr::select(asv, Genus) %>%
    dplyr::left_join(dplyr::select(parts[[1]], asv, clr_1), by = "asv") %>%
    dplyr::left_join(dplyr::select(parts[[2]], asv, clr_0.5), by = "asv") %>%
    dplyr::left_join(dplyr::select(parts[[3]], asv, clr_0.1), by = "asv")
  normalize_wide(wide)
}

normalize_vector <- function(x, mode = imp_norm_mode) {
  if (identical(mode, "none")) return(x)
  gt0 <- is.finite(x) & x > 0
  if (identical(mode, "sum_gt0")) {
    denom <- sum(x[gt0], na.rm = TRUE)
  } else if (identical(mode, "max_gt0")) {
    denom <- max(x[gt0], na.rm = TRUE)
  } else {
    stop("RF_IMP_DIFF_NORM must be sum_gt0, max_gt0, or none, got: ", mode)
  }
  if (!is.finite(denom) || denom <= 0) return(x)
  x / denom
}

normalize_wide <- function(wide) {
  out <- wide
  for (tr in TRANSFORMS) {
    out[[tr]] <- normalize_vector(out[[tr]])
  }
  out
}

norm_note <- function() {
  switch(
    imp_norm_mode,
    sum_gt0 = "importance ÷ Σ(gt0) per transform",
    max_gt0 = "importance ÷ max(gt0) per transform",
    none = "raw importance"
  )
}

pair_diff_ranked <- function(wide, pair_name) {
  p <- PAIRS[[pair_name]]
  a <- p[1]
  b <- p[2]
  wide %>%
    dplyr::filter(
      is.finite(.data[[a]]), is.finite(.data[[b]]),
      .data[[a]] > 0, .data[[b]] > 0
    ) %>%
    dplyr::transmute(
      asv,
      Genus,
      transform_a = a,
      transform_b = b,
      importance_a = .data[[a]],
      importance_b = .data[[b]],
      diff = importance_a - importance_b,
      abs_diff = abs(importance_a - importance_b),
      comparison = PAIR_LABELS[[pair_name]],
      pair = pair_name,
      rank = dplyr::dense_rank(dplyr::desc(importance_a))
    )
}

rank_group_label <- function(rank, width) {
  start <- ((rank - 1L) %/% width) * width + 1L
  end <- start + width - 1L
  sprintf("%d–%d", start, end)
}

add_rank_groups <- function(df, widths) {
  out <- df
  for (w in widths) {
    col <- sprintf("rank_grp_%d", w)
    out[[col]] <- factor(
      vapply(out$rank, rank_group_label, character(1), width = w),
      levels = unique(vapply(
        sort(unique(out$rank)),
        rank_group_label,
        character(1),
        width = w
      ))
    )
  }
  out
}

sym_y_limits <- function(y) {
  lim <- max(abs(y), na.rm = TRUE)
  if (!is.finite(lim) || lim <= 0) lim <- 1
  lim * 1.05
}

pos_y_limits <- function(y) {
  lim <- max(y, na.rm = TRUE)
  if (!is.finite(lim) || lim <= 0) lim <- 1
  lim * 1.05
}

y_scale <- function(y, signed = TRUE) {
  if (signed) {
    lim <- sym_y_limits(y)
    ggplot2::scale_y_continuous(limits = c(-lim, lim))
  } else {
    lim <- pos_y_limits(y)
    ggplot2::scale_y_continuous(limits = c(0, lim))
  }
}

y_label <- function(signed) {
  if (signed) "Normalized importance difference" else "|Normalized importance difference|"
}

taxon_label <- function(genus, asv) {
  g <- as.character(genus)
  bad <- is.na(g) | !nzchar(g) | g %in% c("NA", "unclassified", "Unclassified")
  ifelse(bad, asv, g)
}

add_abs_labels <- function(p, df, x_var, y_var, label_thresh) {
  if (is.null(label_thresh) || !is.finite(label_thresh) || y_var != "abs_diff") {
    return(p)
  }
  lab_df <- df[df[[y_var]] >= label_thresh, , drop = FALSE]
  if (!nrow(lab_df)) return(p)
  lab_df$label <- taxon_label(lab_df$Genus, lab_df$asv)
  p + ggrepel::geom_text_repel(
    data = lab_df,
    ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]], label = label),
    size = 2.6,
    colour = "grey20",
    max.overlaps = Inf,
    min.segment.length = 0,
    box.padding = 0.25,
    point.padding = 0.15,
    segment.size = 0.25,
    inherit.aes = FALSE
  )
}

add_imp_groups <- function(df, n_bins) {
  out <- df
  probs <- seq(0, 1, length.out = n_bins + 1L)
  br <- unique(stats::quantile(df$importance_a, probs = probs, na.rm = TRUE))
  if (length(br) < 3L) {
    out[[sprintf("imp_grp_%d", n_bins)]] <- factor("all")
    return(out)
  }
  labs <- sprintf(
    "%.4g–%.4g",
    br[-length(br)],
    br[-1]
  )
  out[[sprintf("imp_grp_%d", n_bins)]] <- factor(
    cut(df$importance_a, breaks = br, include.lowest = TRUE, labels = labs),
    levels = labs
  )
  out
}

plot_xy_dots <- function(
  df, x_var, y_var, xlab, title, signed = TRUE, ylim = NULL,
  subtitle = "One dot per taxon", label_thresh = NULL
) {
  y_vals <- df[[y_var]]
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]])) +
    ggplot2::geom_point(
      colour = POINT_COLOUR,
      size = 1.1,
      alpha = 0.7
    )
  if (signed) {
    lim <- if (is.null(ylim)) sym_y_limits(y_vals) else ylim
    p <- p + ggplot2::scale_y_continuous(limits = c(-lim, lim))
  } else {
    lim <- if (is.null(ylim)) pos_y_limits(y_vals) else ylim
    p <- p + ggplot2::scale_y_continuous(limits = c(0, lim))
  }
  p <- add_abs_labels(p, df, x_var, y_var, label_thresh)
  p <- p +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = xlab,
      y = y_label(signed)
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 10.5),
      plot.subtitle = ggplot2::element_text(size = 8),
      panel.grid.minor = ggplot2::element_blank()
    )
  if (signed) {
    p <- p + ggplot2::geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.35)
  }
  p
}

plot_rank_dots <- function(
  df, y_var, title, signed = TRUE, ylim = NULL, label_thresh = NULL
) {
  plot_xy_dots(
    df, "rank", y_var,
    xlab = "Importance rank (left transform)",
    title = title,
    signed = signed,
    ylim = ylim,
    subtitle = "One dot per rank",
    label_thresh = label_thresh
  ) + ggplot2::scale_x_continuous(breaks = pretty(df$rank, n = 12))
}

plot_importance_dots <- function(
  df, y_var, title, signed = TRUE, ylim = NULL, label_thresh = NULL
) {
  plot_xy_dots(
    df, "importance_a", y_var,
    xlab = "Normalized importance (left transform)",
    title = title,
    signed = signed,
    ylim = ylim,
    subtitle = "One dot per taxon",
    label_thresh = label_thresh
  )
}

plot_xy_bins <- function(df, grp_col, y_var, xlab, title, signed = TRUE) {
  y_vals <- df[[y_var]]
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[grp_col]], y = .data[[y_var]])) +
    ggplot2::geom_boxplot(
      fill = BOX_FILL,
      colour = "grey40",
      outlier.shape = NA,
      alpha = 0.85,
      width = 0.55,
      linewidth = 0.35
    ) +
    ggplot2::geom_point(
      colour = POINT_COLOUR,
      size = 0.9,
      alpha = 0.45,
      position = ggplot2::position_jitter(width = 0.12, height = 0)
    ) +
    y_scale(y_vals, signed = signed) +
    ggplot2::labs(
      title = title,
      subtitle = "Boxplot with jittered taxa",
      x = xlab,
      y = y_label(signed)
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 10.5),
      plot.subtitle = ggplot2::element_text(size = 8),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.5)
    )
  if (signed) {
    p <- p + ggplot2::geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.35)
  }
  p
}

plot_rank_bins <- function(df, grp_col, y_var, title, signed = TRUE) {
  plot_xy_bins(
    df, grp_col, y_var,
    xlab = "Importance rank group (left transform)",
    title = title,
    signed = signed
  )
}

plot_importance_bins <- function(df, grp_col, y_var, title, signed = TRUE) {
  plot_xy_bins(
    df, grp_col, y_var,
    xlab = "Importance bin (left transform, quantiles)",
    title = title,
    signed = signed
  )
}

build_pair_plot <- function(marker, pair_name, df, by = c("rank", "importance")) {
  by <- match.arg(by)
  n_taxa <- nrow(df)

  if (by == "rank") {
    df <- add_rank_groups(df, widths = c(5L, 10L))
    p1s <- plot_rank_dots(df, "diff", glue("Per rank | signed | n = {n_taxa}"), TRUE)
    p1a <- plot_rank_dots(df, "abs_diff", glue("Per rank | absolute | n = {n_taxa}"), FALSE)
    p5s <- plot_rank_bins(df, "rank_grp_5", "diff", "Rank groups of 5 | signed", TRUE)
    p5a <- plot_rank_bins(df, "rank_grp_5", "abs_diff", "Rank groups of 5 | absolute", FALSE)
    p10s <- plot_rank_bins(df, "rank_grp_10", "diff", "Rank groups of 10 | signed", TRUE)
    p10a <- plot_rank_bins(df, "rank_grp_10", "abs_diff", "Rank groups of 10 | absolute", FALSE)
    main_title <- glue("{marker} | {PAIR_LABELS[[pair_name]]} | importance difference by rank")
  } else {
    df <- add_imp_groups(df, 5L)
    df <- add_imp_groups(df, 10L)
    p1s <- plot_importance_dots(df, "diff", glue("By importance | signed | n = {n_taxa}"), TRUE)
    p1a <- plot_importance_dots(df, "abs_diff", glue("By importance | absolute | n = {n_taxa}"), FALSE)
    p5s <- plot_importance_bins(df, "imp_grp_5", "diff", "Importance quintiles | signed", TRUE)
    p5a <- plot_importance_bins(df, "imp_grp_5", "abs_diff", "Importance quintiles | absolute", FALSE)
    p10s <- plot_importance_bins(df, "imp_grp_10", "diff", "Importance deciles | signed", TRUE)
    p10a <- plot_importance_bins(df, "imp_grp_10", "abs_diff", "Importance deciles | absolute", FALSE)
    main_title <- glue("{marker} | {PAIR_LABELS[[pair_name]]} | importance difference by score")
  }

  patchwork::wrap_plots(
    p1s, p1a,
    p5s, p5a,
    p10s, p10a,
    ncol = 2,
    nrow = 3,
    heights = c(1, 1.1, 1.1)
  ) +
    patchwork::plot_annotation(
      title = main_title,
      subtitle = glue("Taxa with importance > 0 in both transforms; {norm_note()}; rank by left transform"),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 12),
        plot.subtitle = ggplot2::element_text(size = 8.5)
      )
    )
}

build_abs_3panel_plot <- function(marker, pair_dfs, by = c("rank", "importance")) {
  by <- match.arg(by)
  ylim <- pos_y_limits(unlist(lapply(pair_dfs, function(d) d$abs_diff)))

  panels <- lapply(names(PAIRS), function(pair_name) {
    df <- pair_dfs[[pair_name]]
    if (by == "rank") {
      plot_rank_dots(
        df, "abs_diff", PAIR_LABELS[[pair_name]], FALSE, ylim,
        label_thresh = ABS_LABEL_THRESH
      )
    } else {
      plot_importance_dots(
        df, "abs_diff", PAIR_LABELS[[pair_name]], FALSE, ylim,
        label_thresh = ABS_LABEL_THRESH
      )
    }
  })

  label_note <- glue("Labels: genus where |diff| ≥ {ABS_LABEL_THRESH}")
  main_title <- if (by == "rank") {
    glue("{marker} | absolute importance difference by rank | all pc pairs")
  } else {
    glue("{marker} | absolute importance difference by score | all pc pairs")
  }
  sub <- if (by == "rank") {
    glue("One panel per pc pair; rank by left transform; shared y-axis; {norm_note()}; {label_note}")
  } else {
    glue("One panel per pc pair; x = normalized left-transform importance; shared y-axis; {norm_note()}; {label_note}")
  }

  patchwork::wrap_plots(panels, ncol = 3) +
    patchwork::plot_annotation(
      title = main_title,
      subtitle = sub,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 12),
        plot.subtitle = ggplot2::element_text(size = 8.5)
      )
    )
}

for (marker in c("16S", "18S")) {
  wide <- build_wide(marker)
  pair_dfs <- list()

  for (pair_name in names(PAIRS)) {
    df <- pair_diff_ranked(wide, pair_name)
    pair_dfs[[pair_name]] <- df
    marker_dir <- marker_out_dir(marker)

    csv_out <- file.path(
      marker_dir,
      sprintf("rf_imp_%s_%s_diff_by_rank.csv", marker, pair_name)
    )
    utils::write.csv(df, csv_out, row.names = FALSE)

    p <- build_pair_plot(marker, pair_name, df, by = "rank")
    pdf_out <- file.path(
      marker_dir,
      sprintf("rf_imp_%s_%s_diff_by_rank.pdf", marker, pair_name)
    )
    ggplot2::ggsave(pdf_out, p, width = 14, height = 12, useDingbats = FALSE)
    message("Wrote ", pdf_out, " (", PAIR_LABELS[[pair_name]], ", n = ", nrow(df), ")")
    message("Wrote ", csv_out)

    p_imp <- build_pair_plot(marker, pair_name, df, by = "importance")
    pdf_imp <- file.path(
      marker_dir,
      sprintf("rf_imp_%s_%s_diff_by_importance.pdf", marker, pair_name)
    )
    ggplot2::ggsave(pdf_imp, p_imp, width = 14, height = 12, useDingbats = FALSE)
    message("Wrote ", pdf_imp)
  }

  p3 <- build_abs_3panel_plot(marker, pair_dfs, by = "rank")
  pdf3 <- file.path(marker_out_dir(marker), sprintf("rf_imp_%s_abs_diff_3panel.pdf", marker))
  ggplot2::ggsave(pdf3, p3, width = 16, height = 5.5, useDingbats = FALSE)
  message("Wrote ", pdf3)

  p3i <- build_abs_3panel_plot(marker, pair_dfs, by = "importance")
  pdf3i <- file.path(
    marker_out_dir(marker),
    sprintf("rf_imp_%s_abs_diff_3panel_by_importance.pdf", marker)
  )
  ggplot2::ggsave(pdf3i, p3i, width = 16, height = 5.5, useDingbats = FALSE)
  message("Wrote ", pdf3i)
}

message("Done.")
