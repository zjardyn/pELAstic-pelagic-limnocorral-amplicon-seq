# CLR pseudocount importance panels (clr_1, clr_0.5, clr_0.1) side by side with
# genus-coloured ASV trajectories (top-N by mean normalized importance) and rank-mean line.
# Importance normalized per panel: value / sum(importance > 0) within each transform.
#
# Writes: figures/rf_loocv_importance_clr_pc/tracks/{16S,18S}/
#   rf_imp_{marker}_clr_pc_tracks_top{n}.pdf  for n in 10, 20, 30, 40, 50
#   (skips n > available exact3of3 ASVs)
#
# Optional env:
#   RF_IMP_FIG_DIR  — output figure root (default figures/rf_loocv_importance_clr_pc)
#   RF_IMP_ASVS     — exact3of3 (default) | union | panel_gt0
#   RF_IMP_TOP_N    — comma-separated cutoffs (default 10,20,30,40,50)
#
# Run: Rscript R/25_rf_clr_pc_importance_tracks.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(glue)
})

TRANSFORMS <- c("clr_1", "clr_0.5", "clr_0.1")
PANEL_TITLES <- c(
  "clr_1" = "CLR pc = 1",
  "clr_0.5" = "CLR pc = 0.5",
  "clr_0.1" = "CLR pc = 0.1"
)
PANEL_IDX <- stats::setNames(seq_along(TRANSFORMS), TRANSFORMS)
PANEL_WIDTH <- 0.9
TOP_N_DEFAULT <- c(10L, 20L, 30L, 40L, 50L)

rf_dir <- Sys.getenv("RF_OUT_DIR", unset = "output")
fig_dir <- Sys.getenv("RF_IMP_FIG_DIR", unset = "figures/rf_loocv_importance_clr_pc")
tracks_dir <- file.path(fig_dir, "tracks")
asv_mode <- Sys.getenv("RF_IMP_ASVS", unset = "exact3of3")
top_n_vals <- as.integer(strsplit(
  Sys.getenv("RF_IMP_TOP_N", unset = paste(TOP_N_DEFAULT, collapse = ",")),
  ",",
  fixed = TRUE
)[[1]])
top_n_vals <- top_n_vals[is.finite(top_n_vals) & top_n_vals > 0L]
if (!length(top_n_vals)) stop("RF_IMP_TOP_N must list positive integers")

for (mk in c("16S", "18S")) {
  dir.create(file.path(tracks_dir, mk), showWarnings = FALSE, recursive = TRUE)
}

rf_path <- function(marker, transform) {
  file.path(rf_dir, sprintf("rf_%s_9_ws_n9_prev2of9_%s.rds", marker, transform))
}

load_imp_gt0 <- function(marker, transform) {
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
      median_loo_importance = dplyr::coalesce(
        median_loo_importance, full_n_importance, mean_importance
      )
    ) %>%
    dplyr::filter(is.finite(median_loo_importance), median_loo_importance > 0) %>%
    dplyr::mutate(
      median_loo_importance = median_loo_importance / sum(median_loo_importance, na.rm = TRUE)
    ) %>%
    dplyr::arrange(dplyr::desc(median_loo_importance)) %>%
    dplyr::mutate(rank = dplyr::row_number())
}

OKABE_ITO_8 <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#999999"
)
GOLDEN_RATIO <- 0.618033988749895

extended_distinct_palette <- function(n) {
  if (n <= 0L) return(character())
  if (n <= length(OKABE_ITO_8)) return(OKABE_ITO_8[seq_len(n)])
  cols <- character(n)
  cols[seq_len(length(OKABE_ITO_8))] <- OKABE_ITO_8
  i <- seq_len(n - length(OKABE_ITO_8)) + length(OKABE_ITO_8)
  h <- (i * GOLDEN_RATIO) %% 1 * 360
  cols[i] <- grDevices::hcl(
    h = h, c = 62 + (i %% 4L) * 7, l = 48 + (i %% 5L) * 6
  )
  cols
}

genus_colour_map <- function(genus) {
  labs <- sort(unique(as.character(genus)))
  labs <- labs[!is.na(labs) & nzchar(labs)]
  stats::setNames(extended_distinct_palette(length(labs)), labs)
}

genus_label <- function(genus) {
  g <- as.character(genus)
  g[is.na(g) | !nzchar(g) | g %in% c("NA", "unclassified", "Unclassified")] <- "(unassigned)"
  g
}

panel_x <- function(panel_i, rank, n_keep) {
  panel_i + (rank - 1) / max(n_keep, 1L) * PANEL_WIDTH
}

select_asvs <- function(gt0_list) {
  sets <- lapply(gt0_list, function(d) d$asv)
  if (asv_mode == "union") {
    return(sort(unique(unlist(sets))))
  }
  if (asv_mode == "panel_gt0") {
    return(NULL)
  }
  if (asv_mode == "exact3of3") {
    all_asvs <- unique(unlist(sets))
    n_tr <- vapply(all_asvs, function(a) {
      sum(vapply(sets, function(s) a %in% s, logical(1)))
    }, numeric(1))
    return(sort(unname(all_asvs[n_tr == length(TRANSFORMS)])))
  }
  stop("RF_IMP_ASVS must be exact3of3, union, or panel_gt0, got: ", asv_mode)
}

rank_top_asvs <- function(track_df, keep_asvs, n_top) {
  med_imp <- track_df %>%
    dplyr::group_by(asv) %>%
    dplyr::summarise(
      median_loo_importance = median(median_loo_importance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(median_loo_importance))
  if (!is.null(keep_asvs)) {
    med_imp <- dplyr::filter(med_imp, asv %in% keep_asvs)
  }
  head(med_imp$asv, n_top)
}

build_track_plot <- function(marker, n_top) {
  gt0_list <- lapply(TRANSFORMS, function(tr) load_imp_gt0(marker, tr))
  names(gt0_list) <- TRANSFORMS
  n_keep <- vapply(gt0_list, nrow, integer(1))
  keep_asvs <- select_asvs(gt0_list)

  env_rows <- list()
  track_rows <- list()
  for (tr in TRANSFORMS) {
    d <- gt0_list[[tr]]
    pi <- PANEL_IDX[[tr]]
    nk <- n_keep[[tr]]
    d <- d %>%
      dplyr::mutate(
        panel = pi,
        x = panel_x(pi, rank, nk)
      )
    env_rows[[tr]] <- d
    trk <- d
    if (!is.null(keep_asvs)) {
      trk <- dplyr::filter(trk, asv %in% keep_asvs)
    }
    track_rows[[tr]] <- trk
  }
  env_df <- dplyr::bind_rows(env_rows)
  track_df <- dplyr::bind_rows(track_rows)

  top_asvs <- rank_top_asvs(track_df, keep_asvs, n_top)
  track_df <- track_df %>%
    dplyr::filter(asv %in% top_asvs) %>%
    dplyr::mutate(GenusLabel = genus_label(Genus))

  colour_map <- genus_colour_map(track_df$GenusLabel)
  track_df$GenusLabel <- factor(track_df$GenusLabel, levels = names(colour_map))

  vline_x <- mapply(function(pi, nk) {
    panel_x(pi, nk, nk) + 0.04
  }, pi = PANEL_IDX[TRANSFORMS], nk = n_keep[TRANSFORMS])

  panel_centres <- PANEL_IDX[TRANSFORMS] + PANEL_WIDTH / 2

  ggplot2::ggplot() +
    ggplot2::geom_line(
      data = env_df,
      ggplot2::aes(x = x, y = median_loo_importance, group = transform),
      colour = "grey25",
      linewidth = 0.7,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_line(
      data = track_df,
      ggplot2::aes(
        x = x,
        y = median_loo_importance,
        colour = GenusLabel,
        group = asv
      ),
      linewidth = 0.35,
      alpha = 0.55
    ) +
    ggplot2::geom_point(
      data = track_df,
      ggplot2::aes(x = x, y = mean_importance, colour = GenusLabel),
      size = 0.9,
      alpha = 0.75
    ) +
    ggplot2::geom_vline(
      xintercept = c(PANEL_IDX[TRANSFORMS[-1]] - 0.05),
      colour = "grey55",
      linewidth = 0.35
    ) +
    ggplot2::geom_vline(
      xintercept = vline_x,
      colour = "red",
      linetype = "dotted",
      linewidth = 0.45
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      colour = "red",
      linetype = "dotted",
      linewidth = 0.45
    ) +
    ggplot2::scale_x_continuous(
      breaks = panel_centres,
      labels = unname(PANEL_TITLES[TRANSFORMS]),
      limits = c(0.85, max(PANEL_IDX) + PANEL_WIDTH + 0.15),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::scale_colour_manual(
      values = colour_map,
      name = "Genus",
      guide = ggplot2::guide_legend(
        ncol = 1L,
        byrow = FALSE,
        title.position = "top",
        override.aes = list(alpha = 1, linewidth = 0.6, size = 1.2)
      )
    ) +
    ggplot2::labs(
      title = glue(
        "{marker} | CLR pseudocount importance tracks | top {length(top_asvs)} ASVs | {asv_mode}"
      ),
      subtitle = paste(
        "Coloured lines = top ASVs by mean normalized importance across CLR pc;",
        "grey = rank-mean (÷ Σ gt0 per panel); red dotted = last rank > 0 per panel",
        sep = " "
      ),
      x = NULL,
      y = "Normalized permutation importance (÷ Σ gt0 per panel)"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 8.5),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 6),
      legend.text = ggplot2::element_text(size = 4.5),
      legend.key.size = ggplot2::unit(0.22, "cm"),
      legend.key.height = ggplot2::unit(0.22, "cm"),
      legend.key.width = ggplot2::unit(0.35, "cm"),
      legend.spacing.y = ggplot2::unit(0.05, "cm"),
      legend.box.spacing = ggplot2::unit(0.1, "cm"),
      plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5)
    )
}

for (marker in c("16S", "18S")) {
  gt0_list <- lapply(TRANSFORMS, function(tr) load_imp_gt0(marker, tr))
  keep_asvs <- select_asvs(gt0_list)
  n_available <- length(keep_asvs)

  for (n_top in top_n_vals) {
    if (n_top > n_available) {
      message(marker, ": skip top ", n_top, " (only ", n_available, " ASVs)")
      next
    }
    p <- build_track_plot(marker, n_top)
    out <- file.path(
      tracks_dir, marker,
      sprintf("rf_imp_%s_clr_pc_tracks_top%d.pdf", marker, n_top)
    )
    ggplot2::ggsave(out, p, width = 15, height = 9, useDingbats = FALSE)
    message("Wrote ", out)
  }

  if (!(n_available %in% top_n_vals) && n_available > min(top_n_vals)) {
    p <- build_track_plot(marker, n_available)
    out <- file.path(
      tracks_dir, marker,
      sprintf("rf_imp_%s_clr_pc_tracks_top%d.pdf", marker, n_available)
    )
    ggplot2::ggsave(out, p, width = 15, height = 9, useDingbats = FALSE)
    message("Wrote ", out, " (all ", n_available, " ASVs)")
  }
}

message("Done.")
