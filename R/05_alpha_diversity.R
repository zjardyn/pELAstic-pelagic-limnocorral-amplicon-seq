source("R/01_load_files2.R")
source("R/functions.R")
source("R/alpha_diversity_functions.R")

theme_set(theme_bw(base_size = 14))

min_sum_taxa <- 5000L
rarefy_depth <- min_sum_taxa
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

alpha_plot_panel <- function(alpha, marker, metric_col, panel_tag, y_label, row_title, xlab = NULL) {
  y_range <- range(alpha[[metric_col]], na.rm = TRUE)

  p_ws <- alpha %>%
    filter(Location == "WS") %>%
    ggplot(aes(x = sample, y = .data[[metric_col]])) +
    geom_point() +
    facet_grid(
      Location ~ Date, scales = "free_x",
      labeller = labeller(
        Location = function(x) "Wall strip",
        Date = function(x) paste("Week", x)
      )
    ) +
    theme(
      axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
      plot.tag = element_text(face = "bold")
    ) +
    ylim(y_range) +
    labs(tag = panel_tag, x = xlab, y = y_label)

  p_ms <- alpha %>%
    filter(Location == "MS") %>%
    ggplot(aes(x = sample, y = .data[[metric_col]])) +
    geom_point() +
    facet_grid(
      Location ~ Date,
      labeller = labeller(
        Location = function(x) "Microscope slide",
        Date = function(x) paste("Week", x)
      )
    ) +
    theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
    ylim(y_range) +
    labs(x = xlab, y = y_label)

  p_ws + p_ms + p_ms + p_ms +
    plot_layout(design = "AAAB", guides = "collect", axes = "collect") +
    plot_annotation(title = sprintf("%s: %s", marker, row_title))
}

make_alpha_figure <- function(phy, marker, out_file) {
  alpha <- prepare_alpha_data(phy, marker, min_sum_taxa, rarefy_depth) %>%
    arrange_alpha_for_plots()

  p1 <- alpha_plot_panel(
    alpha, marker, "Observed", "A", "Observed ASVs",
    sprintf("Observed ASV richness (rarefied to %d reads)", rarefy_depth)
  )
  p2 <- alpha_plot_panel(
    alpha, marker, "Berger_Parker", "B", "Berger-Parker",
    sprintf("Berger-Parker dominance (rarefied to %d reads)", rarefy_depth)
  )
  p3 <- alpha_plot_panel(
    alpha, marker, "Shannon", "C", "Shannon",
    sprintf("Shannon diversity (rarefied to %d reads)", rarefy_depth),
    xlab = "Sample name"
  )

  p_comb <- wrap_plots(p1, p2, p3, ncol = 1, axes = "collect")
  ggsave(out_file, p_comb, width = 10, height = 12, dpi = 300, scale = 0.9)
  message("Wrote ", out_file)
  invisible(p_comb)
}

message("Alpha diversity exploratory plots (ASV-level, rarefied) | depth = ", rarefy_depth)
make_alpha_figure(phy_16S, "16S", "figures/fig_S3_alpha_diversity_16S.png")
make_alpha_figure(phy_18S, "18S", "figures/fig_S4_alpha_diversity_18S.png")
