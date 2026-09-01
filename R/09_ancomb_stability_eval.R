source("R/functions.R")
source("R/01_load_files2.R")

# Enable markdown rendering globally
theme_set(theme_minimal())

# Week-9 WS prev2of9 trend stability (100 seeds) from R/12
stability_data_16S <- readRDS("output/ancombc_trend_stability_16S_9_ws_prev2of9.rds")$results
stability_data_18S <- readRDS("output/ancombc_trend_stability_18S_9_ws_prev2of9.rds")$results
message("16S trend runs: ", length(stability_data_16S))
message("18S trend runs: ", length(stability_data_18S))

genus_map_16S <- genus_label_map_from_phy(phy_16S)
genus_map_18S <- genus_label_map_from_phy(phy_18S)
rf_gen_16S <- rf_all4_genus_labels("16S")
rf_gen_18S <- rf_all4_genus_labels("18S")
message("RF 4/4 genera for * mark: 16S=", length(rf_gen_16S),
        " 18S=", length(rf_gen_18S))

# Function to extract significant taxa from a single run
extract_significant_taxa <- function(run_result) {
  # Extract trend results
  trend_results <- run_result$res_trend
  
  # Filter for significant taxa (p_val < 0.05 and diff_abn == TRUE)
  significant_taxa <- trend_results %>%
    filter(p_val < 0.05 & diff_abn == TRUE) %>%
    select(taxon, 
           lfc_plastic_levellow, lfc_plastic_levelmedium, lfc_plastic_levelhigh,
           se_plastic_levellow, se_plastic_levelmedium, se_plastic_levelhigh,
           W, p_val, q_val, diff_abn, passed_ss, diff_robust_abn)
  
  return(significant_taxa)
}

# Extract significant taxa from all runs
significant_taxa_list_16S <- map(stability_data_16S, extract_significant_taxa)
significant_taxa_list_18S <- map(stability_data_18S, extract_significant_taxa)

# Count how many runs each taxon appears as significant
all_significant_taxa_16S <- map(significant_taxa_list_16S, ~ .x$taxon) %>%
  unlist() %>%
  table() %>%
  as.data.frame() %>%
  arrange(desc(Freq)) %>%
  rename(taxon = ".", n_runs = "Freq")

# Count how many runs each taxon appears as significant
all_significant_taxa_18S <- map(significant_taxa_list_18S, ~ .x$taxon) %>%
  unlist() %>%
  table() %>%
  as.data.frame() %>%
  arrange(desc(Freq)) %>%
  rename(taxon = ".", n_runs = "Freq")

# Add run ID to each dataset and combine all runs
significant_taxa_with_runs_16S <- map2(significant_taxa_list_16S, 1:length(significant_taxa_list_16S), 
                                  ~ .x %>% mutate(run_id = .y))
all_significant_taxa_data_16S <- bind_rows(significant_taxa_with_runs_16S)

significant_taxa_with_runs_18S <- map2(significant_taxa_list_18S, 1:length(significant_taxa_list_18S), 
                                  ~ .x %>% mutate(run_id = .y))
all_significant_taxa_data_18S <- bind_rows(significant_taxa_with_runs_18S)

# Pivot LFC data to long format
lfc_long_16S <- all_significant_taxa_data_16S %>%
  select(taxon, run_id, lfc_plastic_levellow, lfc_plastic_levelmedium, lfc_plastic_levelhigh, passed_ss) %>%
  pivot_longer(
    cols = c(lfc_plastic_levellow, lfc_plastic_levelmedium, lfc_plastic_levelhigh),
    names_to = "plastic_level",
    values_to = "lfc"
  ) %>%
  mutate(
    plastic_level = case_when(
      plastic_level == "lfc_plastic_levellow" ~ "low",
      plastic_level == "lfc_plastic_levelmedium" ~ "medium", 
      plastic_level == "lfc_plastic_levelhigh" ~ "high"
    )
  )

lfc_long_18S <- all_significant_taxa_data_18S %>%
  select(taxon, run_id, lfc_plastic_levellow, lfc_plastic_levelmedium, lfc_plastic_levelhigh, passed_ss) %>%
  pivot_longer(
    cols = c(lfc_plastic_levellow, lfc_plastic_levelmedium, lfc_plastic_levelhigh),
    names_to = "plastic_level",
    values_to = "lfc"
  ) %>%
  mutate(
    plastic_level = case_when(
      plastic_level == "lfc_plastic_levellow" ~ "low",
      plastic_level == "lfc_plastic_levelmedium" ~ "medium", 
      plastic_level == "lfc_plastic_levelhigh" ~ "high"
    )
  )

# Per-(taxon, level) bootstrap summaries (before stability cut)
bootstrap_stats_16S_all <- lfc_long_16S %>%
  group_by(taxon, plastic_level) %>%
  summarise(
    n_runs = n(),
    lfc = mean(lfc, na.rm = TRUE),
    passed_ss = all(passed_ss, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_runs), plastic_level)

bootstrap_stats_18S_all <- lfc_long_18S %>%
  group_by(taxon, plastic_level) %>%
  summarise(
    n_runs = n(),
    lfc = mean(lfc, na.rm = TRUE),
    passed_ss = all(passed_ss, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_runs), plastic_level)

make_trend_heatmap_combined <- function(stats_16S, stats_18S, min_n_runs, outfile) {
  bootstrap_stats_16S <- stats_16S %>% filter(n_runs >= min_n_runs)
  bootstrap_stats_18S <- stats_18S %>% filter(n_runs >= min_n_runs)

  message(
    if (min_n_runs <= 0L) {
      "Stability cut: all taxa ever trend-significant"
    } else {
      paste0("Stability cut n_runs >= ", min_n_runs, " (", min_n_runs, "%)")
    },
    ": 16S taxa=", n_distinct(bootstrap_stats_16S$taxon),
    ", 18S taxa=", n_distinct(bootstrap_stats_18S$taxon)
  )

  heatmap_data_16S <- bootstrap_stats_16S %>%
    mutate(
      plastic_level = factor(
        plastic_level,
        levels = c("low", "medium", "high"),
        labels = c("Low", "Medium", "High")
      ),
      taxon_lab = label_ancombc_taxa(taxon, genus_map_16S),
      taxon_display = paste0(
        star_rf_overlap_labels(taxon_lab, rf_gen_16S), " (", n_runs, ")"
      )
    )

  heatmap_data_18S <- bootstrap_stats_18S %>%
    mutate(
      plastic_level = factor(
        plastic_level,
        levels = c("low", "medium", "high"),
        labels = c("Low", "Medium", "High")
      ),
      taxon_lab = label_ancombc_taxa(taxon, genus_map_18S),
      taxon_display = paste0(
        star_rf_overlap_labels(taxon_lab, rf_gen_18S), " (", n_runs, ")"
      )
    )
  message(
    "RF* on trend heatmap: 16S=",
    n_distinct(heatmap_data_16S$taxon_lab[
      heatmap_data_16S$taxon_lab %in% rf_gen_16S
    ]),
    " 18S=",
    n_distinct(heatmap_data_18S$taxon_lab[
      heatmap_data_18S$taxon_lab %in% rf_gen_18S
    ])
  )

  make_one_panel <- function(heatmap_data, y_lab, tag) {
    if (nrow(heatmap_data) == 0L) {
      return(
        ggplot() +
          annotate(
            "text", x = 0.5, y = 0.5,
            label = paste0("No taxa with n_runs >= ", min_n_runs),
            size = 4.5
          ) +
          labs(x = "Plastic level", y = y_lab, tag = tag) +
          theme_void() +
          theme(
            axis.title = element_text(size = 16),
            plot.tag = element_text(size = 28, face = "bold"),
            plot.margin = margin(20, 20, 20, 20)
          )
      )
    }

    order_taxa <- order_taxa_plastic_lfc_columns(
      heatmap_data,
      id_col = "taxon_display",
      level_col = "plastic_level",
      value_col = "lfc",
      col_order = c("Low", "Medium", "High"),
      tol = 0.15
    )

    heatmap_data <- heatmap_data %>%
      mutate(taxon_display = factor(taxon_display, levels = order_taxa))

    taxon_display_passed_ss <- heatmap_data %>%
      select(taxon_display, passed_ss) %>%
      distinct() %>%
      deframe()

    ggplot(heatmap_data, aes(x = plastic_level, y = taxon_display, fill = lfc)) +
      geom_tile(color = "white", linewidth = 0.5) +
      labs(x = "Plastic level", y = y_lab, tag = tag) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(
          color = ifelse(
            taxon_display_passed_ss[levels(heatmap_data$taxon_display)],
            "#000000",
            "#8a8a8a"
          ),
          size = 13
        )
      )
  }

  lfc_vals <- c(heatmap_data_16S$lfc, heatmap_data_18S$lfc)
  max_abs_lfc <- if (length(lfc_vals)) max(abs(lfc_vals), na.rm = TRUE) else 1
  if (!is.finite(max_abs_lfc) || max_abs_lfc == 0) max_abs_lfc <- 1

  shared_scale <- scale_fill_gradient2(
    low = "#440154",
    mid = "#F7F7F7",
    high = "#FDE725",
    midpoint = 0,
    limits = c(-max_abs_lfc, max_abs_lfc),
    oob = scales::squish,
    name = "Log fold change",
    guide = guide_colorbar(title.position = "top", title.hjust = 0.5)
  )

  p_heatmap_16S <- make_one_panel(heatmap_data_16S, "Prokaryotic taxa", "A")
  p_heatmap_18S <- make_one_panel(heatmap_data_18S, "Eukaryotic taxa", "B")

  if (nrow(heatmap_data_16S) > 0L) p_heatmap_16S <- p_heatmap_16S + shared_scale
  if (nrow(heatmap_data_18S) > 0L) p_heatmap_18S <- p_heatmap_18S + shared_scale

  layout <- "
  A
  A
  A
  A
  A
  A
  B "

  p_combined <- (p_heatmap_16S / p_heatmap_18S) +
    patchwork::plot_layout(design = layout, guides = "collect", axes = "collect") &
    theme(
      axis.text.x = element_text(size = 13),
      axis.title.y = element_text(size = 16),
      axis.title.x = element_text(size = 16),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      legend.position = "bottom",
      legend.justification = c(0, 0),
      legend.text = element_text(size = 13),
      legend.title = element_text(size = 16),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      plot.tag = element_text(size = 28, face = "bold")
    )

  pdf(
    outfile,
    width = 10 * 0.8,
    height = max(12 * 0.8, 0.28 * (n_distinct(heatmap_data_16S$taxon) +
      n_distinct(heatmap_data_18S$taxon) + 8)),
    family = "Helvetica",
    useDingbats = FALSE
  )
  print(p_combined)
  dev.off()
  message("Wrote ", outfile)
}

# Stability cut: n_runs >= 20 (20% of 100 boots)
make_trend_heatmap_combined(
  bootstrap_stats_16S_all,
  bootstrap_stats_18S_all,
  min_n_runs = 20L,
  outfile = "figures/20pct_fig_5_trend_analysis_heatmap_stable_taxa_combined.pdf"
)

# Stability cut: n_runs >= 50 (50% of 100 boots; primary reproducibility bar)
make_trend_heatmap_combined(
  bootstrap_stats_16S_all,
  bootstrap_stats_18S_all,
  min_n_runs = 50L,
  outfile = "figures/50pct_fig_5_trend_analysis_heatmap_stable_taxa_combined.pdf"
)

# All taxa ever trend-significant in any bootstrap (no stability cut)
make_trend_heatmap_combined(
  bootstrap_stats_16S_all,
  bootstrap_stats_18S_all,
  min_n_runs = 0L,
  outfile = "figures/fig_5_trend_analysis_heatmap_all_taxa_combined.pdf"
)
### Transposed heatmaps - taxa on x-axis, plastic level on y-axis
# Reverse order for 16S transposed plot
# heatmap_data_16S_transposed <- heatmap_data_16S %>%
#   mutate(taxon_display = factor(taxon_display, levels = rev(order_taxa_16S)))

# p_heatmap_16S_transposed <- ggplot(heatmap_data_16S_transposed, aes(x = taxon_display, y = plastic_level, fill = lfc)) +
#   geom_tile(color = "white", linewidth = 0.5) +
#   shared_scale +
#   scale_x_discrete(labels = label_func_16S) +
#   labs(
#     x = "Bacterial & Archeal Taxa", 
#     y = "Plastic Level",
#     tag = "A"
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x = ggtext::element_markdown(
#       size = 10,
#       angle = 45,
#       hjust = 1,
#       vjust = 1
#     ),
#     axis.text.y = element_text(size = 12),
#     plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
#     plot.subtitle = element_text(hjust = 0.5, size = 12),
#     legend.position = "bottom",
#     panel.grid = element_blank(),
#     axis.ticks = element_blank(),
#     plot.tag = element_text(size = 16, face = "bold")
#   )

# p_heatmap_18S_transposed <- ggplot(heatmap_data_18S, aes(x = taxon_display, y = plastic_level, fill = lfc)) +
#   geom_tile(color = "white", linewidth = 0.5) +
#   shared_scale +
#   scale_x_discrete(labels = label_func_18S) +
#   labs(
#     x = "Eukaryotic Taxa", 
#     y = "Plastic Level",
#     tag = "B"
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x = ggtext::element_markdown(
#       size = 10,
#       angle = 45,
#       hjust = 1,
#       vjust = 1
#     ),
#     axis.text.y = element_text(size = 12),
#     plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
#     plot.subtitle = element_text(hjust = 0.5, size = 12),
#     legend.position = "bottom",
#     panel.grid = element_blank(),
#     axis.ticks = element_blank(),
#     plot.tag = element_text(size = 16, face = "bold")
#   )

# # Create transposed combined plot - side by side layout
# layout_transposed <- "AAAAAAB"

# p_combined_transposed <- (p_heatmap_16S_transposed + p_heatmap_18S_transposed) + 
#   patchwork::plot_layout(design = layout_transposed, guides = "collect", axes = "collect") 

# # Save the transposed combined plot
# ggsave("figures/trend_analysis_heatmap_stable_taxa_combined_transposed.png", p_combined_transposed, 
#        width = 14, height = 6, dpi = 300, scale = 0.8)



