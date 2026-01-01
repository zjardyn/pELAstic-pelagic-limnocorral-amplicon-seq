source("R/functions.R")

# Enable markdown rendering globally
theme_set(theme_minimal())

# Read in the stability data
stability_data_16S <- readRDS("output/ancombc_stability_16S_9_ws_trend.rds")
stability_data_18S <- readRDS("output/ancombc_stability_18S_9_ws_trend.rds")

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

# Calculate bootstrap-based statistics
bootstrap_stats_16S <- lfc_long_16S %>%
  group_by(taxon, plastic_level) %>%
  summarise(
    n_runs = n(),
    lfc = mean(lfc, na.rm = TRUE),
    passed_ss = all(passed_ss, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_runs), plastic_level) %>% 
  filter(n_runs > 50)

bootstrap_stats_18S <- lfc_long_18S %>%
  group_by(taxon, plastic_level) %>%
  summarise(
    n_runs = n(),
    lfc = mean(lfc, na.rm = TRUE),
    passed_ss = all(passed_ss, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_runs), plastic_level) %>%
  filter(n_runs > 50)

# Create heatmap data - prepare for plotting
# Create heatmap data - prepare for plotting
# First create labels without asterisks for ordering
heatmap_data_16S <- bootstrap_stats_16S %>%
  mutate(
    plastic_level = factor(plastic_level, levels = c("low", "medium", "high"), 
                          labels = c("Low", "Medium", "High")),
    taxon_display = paste0(taxon, " (", n_runs, ")")
  )

heatmap_data_18S <- bootstrap_stats_18S %>%
  mutate(
    plastic_level = factor(plastic_level, levels = c("low", "medium", "high"), 
                          labels = c("Low", "Medium", "High")),
    taxon_display = paste0(taxon, " (", n_runs, ")")
  )

# Order taxa by summed LFC (descending) for better visualization
order_taxa_16S <- heatmap_data_16S %>%
  group_by(taxon_display) %>%
  summarise(sum_lfc = sum(lfc), .groups = "drop") %>%
  arrange(desc(sum_lfc)) %>%
  pull(taxon_display)

order_taxa_18S <- heatmap_data_18S %>%
  group_by(taxon_display) %>%
  summarise(sum_lfc = sum(lfc), .groups = "drop") %>%
  arrange(desc(sum_lfc)) %>%
  pull(taxon_display)

# Convert to factors for ordering
heatmap_data_16S <- heatmap_data_16S %>% 
  mutate(taxon_display = factor(taxon_display, levels = order_taxa_16S))

heatmap_data_18S <- heatmap_data_18S %>% 
  mutate(taxon_display = factor(taxon_display, levels = order_taxa_18S))

# Create mappings for colors (based on passed_ss)
taxon_display_passed_ss_16S <- heatmap_data_16S %>%
  select(taxon_display, passed_ss) %>%
  distinct() %>%
  deframe()

taxon_display_passed_ss_18S <- heatmap_data_18S %>%
  select(taxon_display, passed_ss) %>%
  distinct() %>%
  deframe()

# Create color scale (yellow to white to purple)
max_abs_lfc <- max(c(abs(heatmap_data_16S$lfc), abs(heatmap_data_18S$lfc)), na.rm = TRUE)
if (!is.finite(max_abs_lfc) || max_abs_lfc == 0) max_abs_lfc <- 1

shared_scale <- scale_fill_gradient2(
  low = "#440154",  # Purple
  mid = "#F7F7F7",  # White
  high = "#FDE725", # Yellow
  midpoint = 0,
  limits = c(-max_abs_lfc, max_abs_lfc),
  oob = scales::squish,
  name = "Log Fold Change",
  guide = guide_colorbar(title.position = "top", title.hjust = 0.5)
)

# Create the heatmap
p_heatmap_16S <- ggplot(heatmap_data_16S, aes(x = plastic_level, y = taxon_display, fill = lfc)) +
  geom_tile(color = "white", linewidth = 0.5) +
  shared_scale +
  labs(
    x = "Plastic Level", 
    y = "Bacterial & Archeal Taxa",
    tag = "A"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(
      color = ifelse(taxon_display_passed_ss_16S[levels(heatmap_data_16S$taxon_display)], "#000000", "#8a8a8a"),
      size = 13
    )
  )

p_heatmap_18S <- ggplot(heatmap_data_18S, aes(x = plastic_level, y = taxon_display, fill = lfc)) +
  geom_tile(color = "white", linewidth = 0.5) +
  shared_scale +
  labs(
    x = "Plastic Level", 
    y = "Eukaryotic Taxa",
    tag = "B"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(
      color = ifelse(taxon_display_passed_ss_18S[levels(heatmap_data_18S$taxon_display)], "#000000", "#8a8a8a"),
      size = 13
    )
  )

# Create layout for combined plots
layout <- "
  A
  A
  A
  A 
  A 
  A
  B "

# Create the combined plot with layout
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

# save as text-editable pdf for inkscape editing, can't italcize the taxa 
pdf(
  "figures/fig_5_trend_analysis_heatmap_stable_taxa_combined.pdf",
  width = 10 * 0.8,
  height = 12 * 0.8,
  family = "Helvetica",
  useDingbats = FALSE
)

print(p_combined)
dev.off()

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



