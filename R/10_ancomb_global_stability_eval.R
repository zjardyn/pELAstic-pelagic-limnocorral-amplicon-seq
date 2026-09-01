source("R/functions.R")
source("R/01_load_files2.R")

# Global-test ANCOM-BC2 heatmap from the dedicated global run (R/12 prev2of9).
# Significance from res_global (raw p < 0.05); LFCs from primary pairwise res.
# Layout: prokaryotes split halfway into two side-by-side panels; eukaryotes below.

theme_set(theme_minimal())

global_16S <- readRDS("output/ancombc_global_16S_9_ws_prev2of9.rds")
global_18S <- readRDS("output/ancombc_global_18S_9_ws_prev2of9.rds")
run_16S <- global_16S$result
run_18S <- global_18S$result
message("Using 16S global seed: ", global_16S$seed)
message("Using 18S global seed: ", global_18S$seed)

genus_map_16S <- genus_label_map_from_phy(phy_16S)
genus_map_18S <- genus_label_map_from_phy(phy_18S)
rf_gen_16S <- rf_all4_genus_labels("16S")
rf_gen_18S <- rf_all4_genus_labels("18S")
message("RF 4/4 genera for * mark: 16S=", length(rf_gen_16S),
        " 18S=", length(rf_gen_18S))

extract_significant_taxa_global <- function(run_result) {
  run_result$res_global %>%
    filter(p_val < 0.05) %>%
    select(taxon, W, p_val, q_val, diff_abn, passed_ss, diff_robust_abn) %>%
    left_join(
      run_result$res %>%
        select(
          taxon,
          lfc_plastic_levellow, lfc_plastic_levelmedium, lfc_plastic_levelhigh
        ),
      by = "taxon"
    )
}

sig_16S <- extract_significant_taxa_global(run_16S)
sig_18S <- extract_significant_taxa_global(run_18S)

message("16S global p < 0.05 taxa: ", nrow(sig_16S))
message("18S global p < 0.05 taxa: ", nrow(sig_18S))

to_heatmap_long <- function(sig_data) {
  if (nrow(sig_data) == 0L) {
    return(tibble(
      taxon = character(),
      plastic_level = factor(levels = c("Low", "Medium", "High")),
      lfc = numeric(),
      passed_ss = logical(),
      taxon_display = character()
    ))
  }
  sig_data %>%
    select(
      taxon, passed_ss,
      lfc_plastic_levellow, lfc_plastic_levelmedium, lfc_plastic_levelhigh
    ) %>%
    pivot_longer(
      cols = c(lfc_plastic_levellow, lfc_plastic_levelmedium, lfc_plastic_levelhigh),
      names_to = "plastic_level",
      values_to = "lfc"
    ) %>%
    mutate(
      plastic_level = factor(
        plastic_level,
        levels = c(
          "lfc_plastic_levellow",
          "lfc_plastic_levelmedium",
          "lfc_plastic_levelhigh"
        ),
        labels = c("Low", "Medium", "High")
      ),
      taxon_display = taxon
    )
}

heatmap_data_16S <- to_heatmap_long(sig_16S) %>%
  mutate(
    taxon_lab = label_ancombc_taxa(taxon, genus_map_16S),
    taxon_display = star_rf_overlap_labels(taxon_lab, rf_gen_16S)
  )
heatmap_data_18S <- to_heatmap_long(sig_18S) %>%
  mutate(
    taxon_lab = label_ancombc_taxa(taxon, genus_map_18S),
    taxon_display = star_rf_overlap_labels(taxon_lab, rf_gen_18S)
  )
message(
  "RF* on global heatmap: 16S=",
  n_distinct(heatmap_data_16S$taxon_lab[
    heatmap_data_16S$taxon_lab %in% rf_gen_16S
  ]),
  " 18S=",
  n_distinct(heatmap_data_18S$taxon_lab[
    heatmap_data_18S$taxon_lab %in% rf_gen_18S
  ])
)

# Order full prokaryote set once, then split at midpoint
order_taxa_16S <- if (nrow(heatmap_data_16S)) {
  order_taxa_plastic_lfc_columns(
    heatmap_data_16S,
    id_col = "taxon_display",
    level_col = "plastic_level",
    value_col = "lfc",
    col_order = c("Low", "Medium", "High"),
    tol = 0.15
  )
} else {
  character()
}

# After ordering (bottom→top), put the upper half on the left and lower half on the right
# (swap relative to a simple first-half | second-half split).
n16 <- length(order_taxa_16S)
mid <- floor(n16 / 2)
taxa_16S_left <- if (n16 > mid) order_taxa_16S[(mid + 1L):n16] else character()
taxa_16S_right <- order_taxa_16S[seq_len(mid)]

heatmap_data_16S_left <- heatmap_data_16S %>%
  filter(taxon_display %in% taxa_16S_left) %>%
  mutate(taxon_display = factor(taxon_display, levels = taxa_16S_left))

heatmap_data_16S_right <- heatmap_data_16S %>%
  filter(taxon_display %in% taxa_16S_right) %>%
  mutate(taxon_display = factor(taxon_display, levels = taxa_16S_right))

message(
  "16S split (swapped): left=", length(taxa_16S_left),
  " right=", length(taxa_16S_right)
)

make_heatmap <- function(heatmap_data, y_lab, tag, reorder = TRUE, show_x_title = FALSE) {
  if (nrow(heatmap_data) == 0L) {
    return(
      ggplot() +
        annotate(
          "text", x = 0.5, y = 0.5,
          label = "No taxa with global p < 0.05",
          size = 4.5
        ) +
        labs(
          x = if (isTRUE(show_x_title)) "Plastic level" else NULL,
          y = y_lab,
          tag = tag
        ) +
        theme_void() +
        theme(
          axis.title = element_text(size = 16),
          plot.tag = element_text(size = 28, face = "bold"),
          plot.margin = margin(20, 20, 20, 20)
        )
    )
  }

  if (isTRUE(reorder)) {
    order_taxa <- order_taxa_plastic_lfc_columns(
      heatmap_data,
      id_col = "taxon_display",
      level_col = "plastic_level",
      value_col = "lfc",
      col_order = c("Low", "Medium", "High"),
      tol = 0.15
    )
    heatmap_data <- heatmap_data %>%
      mutate(taxon_display = factor(as.character(taxon_display), levels = order_taxa))
  }

  taxon_display_passed_ss <- heatmap_data %>%
    select(taxon_display, passed_ss) %>%
    distinct() %>%
    deframe()

  ggplot(heatmap_data, aes(x = plastic_level, y = taxon_display, fill = lfc)) +
    geom_tile(color = "white", linewidth = 0.5) +
    labs(
      x = if (isTRUE(show_x_title)) "Plastic level" else NULL,
      y = y_lab,
      tag = tag
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(
        color = ifelse(
          taxon_display_passed_ss[levels(heatmap_data$taxon_display)],
          "#000000",
          "#8a8a8a"
        ),
        size = 9
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
  guide = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    direction = "horizontal",
    barheight = unit(0.6, "cm"),
    barwidth = unit(8, "cm")
  )
)

p_16S_left <- make_heatmap(
  heatmap_data_16S_left, "Prokaryotic taxa", "A",
  reorder = FALSE, show_x_title = FALSE
) + shared_scale

p_16S_right <- make_heatmap(
  heatmap_data_16S_right, NULL, "B",
  reorder = FALSE, show_x_title = FALSE
) + shared_scale

p_18S <- make_heatmap(
  heatmap_data_18S, "Eukaryotic taxa", "C",
  show_x_title = TRUE
) + shared_scale

n18 <- if (nrow(heatmap_data_18S)) n_distinct(heatmap_data_18S$taxon_display) else 1L
n16_panel <- max(length(taxa_16S_left), length(taxa_16S_right), 1L)

# Top: A | B (equal). Bottom: C same width as A; legend glued to left of right cell (beside C)
layout <- "
AB
CL
"

p_combined <- p_16S_left + p_16S_right + p_18S + patchwork::guide_area() +
  patchwork::plot_layout(
    design = layout,
    widths = c(1, 1),
    heights = c(n16_panel, max(n18, 4L)),
    guides = "collect",
    axes = "keep"
  ) &
  theme(
    axis.text.x = element_text(size = 13),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    # Legend in bottom-right cell, left-aligned beside C (on-page)
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = c(0, 0.5),
    legend.box.just = "left",
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 16),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    plot.tag = element_text(size = 28, face = "bold")
  )

# Landscape A4 (297 × 210 mm) with margins for LaTeX sidewaysfigure
a4_w_in <- 297 / 25.4
a4_h_in <- 210 / 25.4
margin_in <- 0.5

pdf(
  "figures/fig_5_global_analysis_heatmap_stable_taxa_combined.pdf",
  width = a4_w_in - 2 * margin_in,
  height = a4_h_in - 2 * margin_in,
  family = "Helvetica",
  useDingbats = FALSE
)
print(p_combined)
dev.off()

message("Wrote figures/fig_5_global_analysis_heatmap_stable_taxa_combined.pdf")
