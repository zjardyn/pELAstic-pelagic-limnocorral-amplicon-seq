source("R/01_load_files.R")
source("R/functions.R")

# theme_set(theme_linedraw(base_size = 12))

p1 <- pca_plot0(phy_16S, 
    # transform = "hellinger",
    colour = Date, 
    shape = Location,
    # r2_cutoff = 0.017,
    r2_cutoff = 0.019,
    tax_lab_size = 6,
    title = "16S rRNA", 
    italics = TRUE,
    point_size = 5,
    transform = "clr") + 
    scale_shape_manual(name = "Location",
                      values = c("MS" = 21, "WS" = 24),
                      labels = c("MS" = "Microscope slide", "WS" = "Wall strip")) + 
    scale_fill_viridis_d(name = "Week") +
    guides(shape = guide_legend(override.aes = list(color = "black", stroke = 0.5)),
           fill = guide_legend(override.aes = list(color = "black", stroke = 0.5, shape = 21))) +
    labs(tag = "A") +
    plastic_theme

p2 <- pca_plot0(phy_18S, 
    # transform = "hellinger",
    colour = Date, 
    shape = Location,
    # r2_cutoff = 0.03,
    r2_cutoff = 0.035,
    tax_lab_size = 6,
    title = "18S rRNA", 
    italics = TRUE,
    point_size = 5,
    transform = "clr") + 
    scale_shape_manual(name = "Location",
                      values = c("MS" = 21, "WS" = 24),
                      labels = c("MS" = "Microscope slide", "WS" = "Wall strip")) + 
    scale_fill_viridis_d(name = "Week") +
    guides(shape = guide_legend(override.aes = list(color = "black", stroke = 0.5)),
           fill = guide_legend(override.aes = list(color = "black", stroke = 0.5, shape = 21))) +
    labs(tag = "B") +
    plastic_theme

p3 <- (p1 + p2 + patchwork::plot_layout(nrow = 2, guides = "collect", axes = "collect")) +
  patchwork::plot_annotation(
    theme = theme(plot.title = element_text(size = 32, face = "bold"))
  ) &
  plastic_theme &
  theme(legend.position = "bottom", legend.box = "vertical")

# Main: coloured by week
pdf(
  "figures/fig3_2.pdf",
  width = 12 * 0.85,
  height = 20 * 0.85,
  family = "Helvetica",
  useDingbats = FALSE
)
print(p3)
dev.off()

# Week-9 wall strips at the two highest plastic treatments (G = 7071, D = 29240)
hi_week9_ws <- sample_data(phy_16S) %>%
  as_tibble(rownames = "Sample") %>%
  filter(Location == "WS", as.character(Date) == "9", plastic_level == "high") %>%
  pull(Sample)

hi_week9_ws_18S <- sample_data(phy_18S) %>%
  as_tibble(rownames = "Sample") %>%
  filter(Location == "WS", as.character(Date) == "9", plastic_level == "high") %>%
  pull(Sample)

message("Highlighting 16S: ", paste(hi_week9_ws, collapse = ", "))
message("Highlighting 18S: ", paste(hi_week9_ws_18S, collapse = ", "))

# Supplementary: coloured by plastic level; highlight high week-9 wall strips
p1_plastic <- pca_plot0(phy_16S,
    colour = plastic_level,
    shape = Location,
    r2_cutoff = 0.019,
    tax_lab_size = 6,
    title = "16S rRNA",
    italics = TRUE,
    point_size = 5,
    transform = "clr",
    highlight_samples = hi_week9_ws,
    show_vectors = FALSE,
    time_site_ellipses = TRUE) +
    scale_shape_manual(name = "Location",
                      values = c("MS" = 21, "WS" = 24),
                      labels = c("MS" = "Microscope slide", "WS" = "Wall strip")) +
    scale_fill_viridis_d(name = "Plastic level", labels = c("None", "Low", "Medium", "High")) +
    guides(shape = guide_legend(override.aes = list(color = "black", stroke = 0.5)),
           fill = guide_legend(override.aes = list(color = "black", stroke = 0.5, shape = 21))) +
    labs(tag = "A") +
    plastic_theme

p2_plastic <- pca_plot0(phy_18S,
    colour = plastic_level,
    shape = Location,
    r2_cutoff = 0.035,
    tax_lab_size = 6,
    title = "18S rRNA",
    italics = TRUE,
    point_size = 5,
    transform = "clr",
    highlight_samples = hi_week9_ws_18S,
    show_vectors = FALSE,
    time_site_ellipses = TRUE) +
    scale_shape_manual(name = "Location",
                      values = c("MS" = 21, "WS" = 24),
                      labels = c("MS" = "Microscope slide", "WS" = "Wall strip")) +
    scale_fill_viridis_d(name = "Plastic level", labels = c("None", "Low", "Medium", "High")) +
    guides(shape = guide_legend(override.aes = list(color = "black", stroke = 0.5)),
           fill = guide_legend(override.aes = list(color = "black", stroke = 0.5, shape = 21))) +
    labs(tag = "B") +
    plastic_theme

p3_plastic <- (p1_plastic + p2_plastic + patchwork::plot_layout(nrow = 2, guides = "collect", axes = "collect")) +
  patchwork::plot_annotation(
    theme = theme(plot.title = element_text(size = 32, face = "bold"))
  ) &
  plastic_theme &
  theme(legend.position = "bottom", legend.box = "vertical")

pdf(
  "figures/figS3_2.pdf",
  width = 12 * 0.85,
  height = 20 * 0.85,
  family = "Helvetica",
  useDingbats = FALSE
)
print(p3_plastic)
dev.off()
