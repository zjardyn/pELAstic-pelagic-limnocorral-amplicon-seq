source("R/01_load_files.R")
source("R/functions.R")

# theme_set(theme_linedraw(base_size = 12))

p1 <- pca_plot0(phy_16S, 
    # transform = "hellinger",
    colour = Date, 
    shape = Location,
    r2_cutoff = 0.017,
    title = "16S", 
    italics = FALSE,
    point_size = 3) + 
    scale_shape_manual(name = "Location",
                      values = c("MS" = 16, "WS" = 17),
                      labels = c("MS" = "Microscope slide", "WS" = "Wall strip")) + 
    scale_color_viridis_d(name = "Week") +
    labs(tag = "A") +
    plastic_theme

p2 <- pca_plot0(phy_18S, 
    # transform = "hellinger",
    colour = Date, 
    shape = Location,
    r2_cutoff = 0.03,
    title = "18S", 
    italics = FALSE,
    point_size = 3) + 
    scale_shape_manual(name = "Location",
                      values = c("MS" = 16, "WS" = 17),
                      labels = c("MS" = "Microscope slide", "WS" = "Wall strip")) + 
    scale_color_viridis_d(name = "Week") +
    labs(tag = "B") +
    plastic_theme

p3 <- (p1 + p2 + patchwork::plot_layout(ncol = 2, guides = "collect", axes = "collect")) +
  patchwork::plot_annotation(
    title = "PCA of CLR-transformed ASV counts",
    theme = theme(plot.title = element_text(size = 32, face = "bold"))
  ) &
  plastic_theme &
  theme(legend.position = "bottom")

ggsave("figures/fig_2_pca.png", p3, width = 18, height = 10, dpi = 300, scale = 0.8)


# p3 <- pca_plot0(phy_16S, 
#     colour = plastic_conc_log, 
#     shape = Date,
#     r2_cutoff = 0.017,
#     title = "16S PCA", 
#     italics = FALSE,
#     point_size = 3) + 
#     scale_color_viridis_c(name = "log(Plastic concentration)") +
#     scale_shape_discrete(name = "Week")

# ggsave("figures/16S_pca_plastic.png", p3, width = 10, height = 8, dpi = 300, scale = 0.8)

# p4 <- pca_plot0(phy_18S, 
#     colour = plastic_conc_log, 
#     shape = Date,
#     r2_cutoff = 0.03,
#     title = "18S PCA", 
#     italics = FALSE,
#     point_size = 3) + 
#     scale_color_viridis_c(name = "log(Plastic concentration)") +
#     scale_shape_discrete(name = "Week")

# ggsave("figures/18S_pca_plastic.png", p4, width = 10, height = 8, dpi = 300, scale = 0.8)
