source("R/01_load_files.R")
source("R/functions.R")

# Figure 1: experimental design and limnocorral illustration
phy_16S_9 <- phy_16S %>%
    subset_samples(Date == 9) %>%
    subset_samples(Location == "WS")

rn <- rownames(sample_data(phy_16S_9))
meta_data <- sample_data(phy_16S_9) %>% 
    as_tibble()
meta_data <- tibble(sample_id = rn, meta_data)

# for plastic concentration
cutoffs <- c(-Inf, 5, 25, 2000,Inf)
plastic_level <-cut(meta_data$plastic_concentration, breaks = cutoffs, labels =c("none", "low", "medium", "high"))

tb <-tibble(plastic_concentration = meta_data$plastic_concentration, 
        plastic_level = plastic_level,
        Corral = meta_data$CorralLetter
        ) %>%
    distinct() %>%
    arrange(plastic_concentration)

p1 <- ggplot(tb, aes(x = reorder(Corral, plastic_concentration), y = plastic_concentration, fill = plastic_level)) +
  geom_point(size = 7, shape = 21, show.legend = TRUE) +
  geom_text(aes(label = plastic_concentration), vjust = -1, size = 5, fontface = "bold") +
  scale_fill_viridis_d(name = "Plastic Level", labels = c("None", "Low", "Medium", "High")) +
  scale_x_discrete(expand = expansion(add = c(1, 1))) +
  scale_y_log10(
    name = "Particles/L",
    limits = c(0.5, 300000),
    labels = scales::comma
  ) +
  labs(x = "Corral Letter", tag = "B") +
  plastic_theme

p2 <- meta_data %>% 
    select(particles_total_d20, CorralLetter, plastic_concentration, plastic_level) %>%
    ggplot(aes(x = reorder(CorralLetter, plastic_concentration), y = particles_total_d20, fill = plastic_level)) +
    geom_point(size = 7, shape = 21, show.legend = TRUE) +
    scale_fill_viridis_d(name = "Plastic Level", labels = c("None", "Low", "Medium", "High")) +
    scale_x_discrete(expand = expansion(add = c(1, 1))) + 
    ylab(expression("Particles/cm"^2)) + 
    scale_y_continuous(limits = c(0, 30))  + 
    labs(x = "Corral Letter", tag = "C") +
    plastic_theme

# Read in the image
img <- magick::image_read("figures/limnocorral.png")
# Convert to a rasterGrob for ggplot
img_grob <- grid::rasterGrob(img, interpolate = TRUE)

# Create a blank ggplot with the image
p_img <- ggplot() + 
  annotation_custom(img_grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_void() +
  labs(tag = "A") +
  theme(plot.tag = element_text(face = "bold", size = 28))

# Combine all three plots
p3 <- (p_img + p1 + p2 + patchwork::plot_layout(ncol = 3, guides = "collect", axes = "collect")) & 
  theme(legend.position = "bottom")
ggsave("figures/fig_1_plastic_density.png", p3, width = 20, height = 7, dpi = 300, scale = 0.8)

