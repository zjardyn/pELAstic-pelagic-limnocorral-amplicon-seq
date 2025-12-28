source("R/01_load_files.R")
source("R/functions.R")
# library(pairwiseAdonis)
set.seed(123)

phy_16S_9 <- phy_16S %>%
    subset_samples(Date == 9) %>%
    subset_samples(Location == "WS")

phy_18S_9 <- phy_18S %>%
    subset_samples(Date == 9) %>%
    subset_samples(Location == "WS")

# Add PERMANOVA analysis
# First, get the distance matrix and metadata
dist_matrix <- phy_16S_9 %>%
    tax_transform("identity", rank = "Genus") %>%
    dist_calc("aitchison") %>%
    dist_get()

metadata <- sample_data(phy_16S_9) %>%
    as_tibble(rownames = "sample_id")

# Run PERMANOVA using adonis2 from vegan package
# Test the effect of plastic_level on community composition
permanova_result <- adonis2(dist_matrix ~ plastic_level, 
                           data = metadata, 
                           permutations = 9999)
pval <- permanova_result$`Pr(>F)`[1]

# pairwise.adonis2(dist_matrix ~ plastic_level, data = metadata, permutations = 9999)
#Aitchison NMDS, PERMANOVA of plastic level 
p5 <- phy_16S_9 %>%
    tax_transform("identity", rank = "Genus") %>%
    dist_calc("aitchison") %>%
    ord_calc("NMDS") %>%
    ord_plot(color = "particles_total_d20", size = 3, auto_caption = NA) +
    scale_colour_viridis_c(name = expression("Particles/cm"^2)) + 
    ggtitle(glue("16S")) +
    labs(tag = "A") +
    theme_bw() +
    theme(plot.title = element_markdown(), plot.tag = element_text(face = "bold", size = 16))

# ggsave("figures/aitchison_nmds_16S.png", p5, width = 10, height = 8, dpi = 300, scale = 0.7)


p5.1 <- phy_16S_9 %>%
    tax_transform("identity", rank = "Genus") %>%
    dist_calc("aitchison") %>%
    ord_calc("NMDS") %>%
    ord_plot(color = "plastic_level", size = 3, auto_caption = NA) +
    scale_colour_viridis_d(name = "Plastic Level") + 
    ggrepel::geom_text_repel(aes(label = plastic_concentration), size = 3, show.legend = FALSE) +
    ggtitle(glue("16S, PERMANOVA p = {round(pval, 3)}")) +
    labs(tag = "A") +
    theme_bw() +
    theme(plot.title = element_markdown(), plot.tag = element_text(face = "bold", size = 16))

# ggsave("figures/aitchison_nmds_16S_plastic_level.png", p5.1, width = 10, height = 8, dpi = 300, scale = 0.7)

# 18S
# Add PERMANOVA analysis
# First, get the distance matrix and metadata
dist_matrix <- phy_18S_9 %>%
    tax_transform("identity", rank = "Genus") %>%
    dist_calc("aitchison") %>%
    dist_get()

metadata <- sample_data(phy_18S_9) %>%
    as_tibble(rownames = "sample_id")

# Run PERMANOVA using adonis2 from vegan package
# Test the effect of plastic_level on community composition
permanova_result <- adonis2(dist_matrix ~ plastic_level, 
                           data = metadata, 
                           permutations = 9999)

pval <- permanova_result$`Pr(>F)`[1]

# pairwise.adonis2(dist_matrix ~ plastic_level, data = metadata, permutations = 9999)

p6 <- phy_18S_9 %>%
    tax_transform("identity", rank = "Genus") %>%
    dist_calc("aitchison") %>%
    ord_calc("NMDS") %>%
    ord_plot(color = "particles_total_d20", size = 3, auto_caption = NA) +
    scale_colour_viridis_c(name = expression("Particles/cm"^2)) + 
    ggtitle(glue("18S")) +
    labs(tag = "B") +
    theme_bw() +
    theme(plot.title = element_markdown(), plot.tag = element_text(face = "bold", size = 16))

# ggsave("figures/aitchison_nmds_18S.png", p6, width = 10, height = 8, dpi = 300, scale = 0.7)

p6.1 <- phy_18S_9 %>%
    tax_transform("identity", rank = "Genus") %>%
    dist_calc("aitchison") %>%
    ord_calc("NMDS") %>%
    ord_plot(color = "plastic_level", size = 3, auto_caption = NA) +
    scale_colour_viridis_d(name = "Plastic Level") + 
    ggrepel::geom_text_repel(aes(label = plastic_concentration), size = 3, show.legend = FALSE) +
    ggtitle(glue("18S, PERMANOVA p = {round(pval, 3)}")) +
    labs(tag = "B") +
    theme_bw() +
    theme(plot.title = element_markdown(), plot.tag = element_text(face = "bold", size = 16))

# ggsave("figures/aitchison_nmds_18S_plastic_level.png", p6.1, width = 10, height = 8, dpi = 300, scale = 0.7)


# Combine the two NMDS plots side by side with a shared title and consistent title styling
p_comb <- (p5.1 + p6.1) + 
  patchwork::plot_layout(ncol = 2, guides = "collect", axes = "collect") + 
  plot_annotation(title = "NMDS of ASV Aitchison distances", theme = theme(plot.title = element_text(face = "bold", size = 16))) & plastic_theme
ggsave("figures/fig_4_aitchison_nmds_16S_18S.png", p_comb & theme(legend.position = "bottom"), width = 18, height = 10, dpi = 300, scale = 0.8) 

# p_comb <- (p5 + p6) + 
#   patchwork::plot_layout(ncol = 2, guides = "collect", axes = "collect") + 
#   plot_annotation(title = "NMDS of ASV Aitchison distances", theme = theme(plot.title = element_text(face = "bold", size = 16))) & plastic_theme
# ggsave("figures/fig_S5_aitchison_nmds_16S_18S.png", p_comb & theme(legend.position = "bottom"), width = 18, height = 10, dpi = 300, scale = 0.8) 
