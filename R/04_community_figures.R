source("R/01_load_files.R")
source("R/functions.R")

theme_set(theme_bw(base_size = 14))

p1 <- plot_stacked_barchart(gene = "16S", taxa_level = "Phylum", n_taxa = 13, italics = FALSE, tag = "A") 
p2 <- plot_stacked_barchart(gene = "18S", taxa_level = "Family", n_taxa = 13, italics = FALSE, tag = "B") 

p3 <- patchwork::wrap_plots(p1, p2, nrow = 2)
ggsave("figures/fig_3_stacked_barcharts.png", p3, width = 13, height = 15, dpi = 300, scale = 0.8)

# Bubble plots of specific lineages

phy_16S_cyano <- phy_16S %>% 
    subset_taxa(Phylum == "Cyanobacteriota")

phy_16S_proteo <- phy_16S %>% 
    subset_taxa(Phylum == "Pseudomonadota")
   
taxa_level <- "Genus"
n_taxa <- 31 
italics <- TRUE
rel_abund <- rel_abund_phy(get(glue("phy_16S_cyano")), 
                               meta_data = TRUE, 
                               taxa_level = taxa_level, 
                               var = "sample_id") %>% 
            {if(italics) {
                # Apply italics but remove them from Family taxa
                d <- taxon_italics(.)
                d %>% mutate(taxon = ifelse(grepl("Family", taxon), 
                                          gsub("^\\*|\\*$", "", taxon), 
                                          taxon))
            } else .} %>%
            pool_taxa(n_taxa = n_taxa, keep_metadata = TRUE) %>%
            filter(taxon != "Other") %>%
            arrange_taxa() %>%
            mutate(sample_id = factor(sample_id)) %>%
            mutate(sample_id = fct_reorder(sample_id, plastic_concentration, .desc = F)) %>%
            mutate(Date = factor(Date, levels = c(3, 6, 9)))

# Get max and min rel_abund for consistent scale between plots
rel_1 <- rel_abund %>% filter(Location == "WS")
rel_2 <- rel_abund %>% filter(Location == "MS")

rel_abund_min <- min(rel_abund$rel_abund, na.rm = TRUE)
rel_abund_max <- max(rel_abund$rel_abund, na.rm = TRUE)

p1 <- rel_1 %>% 
    bubbler::bubble_plot(x = "sample_id", italics = italics) + 
    scale_fill_continuous(limits = c(rel_abund_min, rel_abund_max), name = "Relative abundance") +
    scale_size_continuous(range = c(0, 6), limits = c(rel_abund_min, rel_abund_max), name = "Relative abundance") +
    facet_wrap(~ Date + Location, scales = "free_x", 
        labeller = labeller(
            Date = function(x) paste0("Week ", x),
            Location = function(x) ifelse(x == "WS", "Wall strip", ifelse(x == "MS", "Microscope slide", x))
        )
    ) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    labs(y = "Cyanobacteriota" , x = "Sample name")

p2 <- rel_2 %>% 
    bubbler::bubble_plot(x = "sample_id", italics = italics) + 
    scale_fill_continuous(limits = c(rel_abund_min, rel_abund_max), name = "Relative abundance") +
    scale_size_continuous(range = c(0, 6), limits = c(rel_abund_min, rel_abund_max), name = "Relative abundance") +
    facet_wrap(~ Date + Location, scales = "free_x", 
        labeller = labeller(
            Date = function(x) paste0("Week ", x),
            Location = function(x) ifelse(x == "WS", "Wall strip", ifelse(x == "MS", "Microscope slide", x))
        )
    ) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    labs(y = "Cyanobacteriota" , x = "Sample name")


layout <- "AAAB"
p3 <- p1 + p2 + patchwork::plot_layout(guides = "collect", axes = "collect", design = layout)
ggsave("figures/fig_S2_bubble_plots_cyanobacteriota.png", p3, width = 10, height = 8, dpi = 300, scale = 1)

taxa_level <- "Genus"
n_taxa <- 31 
italics <- TRUE
rel_abund <- rel_abund_phy(get(glue("phy_16S_proteo")), 
                               meta_data = TRUE, 
                               taxa_level = taxa_level, 
                               var = "sample_id") %>% 
            {if(italics) {
                # Apply italics but remove them from Family taxa
                d <- taxon_italics(.)
                d %>% mutate(taxon = ifelse(grepl("Family", taxon), 
                                          gsub("^\\*|\\*$", "", taxon), 
                                          taxon))
            } else .} %>%
            pool_taxa(n_taxa = n_taxa, keep_metadata = TRUE) %>%
            filter(taxon != "Other") %>%
            arrange_taxa() %>%
            mutate(sample_id = factor(sample_id)) %>%
            mutate(sample_id = fct_reorder(sample_id, plastic_concentration, .desc = F)) %>%
            mutate(Date = factor(Date, levels = c(3, 6, 9)))

# Get max and min rel_abund for consistent scale between plots
rel_1 <- rel_abund %>% filter(Location == "WS")
rel_2 <- rel_abund %>% filter(Location == "MS")

rel_abund_min <- min(rel_abund$rel_abund, na.rm = TRUE)
rel_abund_max <- max(rel_abund$rel_abund, na.rm = TRUE)

p1 <- rel_1 %>% 
    bubbler::bubble_plot(x = "sample_id", italics = italics) + 
    scale_fill_continuous(limits = c(rel_abund_min, rel_abund_max), name = "Relative abundance") +
    scale_size_continuous(range = c(0, 6), limits = c(rel_abund_min, rel_abund_max), name = "Relative abundance") +
    facet_wrap(~ Date + Location, scales = "free_x", 
            labeller = labeller(
                Date = function(x) paste0("Week ", x),
                Location = function(x) ifelse(x == "WS", "Wall strip", ifelse(x == "MS", "Microscope slide", x))
            )) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    labs(y = "Pseudomonadota" , x = "Sample name")

p2 <- rel_2 %>% 
    bubbler::bubble_plot(x = "sample_id", italics = italics) + 
    scale_fill_continuous(limits = c(rel_abund_min, rel_abund_max), name = "Relative abundance") +
    scale_size_continuous(range = c(0, 6), limits = c(rel_abund_min, rel_abund_max), name = "Relative abundance") +
    facet_wrap(~ Date + Location, scales = "free_x", 
            labeller = labeller(
                Date = function(x) paste0("Week ", x),
                Location = function(x) ifelse(x == "WS", "Wall strip", ifelse(x == "MS", "Microscope slide", x))
            )) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    labs(y = "Pseudomonadota" , x = "Sample name")


layout <- "AAAB"
p3 <- p1 + p2 + patchwork::plot_layout(guides = "collect", axes = "collect", design = layout)
ggsave("figures/fig_S1_bubble_plots_pseudomonadota.png", p3, width = 10, height = 8, dpi = 300, scale = 1)

# Save phyloseq objects as RDS files
# saveRDS(phy_16S_cyano, "out/phy_particles_16S_cyano.RDS")
# saveRDS(phy_16S_proteo, "out/phy_particles_16S_proteo.RDS")
