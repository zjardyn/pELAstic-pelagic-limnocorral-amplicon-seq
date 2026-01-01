source("R/01_load_files.R")
source("R/functions.R")

rf_16S <-readRDS("output/rf_runs_1000_16S.rds")
rf_18S <-readRDS("output/rf_runs_1000_18S.rds")

phy_16S_9_ws <- phy_16S %>%
    subset_samples(Date == 9) %>%
    subset_samples(Location == "WS") %>%
    subset_samples(!is.na(particles_total_d20)) 

phy_18S_9_ws <- phy_18S %>%
    subset_samples(Date == 9) %>%
    subset_samples(Location == "WS") %>%
    subset_samples(!is.na(particles_total_d20))

top_n <- 20
df <- rf_16S %>%
        slice_max(mean_importance, n = top_n) %>%
        mutate(base = ifelse(is.na(Genus) | Genus == "", asv, Genus),
               base = make.unique(base),
               label = paste0(base, " (", round(topk_freq, 2), ")")) %>%
         mutate(ymin = mean_importance - sd_importance,
         ymax = mean_importance + sd_importance)

p1 <- ggplot(df, aes(x = reorder(label, mean_importance), y = mean_importance)) +
  geom_col(fill = "grey40") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2, color = "black") +
  coord_flip() +
  labs(x = "Bacterial and Archaeal Taxa", y = "Mean permutation importance ± SD", tag = "A") +
  theme(
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    plot.tag = element_text(size = 28, face = "bold")
  )

df <- rf_18S %>%
        slice_max(mean_importance, n = top_n) %>%
        mutate(base = ifelse(is.na(Genus) | Genus == "", asv, Genus),
               base = make.unique(base),
               label = paste0(base, " (", round(topk_freq, 2), ")")) %>%
         mutate(ymin = mean_importance - sd_importance,
         ymax = mean_importance + sd_importance)

p2 <- ggplot(df, aes(x = reorder(label, mean_importance), y = mean_importance)) +
  geom_col(fill = "grey40") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2, color = "black") +
  coord_flip() +
  labs(x = "Eukaryotic Taxa", y = "Mean permutation importance ± SD", tag = "C") +
  theme(
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    plot.tag = element_text(size = 28, face = "bold")
  )

### Bubble plots for top 20 taxa - 16S
rf_taxa_16S <- rf_16S %>%
    pull(asv)

# make Hellinger-transformed relative abundances (sqrt of proportions)
otu_tbl_16S <- otu_table(phy_16S_9_ws) %>%
    as.data.frame() %>%
    rownames_to_column("asv") %>%
    as_tibble() %>%
    pivot_longer(-asv, names_to = "SampleID", values_to = "abundance") %>%
    group_by(SampleID) %>%
    mutate(rel_abund = sqrt(abundance / sum(abundance))) %>%
    ungroup() %>%
    select(asv, SampleID, rel_abund) %>%
    filter(asv %in% rf_taxa_16S)  %>% 
    pivot_wider(names_from = "SampleID", values_from = "rel_abund") 

# order rows by rf_16S importance (descending)
asv_order_16S <- rf_16S %>% arrange(desc(mean_importance)) %>% pull(asv)
otu_tbl_16S <- otu_tbl_16S %>%
    mutate(asv = factor(asv, levels = asv_order_16S)) %>%
    arrange(asv)

abund_tb_16S <- otu_tbl_16S %>% 
    pivot_longer(-asv, names_to = "SampleID", values_to = "rel_abund") %>% 
    inner_join(rf_16S, by = "asv") %>%
    arrange(desc(mean_importance)) 

top_n <- 20
top_asv_16S <- rf_16S %>%
        slice_max(mean_importance, n = top_n)  %>% pull(asv)

abund_tb_16S <- abund_tb_16S %>%
    filter(asv %in% top_asv_16S) %>%
    mutate(asv = factor(asv, levels = top_asv_16S),
           Genus = fct_reorder(Genus, mean_importance)) 

meta_16S <- sample_data(phy_16S_9_ws) 
rn <- rownames(meta_16S)

meta_16S_tb <- meta_16S %>% as_tibble() %>% bind_cols(SampleID = rn)

abund_tb_16S_joined <- abund_tb_16S %>%
    inner_join(meta_16S_tb, by = "SampleID") %>%
    mutate(Corral = fct_reorder(CorralLetter, plastic_concentration)) %>%
    mutate(Genus_label = paste0(Genus, " (" , round(mean_importance, 2), ") " , " (", round(topk_freq, 2), ")")) %>%
    mutate(Genus_label = fct_reorder(Genus_label, mean_importance))

### Bubble plots for top 20 taxa - 18S
rf_taxa_18S <- rf_18S %>%
    pull(asv)

# make Hellinger-transformed relative abundances (sqrt of proportions)
otu_tbl_18S <- otu_table(phy_18S_9_ws) %>%
    as.data.frame() %>%
    rownames_to_column("asv") %>%
    as_tibble() %>%
    pivot_longer(-asv, names_to = "SampleID", values_to = "abundance") %>%
    group_by(SampleID) %>%
    mutate(rel_abund = sqrt(abundance / sum(abundance))) %>%
    ungroup() %>%
    select(asv, SampleID, rel_abund) %>%
    filter(asv %in% rf_taxa_18S)  %>% 
    pivot_wider(names_from = "SampleID", values_from = "rel_abund") 

# order rows by rf_18S importance (descending)
asv_order_18S <- rf_18S %>% arrange(desc(mean_importance)) %>% pull(asv)
otu_tbl_18S <- otu_tbl_18S %>%
    mutate(asv = factor(asv, levels = asv_order_18S)) %>%
    arrange(asv)

abund_tb_18S <- otu_tbl_18S %>% 
    pivot_longer(-asv, names_to = "SampleID", values_to = "rel_abund") %>% 
    inner_join(rf_18S, by = "asv") %>%
    arrange(desc(mean_importance)) 

top_n <- 20
top_asv_18S <- rf_18S %>%
        slice_max(mean_importance, n = top_n)  %>% pull(asv)

abund_tb_18S <- abund_tb_18S %>%
    filter(asv %in% top_asv_18S) %>%
    mutate(asv = factor(asv, levels = top_asv_18S),
           Genus = fct_reorder(Genus, mean_importance)) 

meta_18S <- sample_data(phy_18S_9_ws) 
rn <- rownames(meta_18S)

meta_18S_tb <- meta_18S %>% as_tibble() %>% bind_cols(SampleID = rn)

abund_tb_18S_joined <- abund_tb_18S %>%
    inner_join(meta_18S_tb, by = "SampleID") %>%
    mutate(Corral = fct_reorder(CorralLetter, plastic_concentration)) %>%
    mutate(Genus_label = paste0(Genus, " (" , round(mean_importance, 2), ") " , " (", round(topk_freq, 2), ")")) %>%
    mutate(Genus_label = fct_reorder(Genus_label, mean_importance))

# Calculate maximum Hellinger-transformed relative abundance across both datasets for consistent scaling
# Hellinger transformation (sqrt of proportions) provides scale invariance
max_hellinger_rel_abund <- max(
  max(abund_tb_16S_joined$rel_abund, na.rm = TRUE),
  max(abund_tb_18S_joined$rel_abund, na.rm = TRUE),
  na.rm = TRUE
)

p3 <- ggplot(abund_tb_16S_joined, aes(x = Corral, y = Genus_label, size = rel_abund)) +
geom_point(alpha = 0.8) +
  scale_size_continuous(range = c(0, 8), limits = c(0, max_hellinger_rel_abund), name = "Hellinger abundance") +
  labs(x = "Corral Letter", y = NULL, tag = "B")

p4 <- ggplot(abund_tb_18S_joined, aes(x = Corral, y = Genus_label, size = rel_abund)) +
geom_point(alpha = 0.8) +
  scale_size_continuous(range = c(0, 8), limits = c(0, max_hellinger_rel_abund), name = "Hellinger abundance") +
  labs(x = "Corral Letter", y = NULL, tag = "D")

# Create four-panel layout in columns: A (16S importance) / C (16S bubbles) | B (18S importance) / D (18S bubbles)
layout <- "
    AB
    CD
"

# Apply y-axis removal and text sizes to p3 and p4 before combining
p3 <- p3 + theme(
    axis.text.y = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 13),
    axis.title.x = element_text(size = 16),
    plot.tag = element_text(size = 28, face = "bold"),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 16)
  )
p4 <- p4 + theme(
    axis.text.y = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 13),
    axis.title.x = element_text(size = 16),
    plot.tag = element_text(size = 28, face = "bold"),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 16)
  )

p_combined <- (p1 + p3 + p2 + p4) + 
  patchwork::plot_layout(design = layout, guides = "collect", axes = "collect") &
  theme(
    legend.position = "bottom", 
    legend.justification = "left",
    legend.box.just = "left",
    legend.margin = margin(t = 0, r = 0, b = 0, l = -20, unit = "pt")
  )

pdf(
  "figures/fig_6_random_forests_combined.pdf",
  width = 11 * 0.8,
  height = 13 * 0.8,
  family = "Helvetica",
  useDingbats = FALSE
)

print(p_combined)
dev.off()