source("R/01_load_files.R")
source("R/functions.R")
library(tidyverse)
library(tidyr)
library(scales)


rf_16S <-readRDS("out/rf_runs_1000_16S.rds")
rf_18S <-readRDS("out/rf_runs_1000_18S.rds")

# nrow(rf_16S)
# nrow(rf_18S)

# if the variable "total particles/m2" is used, remove NAs
phy_16S_9_ws <- phy_16S %>%
    subset_samples(Date == 9) %>%
    subset_samples(Location == "WS") %>%
    subset_samples(!is.na(particles_total_d20)) %>%
    ps_mutate(plastic_retained = if_else(particles_total_d20 > 10,  "retained", "not_retained")) 

phy_18S_9_ws <- phy_18S %>%
    subset_samples(Date == 9) %>%
    subset_samples(Location == "WS") %>%
    subset_samples(!is.na(particles_total_d20)) %>%
    ps_mutate(plastic_retained = if_else(particles_total_d20 > 10,  "retained", "not_retained"))

# rf_16S %>% as_tibble() %>%
#     filter(topk_freq > 0.5) 
# rf_18S %>% as_tibble() %>% 
#     filter(topk_freq > 0.5)  %>%
#     print(n = 50)

top_n <- 20
df <- rf_16S %>%
        # filter(topk_freq > 0.5) %>%
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
  labs(x = "Genus", y = "Mean permutation importance ± SD",
       title = paste0("Top ", top_n, " taxa (16S)"), tag = "A") +
  theme_minimal() +
  theme(plot.tag = element_text(face = "bold"))

df <- rf_18S %>%
        # filter(topk_freq > 0.5) %>%
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
  labs(x = "Genus", y = "Mean permutation importance ± SD",
       title = paste0("Top ", top_n, " taxa (18S)"), tag = "A") +
  theme_minimal() +
  theme(plot.tag = element_text(face = "bold"))

### Bubble plots for top 20 taxa - 16S
rf_taxa_16S <- rf_16S %>%
    pull(asv)

# have to somehow make abundances, here I use rclr scaled to 0-1 but it might not be right.
otu_tbl_16S <- otu_table(phy_16S_9_ws) %>%
    as.data.frame() %>%
    rownames_to_column("asv") %>%
    filter(asv %in% rf_taxa_16S)  %>% 
    as_tibble() %>%
    pivot_longer(-asv, names_to = "SampleID", values_to = "abundance") %>%
    group_by(SampleID) %>%
    mutate(rel_abund = abundance / sum(abundance)) %>%
    ungroup() %>%
    select(asv, SampleID, rel_abund) %>%
    pivot_wider(names_from = "SampleID", values_from = "rel_abund") 
    # column_to_rownames("asv") %>% 
    # as.matrix() %>%
    # decostand("rclr", MARGIN = 1) %>%
    # as.data.frame() %>%
    # rownames_to_column("asv") %>%
    # as_tibble() %>%
    # pivot_longer(-asv, names_to = "SampleID", values_to = "rclr_abundance") %>%
    # group_by(SampleID) %>%
    # ungroup() %>%
    # pivot_wider(names_from = "SampleID", values_from = "rclr_abundance") %>%
    # column_to_rownames("asv") %>% 
    # as.matrix()

# normalize per-sample so each column sums to 1 (back-transform rclr to positive first)
# otu_tbl_16S <- exp(otu_tbl_16S)
# otu_tbl_16S <- sweep(otu_tbl_16S, 2, colSums(otu_tbl_16S), "/") %>%
#     as.data.frame() %>%
#     rownames_to_column("asv") %>%  as_tibble()

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
        # filter(topk_freq > 0.5) %>%
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
    mutate(SampleID = fct_reorder(SampleID, plastic_concentration)) %>%
    mutate(Genus_label = paste0(Genus, " (" , round(mean_importance, 2), ") " , " (", round(topk_freq, 2), ")")) %>%
    mutate(Genus_label = fct_reorder(Genus_label, mean_importance))

p3 <- ggplot(abund_tb_16S_joined, aes(x = SampleID, y = Genus_label, size = rel_abund)) +
geom_point(alpha = 0.8) +
  scale_size(range = c(0, 8), name = "Relative abundance") +
#   scale_fill_viridis_c(name = "RF importance") +
  labs(x = "Sample", y = "Taxon", title = "Top taxa bubble chart (16S)", tag = "B") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      plot.tag = element_text(face = "bold"))

### Bubble plots for top 20 taxa - 18S
rf_taxa_18S <- rf_18S %>%
    pull(asv)

# have to somehow make abundances, here I use rclr scaled to 0-1 but it might not be right.
otu_tbl_18S <- otu_table(phy_18S_9_ws) %>%
    as.data.frame() %>%
    rownames_to_column("asv") %>%
    filter(asv %in% rf_taxa_18S)  %>% 
    as_tibble() %>%
    pivot_longer(-asv, names_to = "SampleID", values_to = "abundance") %>%
    group_by(SampleID) %>%
    mutate(rel_abund = abundance / sum(abundance)) %>%
    ungroup() %>%
    select(asv, SampleID, rel_abund) %>%
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
        # filter(topk_freq > 0.5) %>%
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
    mutate(SampleID = fct_reorder(SampleID, plastic_concentration)) %>%
    mutate(Genus_label = paste0(Genus, " (" , round(mean_importance, 2), ") " , " (", round(topk_freq, 2), ")")) %>%
    mutate(Genus_label = fct_reorder(Genus_label, mean_importance))

p4 <- ggplot(abund_tb_18S_joined, aes(x = SampleID, y = Genus_label, size = rel_abund)) +
geom_point(alpha = 0.8) +
  scale_size(range = c(0, 8), name = "Relative abundance") +
#   scale_fill_viridis_c(name = "RF importance") +
  labs(x = "Sample", y = "Taxon", title = "Top taxa bubble chart (18S)", tag = "B") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      plot.tag = element_text(face = "bold"))

library(patchwork)
p5 <- p1 + p3
ggsave("figures/fig_S3_random_forests_16S.png", p5, width = 10, height = 8, dpi = 300, scale = 0.8)
p6 <- p2 + p4
ggsave("figures/fig_S4_random_forests_18S.png", p6, width = 10, height = 8, dpi = 300, scale = 0.8)