source("R/01_load_files.R")
source("R/functions.R")

# Supplementary RF figure: all taxa (same layout as fig_6 top-20)
rf_16S <- readRDS("output/rf_runs_1000_16S.rds")
rf_18S <- readRDS("output/rf_runs_1000_18S.rds")

phy_16S_9_ws <- phy_16S %>%
  subset_samples(Date == 9) %>%
  subset_samples(Location == "WS") %>%
  subset_samples(!is.na(particles_total_d20))

phy_18S_9_ws <- phy_18S %>%
  subset_samples(Date == 9) %>%
  subset_samples(Location == "WS") %>%
  subset_samples(!is.na(particles_total_d20))

prep_importance_df <- function(rf) {
  rf %>%
    arrange(desc(mean_importance)) %>%
    mutate(base = ifelse(is.na(Genus) | Genus == "", asv, Genus)) %>%
    group_by(base) %>%
    mutate(
      base_count = n(),
      base = paste0(base, " ", row_number())
    ) %>%
    ungroup() %>%
    mutate(base = ifelse(base_count == 1, sub(" 1$", "", base), base)) %>%
    mutate(label = paste0(base, " (", round(topk_freq, 2), ")")) %>%
    mutate(
      ymin = mean_importance - sd_importance,
      ymax = mean_importance + sd_importance
    )
}

make_importance_plot <- function(df, y_lab, tag, y_text_size = 5) {
  ggplot(df, aes(x = reorder(label, mean_importance), y = mean_importance)) +
    geom_col(fill = "grey40") +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2, color = "black") +
    coord_flip() +
    labs(x = y_lab, y = "Mean permutation importance ± SD", tag = tag) +
    theme(
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = y_text_size),
      axis.title.y = element_text(size = 14),
      axis.title.x = element_text(size = 14),
      plot.tag = element_text(size = 28, face = "bold")
    )
}

hellinger_long <- function(phy, rf) {
  asv_keep <- rf$asv
  otu_tbl <- otu_table(phy) %>%
    as.data.frame() %>%
    rownames_to_column("asv") %>%
    as_tibble() %>%
    pivot_longer(-asv, names_to = "SampleID", values_to = "abundance") %>%
    group_by(SampleID) %>%
    mutate(rel_abund = sqrt(abundance / sum(abundance))) %>%
    ungroup() %>%
    select(asv, SampleID, rel_abund) %>%
    filter(asv %in% asv_keep)

  meta <- sample_data(phy)
  meta_tb <- as_tibble(meta) %>% bind_cols(SampleID = rownames(meta))

  # Unique y labels (same scheme as importance bars)
  lab_map <- prep_importance_df(rf) %>%
    select(asv, Genus_label = label, mean_importance, topk_freq)

  otu_tbl %>%
    inner_join(rf, by = "asv") %>%
    inner_join(lab_map, by = "asv") %>%
    inner_join(meta_tb, by = "SampleID") %>%
    mutate(Corral = fct_reorder(CorralLetter, plastic_concentration)) %>%
    mutate(
      Genus_label = factor(Genus_label, levels = lab_map$Genus_label[order(lab_map$mean_importance)])
    )
}

df_16S <- prep_importance_df(rf_16S)
df_18S <- prep_importance_df(rf_18S)
n16 <- nrow(df_16S)
n18 <- nrow(df_18S)
message("All RF taxa: 16S=", n16, " 18S=", n18)

p1 <- make_importance_plot(df_16S, "Prokaryotic taxa", "A", y_text_size = 4)
p2 <- make_importance_plot(df_18S, "Eukaryotic taxa", "C", y_text_size = 5)

abund_tb_16S_joined <- hellinger_long(phy_16S_9_ws, rf_16S)
abund_tb_18S_joined <- hellinger_long(phy_18S_9_ws, rf_18S)

max_hellinger_rel_abund <- max(
  max(abund_tb_16S_joined$rel_abund, na.rm = TRUE),
  max(abund_tb_18S_joined$rel_abund, na.rm = TRUE),
  na.rm = TRUE
)

p3 <- ggplot(abund_tb_16S_joined, aes(x = Corral, y = Genus_label, size = rel_abund)) +
  geom_point(alpha = 0.8) +
  scale_size_continuous(
    range = c(0, 5),
    limits = c(0, max_hellinger_rel_abund),
    name = "Hellinger abundance"
  ) +
  labs(x = "Limnocorral ID", y = NULL, tag = "B")

p4 <- ggplot(abund_tb_18S_joined, aes(x = Corral, y = Genus_label, size = rel_abund)) +
  geom_point(alpha = 0.8) +
  scale_size_continuous(
    range = c(0, 5),
    limits = c(0, max_hellinger_rel_abund),
    name = "Hellinger abundance"
  ) +
  labs(x = "Limnocorral ID", y = NULL, tag = "D")

p3 <- p3 + theme(
  axis.text.y = element_blank(),
  axis.title.y = element_blank(),
  axis.text.x = element_blank(),
  axis.title.x = element_text(size = 14),
  plot.tag = element_text(size = 28, face = "bold"),
  legend.text = element_text(size = 11),
  legend.title = element_text(size = 13)
)
p4 <- p4 + theme(
  axis.text.y = element_blank(),
  axis.title.y = element_blank(),
  axis.text.x = element_text(size = 11),
  axis.title.x = element_text(size = 14),
  plot.tag = element_text(size = 28, face = "bold"),
  legend.text = element_text(size = 11),
  legend.title = element_text(size = 13)
)

layout <- "
    AB
    CD
"

p_combined <- (p1 + p3 + p2 + p4) +
  patchwork::plot_layout(
    design = layout,
    heights = c(n16, n18),
    guides = "collect",
    axes = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.justification = "left",
    legend.box.just = "left",
    legend.margin = margin(t = 0, r = 0, b = 0, l = -20, unit = "pt")
  )

outfile <- "figures/fig_S_random_forests_all_taxa_combined.pdf"
pdf(
  outfile,
  width = 11 * 0.8,
  height = max(20, 0.09 * (n16 + n18) + 3),
  family = "Helvetica",
  useDingbats = FALSE
)
print(p_combined)
dev.off()
message("Wrote ", outfile, " (height scales with ", n16 + n18, " taxa)")
