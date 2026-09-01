# Aitchison NMDS + PCoA on the full dataset (same loading / taxonomy as full-dataset CLR PCA)
# Figure: 2x2 — left time colour, right plastic colour; rows 16S / 18S
# Captions: left = time+site PERMANOVA; right = plastic PERMANOVA; stress far right (NMDS only)
source("R/01_load_files.R")
source("R/functions.R")

# Week-9 wall-strip panels (C/F) — keep in sync with R/07_aitchison_nmds.R
NMDS_SEED_W9 <- 456L
PERMANOVA_SEED_W9 <- 3L

phy_16S <- ms_as_week9(phy_16S)
phy_18S <- ms_as_week9(phy_18S)
message("MS samples coded as Date = 9 (same time as week-9 WS)")

ord_16S <- aitchison_ordinations(phy_16S)
ord_18S <- aitchison_ordinations(phy_18S)

message(sprintf(
  "16S Aitchison PCoA var = %.1f%% / %.1f%%; NMDS stress = %.3f",
  ord_16S$pcoa_var[1], ord_16S$pcoa_var[2], ord_16S$stress
))
message(sprintf(
  "18S Aitchison PCoA var = %.1f%% / %.1f%%; NMDS stress = %.3f",
  ord_18S$pcoa_var[1], ord_18S$pcoa_var[2], ord_18S$stress
))

# Separate PERMANOVAs for captions (shared across NMDS and PCoA — same distance)
ts_16S <- permanova_time_site(ord_16S$dist, ord_16S$phy)
pl_16S <- permanova_plastic(ord_16S$dist, ord_16S$phy)
ts_18S <- permanova_time_site(ord_18S$dist, ord_18S$phy)
pl_18S <- permanova_plastic(ord_18S$dist, ord_18S$phy)

full_16S <- permanova_time_site_plastic(ord_16S$dist, ord_16S$phy)
full_18S <- permanova_time_site_plastic(ord_18S$dist, ord_18S$phy)

dir.create("results", showWarnings = FALSE)

# --- Week-9 wall-strip subset; Aitchison NMDS + plastic PERMANOVA ---
phy_16S_w9 <- subset_week9_ws(phy_16S, "16S")
phy_18S_w9 <- subset_week9_ws(phy_18S, "18S")
ord_16S_w9 <- aitchison_ordinations(phy_16S_w9, nmds_seed = NMDS_SEED_W9)
ord_18S_w9 <- aitchison_ordinations(phy_18S_w9, nmds_seed = NMDS_SEED_W9)
pl_16S_w9 <- permanova_plastic(ord_16S_w9$dist, ord_16S_w9$phy, seed = PERMANOVA_SEED_W9)
pl_18S_w9 <- permanova_plastic(ord_18S_w9$dist, ord_18S_w9$phy, seed = PERMANOVA_SEED_W9)

message(sprintf("16S week-9 WS Aitchison NMDS stress = %.3f", ord_16S_w9$stress))
message(sprintf("18S week-9 WS Aitchison NMDS stress = %.3f", ord_18S_w9$stress))

write_csv(
  bind_rows(
    format_permanova(ts_16S, "16S Aitchison ~ Location + Date") %>%
      mutate(amplicon = "16S", distance = "aitchison", model = "time_site"),
    format_permanova(pl_16S, "16S Aitchison ~ plastic_level") %>%
      mutate(amplicon = "16S", distance = "aitchison", model = "plastic"),
    format_permanova(ts_18S, "18S Aitchison ~ Location + Date") %>%
      mutate(amplicon = "18S", distance = "aitchison", model = "time_site"),
    format_permanova(pl_18S, "18S Aitchison ~ plastic_level") %>%
      mutate(amplicon = "18S", distance = "aitchison", model = "plastic"),
    format_permanova(full_16S, "16S Aitchison ~ Location + Date + plastic_level") %>%
      mutate(amplicon = "16S", distance = "aitchison", model = "full"),
    format_permanova(full_18S, "18S Aitchison ~ Location + Date + plastic_level") %>%
      mutate(amplicon = "18S", distance = "aitchison", model = "full"),
    format_permanova(pl_16S_w9, "16S Aitchison week-9 WS ~ plastic_level") %>%
      mutate(amplicon = "16S", distance = "aitchison", model = "plastic_week9_ws"),
    format_permanova(pl_18S_w9, "18S Aitchison week-9 WS ~ plastic_level") %>%
      mutate(amplicon = "18S", distance = "aitchison", model = "plastic_week9_ws")
  ),
  "results/permanova_aitchison_full.csv"
)
message("Wrote results/permanova_aitchison_full.csv")

# --- Aitchison NMDS 2x3: time | plastic (full) | plastic (week-9 WS) ---
p_nmds <- combine_time_plastic_w9_6panel(
  plot_ord_points(
    ord_16S$nmds_scores, x = NMDS1, y = NMDS2, colour = Date,
    title = "16S rRNA (Aitchison NMDS)", tag = "A",
    fill_name = "Week"
  ),
  plot_ord_points(
    ord_16S$nmds_scores, x = NMDS1, y = NMDS2, colour = plastic_level,
    title = "16S rRNA (Aitchison NMDS)", tag = "B",
    fill_name = "Plastic level",
    fill_labels = c("None", "Low", "Medium", "High")
  ),
  plot_ord_points(
    add_retention_label(ord_16S_w9$nmds_scores),
    x = NMDS1, y = NMDS2, colour = plastic_level,
    title = "16S week-9 wall strip", tag = "C",
    fill_name = "Plastic level",
    fill_labels = c("None", "Low", "Medium", "High"),
    sample_label = "retention_label"
  ),
  plot_ord_points(
    ord_18S$nmds_scores, x = NMDS1, y = NMDS2, colour = Date,
    title = "18S rRNA (Aitchison NMDS)", tag = "D",
    fill_name = "Week"
  ),
  plot_ord_points(
    ord_18S$nmds_scores, x = NMDS1, y = NMDS2, colour = plastic_level,
    title = "18S rRNA (Aitchison NMDS)", tag = "E",
    fill_name = "Plastic level",
    fill_labels = c("None", "Low", "Medium", "High")
  ),
  plot_ord_points(
    add_retention_label(ord_18S_w9$nmds_scores),
    x = NMDS1, y = NMDS2, colour = plastic_level,
    title = "18S week-9 wall strip", tag = "F",
    fill_name = "Plastic level",
    fill_labels = c("None", "Low", "Medium", "High"),
    sample_label = "retention_label"
  ),
  caption_16s_time = caption_time_site(ts_16S),
  caption_16s_plastic = caption_plastic(pl_16S),
  caption_16s_w9 = caption_plastic_w9(pl_16S_w9),
  caption_18s_time = caption_time_site(ts_18S),
  caption_18s_plastic = caption_plastic(pl_18S),
  caption_18s_w9 = caption_plastic_w9(pl_18S_w9),
  stress_16s = ord_16S$stress,
  stress_16s_w9 = ord_16S_w9$stress,
  stress_18s = ord_18S$stress,
  stress_18s_w9 = ord_18S_w9$stress
)
save_ord_pdf(p_nmds, "figures/fig_S_aitchison_nmds_full.pdf", width = 24 * 0.85, height = 16 * 0.85)

# --- Aitchison PCoA 2x2 ---
p_pcoa <- combine_time_plastic_4panel(
  plot_ord_points(
    ord_16S$pcoa_scores, x = PCoA1, y = PCoA2, colour = Date,
    title = "16S rRNA (Aitchison PCoA)", tag = "A",
    fill_name = "Week",
    xlab = glue("PCoA1 [{round(ord_16S$pcoa_var[1], 1)}%]"),
    ylab = glue("PCoA2 [{round(ord_16S$pcoa_var[2], 1)}%]")
  ),
  plot_ord_points(
    ord_16S$pcoa_scores, x = PCoA1, y = PCoA2, colour = plastic_level,
    title = "16S rRNA (Aitchison PCoA)", tag = "B",
    fill_name = "Plastic level",
    fill_labels = c("None", "Low", "Medium", "High"),
    xlab = glue("PCoA1 [{round(ord_16S$pcoa_var[1], 1)}%]"),
    ylab = glue("PCoA2 [{round(ord_16S$pcoa_var[2], 1)}%]")
  ),
  plot_ord_points(
    ord_18S$pcoa_scores, x = PCoA1, y = PCoA2, colour = Date,
    title = "18S rRNA (Aitchison PCoA)", tag = "C",
    fill_name = "Week",
    xlab = glue("PCoA1 [{round(ord_18S$pcoa_var[1], 1)}%]"),
    ylab = glue("PCoA2 [{round(ord_18S$pcoa_var[2], 1)}%]")
  ),
  plot_ord_points(
    ord_18S$pcoa_scores, x = PCoA1, y = PCoA2, colour = plastic_level,
    title = "18S rRNA (Aitchison PCoA)", tag = "D",
    fill_name = "Plastic level",
    fill_labels = c("None", "Low", "Medium", "High"),
    xlab = glue("PCoA1 [{round(ord_18S$pcoa_var[1], 1)}%]"),
    ylab = glue("PCoA2 [{round(ord_18S$pcoa_var[2], 1)}%]")
  ),
  caption_16s_left = caption_time_site(ts_16S),
  caption_16s_right = caption_plastic(pl_16S),
  caption_18s_left = caption_time_site(ts_18S),
  caption_18s_right = caption_plastic(pl_18S)
)
save_ord_pdf(p_pcoa, "figures/fig_S_aitchison_pcoa_full.pdf", width = 18 * 0.85, height = 18 * 0.85)

# --- Same 2x2 PCoA with CLR-PCA native loadings on plastic panels (B, D) ---
# Top 10 genera by max(r2_PC1, r2_PC2), same r2 definition as fig_S2_pca.
taxa_16S <- pcoa_pca_loadings(
  ord_16S$otu, ord_16S$pcoa_scores, n_top = 10L, arrow_scale = 0.7
)
taxa_18S <- pcoa_pca_loadings(
  ord_18S$otu, ord_18S$pcoa_scores, n_top = 10L, arrow_scale = 0.7
)
message(sprintf(
  "CLR-PCA loadings on PCoA (top 10 by max r2): 16S n=%d (threshold max_r2=%.4f) | 18S n=%d (threshold max_r2=%.4f)",
  nrow(taxa_16S), min(taxa_16S$max_r2),
  nrow(taxa_18S), min(taxa_18S$max_r2)
))
thresh_16S <- min(taxa_16S$max_r2)
thresh_18S <- min(taxa_18S$max_r2)
write_csv(
  bind_rows(
    taxa_16S %>% mutate(amplicon = "16S", n_top = 10L, r2_threshold = thresh_16S),
    taxa_18S %>% mutate(amplicon = "18S", n_top = 10L, r2_threshold = thresh_18S)
  ),
  "results/aitchison_pcoa_pca_loadings_top10.csv"
)

p_pcoa_taxa <- combine_time_plastic_4panel(
  plot_ord_points(
    ord_16S$pcoa_scores, x = PCoA1, y = PCoA2, colour = Date,
    title = "16S rRNA (Aitchison PCoA)", tag = "A",
    fill_name = "Week",
    xlab = glue("PCoA1 [{round(ord_16S$pcoa_var[1], 1)}%]"),
    ylab = glue("PCoA2 [{round(ord_16S$pcoa_var[2], 1)}%]")
  ),
  plot_ord_points(
    ord_16S$pcoa_scores, x = PCoA1, y = PCoA2, colour = plastic_level,
    title = "16S rRNA (Aitchison PCoA)", tag = "B",
    fill_name = "Plastic level",
    fill_labels = c("None", "Low", "Medium", "High"),
    xlab = glue("PCoA1 [{round(ord_16S$pcoa_var[1], 1)}%]"),
    ylab = glue("PCoA2 [{round(ord_16S$pcoa_var[2], 1)}%]"),
    taxa = taxa_16S,
    taxa_size = 2.4
  ),
  plot_ord_points(
    ord_18S$pcoa_scores, x = PCoA1, y = PCoA2, colour = Date,
    title = "18S rRNA (Aitchison PCoA)", tag = "C",
    fill_name = "Week",
    xlab = glue("PCoA1 [{round(ord_18S$pcoa_var[1], 1)}%]"),
    ylab = glue("PCoA2 [{round(ord_18S$pcoa_var[2], 1)}%]")
  ),
  plot_ord_points(
    ord_18S$pcoa_scores, x = PCoA1, y = PCoA2, colour = plastic_level,
    title = "18S rRNA (Aitchison PCoA)", tag = "D",
    fill_name = "Plastic level",
    fill_labels = c("None", "Low", "Medium", "High"),
    xlab = glue("PCoA1 [{round(ord_18S$pcoa_var[1], 1)}%]"),
    ylab = glue("PCoA2 [{round(ord_18S$pcoa_var[2], 1)}%]"),
    taxa = taxa_18S,
    taxa_size = 2.4
  ),
  caption_16s_left = caption_time_site(ts_16S),
  caption_16s_right = paste0(
    caption_plastic(pl_16S),
    sprintf(
      "; arrows: top 10 CLR-PCA loadings by max(r2_PC1, r2_PC2) (threshold = %.3f)",
      thresh_16S
    )
  ),
  caption_18s_left = caption_time_site(ts_18S),
  caption_18s_right = paste0(
    caption_plastic(pl_18S),
    sprintf(
      "; arrows: top 10 CLR-PCA loadings by max(r2_PC1, r2_PC2) (threshold = %.3f)",
      thresh_18S
    )
  )
)
save_ord_pdf(
  p_pcoa_taxa,
  "figures/fig_S_aitchison_pcoa_full_taxa.pdf",
  width = 18 * 0.85,
  height = 18 * 0.85
)

message(
  "Wrote figures/fig_S_aitchison_nmds_full.pdf (2x3 with week-9 WS), ",
  "fig_S_aitchison_pcoa_full.pdf (2x2), and ",
  "fig_S_aitchison_pcoa_full_taxa.pdf (2x2 + top-10 CLR-PCA loadings on B/D)"
)
