# Aitchison NMDS + PCoA on the full dataset (same loading / taxonomy as full-dataset CLR PCA)
# Figure: 2x2 — left time colour, right plastic colour; rows 16S / 18S
# Captions: left = time+site PERMANOVA; right = plastic PERMANOVA; stress far right (NMDS only)
source("R/01_load_files.R")
source("R/functions.R")

# MS samples are week-9 collection; treat as Date = 9 (not the load-file Date = 10 code)
ms_as_week9 <- function(phy) {
  phy %>%
    ps_mutate(Date = ifelse(as.character(Location) == "MS", 9, as.numeric(as.character(Date)))) %>%
    ps_mutate(Date = factor(Date, levels = c(3, 6, 9)))
}
phy_16S <- ms_as_week9(phy_16S)
phy_18S <- ms_as_week9(phy_18S)
message("MS samples coded as Date = 9 (same time as week-9 WS)")

aitchison_ordinations <- function(phy, nmds_seed = 123L) {
  g <- genus_otu(phy, tax_level = "Genus")
  d <- aitchison_dist(g$otu)
  nmds <- nmds_from_dist(d, seed = nmds_seed)
  pcoa <- pcoa_from_dist(d, k = 2L)
  list(
    dist = d,
    phy = g$phy,
    nmds_scores = ord_scores_df(nmds$points, g$phy, c("NMDS1", "NMDS2")),
    stress = nmds$stress,
    pcoa_scores = ord_scores_df(pcoa$points, g$phy, c("PCoA1", "PCoA2")),
    pcoa_var = pcoa$var_explained
  )
}

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
subset_week9_ws <- function(phy, label) {
  keep <- as.character(sample_data(phy)$Location) == "WS" &
    as.character(sample_data(phy)$Date) == "9"
  n_keep <- sum(keep)
  message(label, " week-9 WS: keeping ", n_keep, " / ", nsamples(phy), " samples")
  if (n_keep < 3L) {
    stop(label, ": fewer than 3 week-9 wall-strip samples")
  }
  phy <- prune_samples(keep, phy)
  prune_taxa(taxa_sums(phy) > 0, phy)
}

phy_16S_w9 <- subset_week9_ws(phy_16S, "16S")
phy_18S_w9 <- subset_week9_ws(phy_18S, "18S")
ord_16S_w9 <- aitchison_ordinations(phy_16S_w9, nmds_seed = 456L)
ord_18S_w9 <- aitchison_ordinations(phy_18S_w9, nmds_seed = 456L)
pl_16S_w9 <- permanova_plastic(ord_16S_w9$dist, ord_16S_w9$phy, seed = 3L)
pl_18S_w9 <- permanova_plastic(ord_18S_w9$dist, ord_18S_w9$phy, seed = 3L)

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
    ord_16S_w9$nmds_scores, x = NMDS1, y = NMDS2, colour = plastic_level,
    title = "16S week-9 wall strip", tag = "C",
    fill_name = "Plastic level",
    fill_labels = c("None", "Low", "Medium", "High"),
    sample_label = "CorralLetter"
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
    ord_18S_w9$nmds_scores, x = NMDS1, y = NMDS2, colour = plastic_level,
    title = "18S week-9 wall strip", tag = "F",
    fill_name = "Plastic level",
    fill_labels = c("None", "Low", "Medium", "High"),
    sample_label = "CorralLetter"
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

message(
  "Wrote figures/fig_S_aitchison_nmds_full.pdf (2x3 with week-9 WS) and ",
  "fig_S_aitchison_pcoa_full.pdf (2x2)"
)
