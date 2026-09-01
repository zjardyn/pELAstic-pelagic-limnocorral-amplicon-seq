# Figure 4: Aitchison NMDS on week-9 wall-strip samples (16S + 18S).
# Same pipeline as Fig S panels C/F (R/13_aitchison_nmds_full.R):
#   genus_otu → aitchison_dist → metaMDS; plot_ord_points; permanova_plastic.
#
# Run: Rscript R/07_aitchison_nmds.R

source("R/01_load_files.R")
source("R/functions.R")

NMDS_SEED_W9 <- 456L
PERMANOVA_SEED_W9 <- 3L

phy_16S <- ms_as_week9(phy_16S)
phy_18S <- ms_as_week9(phy_18S)

phy_16S_w9 <- subset_week9_ws(phy_16S, "16S")
phy_18S_w9 <- subset_week9_ws(phy_18S, "18S")

ord_16S_w9 <- aitchison_ordinations(phy_16S_w9, nmds_seed = NMDS_SEED_W9)
ord_18S_w9 <- aitchison_ordinations(phy_18S_w9, nmds_seed = NMDS_SEED_W9)

pl_16S_w9 <- permanova_plastic(
  ord_16S_w9$dist, ord_16S_w9$phy, seed = PERMANOVA_SEED_W9
)
pl_18S_w9 <- permanova_plastic(
  ord_18S_w9$dist, ord_18S_w9$phy, seed = PERMANOVA_SEED_W9
)

message(sprintf("16S week-9 WS NMDS stress = %.3f", ord_16S_w9$stress))
message(sprintf("18S week-9 WS NMDS stress = %.3f", ord_18S_w9$stress))
format_permanova(pl_16S_w9, "16S Aitchison week-9 WS ~ plastic_level")
format_permanova(pl_18S_w9, "18S Aitchison week-9 WS ~ plastic_level")

p_16S <- plot_ord_points(
  add_retention_label(ord_16S_w9$nmds_scores),
  x = NMDS1,
  y = NMDS2,
  colour = plastic_level,
  title = "16S week-9 wall strip",
  tag = "A",
  fill_name = "Plastic level",
  fill_labels = c("None", "Low", "Medium", "High"),
  sample_label = "retention_label"
)

p_18S <- plot_ord_points(
  add_retention_label(ord_18S_w9$nmds_scores),
  x = NMDS1,
  y = NMDS2,
  colour = plastic_level,
  title = "18S week-9 wall strip",
  tag = "B",
  fill_name = "Plastic level",
  fill_labels = c("None", "Low", "Medium", "High"),
  sample_label = "retention_label"
)

p_comb <- combine_plastic_w9_2panel(
  p_16S,
  p_18S,
  caption_16s = caption_plastic_w9(pl_16S_w9),
  caption_18s = caption_plastic_w9(pl_18S_w9),
  stress_16s = ord_16S_w9$stress,
  stress_18s = ord_18S_w9$stress
)

dir.create("figures", showWarnings = FALSE, recursive = TRUE)
ggsave(
  "figures/fig_4_aitchison_nmds_16S_18S.png",
  p_comb,
  width = 14,
  height = 8,
  dpi = 300,
  bg = "white"
)
save_ord_pdf(
  p_comb,
  "figures/fig_4_aitchison_nmds_16S_18S.pdf",
  width = 14 * 0.85,
  height = 8 * 0.85
)

message("Wrote figures/fig_4_aitchison_nmds_16S_18S.png")
message("Wrote figures/fig_4_aitchison_nmds_16S_18S.pdf")
