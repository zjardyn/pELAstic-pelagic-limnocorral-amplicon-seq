# ALDEx2 exploratory figures — week-9 WS prev2of9 (genus level).
#
# Outcomes: plastic_level (KW; primary kw.eBH) and log10(particles_total_d20) retention
# (aldex.corr; Spearman designated primary — Pearson/Kendall also in RDS).
# No taxa pass BH <= 0.10 at n = 9; plots show effect size and nominal evidence.
#
#   Rscript R/17_aldex_figures_week9_ws_prev2of9.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
})

theme_set(theme_minimal(base_size = 12))

out_dir <- Sys.getenv("ALDEX_FIG_DIR", unset = "")
if (!nzchar(out_dir)) {
  tax_default <- toupper(Sys.getenv("ALDEX_TAX_LEVEL", unset = "Genus"))
  out_dir <- if (tax_default %in% c("ASV", "A")) "figures/aldex_asv" else "figures/aldex_genus"
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

tax_level_raw <- trimws(Sys.getenv("ALDEX_TAX_LEVEL", unset = "Genus"))
tag_suffix <- if (toupper(tax_level_raw) %in% c("GENUS", "G")) "_genus" else ""
aldex_rds <- function(marker, arm) {
  file.path(
    Sys.getenv("ALDEX_OUT_DIR", unset = "output"),
    sprintf("aldex_%s_9_ws_prev2of9%s_%s.rds", marker, tag_suffix, arm)
  )
}

top_n_bars <- as.integer(Sys.getenv("ALDEX_FIG_TOP_N", unset = "20"))
label_n_volcano <- as.integer(Sys.getenv("ALDEX_FIG_LABEL_N", unset = "8"))

arm_files <- tibble::tribble(
  ~marker, ~outcome, ~arm,
  "16S", "Plastic level", "plastic_level_kw",
  "16S", "MP retention", "retention_log10_corr",
  "18S", "Plastic level", "plastic_level_kw",
  "18S", "MP retention", "retention_log10_corr"
) %>%
  mutate(path = purrr::map2_chr(marker, arm, aldex_rds))

taxon_label <- function(genus, family, asv) {
  g <- as.character(genus)
  f <- as.character(family)
  out <- ifelse(!is.na(g) & nzchar(g), g, NA_character_)
  out <- ifelse(is.na(out) & !is.na(f) & nzchar(f), f, out)
  ifelse(is.na(out) | !nzchar(out), substr(as.character(asv), 1L, 8L), out)
}

load_arm <- function(path, marker, outcome) {
  if (!file.exists(path)) {
    stop("Missing ALDEx RDS: ", path, " — run R/15_aldex_week9_ws_prev2of9.R first")
  }
  x <- readRDS(path)
  ct <- x$correlation_table
  if (!is.null(ct) && nrow(ct) > 0 &&
      all(c("spearman.erho", "spearman.ep", "spearman.eBH") %in% names(ct))) {
    ct %>%
      mutate(
        marker = marker,
        outcome = outcome,
        arm_type = "corr",
        taxon = taxon_label(Genus, Family, asv),
        rho = spearman.erho,
        p = spearman.ep,
        eBH = spearman.eBH,
        neglog10p = -log10(pmax(p, .Machine$double.xmin)),
        direction = if_else(rho >= 0, "Positive", "Negative"),
        panel = paste(marker, outcome, sep = " | ")
      ) %>%
      select(marker, outcome, arm_type, panel, asv, taxon, rho, p, eBH, neglog10p, direction)
  } else {
    kt <- x$kw_table
    stopifnot(all(c("kw.ep", "kw.eBH") %in% names(kt)))
    kt %>%
      mutate(
        marker = marker,
        outcome = outcome,
        arm_type = "kw",
        taxon = taxon_label(Genus, Family, asv),
        rho = NA_real_,
        p = kw.ep,
        eBH = kw.eBH,
        neglog10p = -log10(pmax(p, .Machine$double.xmin)),
        direction = "KW",
        panel = paste(marker, outcome, sep = " | ")
      ) %>%
      select(marker, outcome, arm_type, panel, asv, taxon, rho, p, eBH, neglog10p, direction)
  }
}

save_plot <- function(plot, stem, width, height) {
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  png_path <- file.path(out_dir, paste0(stem, ".png"))
  ggsave(pdf_path, plot, width = width, height = height, useDingbats = FALSE)
  ggsave(png_path, plot, width = width, height = height, dpi = 200)
  message("Wrote ", pdf_path, " and ", png_path)
}

df <- pmap_dfr(arm_files[c("path", "marker", "outcome")], load_arm)
df_corr <- df %>% filter(arm_type == "corr")
df_kw <- df %>% filter(arm_type == "kw")

message(
  "Loaded genus ALDEx arms: ",
  paste(unique(df$panel), collapse = "; "),
  " | n_taxa=", nrow(df),
  " | min eBH=", signif(min(df$eBH, na.rm = TRUE), 3)
)

ymax <- max(df_corr$neglog10p, na.rm = TRUE)
ymax <- max(ceiling(ymax * 10) / 10, 1.5)
rho_lim <- range(df_corr$rho, na.rm = TRUE)

panel_tag <- function(marker, outcome) {
  dplyr::case_when(
    marker == "16S" & outcome == "Plastic level" ~ "A",
    marker == "16S" & outcome == "MP retention" ~ "B",
    marker == "18S" & outcome == "Plastic level" ~ "C",
    marker == "18S" & outcome == "MP retention" ~ "D",
    TRUE ~ ""
  )
}

make_volcano <- function(d, tag = NULL) {
  lab_df <- d %>%
    arrange(p, desc(abs(rho))) %>%
    slice_head(n = label_n_volcano)

  p <- ggplot(d, aes(x = rho, y = neglog10p, color = direction)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey55", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_point(alpha = 0.75, size = 2) +
    ggrepel::geom_text_repel(
      data = lab_df,
      aes(label = taxon),
      size = 2.8,
      max.overlaps = 25,
      min.segment.length = 0,
      box.padding = 0.3,
      show.legend = FALSE
    ) +
    scale_color_manual(values = c(Positive = "#B2182B", Negative = "#2166AC")) +
    coord_cartesian(ylim = c(0, ymax), xlim = rho_lim) +
    labs(
      x = expression("Spearman" ~ rho),
      y = expression(-log[10](p)),
      color = NULL,
      title = d$panel[[1]],
      subtitle = expression(log[10](particles[total]~d20))
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  if (!is.null(tag)) {
    p <- p + labs(tag = tag) +
      theme(plot.tag = element_text(face = "bold", size = 13))
  }
  p
}

make_top_bars <- function(d, tag = NULL) {
  top <- d %>%
    mutate(abs_rho = abs(rho)) %>%
    arrange(desc(abs_rho)) %>%
    slice_head(n = top_n_bars) %>%
    mutate(taxon = factor(taxon, levels = rev(unique(taxon))))

  p <- ggplot(top, aes(x = abs_rho, y = taxon, fill = direction)) +
    geom_col(width = 0.75) +
    geom_text(
      aes(label = sprintf("rho=%0.2f  p=%0.3g", rho, p)),
      hjust = -0.02,
      size = 2.5,
      color = "grey20"
    ) +
    scale_fill_manual(values = c(Positive = "#B2182B", Negative = "#2166AC")) +
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.32)),
      limits = c(0, max(1, max(top$abs_rho, na.rm = TRUE) * 1.05))
    ) +
    labs(
      x = expression("|" * rho * "|"),
      y = NULL,
      fill = NULL,
      title = d$panel[[1]],
      subtitle = expression(log[10](particles[total]~d20))
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40"),
      legend.position = "bottom",
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
  if (!is.null(tag)) {
    p <- p + labs(tag = tag) +
      theme(plot.tag = element_text(face = "bold", size = 13))
  }
  p
}

make_top_kw_bars <- function(d, tag = NULL) {
  top <- d %>%
    arrange(p) %>%
    slice_head(n = top_n_bars) %>%
    mutate(taxon = factor(taxon, levels = rev(unique(taxon))))

  p <- ggplot(top, aes(x = neglog10p, y = taxon)) +
    geom_col(width = 0.75, fill = "#4D4D4D") +
    geom_text(
      aes(label = sprintf("p=%0.3g", p)),
      hjust = -0.02,
      size = 2.5,
      color = "grey20"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.32))) +
    labs(
      x = expression(-log[10](p)),
      y = NULL,
      title = d$panel[[1]],
      subtitle = "plastic_level (KW)"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
  if (!is.null(tag)) {
    p <- p + labs(tag = tag) +
      theme(plot.tag = element_text(face = "bold", size = 13))
  }
  p
}

panel_for <- function(marker, outcome, type = c("kw", "corr", "volcano")) {
  type <- match.arg(type)
  d <- df %>% filter(marker == !!marker, outcome == !!outcome)
  tag <- panel_tag(marker, outcome)
  if (type == "kw") {
    make_top_kw_bars(d, tag = tag)
  } else if (type == "corr") {
    make_top_bars(d, tag = tag)
  } else {
    make_volcano(d, tag = tag)
  }
}

# --- combined 2x2 summary: rows = marker, cols = outcome ---
p_summary <- (
  panel_for("16S", "Plastic level", "kw") | panel_for("16S", "MP retention", "volcano")
) / (
  panel_for("18S", "Plastic level", "kw") | panel_for("18S", "MP retention", "volcano")
) +
  plot_annotation(
    title = "ALDEx2 genus | week-9 WS, prev2of9 (n = 9)",
    subtitle = paste0(
      "Left: plastic_level (KW; kw.eBH). Right: retention Spearman (designated primary) vs log10(particles_total_d20). ",
      "Dashed: nominal p = 0.05. No taxa pass BH <= 0.10."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey30")
    )
  )

save_plot(p_summary, "fig_aldex_genus_4panel_summary", width = 13, height = 12)

# --- individual panel PDFs ---
volcano_panels <- df_corr %>%
  group_split(marker, outcome, .keep = TRUE) %>%
  lapply(make_volcano)

p_volcano <- wrap_plots(volcano_panels, ncol = 2, guides = "collect") +
  plot_annotation(
    title = "ALDEx2 genus Spearman (primary) vs log10 MP retention (week-9 WS, prev2of9)",
    subtitle = "No taxa pass BH <= 0.10 (n = 9). Dashed line: nominal p = 0.05.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey30")
    )
  ) &
  theme(legend.position = "bottom")

save_plot(p_volcano, "fig_aldex_genus_volcano_mp_retention", width = 11, height = 9)

bar_panels <- df_corr %>%
  group_split(marker, outcome, .keep = TRUE) %>%
  lapply(make_top_bars)

p_bars <- wrap_plots(bar_panels, ncol = 2, guides = "collect") +
  plot_annotation(
    title = sprintf("ALDEx2 genus top-%d |Spearman rho| vs MP retention", top_n_bars),
    subtitle = "Exploratory ranking only; none pass BH <= 0.10.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey30")
    )
  ) &
  theme(legend.position = "bottom")

save_plot(p_bars, "fig_aldex_genus_top_abs_rho_mp_retention", width = 12, height = 11)

kw_bar_panels <- df_kw %>%
  group_split(marker, outcome, .keep = TRUE) %>%
  lapply(make_top_kw_bars)

p_kw <- wrap_plots(kw_bar_panels, ncol = 2) +
  plot_annotation(
    title = sprintf("ALDEx2 genus top-%d KW p (plastic level)", top_n_bars),
    subtitle = "Week-9 WS, prev2of9. Exploratory only.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey30")
    )
  )

save_plot(p_kw, "fig_aldex_genus_top_kw_plastic_level", width = 12, height = 11)

# --- bridge: genera in top-N of both arms (per marker) ---
bridge_n <- min(20L, top_n_bars)
bridge_tbl <- df %>%
  group_by(marker) %>%
  group_modify(function(d, key) {
    kw_top <- d %>%
      filter(arm_type == "kw") %>%
      arrange(p) %>%
      slice_head(n = bridge_n) %>%
      pull(taxon)
    ret_top <- d %>%
      filter(arm_type == "corr") %>%
      mutate(abs_rho = abs(rho)) %>%
      arrange(desc(abs_rho)) %>%
      slice_head(n = bridge_n) %>%
      pull(taxon)
    tibble(
      taxon = union(kw_top, ret_top),
      in_plastic_kw_top = taxon %in% kw_top,
      in_retention_top = taxon %in% ret_top
    ) %>%
      mutate(bridge = in_plastic_kw_top & in_retention_top)
  }) %>%
  ungroup()

bridge_counts <- bridge_tbl %>%
  summarise(
    n_kw_top = sum(in_plastic_kw_top),
    n_ret_top = sum(in_retention_top),
    n_both = sum(bridge),
    .by = marker
  )
message("Bridge top-", bridge_n, " overlap: ", paste(
  sprintf("%s both=%d", bridge_counts$marker, bridge_counts$n_both),
  collapse = "; "
))

if (any(bridge_tbl$bridge)) {
  p_bridge <- bridge_tbl %>%
    filter(bridge) %>%
    ggplot(aes(x = marker, y = taxon)) +
    geom_point(size = 3, color = "#2166AC") +
    labs(
      x = NULL,
      y = NULL,
      title = sprintf("ALDEx2 bridge: genera in top-%d plastic KW and top-%d |rho| retention", bridge_n, bridge_n),
      subtitle = "Week-9 WS genus, prev2of9"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey35")
    )
  save_plot(p_bridge, "fig_aldex_genus_bridge_overlap", width = 7, height = max(4, 0.25 * sum(bridge_tbl$bridge) + 2))
} else {
  message("No genera in top-", bridge_n, " of both arms for either marker (bridge plot skipped).")
}

write.csv(
  df %>% arrange(marker, outcome, p),
  file.path(out_dir, "aldex_genus_results_long.csv"),
  row.names = FALSE
)
message("Wrote ", file.path(out_dir, "aldex_genus_results_long.csv"))
message("Done.")
