# ALDEx2 exploratory figures — week-9 WS prev2of9 (genus level).
#
# No taxa pass spearman/kw eBH ≤ 0.10 (n = 9). Plots show effect size and
# nominal evidence: volcano (ρ vs −log10 p) and top-|ρ| bars.
#
#   Rscript R/17_aldex_figures_week9_ws_prev2of9.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
})

theme_set(theme_minimal(base_size = 12))

out_dir <- "figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

top_n_bars <- 20L
label_n_volcano <- 8L

arm_files <- tibble::tribble(
  ~marker, ~outcome, ~path,
  "16S", "Plastic level", "output/aldex_16S_9_ws_prev2of9_genus_plastic_level_kw.rds",
  "16S", "MP retention", "output/aldex_16S_9_ws_prev2of9_genus_mp_retention_corr_kw.rds",
  "18S", "Plastic level", "output/aldex_18S_9_ws_prev2of9_genus_plastic_level_kw.rds",
  "18S", "MP retention", "output/aldex_18S_9_ws_prev2of9_genus_mp_retention_corr_kw.rds"
)

taxon_label <- function(genus, family, asv) {
  g <- as.character(genus)
  f <- as.character(family)
  out <- ifelse(!is.na(g) & nzchar(g), g, NA_character_)
  out <- ifelse(is.na(out) & !is.na(f) & nzchar(f), f, out)
  ifelse(is.na(out) | !nzchar(out), substr(as.character(asv), 1L, 8L), out)
}

load_arm <- function(path, marker, outcome) {
  x <- readRDS(path)
  ct <- x$correlation_table
  if (!is.null(ct) && nrow(ct) > 0 && all(c("spearman.erho", "spearman.ep", "spearman.eBH") %in% names(ct))) {
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
        panel = paste(marker, outcome, sep = " · ")
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
        panel = paste(marker, outcome, sep = " · ")
      ) %>%
      select(marker, outcome, arm_type, panel, asv, taxon, rho, p, eBH, neglog10p, direction)
  }
}

df <- purrr::pmap_dfr(arm_files, load_arm)
df_corr <- df %>% dplyr::filter(arm_type == "corr")
message(
  "Loaded genus ALDEx arms: ",
  paste(unique(df$panel), collapse = "; "),
  " | n_taxa=", nrow(df),
  " | min eBH=", signif(min(df$eBH, na.rm = TRUE), 3)
)

# Shared volcano y-limit across correlation panels
ymax <- max(df_corr$neglog10p, na.rm = TRUE)
ymax <- max(ceiling(ymax * 10) / 10, 1.5)

make_volcano <- function(d) {
  lab_df <- d %>%
    arrange(p, desc(abs(rho))) %>%
    slice_head(n = label_n_volcano)

  ggplot(d, aes(x = rho, y = neglog10p, color = direction)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey55", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_point(alpha = 0.7, size = 1.8) +
    ggrepel::geom_text_repel(
      data = lab_df,
      aes(label = taxon),
      size = 2.8,
      max.overlaps = 20,
      min.segment.length = 0,
      box.padding = 0.25,
      show.legend = FALSE
    ) +
    scale_color_manual(values = c(Positive = "#B2182B", Negative = "#2166AC")) +
    coord_cartesian(ylim = c(0, ymax), xlim = range(df_corr$rho, na.rm = TRUE)) +
    labs(
      x = expression("Spearman" ~ rho ~ "(ALDEx CLR)"),
      y = expression(-log[10](p)),
      color = NULL,
      title = d$panel[[1]]
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
}

volcano_panels <- df_corr %>%
  group_split(marker, outcome, .keep = TRUE) %>%
  lapply(make_volcano)

p_volcano <- wrap_plots(volcano_panels, ncol = 2, guides = "collect") +
  plot_annotation(
    title = "ALDEx2 genus Spearman vs MP retention (week-9 WS, prev2of9)",
    subtitle = "No taxa pass BH ≤ 0.10 (n = 9). Dashed line: nominal p = 0.05. Top labels: lowest p.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey30")
    )
  ) &
  theme(legend.position = "bottom")

volcano_pdf <- file.path(out_dir, "fig_aldex_genus_volcano_mp_retention.pdf")
ggsave(volcano_pdf, p_volcano, width = 11, height = 9)
message("Wrote ", volcano_pdf)

make_top_bars <- function(d) {
  top <- d %>%
    mutate(abs_rho = abs(rho)) %>%
    arrange(desc(abs_rho)) %>%
    slice_head(n = top_n_bars) %>%
    mutate(taxon = factor(taxon, levels = rev(unique(taxon))))

  ggplot(top, aes(x = abs_rho, y = taxon, fill = direction)) +
    geom_col(width = 0.75) +
    geom_text(
      aes(label = sprintf("ρ=%0.2f  p=%0.3g", rho, p)),
      hjust = -0.05,
      size = 2.6,
      color = "grey20"
    ) +
    scale_fill_manual(values = c(Positive = "#B2182B", Negative = "#2166AC")) +
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.28)),
      limits = c(0, max(1, max(top$abs_rho, na.rm = TRUE)))
    ) +
    labs(
      x = expression("|" * rho * "|"),
      y = NULL,
      fill = NULL,
      title = d$panel[[1]]
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      legend.position = "bottom",
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
}

make_top_kw_bars <- function(d) {
  top <- d %>%
    arrange(p) %>%
    slice_head(n = top_n_bars) %>%
    mutate(taxon = factor(taxon, levels = rev(unique(taxon))))

  ggplot(top, aes(x = neglog10p, y = taxon)) +
    geom_col(width = 0.75, fill = "#4D4D4D") +
    geom_text(
      aes(label = sprintf("p=%0.3g", p)),
      hjust = -0.05,
      size = 2.6,
      color = "grey20"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.28))) +
    labs(
      x = expression(-log[10](p) ~ "(KW)"),
      y = NULL,
      title = d$panel[[1]]
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
}

bar_panels <- df_corr %>%
  group_split(marker, outcome, .keep = TRUE) %>%
  lapply(make_top_bars)

kw_bar_panels <- df %>%
  dplyr::filter(arm_type == "kw") %>%
  group_split(marker, outcome, .keep = TRUE) %>%
  lapply(make_top_kw_bars)

p_bars <- wrap_plots(bar_panels, ncol = 2, guides = "collect") +
  plot_annotation(
    title = sprintf(
      "ALDEx2 genus top-%d |Spearman ρ| (week-9 WS, prev2of9)",
      top_n_bars
    ),
    subtitle = "Exploratory ranking only — none pass BH ≤ 0.10.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey30")
    )
  ) &
  theme(legend.position = "bottom")

bars_pdf <- file.path(out_dir, "fig_aldex_genus_top_abs_rho_mp_retention.pdf")
ggsave(bars_pdf, p_bars, width = 12, height = 11)
message("Wrote ", bars_pdf)

if (length(kw_bar_panels)) {
  p_kw <- wrap_plots(kw_bar_panels, ncol = 2) +
    plot_annotation(
      title = sprintf(
        "ALDEx2 genus top-%d KW p (plastic level; week-9 WS, prev2of9)",
        top_n_bars
      ),
      theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
    )
  kw_pdf <- file.path(out_dir, "fig_aldex_genus_top_kw_plastic_level.pdf")
  ggsave(kw_pdf, p_kw, width = 12, height = 11)
  message("Wrote ", kw_pdf)
}
