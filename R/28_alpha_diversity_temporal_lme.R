# Alpha diversity LME — ASV-level rarefied metrics.
#
# A) WS temporal: metric ~ Week + (1|Corral), varIdent(~1|Week)
# B) Endpoint: metric ~ Endpoint + (1|Corral), varIdent(~1|Endpoint)
# C) WS plastic: metric ~ Week * log10_plastic + (1|Corral), varIdent(~1|Week)
# D) Endpoint plastic: metric ~ Endpoint * log10_plastic + (1|Corral), varIdent(~1|Endpoint)
# E) Week-9 WS retention: lm(metric ~ log10_retention), n = 9 corrals
#
# BH adjustment across the three metrics within each marker per analysis arm.
# Berger-Parker: model on logit scale; plots show raw Berger-Parker.
#
#   Rscript R/28_alpha_diversity_temporal_lme.R

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(nlme)
  library(emmeans)
})

if (!requireNamespace("ggsignif", quietly = TRUE)) {
  install.packages("ggsignif", repos = "https://cloud.r-project.org", quiet = TRUE)
}
suppressPackageStartupMessages(library(ggsignif))

source("R/01_load_files2.R")
source("R/alpha_diversity_functions.R")

min_sum_taxa <- as.integer(Sys.getenv("ALPHA_MIN_DEPTH", unset = "5000"))
rarefy_depth <- min_sum_taxa
out_main <- Sys.getenv("ALPHA_FIG_DIR", unset = "figures/alpha_diversity_lme")
out_plastic <- Sys.getenv("ALPHA_PLASTIC_FIG_DIR", unset = "figures/alpha_diversity_plastic_lme")
dir.create(out_main, showWarnings = FALSE, recursive = TRUE)
dir.create(out_plastic, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_main, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_plastic, "tables"), showWarnings = FALSE, recursive = TRUE)

theme_set(theme_bw(base_size = 12))

fit_lme <- function(form, data, var_group = NULL) {
  args <- list(
    fixed = form,
    random = ~ 1 | Corral,
    data = data,
    method = "REML",
    na.action = na.omit
  )
  if (!is.null(var_group) && var_group %in% names(data)) {
    args$weights <- nlme::varIdent(form = stats::as.formula(paste("~ 1 |", var_group)))
  }
  do.call(nlme::lme, args)
}

fit_ws_temporal <- function(data, response, label) {
  df <- data %>%
    filter(Location == "WS", Week %in% c("3", "6", "9")) %>%
    mutate(Week = factor(Week, levels = c("3", "6", "9")))
  fit <- fit_lme(as.formula(paste(response, "~ Week")), df, "Week")
  list(label = label, response = response, fit = fit, anova = anova(fit),
       emmeans = emmeans(fit, ~ Week),
       pairwise = contrast(emmeans(fit, ~ Week), "pairwise", adjust = "BH"),
       data = df, term = "Week")
}

fit_endpoint <- function(data, response, label) {
  df <- data %>%
    filter((Location == "WS" & Week == "9") | (Location == "MS" & Week == "10"))
  fit <- fit_lme(as.formula(paste(response, "~ Endpoint")), df, "Endpoint")
  list(label = label, response = response, fit = fit, anova = anova(fit),
       emmeans = emmeans(fit, ~ Endpoint),
       pairwise = contrast(emmeans(fit, ~ Endpoint), "pairwise", adjust = "BH"),
       data = df, term = "Endpoint")
}

fit_ws_plastic <- function(data, response, label) {
  df <- data %>%
    filter(Location == "WS", Week %in% c("3", "6", "9")) %>%
    mutate(Week = factor(Week, levels = c("3", "6", "9")))
  fit <- fit_lme(as.formula(paste(response, "~ Week * log10_plastic")), df, "Week")
  slopes <- emtrends(fit, ~ Week, var = "log10_plastic")
  list(label = label, response = response, fit = fit, anova = anova(fit),
       slopes = slopes, slope_test = test(slopes, adjust = "BH"),
       slope_pairs = pairs(slopes, adjust = "BH"), data = df, term = "Week:log10_plastic")
}

fit_endpoint_plastic <- function(data, response, label) {
  df <- data %>%
    filter((Location == "WS" & Week == "9") | (Location == "MS" & Week == "10"))
  fit <- fit_lme(as.formula(paste(response, "~ Endpoint * log10_plastic")), df, "Endpoint")
  slopes <- emtrends(fit, ~ Endpoint, var = "log10_plastic")
  list(label = label, response = response, fit = fit, anova = anova(fit),
       slopes = slopes, slope_test = test(slopes, adjust = "BH"),
       slope_pairs = pairs(slopes, adjust = "BH"), data = df, term = "Endpoint:log10_plastic")
}

fit_retention_ws9 <- function(data, response, label) {
  df <- data %>%
    filter(Location == "WS", Week == "9", !is.na(particles_total_d20)) %>%
    mutate(log10_retention = log10(particles_total_d20))
  fit <- stats::lm(as.formula(paste(response, "~ log10_retention")), data = df)
  list(label = label, response = response, fit = fit, anova = anova(fit),
       data = df, term = "log10_retention")
}

plastic_level_scale <- function(name = "Plastic level") {
  ggplot2::scale_color_viridis_d(
    name = name,
    labels = c("None", "Low", "Medium", "High"),
    guide = ggplot2::guide_legend(override.aes = list(size = 3))
  )
}

plot_week_temporal <- function(res, m, p_bh = NA, col_title = NULL) {
  df <- res$data
  cap <- format_p_line("Week effect (BH)", res$p_omnibus, p_bh)
  if (m$logit_model) cap <- paste0(cap, "\n(model: logit Berger-Parker)")

  p <- ggplot(df, aes(Week, .data[[m$plot_col]], color = Week)) +
    geom_jitter(width = 0.12, size = 2.5, alpha = 0.85) +
    stat_summary(fun = mean, geom = "point", size = 4, color = "black") +
    stat_summary(fun = mean, geom = "line", aes(group = 1), color = "black", linewidth = 0.5) +
    labs(subtitle = col_title, x = "Week", y = m$ylab, caption = cap) +
    theme(legend.position = "none", plot.caption = element_text(hjust = 0, size = 8, lineheight = 1.1))

  add_pairwise_signif(p, df, res$pairwise, m$plot_col)
}

plot_endpoint <- function(res, m, p_bh = NA, col_title = NULL) {
  df <- res$data
  cap <- format_p_line("Endpoint effect (BH)", res$p_omnibus, p_bh)
  if (m$logit_model) cap <- paste0(cap, "\n(model: logit Berger-Parker)")

  p <- ggplot(df, aes(Endpoint, .data[[m$plot_col]], color = Endpoint)) +
    geom_line(aes(group = Corral), color = "grey70", linewidth = 0.4) +
    geom_point(aes(fill = Endpoint), size = 3, shape = 21, color = "black", stroke = 0.3) +
    labs(subtitle = col_title, x = NULL, y = m$ylab, caption = cap) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 18, hjust = 1),
          plot.caption = element_text(hjust = 0, size = 8, lineheight = 1.1))

  add_pairwise_signif(p, df, res$pairwise, m$plot_col, step_frac = 0.12)
}

plot_plastic_facets <- function(res, m, facet_var, col_title = NULL, p_int_bh = NA, show_legend = FALSE) {
  df <- res$data
  p_int <- res$p_omnibus
  st <- as.data.frame(res$slope_test)
  pcol <- intersect(c("p.value", "P.value"), names(st))[1]
  slope_txt <- paste(
    paste0(st[[facet_var]], ": p = ", format.pval(st[[pcol]], digits = 3, eps = 0.001)),
    collapse = "\n"
  )
  cap <- paste0(
    format_p_line("Interaction", p_int, p_int_bh), "\n",
    "Within-group plastic slopes (BH):\n", slope_txt
  )
  if (m$logit_model) cap <- paste0(cap, "\n(model: logit Berger-Parker)")

  p <- ggplot(df, aes(log10_plastic, .data[[m$plot_col]], color = plastic_level)) +
    geom_point(size = 2.5, alpha = 0.9) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5,
                inherit.aes = FALSE,
                aes(x = log10_plastic, y = .data[[m$plot_col]])) +
    plastic_level_scale() +
    facet_wrap(stats::as.formula(paste("~", facet_var)), nrow = 1, scales = "free_y") +
    labs(
      subtitle = col_title,
      x = expression(log[10](nominal~concentration~+~1)), y = m$ylab, caption = cap
    ) +
    theme(
      plot.caption = element_text(hjust = 0, size = 7.5, lineheight = 1.05),
      strip.text = element_text(size = 8),
      legend.position = if (show_legend) "bottom" else "none"
    )
  p
}

plot_retention <- function(res, m, p_bh = NA, show_legend = FALSE) {
  df <- res$data
  p <- res$anova["log10_retention", "Pr(>F)"]
  cap <- format_p_line("log10(retention) slope", p, p_bh)
  if (m$logit_model) cap <- paste0(cap, "\n(model: logit Berger-Parker)")

  ggplot(df, aes(log10_retention, .data[[m$plot_col]], color = plastic_level, label = Corral)) +
    geom_point(size = 3) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.6) +
    plastic_level_scale() +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 12, show.legend = FALSE) +
    labs(subtitle = "Week-9 wall strip | n = 9 corrals",
         x = expression(log[10](particles[total]~d20)), y = m$ylab, caption = cap) +
    theme(
      plot.caption = element_text(hjust = 0, size = 8, lineheight = 1.1),
      legend.position = if (show_legend) "bottom" else "none"
    )
}

metric_row_grid <- function(metrics, left_fn, right_fn, col_titles = c("Weeks", "Location")) {
  rows <- Map(function(m, i) {
    col_title <- if (i == 1L) col_titles else c(NULL, NULL)
    left_fn(m, col_title[1]) + right_fn(m, col_title[2])
  }, metrics, seq_along(metrics))
  wrap_plots(rows, ncol = 1)
}

run_marker <- function(phy, marker) {
  metrics <- alpha_metrics()
  alpha <- prepare_alpha_data(phy, marker, min_sum_taxa, rarefy_depth)

  fit_all <- function(fit_fn, term) {
    res <- setNames(lapply(metrics, function(m) fit_fn(alpha, m$response, m$label)),
                    vapply(metrics, `[[`, "", "id"))
    attach_bh(res, term)
  }

  ws <- fit_all(fit_ws_temporal, "Week")
  ep <- fit_all(fit_endpoint, "Endpoint")
  ws_pl <- fit_all(fit_ws_plastic, "Week:log10_plastic")
  ep_pl <- fit_all(fit_endpoint_plastic, "Endpoint:log10_plastic")
  ret <- setNames(lapply(metrics, function(m) fit_retention_ws9(alpha, m$response, m$label)),
                  vapply(metrics, `[[`, "", "id"))
  ret <- attach_bh(ret, "log10_retention")

  p_main <- metric_row_grid(
    metrics,
    left_fn = function(m, col_title) {
      plot_week_temporal(ws[[m$id]], m, ws[[m$id]]$p_omnibus_bh, col_title)
    },
    right_fn = function(m, col_title) {
      plot_endpoint(ep[[m$id]], m, ep[[m$id]]$p_omnibus_bh, col_title)
    }
  ) +
    plot_annotation(
      title = sprintf("%s ASV alpha diversity (rarefied to %d reads)", marker, rarefy_depth),
      theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14))
    )

  ggsave(file.path(out_main, sprintf("fig_alpha_diversity_lme_%s.pdf", marker)),
         p_main, width = 11, height = 12, useDingbats = FALSE)
  ggsave(file.path(out_main, sprintf("fig_alpha_diversity_lme_%s.png", marker)),
         p_main, width = 11, height = 12, dpi = 200)
  message("Wrote ", file.path(out_main, sprintf("fig_alpha_diversity_lme_%s.pdf", marker)))

  p_plastic_rows <- metric_row_grid(
    metrics,
    col_titles = c("Weeks", "Location"),
    left_fn = function(m, col_title) {
      plot_plastic_facets(ws_pl[[m$id]], m, "Week", col_title, ws_pl[[m$id]]$p_omnibus_bh)
    },
    right_fn = function(m, col_title) {
      plot_plastic_facets(
        ep_pl[[m$id]], m, "Endpoint", col_title, ep_pl[[m$id]]$p_omnibus_bh,
        show_legend = identical(m$id, "Shannon")
      )
    }
  )

  p_ret <- wrap_plots(Map(function(m) {
    plot_retention(ret[[m$id]], m, ret[[m$id]]$p_omnibus_bh, show_legend = FALSE)
  }, metrics), ncol = 3)

  p_plastic <- (p_plastic_rows / p_ret) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  p_plastic <- p_plastic +
    plot_annotation(
      title = sprintf("%s ASV alpha diversity — plastic loading models", marker),
      theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14))
    )

  ggsave(file.path(out_plastic, sprintf("fig_alpha_plastic_lme_%s.pdf", marker)),
         p_plastic, width = 15, height = 13, useDingbats = FALSE)
  ggsave(file.path(out_plastic, sprintf("fig_alpha_plastic_lme_%s.png", marker)),
         p_plastic, width = 15, height = 13, dpi = 200)
  message("Wrote ", file.path(out_plastic, sprintf("fig_alpha_plastic_lme_%s.pdf", marker)))

  results <- list(marker = marker, alpha_data = alpha,
                  ws_temporal = ws, endpoint = ep,
                  ws_plastic = ws_pl, endpoint_plastic = ep_pl, retention_ws9 = ret)
  saveRDS(results, file.path(out_main, "tables", sprintf("alpha_lme_%s.rds", marker)))
  saveRDS(results, file.path(out_plastic, "tables", sprintf("alpha_plastic_lme_%s.rds", marker)))
  write.csv(alpha, file.path(out_main, "tables", sprintf("alpha_data_%s.csv", marker)), row.names = FALSE)
  invisible(results)
}

message("Alpha diversity LME (ASV, rarefied) | depth = ", rarefy_depth)
run_marker(phy_16S, "16S")
run_marker(phy_18S, "18S")
message("Done.")
