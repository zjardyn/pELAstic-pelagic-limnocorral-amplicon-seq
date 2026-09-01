# ALDEx2 GLM / sample-level figures — week-9 WS prev2of9 (ASV or Genus).
#
# 1) GLM MW + MA plots per plastic_level contrast vs none (aldex.glm.plot)
# 2) Sample-level posterior-median CLR by loading group (all samples, jitter)
# 3) Retention Spearman summary + sample-level CLR vs log10(particles_total_d20)
#
# Output: figures/aldex_glm/{16S,18S}/
#
#   Rscript R/18_aldex_glm_sample_figures_week9_ws_prev2of9.R

suppressPackageStartupMessages({
  library(phyloseq)
  library(ALDEx2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

source("R/01_load_files2.R")

mc_samples <- as.integer(Sys.getenv("ALDEX_MC_SAMPLES", unset = "128"))
min_reads <- as.integer(Sys.getenv("ALDEX_MIN_READS", unset = "3"))
min_samples <- as.integer(Sys.getenv("ALDEX_MIN_SAMPLES", unset = "2"))
tax_level_raw <- trimws(Sys.getenv("ALDEX_TAX_LEVEL", unset = "Genus"))
tax_level <- if (toupper(tax_level_raw) %in% c("GENUS", "G")) "Genus" else "ASV"
tag_suffix <- if (identical(tax_level, "Genus")) "_genus" else ""
unit_label <- tolower(if (identical(tax_level, "Genus")) "genus" else "ASV")
out_root <- Sys.getenv("ALDEX_GLM_FIG_DIR", unset = "")
if (!nzchar(out_root)) {
  out_root <- if (identical(tax_level, "Genus")) "figures/aldex_glm" else "figures/aldex_glm_asv"
}
aldex_out <- Sys.getenv("ALDEX_OUT_DIR", unset = "output")
n_group_genera <- as.integer(Sys.getenv("ALDEX_GLM_TOP_N", unset = "6"))
n_retention_genera <- as.integer(Sys.getenv("ALDEX_RET_TOP_N", unset = "6"))
glm_cutoff <- as.numeric(Sys.getenv("ALDEX_GLM_CUTOFF", unset = "0.05"))

plastic_levels <- c("none", "low", "medium", "high")
plastic_colors <- c(
  none = "#4D4D4D",
  low = "#92C5DE",
  medium = "#F4A582",
  high = "#B2182B"
)

dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

subset_week9_ws <- function(phy) {
  phy %>%
    phyloseq::subset_samples(as.character(Date) == "9") %>%
    phyloseq::subset_samples(as.character(Location) == "WS")
}

aggregate_to_genus <- function(phy) {
  phyloseq::tax_glom(phy, taxrank = "Genus", NArm = FALSE)
}

filter_prevalence_counts <- function(phy, min_reads = 3L, min_samples = 2L) {
  otu <- as(phyloseq::otu_table(phy), "matrix")
  if (!phyloseq::taxa_are_rows(phy)) {
    otu <- t(otu)
  }
  keep <- rowSums(otu >= min_reads) >= min_samples
  otu[keep, , drop = FALSE]
}

tax_annot <- function(phy, taxon_ids) {
  tax <- as.data.frame(phyloseq::tax_table(phy)) %>%
    tibble::rownames_to_column("asv")
  if (!"Genus" %in% names(tax)) {
    tax$Genus <- NA_character_
  }
  tax %>%
    dplyr::filter(asv %in% taxon_ids) %>%
    dplyr::select(asv, Genus, dplyr::any_of("Family"))
}

taxon_label <- function(genus, family, asv) {
  g <- as.character(genus)
  f <- as.character(family)
  out <- ifelse(!is.na(g) & nzchar(g), g, NA_character_)
  out <- ifelse(is.na(out) & !is.na(f) & nzchar(f), f, out)
  ifelse(is.na(out) | !nzchar(out), substr(as.character(asv), 1L, 10L), out)
}

prep_marker <- function(phy, label) {
  phy_sub <- subset_week9_ws(phy)
  phy_work <- if (identical(tax_level, "Genus")) {
    aggregate_to_genus(phy_sub)
  } else {
    phy_sub
  }
  counts <- filter_prevalence_counts(
    phy_work,
    min_reads = min_reads,
    min_samples = min_samples
  )
  sd <- as(phyloseq::sample_data(phy_work), "data.frame")
  sd$sample_id <- rownames(sd)
  sample_ids <- colnames(counts)

  loading_group <- factor(
    as.character(sd[sample_ids, "plastic_level"]),
    levels = plastic_levels
  )
  names(loading_group) <- sample_ids

  retention_raw <- setNames(
    as.numeric(sd[sample_ids, "particles_total_d20"]),
    sample_ids
  )
  if (anyNA(retention_raw)) {
    stop(label, ": NA in particles_total_d20")
  }
  retention_log10 <- log10(retention_raw)

  annot <- tax_annot(phy_work, rownames(counts)) %>%
    mutate(feature = taxon_label(Genus, Family, asv))

  group_model <- stats::model.matrix(~ loading_group)
  rownames(group_model) <- sample_ids

  message(sprintf(
    "%s %s prev2of9 (>= %d reads in >= %d/9): n=%d samples | %d taxa",
    label, unit_label, min_reads, min_samples, length(sample_ids), nrow(counts)
  ))

  list(
    marker = label,
    counts = counts,
    sample_ids = sample_ids,
    loading_group = loading_group,
    retention_raw = retention_raw,
    retention_log10 = retention_log10,
    annot = annot,
    group_model = group_model,
    phy = phy_work
  )
}

posterior_median_clr <- function(x_clr, sample_ids) {
  mc <- ALDEx2::getMonteCarloInstances(x_clr)
  sample_clr <- do.call(
    cbind,
    lapply(mc, function(z) apply(z, 1, stats::median))
  )
  colnames(sample_clr) <- sample_ids
  sample_clr
}

feature_display <- function(taxa_ids, annot) {
  lab <- annot$feature[match(taxa_ids, annot$asv)]
  ifelse(is.na(lab), taxa_ids, lab)
}

save_gg <- function(plot, path_no_ext, width, height) {
  ggsave(paste0(path_no_ext, ".pdf"), plot, width = width, height = height, useDingbats = FALSE)
  ggsave(paste0(path_no_ext, ".png"), plot, width = width, height = height, dpi = 200)
  message("Wrote ", path_no_ext, ".pdf and .png")
}

run_glm_contrast_plots <- function(x_group, group_results, group_effect, contrasts, out_dir, marker) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  for (ct in contrasts) {
    ct_safe <- gsub("[^A-Za-z0-9]+", "_", ct)
    pdf_path <- file.path(out_dir, sprintf("ALDEx2_%s_%s_plots.pdf", marker, ct_safe))
    grDevices::pdf(pdf_path, width = 10, height = 5)
    graphics::par(mfrow = c(1, 2))
    ALDEx2::aldex.glm.plot(
      group_results,
      eff = group_effect,
      contrast = ct,
      type = "MW",
      test = "fdr",
      cutoff.pval = glm_cutoff
    )
    ALDEx2::aldex.glm.plot(
      group_results,
      eff = group_effect,
      contrast = ct,
      type = "MA",
      test = "fdr",
      cutoff.pval = glm_cutoff
    )
    grDevices::dev.off()
    message("Wrote ", pdf_path)
  }
}

top_glm_taxa <- function(group_results, n = 6L) {
  p_cols <- grep(":pval\\.padj$", colnames(group_results), value = TRUE)
  p_cols <- setdiff(p_cols, "(Intercept):pval.padj")
  if (!length(p_cols)) {
    stop("No contrast pval.padj columns in aldex.glm output")
  }
  min_p <- apply(group_results[, p_cols, drop = FALSE], 1, min, na.rm = TRUE)
  head(names(sort(min_p)), n)
}

plot_group_clr_samples <- function(
    sample_clr,
    taxa_ids,
    loading_group,
    sample_ids,
    feature_labels,
    marker
) {
  group_plot_data <- as.data.frame(t(sample_clr[taxa_ids, , drop = FALSE])) %>%
    tibble::rownames_to_column("sample") %>%
    pivot_longer(-sample, names_to = "asv", values_to = "clr") %>%
    mutate(
      genus = feature_labels[asv],
      group = loading_group[sample]
    )

  ggplot(group_plot_data, aes(group, clr, color = group)) +
    geom_jitter(width = 0.08, height = 0, size = 3) +
    stat_summary(
      fun = stats::median,
      geom = "crossbar",
      width = 0.45,
      color = "black",
      linewidth = 0.4
    ) +
    facet_wrap(~ genus, scales = "free_y", ncol = 3) +
    scale_color_manual(values = plastic_colors, drop = FALSE) +
    labs(
      x = "Loading group (plastic_level)",
      y = "Posterior median CLR abundance",
      title = sprintf("%s | ALDEx2 %s CLR by plastic level", marker, unit_label),
      subtitle = "All samples shown (jitter); black bar = group median"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey35"),
      strip.text = element_text(face = "bold")
    )
}

plot_retention_corr <- function(corr_df, marker) {
  corr_plot <- corr_df %>%
    mutate(
      neglog10_BH = -log10(pmax(spearman.eBH, .Machine$double.xmin)),
      sig = spearman.eBH < 0.05
    )

  label_data <- corr_plot %>%
    arrange(spearman.eBH, desc(abs(spearman.erho))) %>%
    slice_head(n = 10)

  ggplot(corr_plot, aes(spearman.erho, neglog10_BH)) +
    geom_point(aes(color = sig), alpha = 0.7, size = 2) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
    geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey60") +
    ggrepel::geom_text_repel(
      data = label_data,
      aes(label = feature),
      size = 3,
      max.overlaps = 15,
      min.segment.length = 0
    ) +
    scale_color_manual(
      values = c(`FALSE` = "grey60", `TRUE` = "red"),
      labels = c("FALSE" = "BH >= 0.05", "TRUE" = "BH < 0.05"),
      name = NULL
    ) +
    labs(
      x = "Expected Spearman correlation",
      y = expression(-log[10]("BH-adjusted p")),
      title = sprintf("%s | ALDEx2 %s vs MP retention", marker, unit_label),
      subtitle = "particles_total_d20 (week-9 WS, prev2of9)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey35")
    )
}

plot_retention_clr_samples <- function(
    sample_clr,
    taxa_ids,
    loading_group,
    retention_log10,
    sample_ids,
    feature_labels,
    marker
) {
  retention_plot_data <- as.data.frame(t(sample_clr[taxa_ids, , drop = FALSE])) %>%
    tibble::rownames_to_column("sample") %>%
    pivot_longer(-sample, names_to = "asv", values_to = "clr") %>%
    mutate(
      genus = feature_labels[asv],
      retention_log10 = retention_log10[sample],
      group = loading_group[sample]
    )

  ggplot(retention_plot_data, aes(retention_log10, clr, color = group, label = sample)) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(
      aes(label = sample),
      size = 2.5,
      max.overlaps = 20,
      min.segment.length = 0,
      box.padding = 0.2,
      show.legend = FALSE
    ) +
    facet_wrap(~ genus, scales = "free_y", ncol = 3) +
    scale_color_manual(values = plastic_colors, name = "Loading group") +
    labs(
      x = expression(log[10]("retained microplastic density")),
      y = "Posterior median CLR abundance",
      title = sprintf("%s | %s CLR vs MP retention", marker, unit_label),
      subtitle = sprintf(
        "Each point is one sample; top %s by |Spearman rho|",
        unit_label
      )
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey35"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

run_marker <- function(phy, label) {
  prep <- prep_marker(phy, label)
  out_dir <- file.path(out_root, label)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  set.seed(2026L)
  x_group <- ALDEx2::aldex.clr(
    prep$counts,
    conds = prep$group_model,
    mc.samples = mc_samples,
    denom = "all",
    gamma = NULL,
    verbose = FALSE
  )
  group_results <- ALDEx2::aldex.glm(x_group, verbose = FALSE)
  group_effect <- ALDEx2::aldex.glm.effect(
    x_group,
    CI = TRUE,
    useMC = TRUE,
    verbose = FALSE
  )

  contrasts <- colnames(prep$group_model)[-1]
  run_glm_contrast_plots(
    x_group,
    group_results,
    group_effect,
    contrasts,
    file.path(out_dir, "glm_contrasts"),
    label
  )

  sample_clr <- posterior_median_clr(x_group, prep$sample_ids)
  rownames(sample_clr) <- rownames(prep$counts)
  feature_labels <- setNames(prep$annot$feature, prep$annot$asv)

  glm_top <- top_glm_taxa(group_results, n = n_group_genera)
  p_group <- plot_group_clr_samples(
    sample_clr,
    glm_top,
    prep$loading_group,
    prep$sample_ids,
    feature_labels,
    label
  )
  save_gg(
    p_group,
    file.path(out_dir, "sample_clr_by_plastic_level_top_glm"),
    width = 11,
    height = max(7, 2.2 * ceiling(n_group_genera / 3))
  )

  ret_path <- file.path(
    aldex_out,
    sprintf("aldex_%s_9_ws_prev2of9%s_mp_retention_corr_kw.rds", label, tag_suffix)
  )
  if (!file.exists(ret_path)) {
    warning("Missing retention RDS: ", ret_path, " — run R/15 first")
    return(invisible(NULL))
  }

  ret_obj <- readRDS(ret_path)
  corr_raw <- ret_obj$raw$correlation
  corr_df <- corr_raw %>%
    tibble::rownames_to_column("asv") %>%
    left_join(prep$annot, by = "asv") %>%
    mutate(feature = taxon_label(Genus, Family, asv))

  p_corr <- plot_retention_corr(corr_df, label)
  save_gg(
    p_corr,
    file.path(out_dir, "retention_spearman_summary"),
    width = 8,
    height = 6
  )

  ret_top <- corr_df %>%
    arrange(spearman.eBH, desc(abs(spearman.erho))) %>%
    slice_head(n = n_retention_genera) %>%
    pull(asv)

  p_ret <- plot_retention_clr_samples(
    sample_clr,
    ret_top,
    prep$loading_group,
    prep$retention_log10,
    prep$sample_ids,
    feature_labels,
    label
  )
  save_gg(
    p_ret,
    file.path(out_dir, "sample_clr_vs_retention_top_spearman"),
    width = 11,
    height = max(7, 2.2 * ceiling(n_retention_genera / 3))
  )

  invisible(list(
    group_results = group_results,
    group_effect = group_effect,
    sample_clr = sample_clr
  ))
}

message(sprintf(
  "ALDEx2 GLM/sample figures | tax=%s | mc.samples=%d | filter=>=%d in >=%d/9 | out=%s",
  tax_level, mc_samples, min_reads, min_samples, out_root
))

run_marker(phy_16S, "16S")
run_marker(phy_18S, "18S")

message("Done.")
