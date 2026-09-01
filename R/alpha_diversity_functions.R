# Shared alpha-diversity prep and LME helpers (ASV-level, rarefied).

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
})

alpha_metrics <- function() {
  list(
    list(
      id = "Observed",
      response = "Observed",
      plot_col = "Observed",
      label = "Observed ASV richness",
      ylab = "Observed ASVs",
      logit_model = FALSE
    ),
    list(
      id = "Berger_Parker",
      response = "Berger_Parker_logit",
      plot_col = "Berger_Parker",
      label = "Berger-Parker dominance",
      ylab = "Berger-Parker",
      logit_model = TRUE
    ),
    list(
      id = "Shannon",
      response = "Shannon",
      plot_col = "Shannon",
      label = "Shannon diversity",
      ylab = "Shannon",
      logit_model = FALSE
    )
  )
}

collapse_technical_replicates <- function(phy) {
  otu <- as(otu_table(phy), "matrix")
  if (!taxa_are_rows(phy)) otu <- t(otu)
  for (pair in list(c("WS_E_6", "WS_E_6-2"), c("WS_I_3", "WS_I_3-2"))) {
    if (all(pair %in% colnames(otu))) {
      otu[, pair[1]] <- otu[, pair[1]] + otu[, pair[2]]
      otu <- otu[, setdiff(colnames(otu), pair[2]), drop = FALSE]
    }
  }
  keep <- intersect(sample_names(phy), colnames(otu))
  phy <- prune_samples(keep, phy)
  otu_table(phy) <- otu_table(otu[, sample_names(phy), drop = FALSE], taxa_are_rows = TRUE)
  phy
}

prune_samples_min_depth <- function(phy, label, min_reads) {
  if (min_reads <= 0L) return(phy)
  ss <- sample_sums(phy)
  keep <- unname(ss) >= min_reads
  if (any(!keep)) {
    message(label, ": dropping ", sum(!keep), " sample(s) below ", min_reads, ": ",
            paste(names(ss)[!keep], collapse = ", "))
  }
  phy <- prune_samples(keep, phy)
  prune_taxa(taxa_sums(phy) > 0, phy)
}

berger_parker_from_phy <- function(phy) {
  otu <- as(otu_table(phy), "matrix")
  if (!taxa_are_rows(phy)) otu <- t(otu)
  apply(otu, 2, function(x) max(x) / sum(x))
}

prepare_alpha_data <- function(phy, label, min_sum_taxa = 5000L, rarefy_depth = 5000L) {
  phy <- collapse_technical_replicates(phy)
  phy <- prune_samples_min_depth(phy, label, min_sum_taxa)
  if (rarefy_depth > 0L) {
    m <- min(sample_sums(phy))
    if (m < rarefy_depth) {
      stop(label, ": min library size (", m, ") < rarefy depth (", rarefy_depth, ")")
    }
    phy <- rarefy_even_depth(
      phy,
      sample.size = rarefy_depth,
      rngseed = 1L,
      replace = FALSE,
      trimOTUs = TRUE,
      verbose = FALSE
    )
  }

  phy_asv <- phyloseq::prune_taxa(taxa_sums(phy) > 0, phy)

  bp <- tibble(
    sample = sample_names(phy_asv),
    Berger_Parker = unname(berger_parker_from_phy(phy_asv))
  )
  meta <- as(sample_data(phy_asv), "data.frame") %>% tibble::rownames_to_column("sample")

  alpha <- phy_asv %>%
    estimate_richness(measures = c("Observed", "Shannon")) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample") %>%
    mutate(sample = gsub("\\.", "-", sample)) %>%
    left_join(meta %>% mutate(sample = gsub("\\.", "-", sample)), by = "sample") %>%
    left_join(bp, by = "sample") %>%
    mutate(
      Week = factor(as.character(Date), levels = c("3", "6", "9", "10")),
      Corral = factor(CorralLetter),
      log10_plastic = log10(plastic_concentration + 1),
      Berger_Parker_logit = qlogis(pmin(pmax(Berger_Parker, 1e-6), 1 - 1e-6)),
      Endpoint = factor(
        ifelse(Location == "WS", "Wall strip (week 9)", "Microscope slide (week 10)"),
        levels = c("Wall strip (week 9)", "Microscope slide (week 10)")
      ),
      plastic_level = factor(
        as.character(plastic_level),
        levels = c("none", "low", "medium", "high")
      )
    )

  alpha
}

# Order samples by plastic concentration within each location (consistent across time facets).
arrange_alpha_for_plots <- function(alpha) {
  alpha %>%
    mutate(
      sample = as.character(sample),
      Location = factor(Location, levels = c("WS", "MS"))
    ) %>%
    group_by(Location) %>%
    mutate(sample = forcats::fct_reorder(sample, plastic_concentration, .desc = FALSE)) %>%
    ungroup() %>%
    group_by(Location) %>%
    mutate(sample = factor(sample, levels = unique(sample[order(plastic_concentration)]))) %>%
    ungroup()
}

anova_p <- function(aov_tab, term) {
  if (is.null(aov_tab) || !term %in% rownames(aov_tab)) return(NA_real_)
  if ("p-value" %in% colnames(aov_tab)) return(aov_tab[term, "p-value"])
  if ("Pr(>F)" %in% colnames(aov_tab)) return(aov_tab[term, "Pr(>F)"])
  NA_real_
}

attach_bh <- function(results, term) {
  ids <- names(results)
  p <- vapply(results, function(r) anova_p(r$anova, term), numeric(1))
  names(p) <- ids
  p_bh <- bh_adjust_named(p)
  stats::setNames(lapply(ids, function(id) {
    r <- results[[id]]
    r$p_omnibus <- p[[id]]
    r$p_omnibus_bh <- p_bh[[id]]
    r
  }), ids)
}

bh_adjust_named <- function(p_named) {
  p <- unname(p_named)
  names(p) <- names(p_named)
  stats::p.adjust(p, method = "BH")
}

format_p_line <- function(label, p, p_bh = NULL) {
  if (is.na(p)) return(paste0(label, ": NA"))
  line <- sprintf("%s: p = %s", label, format.pval(p, digits = 3, eps = 0.001))
  if (!is.null(p_bh) && !is.na(p_bh)) {
    line <- paste0(line, " (BH = ", format.pval(p_bh, digits = 3, eps = 0.001), ")")
  }
  line
}

pairwise_caption <- function(pw, max_lines = 3L) {
  d <- as.data.frame(pw)
  if (!nrow(d)) return("")
  pcol <- intersect(c("p.value", "P.value"), names(d))[1]
  lines <- paste0(d$contrast, ": p = ", format.pval(d[[pcol]], digits = 3, eps = 0.001))
  paste(c("Pairwise (BH):", head(lines, max_lines)), collapse = "\n")
}

format_signif_p <- function(p) {
  if (is.na(p)) return("NA")
  stars <- if (p < 0.001) {
    "***"
  } else if (p < 0.01) {
    "**"
  } else if (p < 0.05) {
    "*"
  } else {
    ""
  }
  paste0(format(p, digits = 4, scientific = FALSE), stars)
}

parse_emmeans_contrast <- function(contrast_str) {
  parts <- trimws(strsplit(contrast_str, " - ", fixed = TRUE)[[1]])
  if (length(parts) != 2L) return(NULL)
  gsub("^Week", "", parts)
}

pairwise_comparisons <- function(pw) {
  d <- as.data.frame(pw)
  if (!nrow(d)) {
    return(list(comparisons = list(), annotations = character()))
  }
  pcol <- intersect(c("p.value", "P.value"), names(d))[1]
  comps <- lapply(d$contrast, parse_emmeans_contrast)
  keep <- vapply(comps, function(x) !is.null(x), logical(1))
  list(
    comparisons = comps[keep],
    annotations = vapply(d[[pcol]][keep], format_signif_p, character(1))
  )
}

add_pairwise_signif <- function(p, df, pw, y_col, step_frac = 0.1) {
  if (!requireNamespace("ggsignif", quietly = TRUE)) {
    return(p)
  }
  pk <- pairwise_comparisons(pw)
  n <- length(pk$comparisons)
  if (!n) return(p)

  y_vals <- df[[y_col]]
  y_max <- max(y_vals, na.rm = TRUE)
  y_span <- diff(range(y_vals, na.rm = TRUE))
  if (!is.finite(y_span) || y_span == 0) {
    y_span <- max(abs(y_max), 1) * 0.1
  }
  step <- step_frac * y_span
  y_pos <- y_max + step * seq_len(n)

  p +
    ggsignif::geom_signif(
      comparisons = pk$comparisons,
      annotations = pk$annotations,
      y_position = y_pos,
      tip_length = 0.02,
      vjust = -0.1,
      textsize = 3,
      na.rm = TRUE
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.08 + 0.1 * n))) +
    ggplot2::theme(plot.margin = ggplot2::margin(5.5, 5.5, 15, 5.5))
}
