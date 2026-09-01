# Family ASV line plots — week-9 wall-strip RF ASV subsets.
#
# Same samples / prevalence as R/14: Date==9, Location==WS, n=9;
#   >=3 reads in >=2 of 9 samples.
# Subsets from LOO mean importance across transforms:
#   16S: mean_loo_importance > 0.001 in >=3 of 4 (clr_0.1, clr_0.5, clr_1, rclr)
#        — excludes rclr_optspace
#   18S: mean_loo_importance > 0.001 in all 5 transforms (incl. rclr_optspace)
#
# Values (PDF sets):
#   CLR       — vegan CLR with pseudocount=1 (matches RF `clr_1`)
#   rCLR      — vegan rCLR with impute=FALSE (matches RF `rclr`; no optspace)
#   clr_sum4  — mean of clr_0.1 + clr_0.5 + clr_1 + rclr
#               (elementwise sum / 4; each sample still sums ~0 across taxa)
#
# Facets: one panel per Family; straight lines = ASVs;
#   loess ± SE = family mean only when n_ASV >= 2.
#   n=1 panels titled by lowest available taxon rank (Genus > Family > ...).
# Pages grouped by family size (n ASVs).
# X-axis: samples ordered by increasing plastic_concentration, evenly spaced.
# PDF only (no PNG).
#
# Run: Rscript R/16_rclr_family_lineplots.R

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(glue)
  library(patchwork)
  library(vegan)
})

out_dir <- Sys.getenv(
  "CLR_FAMILY_OUT_DIR",
  unset = "figures/rf_loocv_importance"
)
rf_dir <- Sys.getenv("RF_OUT_DIR", unset = "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

CLR_PSEUDOCOUNT <- 1
TRANSFORM_16S <- c("clr_0.1", "clr_0.5", "clr_1", "rclr")
TRANSFORM_18S <- c("clr_0.1", "clr_0.5", "clr_1", "rclr", "rclr_optspace")
IMP_THRESH <- 0.001
LOESS_SPAN <- 0.85

rf_rds_path <- function(marker, transform) {
  file.path(
    rf_dir,
    sprintf("rf_%s_9_ws_n9_prev2of9_%s.rds", marker, transform)
  )
}

load_importance <- function(marker, transforms) {
  rows <- lapply(transforms, function(tr) {
    path <- rf_rds_path(marker, tr)
    if (!file.exists(path)) {
      stop("Missing RF RDS: ", path)
    }
    obj <- readRDS(path)
    imp <- obj$result$importance
    if (!"mean_loo_importance" %in% names(imp)) {
      stop("mean_loo_importance missing in ", path)
    }
    imp %>%
      dplyr::select(asv, mean_loo_importance) %>%
      dplyr::mutate(transform = tr)
  })
  dplyr::bind_rows(rows)
}

select_asvs <- function(imp_long, min_transforms) {
  imp_long %>%
    dplyr::filter(mean_loo_importance > IMP_THRESH) %>%
    dplyr::count(asv, name = "n_transforms_pass") %>%
    dplyr::filter(n_transforms_pass >= min_transforms) %>%
    dplyr::arrange(asv)
}

load_filtered_phy <- function(marker) {
  path <- rf_rds_path(marker, "clr_1")
  if (!file.exists(path)) {
    path <- rf_rds_path(marker, "rclr")
  }
  if (!file.exists(path)) {
    stop("Missing RF RDS for filtered phy: ", path)
  }
  obj <- readRDS(path)
  phy <- obj$filter$phy
  if (is.null(phy) || !inherits(phy, "phyloseq")) {
    stop("filter$phy missing/invalid in ", path)
  }
  phy
}

clean_taxon <- function(x) {
  x <- as.character(x)
  bad <- is.na(x) | !nzchar(x) | x %in% c("NA", "NaN", "unclassified", "Unclassified")
  x[bad] <- NA_character_
  x
}

family_label <- function(family, genus) {
  fam <- clean_taxon(family)
  gen <- clean_taxon(genus)
  out <- fam
  miss <- is.na(out)
  out[miss] <- gen[miss]
  out[is.na(out)] <- "(unassigned)"
  out
}

asv_display_label <- function(asv, genus, prefix_n = 8L) {
  short <- substr(as.character(asv), 1L, prefix_n)
  gen <- clean_taxon(genus)
  ifelse(is.na(gen), short, paste0(gen, " (", short, ")"))
}

otu_samples_x_taxa <- function(phy) {
  X <- as(phyloseq::otu_table(phy), "matrix")
  if (phyloseq::taxa_are_rows(phy)) {
    X <- t(X)
  }
  X <- X[, colSums(X, na.rm = TRUE) > 0, drop = FALSE]
  storage.mode(X) <- "double"
  X
}

finish_Z <- function(Z) {
  Z[!is.finite(Z)] <- 0
  Z <- as.matrix(Z)
  bad_col <- vapply(seq_len(ncol(Z)), function(j) any(!is.finite(Z[, j])), logical(1))
  if (any(bad_col)) {
    Z <- Z[, !bad_col, drop = FALSE]
  }
  Z
}

# Single-method transforms matching R/14
compute_one_transform <- function(X, transform_id) {
  if (startsWith(transform_id, "clr_")) {
    pc <- as.numeric(sub("^clr_", "", transform_id))
    Z <- vegan::decostand(X, method = "clr", MARGIN = 1L, pseudocount = pc)
  } else if (identical(transform_id, "rclr")) {
    Z <- vegan::decostand(X, method = "rclr", MARGIN = 1L, impute = FALSE)
  } else {
    stop("Unsupported transform_id: ", transform_id)
  }
  finish_Z(Z)
}

# method: "clr_1", "rclr", or "clr_sum4" (clr_0.1+clr_0.5+clr_1+rclr)
compute_transform <- function(phy, method = c("clr_1", "rclr", "clr_sum4")) {
  method <- match.arg(method)
  X <- otu_samples_x_taxa(phy)
  if (identical(method, "clr_sum4")) {
    ids <- c("clr_0.1", "clr_0.5", "clr_1", "rclr")
    mats <- lapply(ids, function(id) compute_one_transform(X, id))
    # Align columns to shared taxa (intersection)
    common <- Reduce(intersect, lapply(mats, colnames))
    if (!length(common)) {
      stop("No shared taxa across clr_sum4 components")
    }
    Z <- Reduce(`+`, lapply(mats, function(m) m[, common, drop = FALSE])) / length(ids)
    return(finish_Z(Z))
  }
  if (identical(method, "clr_1")) {
    return(compute_one_transform(X, "clr_1"))
  }
  compute_one_transform(X, "rclr")
}

# Lowest available rank for n=1 panel titles (Genus preferred)
lowest_taxon_label <- function(asv, Kingdom = NA, Phylum = NA, Class = NA,
                               Order = NA, Family = NA, Genus = NA) {
  ranks <- list(
    clean_taxon(Genus),
    clean_taxon(Family),
    clean_taxon(Order),
    clean_taxon(Class),
    clean_taxon(Phylum),
    clean_taxon(Kingdom)
  )
  for (r in ranks) {
    if (!is.na(r) && nzchar(r)) {
      return(r)
    }
  }
  substr(as.character(asv), 1L, 8L)
}

transform_asv_long <- function(phy, asv_ids, method = c("clr_1", "rclr", "clr_sum4")) {
  method <- match.arg(method)
  Z <- compute_transform(phy, method = method)
  asv_ids <- intersect(asv_ids, colnames(Z))
  if (!length(asv_ids)) {
    stop("No selected ASVs present after ", method)
  }
  Z_sub <- Z[, asv_ids, drop = FALSE]

  tax <- as.data.frame(phyloseq::tax_table(phy))
  for (rnk in c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")) {
    if (!rnk %in% names(tax)) tax[[rnk]] <- NA_character_
  }
  tax <- tax %>%
    tibble::rownames_to_column("asv") %>%
    dplyr::filter(asv %in% asv_ids) %>%
    dplyr::mutate(
      FamilyLabel = family_label(Family, Genus),
      asv_label = asv_display_label(asv, Genus),
      LowestLabel = mapply(
        lowest_taxon_label,
        asv, Kingdom, Phylum, Class, Order, Family, Genus,
        USE.NAMES = FALSE
      )
    ) %>%
    dplyr::group_by(FamilyLabel, asv_label) %>%
    dplyr::mutate(
      asv_label = dplyr::if_else(
        dplyr::n() > 1L,
        paste0(asv_label, ".", dplyr::row_number()),
        asv_label
      )
    ) %>%
    dplyr::ungroup()

  as.data.frame(Z_sub) %>%
    tibble::rownames_to_column("sample_id") %>%
    tidyr::pivot_longer(-sample_id, names_to = "asv", values_to = "value") %>%
    dplyr::left_join(
      tax %>% dplyr::select(asv, FamilyLabel, asv_label, LowestLabel),
      by = "asv"
    )
}

ylab_for_method <- function(method) {
  if (identical(method, "clr_1")) {
    sprintf("CLR (pseudocount=%g)", CLR_PSEUDOCOUNT)
  } else if (identical(method, "rclr")) {
    "rCLR (impute=FALSE)"
  } else {
    "mean(clr_0.1, clr_0.5, clr_1, rclr)"
  }
}

title_tag_for_method <- function(method) {
  if (identical(method, "clr_1")) {
    sprintf("CLR (pc=%g)", CLR_PSEUDOCOUNT)
  } else if (identical(method, "rclr")) {
    "rCLR (impute=FALSE)"
  } else {
    "CLR mean4 ((0.1+0.5+1+rclr)/4)"
  }
}

sample_order_tbl <- function(phy) {
  sd <- as(phyloseq::sample_data(phy), "data.frame")
  sd %>%
    tibble::rownames_to_column("sample_id") %>%
    dplyr::mutate(
      plastic_concentration = as.numeric(plastic_concentration),
      particles_total_d20 = as.numeric(particles_total_d20)
    ) %>%
    dplyr::arrange(plastic_concentration, sample_id) %>%
    dplyr::mutate(
      x_pos = dplyr::row_number(),
      tick_label = sprintf("%s\n(%.0f)", sample_id, plastic_concentration)
    ) %>%
    dplyr::select(
      sample_id, x_pos, tick_label,
      plastic_concentration, particles_total_d20
    )
}

plot_one_family <- function(df_fam, order_tbl, method,
                            show_x_labels = TRUE, compact = FALSE) {
  n_asv <- dplyr::n_distinct(df_fam$asv)
  # n=1: title by lowest taxon rank; n>=2: Family (+ count)
  if (n_asv == 1L) {
    panel_title <- as.character(df_fam$LowestLabel[[1]])
  } else {
    fam <- as.character(df_fam$FamilyLabel[[1]])
    panel_title <- sprintf("%s (n=%d)", fam, n_asv)
  }
  df_fam <- df_fam %>%
    dplyr::mutate(asv_label = factor(asv_label, levels = sort(unique(asv_label))))

  df_mean <- df_fam %>%
    dplyr::group_by(x_pos, sample_id) %>%
    dplyr::summarise(family_mean = mean(value), .groups = "drop")

  legend_size <- if (compact) {
    if (n_asv <= 6L) 4.5 else 3.8
  } else if (n_asv <= 6L) {
    6.5
  } else if (n_asv <= 12L) {
    5.5
  } else {
    4.5
  }
  ncol_leg <- if (n_asv <= 4L) 1L else if (n_asv <= 10L) 2L else 3L

  p <- ggplot2::ggplot(
    df_fam,
    ggplot2::aes(
      x = x_pos,
      y = value,
      colour = asv_label,
      group = asv_label
    )
  ) +
    ggplot2::geom_line(
      linewidth = if (compact) 0.45 else 0.55,
      alpha = 0.85
    ) +
    ggplot2::geom_point(size = if (compact) 1.0 else 1.4, alpha = 0.45)

  if (n_asv > 1L) {
    p <- p +
      ggplot2::geom_smooth(
        data = df_mean,
        ggplot2::aes(x = x_pos, y = family_mean),
        inherit.aes = FALSE,
        method = "loess",
        formula = y ~ x,
        span = LOESS_SPAN,
        se = TRUE,
        colour = "black",
        fill = "grey35",
        linewidth = if (compact) 0.9 else 1.1,
        alpha = 0.25,
        method.args = list(degree = 1L)
      )
  }

  p <- p +
    ggplot2::scale_x_continuous(
      breaks = order_tbl$x_pos,
      labels = if (show_x_labels) order_tbl$tick_label else order_tbl$sample_id
    ) +
    ggplot2::scale_colour_discrete(name = NULL) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(ncol = ncol_leg, byrow = TRUE)
    ) +
    ggplot2::labs(
      title = panel_title,
      x = NULL,
      y = ylab_for_method(method)
    ) +
    ggplot2::theme_bw(base_size = if (compact) 8 else 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = if (compact) 7.5 else 9),
      axis.text.x = ggplot2::element_text(size = if (compact) 4.5 else 5.5, lineheight = 0.85),
      legend.text = ggplot2::element_text(size = legend_size),
      legend.key.height = ggplot2::unit(0.3, "lines"),
      legend.key.width = ggplot2::unit(1.0, "lines"),
      legend.position = "bottom",
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (n_asv == 1L) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  p
}

arrange_pages_by_n <- function(plot_df, order_tbl, title, method,
                               ncol = 4L, nrow = 3L) {
  fam_n <- plot_df %>%
    dplyr::distinct(FamilyLabel, asv) %>%
    dplyr::count(FamilyLabel, name = "n_asv") %>%
    dplyr::arrange(n_asv, FamilyLabel)

  sizes <- sort(unique(fam_n$n_asv))
  page_plots <- list()
  page_meta <- list()

  for (n_size in sizes) {
    fams <- fam_n$FamilyLabel[fam_n$n_asv == n_size]
    per_page <- ncol * nrow
    n_pages_size <- max(1L, as.integer(ceiling(length(fams) / per_page)))

    for (pg in seq_len(n_pages_size)) {
      idx <- ((pg - 1L) * per_page + 1L):min(pg * per_page, length(fams))
      fams_pg <- fams[idx]
      panels <- lapply(fams_pg, function(fam) {
        df_fam <- dplyr::filter(plot_df, FamilyLabel == fam)
        plot_one_family(df_fam, order_tbl, method = method, show_x_labels = TRUE)
      })
      while (length(panels) < per_page) {
        panels[[length(panels) + 1L]] <- patchwork::plot_spacer()
      }
      grid <- patchwork::wrap_plots(panels, ncol = ncol, nrow = nrow)
      page_i <- length(page_plots) + 1L
      loess_note <- if (n_size >= 2L) {
        "black = family mean loess ± SE (n≥2 only)"
      } else {
        "n=1: ASV line only (no loess)"
      }
      page_plots[[page_i]] <- grid +
        patchwork::plot_annotation(
          title = title,
          subtitle = sprintf(
            "Family size n=%d | page %d/%d within n=%d (%d families) | %s | x: plastic_concentration (even)",
            n_size, pg, n_pages_size, n_size, length(fams), loess_note
          ),
          theme = ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 12),
            plot.subtitle = ggplot2::element_text(size = 9)
          )
        )
      page_meta[[page_i]] <- list(n_asv = n_size, page_within = pg, n_families = length(fams_pg))
    }
  }

  list(
    plots = page_plots,
    n_families = nrow(fam_n),
    n_pages = length(page_plots),
    family_sizes = fam_n,
    sizes = sizes,
    page_meta = page_meta
  )
}

save_multipage_pdf <- function(plots, path, width = 16, height = 12) {
  grDevices::pdf(path, width = width, height = height, useDingbats = FALSE)
  for (p in plots) print(p)
  grDevices::dev.off()
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

run_marker <- function(marker, transforms, min_transforms, method, pdf_name,
                       ncol_pdf = 4L, nrow_pdf = 3L) {
  method <- match.arg(method, c("clr_1", "rclr", "clr_sum4"))
  message(glue(
    "=== {marker} | {method} | min_transforms={min_transforms} of {length(transforms)} ==="
  ))
  imp <- load_importance(marker, transforms)
  sel <- select_asvs(imp, min_transforms)
  phy <- load_filtered_phy(marker)
  order_tbl <- sample_order_tbl(phy)

  src <- if (identical(method, "clr_1")) {
    sprintf(
      "vegan::decostand(clr, pseudocount=%g) on filter$phy from %s",
      CLR_PSEUDOCOUNT,
      basename(rf_rds_path(marker, "clr_1"))
    )
  } else if (identical(method, "rclr")) {
    sprintf(
      "vegan::decostand(rclr, impute=FALSE) on filter$phy from %s",
      basename(rf_rds_path(marker, "clr_1"))
    )
  } else {
    sprintf(
      "mean(clr_0.1 + clr_0.5 + clr_1 + rclr[impute=FALSE]) on filter$phy from %s",
      basename(rf_rds_path(marker, "clr_1"))
    )
  }
  message("  Transform source: ", src)

  asv_long <- transform_asv_long(phy, sel$asv, method = method)
  # Sanity: sample sums of selected ASVs need not be 0 (subset), but full Z does.
  # Report mean |sample sum| on full transform matrix for selected taxa only.
  Z_chk <- compute_transform(phy, method = method)
  samp_sums <- rowSums(Z_chk[, intersect(sel$asv, colnames(Z_chk)), drop = FALSE])
  message(sprintf(
    "  Selected-taxa sample sums: mean=%.3g | max|abs|=%.3g (subset need not be 0)",
    mean(samp_sums), max(abs(samp_sums))
  ))
  full_sums <- rowSums(Z_chk)
  message(sprintf(
    "  Full-matrix sample sums: mean=%.3g | max|abs|=%.3g (expect ~0)",
    mean(full_sums), max(abs(full_sums))
  ))

  plot_df <- asv_long %>% dplyr::inner_join(order_tbl, by = "sample_id")

  n_asv <- nrow(sel)
  n_fam <- dplyr::n_distinct(plot_df$FamilyLabel)
  fam_n <- plot_df %>%
    dplyr::distinct(FamilyLabel, asv) %>%
    dplyr::count(FamilyLabel, name = "n_asv") %>%
    dplyr::arrange(n_asv, FamilyLabel)

  title <- glue(
    "{marker} {title_tag_for_method(method)} by Family | mean_loo_imp > {IMP_THRESH} in {min_transforms}/{length(transforms)} | n_ASV={n_asv} | n_Family={n_fam}"
  )
  pages <- arrange_pages_by_n(
    plot_df, order_tbl, title = title, method = method,
    ncol = ncol_pdf, nrow = nrow_pdf
  )

  marker_dir <- file.path(out_dir, marker)
  dir.create(marker_dir, showWarnings = FALSE, recursive = TRUE)
  pdf_path <- file.path(marker_dir, pdf_name)
  pdf_out <- save_multipage_pdf(pages$plots, pdf_path)

  message(glue(
    "{marker}/{method}: n_ASV={n_asv} | n_Family={n_fam} | pdf pages={pages$n_pages}"
  ))
  message("  PDF: ", pdf_out)

  list(
    marker = marker,
    method = method,
    n_asv = n_asv,
    n_family = n_fam,
    pdf = pdf_out,
    n_pages = pages$n_pages
  )
}

# ---- CLR (pc=1) ----
res_16s_clr <- run_marker(
  marker = "16S",
  transforms = TRANSFORM_16S,
  min_transforms = 3L,
  method = "clr_1",
  pdf_name = "clr1_family_16S_ge3of4_no_optspace.pdf",
  ncol_pdf = 4L,
  nrow_pdf = 3L
)

res_18s_clr <- run_marker(
  marker = "18S",
  transforms = TRANSFORM_18S,
  min_transforms = 5L,
  method = "clr_1",
  pdf_name = "clr1_family_18S_all5_gt0.001.pdf",
  ncol_pdf = 3L,
  nrow_pdf = 3L
)

# ---- rCLR (impute=FALSE; no optspace) ----
res_16s_rclr <- run_marker(
  marker = "16S",
  transforms = TRANSFORM_16S,
  min_transforms = 3L,
  method = "rclr",
  pdf_name = "rclr_family_16S_ge3of4_no_optspace.pdf",
  ncol_pdf = 4L,
  nrow_pdf = 3L
)

res_18s_rclr <- run_marker(
  marker = "18S",
  transforms = TRANSFORM_18S,
  min_transforms = 5L,
  method = "rclr",
  pdf_name = "rclr_family_18S_all5_gt0.001.pdf",
  ncol_pdf = 3L,
  nrow_pdf = 3L
)

# ---- Mean of 4: (clr_0.1 + clr_0.5 + clr_1 + rclr) / 4 ----
res_16s_sum4 <- run_marker(
  marker = "16S",
  transforms = TRANSFORM_16S,
  min_transforms = 3L,
  method = "clr_sum4",
  pdf_name = "clr_mean4_family_16S_ge3of4_no_optspace.pdf",
  ncol_pdf = 4L,
  nrow_pdf = 3L
)

res_18s_sum4 <- run_marker(
  marker = "18S",
  transforms = TRANSFORM_18S,
  min_transforms = 5L,
  method = "clr_sum4",
  pdf_name = "clr_mean4_family_18S_all5_gt0.001.pdf",
  ncol_pdf = 3L,
  nrow_pdf = 3L
)

message("Done.")
message(glue("CLR 16S: {res_16s_clr$pdf}"))
message(glue("CLR 18S: {res_18s_clr$pdf}"))
message(glue("rCLR 16S: {res_16s_rclr$pdf}"))
message(glue("rCLR 18S: {res_18s_rclr$pdf}"))
message(glue("mean4 16S: {res_16s_sum4$pdf}"))
message(glue("mean4 18S: {res_18s_sum4$pdf}"))
