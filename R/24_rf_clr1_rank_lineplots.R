# CLR (pc=1) rank line plots — genus-coloured ASV loess lines + points.
#
# Faceted by taxonomic rank (Phylum → Genus). One loess line per ASV;
# colour = Genus (global stable palette); black loess ± SE = rank-group mean
# (drawn in background). Pages grouped by panel size (n ASVs).
#
# ASV subset: exact k-of-3 full-n importance > 0 for clr_1, clr_0.5, clr_0.1.
# Y-axis: vegan CLR pseudocount=1. Week-9 WS, n=9, prev2of9.
#
# Writes: figures/rf_loocv_importance_clr_pc/ranks/{16S,18S}/
# Override: RF_RANK_FIG_DIR, RF_RANK_EXACT_N (default = all 3), RF_OUT_DIR.
#
# Run: Rscript R/24_rf_clr1_rank_lineplots.R

lines <- readLines("R/16_rclr_family_lineplots.R")
cut <- which(grepl("^# ---- CLR \\(pc=1\\)", lines))[1]
if (is.na(cut) || cut < 2L) {
  stop("Could not find CLR run block in R/16_rclr_family_lineplots.R")
}
eval(parse(text = lines[seq_len(cut - 1L)]), envir = environment())

RANKS <- c("Phylum", "Class", "Order", "Family", "Genus")
TRANSFORM_SUITE <- c("clr_1", "clr_0.5", "clr_0.1")
N_SUITE <- length(TRANSFORM_SUITE)
N_EXACT <- as.integer(Sys.getenv("RF_RANK_EXACT_N", unset = as.character(N_SUITE)))
if (!is.finite(N_EXACT) || N_EXACT < 1L || N_EXACT > N_SUITE) {
  stop("RF_RANK_EXACT_N must be in 1..", N_SUITE, ", got: ", N_EXACT)
}
SUITE_TAG <- sprintf("exact%dof%d", N_EXACT, N_SUITE)
SUITE_FIG_DIR <- Sys.getenv(
  "RF_RANK_FIG_DIR",
  unset = "figures/rf_loocv_importance_clr_pc"
)
RANK_FALLBACK <- c(
  Phylum = "Class",
  Class = "Order",
  Order = "Family",
  Family = "Genus",
  Genus = "Family"
)

out_dirs <- file.path(SUITE_FIG_DIR, "ranks")
for (mk in c("16S", "18S")) {
  dir.create(file.path(out_dirs, mk), showWarnings = FALSE, recursive = TRUE)
}

pos_asvs_one_transform <- function(marker, transform) {
  path <- rf_rds_path(marker, transform)
  if (!file.exists(path)) {
    stop("Missing RF RDS: ", path)
  }
  obj <- readRDS(path)
  n_trees <- as.integer(obj$design$num_trees)
  if (!identical(n_trees, 10000L)) {
    warning(glue("{marker} {transform} RDS has num_trees={n_trees}, expected 10000"))
  }
  imp <- obj$result$importance
  vi <- dplyr::coalesce(
    imp$median_loo_importance, imp$full_n_importance, imp$mean_importance
  )
  imp$asv[is.finite(vi) & vi > 0]
}

select_asvs_exact_suite <- function(marker) {
  sets <- lapply(TRANSFORM_SUITE, function(tr) pos_asvs_one_transform(marker, tr))
  names(sets) <- TRANSFORM_SUITE
  all_asvs <- unique(unlist(sets))
  n_tr <- vapply(all_asvs, function(a) {
    sum(vapply(sets, function(s) a %in% s, logical(1)))
  }, numeric(1))
  asvs <- sort(unname(all_asvs[n_tr == N_EXACT]))

  csv_path <- file.path(
    SUITE_FIG_DIR,
    marker,
    sprintf("intersection_fulln_gt0_%s_%s.csv", SUITE_TAG, marker)
  )
  tax <- readRDS(rf_rds_path(marker, "clr_1"))$result$importance %>%
    dplyr::select(asv, dplyr::any_of("Genus"))
  wide <- tax %>% dplyr::filter(asv %in% asvs)
  for (tr in TRANSFORM_SUITE) {
    imp <- readRDS(rf_rds_path(marker, tr))$result$importance
    vi <- if ("full_n_importance" %in% names(imp)) {
      dplyr::coalesce(imp$full_n_importance, imp$mean_importance)
    } else {
      imp$mean_importance
    }
    tmp <- tibble::tibble(asv = imp$asv, !!tr := vi)
    wide <- wide %>% dplyr::left_join(tmp, by = "asv")
  }
  which_tr <- vapply(asvs, function(a) {
    paste(TRANSFORM_SUITE[vapply(sets, function(s) a %in% s, logical(1))], collapse = ";")
  }, character(1))
  wide <- wide %>%
    dplyr::mutate(
      n_transforms = N_EXACT,
      transforms = which_tr[match(asv, asvs)]
    ) %>%
    dplyr::arrange(Genus, asv)
  utils::write.csv(wide, csv_path, row.names = FALSE)
  message(glue("{marker}: {SUITE_TAG} n_ASV={length(asvs)} -> {csv_path}"))
  tibble::tibble(asv = asvs)
}

rank_panel_label <- function(primary, fallback) {
  family_label(primary, fallback)
}

transform_asv_long_rank <- function(phy, asv_ids, rank_name) {
  Z <- compute_transform(phy, method = "clr_1")
  asv_ids <- intersect(asv_ids, colnames(Z))
  if (!length(asv_ids)) {
    stop("No selected ASVs present after clr_1")
  }
  Z_sub <- Z[, asv_ids, drop = FALSE]

  tax <- as.data.frame(phyloseq::tax_table(phy))
  for (rnk in c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")) {
    if (!rnk %in% names(tax)) tax[[rnk]] <- NA_character_
  }
  fb <- RANK_FALLBACK[[rank_name]]
  tax <- tax %>%
    tibble::rownames_to_column("asv") %>%
    dplyr::filter(asv %in% asv_ids) %>%
    dplyr::mutate(
      FamilyLabel = rank_panel_label(.data[[rank_name]], .data[[fb]]),
      GenusLabel = rank_panel_label(Genus, Family),
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
      tax %>% dplyr::select(asv, FamilyLabel, GenusLabel, asv_label, LowestLabel),
      by = "asv"
    )
}

OKABE_ITO_8 <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#999999"
)
GOLDEN_RATIO <- 0.618033988749895

extended_distinct_palette <- function(n) {
  if (n <= 0L) return(character())
  if (n <= length(OKABE_ITO_8)) return(OKABE_ITO_8[seq_len(n)])
  cols <- character(n)
  cols[seq_len(length(OKABE_ITO_8))] <- OKABE_ITO_8
  i <- seq_len(n - length(OKABE_ITO_8)) + length(OKABE_ITO_8)
  h <- (i * GOLDEN_RATIO) %% 1 * 360
  c <- 62 + (i %% 4L) * 7
  l <- 48 + (i %% 5L) * 6
  cols[i] <- grDevices::hcl(h = h, c = c, l = l)
  cols
}

distinct_colour_map <- function(labels) {
  labs <- sort(unique(as.character(labels)))
  labs <- labs[!is.na(labs) & nzchar(labs)]
  if (!length(labs)) return(character())
  stats::setNames(extended_distinct_palette(length(labs)), labs)
}

subset_colour_map <- function(colour_map, labels) {
  labs <- sort(unique(as.character(labels)))
  labs <- labs[!is.na(labs) & nzchar(labs)]
  out <- colour_map[labs]
  if (any(is.na(out))) {
    stop("Colour map missing labels: ", paste(labs[is.na(out)], collapse = ", "))
  }
  out
}

build_marker_genus_colour_map <- function(phy, asv_ids) {
  genus_long <- transform_asv_long_rank(phy, asv_ids, rank_name = "Genus")
  distinct_colour_map(genus_long$GenusLabel)
}

plot_one_rank_panel <- function(df_grp, order_tbl, rank_name,
                                colour_map,
                                show_x_labels = TRUE, compact = FALSE) {
  n_asv <- dplyr::n_distinct(df_grp$asv)
  if (n_asv == 1L) {
    panel_title <- as.character(df_grp$LowestLabel[[1]])
  } else {
    lab <- as.character(df_grp$FamilyLabel[[1]])
    panel_title <- sprintf("%s (n=%d)", lab, n_asv)
  }
  all_labs <- names(colour_map)
  legend_ncol <- max(1L, as.integer(ceiling(length(all_labs) / 8L)))

  df_grp <- df_grp %>%
    dplyr::mutate(
      GenusLabel = factor(as.character(GenusLabel), levels = all_labs),
      asv_label = factor(asv_label, levels = sort(unique(asv_label)))
    )

  df_mean <- df_grp %>%
    dplyr::group_by(x_pos, sample_id) %>%
    dplyr::summarise(rank_mean = mean(value), .groups = "drop")

  p <- ggplot2::ggplot()

  if (n_asv > 1L) {
    p <- p +
      ggplot2::geom_smooth(
        data = df_mean,
        ggplot2::aes(x = x_pos, y = rank_mean),
        method = "loess",
        formula = y ~ x,
        span = LOESS_SPAN,
        se = TRUE,
        colour = "black",
        fill = "grey35",
        linewidth = if (compact) 0.9 else 1.1,
        alpha = 0.2,
        method.args = list(degree = 1L)
      )
  }

  p <- p +
    ggplot2::geom_smooth(
      data = df_grp,
      ggplot2::aes(
        x = x_pos,
        y = value,
        colour = GenusLabel,
        group = asv_label
      ),
      method = "loess",
      formula = y ~ x,
      span = LOESS_SPAN,
      se = FALSE,
      linewidth = if (compact) 0.55 else 0.7,
      alpha = 0.85,
      method.args = list(degree = 1L)
    ) +
    ggplot2::geom_point(
      data = df_grp,
      ggplot2::aes(x = x_pos, y = value, colour = GenusLabel),
      size = if (compact) 1.0 else 1.3,
      alpha = 0.55,
      position = ggplot2::position_jitter(width = 0.18, height = 0.04, seed = 1L)
    ) +
    ggplot2::scale_x_continuous(
      breaks = order_tbl$x_pos,
      labels = if (show_x_labels) order_tbl$tick_label else order_tbl$sample_id
    ) +
    ggplot2::scale_colour_manual(
      values = colour_map,
      drop = FALSE,
      name = "Genus",
      guide = ggplot2::guide_legend(
        ncol = legend_ncol,
        override.aes = list(alpha = 1, linewidth = 1.2, size = 2)
      )
    ) +
    ggplot2::labs(
      title = panel_title,
      x = NULL,
      y = ylab_for_method("clr_1")
    ) +
    ggplot2::theme_bw(base_size = if (compact) 8 else 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = if (compact) 7.5 else 9),
      axis.text.x = ggplot2::element_text(size = if (compact) 4.5 else 5.5, lineheight = 0.85),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = if (compact) 6.5 else 8, face = "bold"),
      legend.text = ggplot2::element_text(size = if (compact) 5.5 else 6.5),
      legend.key.size = ggplot2::unit(if (compact) 0.28 else 0.35, "cm"),
      panel.grid.minor = ggplot2::element_blank()
    )
  p
}

arrange_pages_rank <- function(plot_df, order_tbl, title, rank_name,
                               genus_colour_map,
                               ncol = 4L, nrow = 3L) {
  grp_n <- plot_df %>%
    dplyr::distinct(FamilyLabel, asv) %>%
    dplyr::count(FamilyLabel, name = "n_asv") %>%
    dplyr::arrange(n_asv, FamilyLabel)

  sizes <- sort(unique(grp_n$n_asv))
  page_plots <- list()

  for (n_size in sizes) {
    grps <- grp_n$FamilyLabel[grp_n$n_asv == n_size]
    per_page <- ncol * nrow
    n_pages_size <- max(1L, as.integer(ceiling(length(grps) / per_page)))

    for (pg in seq_len(n_pages_size)) {
      idx <- ((pg - 1L) * per_page + 1L):min(pg * per_page, length(grps))
      grps_pg <- grps[idx]
      panels <- lapply(grps_pg, function(grp) {
        df_grp <- dplyr::filter(plot_df, FamilyLabel == grp)
        plot_one_rank_panel(
          df_grp, order_tbl,
          rank_name = rank_name,
          colour_map = subset_colour_map(genus_colour_map, df_grp$GenusLabel),
          show_x_labels = TRUE
        )
      })
      while (length(panels) < per_page) {
        panels[[length(panels) + 1L]] <- patchwork::plot_spacer()
      }
      grid <- patchwork::wrap_plots(panels, ncol = ncol, nrow = nrow) +
        patchwork::plot_layout(guides = "keep")
      loess_note <- if (n_size >= 2L) {
        sprintf("background black = %s mean loess ± SE", tolower(rank_name))
      } else {
        "n=1: ASV loess + points only"
      }
      page_plots[[length(page_plots) + 1L]] <- grid +
        patchwork::plot_annotation(
          title = title,
          subtitle = sprintf(
            "%s size n=%d | page %d/%d | colour = Genus | one loess line per ASV + points | %s | x: plastic_concentration (even)",
            rank_name, n_size, pg, n_pages_size, loess_note
          ),
          theme = ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 12),
            plot.subtitle = ggplot2::element_text(size = 9)
          )
        )
    }
  }

  list(plots = page_plots, n_groups = nrow(grp_n), n_pages = length(page_plots))
}

run_marker_rank <- function(marker, rank_name, pdf_name, sel,
                            genus_colour_map,
                            ncol_pdf = 4L, nrow_pdf = 3L) {
  message(glue(
    "=== {marker} | clr_1 | facet={rank_name} | {SUITE_TAG} n={nrow(sel)} ==="
  ))
  if (!nrow(sel)) {
    stop(glue("No ASVs in {SUITE_TAG} set for {marker}"))
  }
  phy <- load_filtered_phy(marker)
  order_tbl <- sample_order_tbl(phy)
  asv_long <- transform_asv_long_rank(phy, sel$asv, rank_name = rank_name)
  plot_df <- asv_long %>% dplyr::inner_join(order_tbl, by = "sample_id")

  n_asv <- nrow(sel)
  n_grp <- dplyr::n_distinct(plot_df$FamilyLabel)
  title <- glue(
    "{marker} CLR (pc=1) by {rank_name} | {SUITE_TAG} | n_ASV={n_asv} | n_{rank_name}={n_grp}"
  )
  pages <- arrange_pages_rank(
    plot_df, order_tbl, title = title, rank_name = rank_name,
    genus_colour_map = genus_colour_map,
    ncol = ncol_pdf, nrow = nrow_pdf
  )

  pdf_path <- file.path(out_dirs, marker, pdf_name)
  pdf_out <- save_multipage_pdf(pages$plots, pdf_path)
  message("  PDF: ", pdf_out)
  message(glue(
    "{marker}/{rank_name}: n_ASV={n_asv} | n_{rank_name}={n_grp} | pages={pages$n_pages}"
  ))
  invisible(list(
    n_asv = n_asv,
    n_group = n_grp,
    n_pages = pages$n_pages,
    pdf = pdf_out
  ))
}

ncol_for <- function(marker) if (identical(marker, "18S")) 3L else 4L

results <- list()
for (marker in c("16S", "18S")) {
  sel <- select_asvs_exact_suite(marker)
  phy <- load_filtered_phy(marker)
  genus_colour_map <- build_marker_genus_colour_map(phy, sel$asv)
  message(glue("{marker}: global genus palette — {length(genus_colour_map)} genera"))
  for (rnk in RANKS) {
    key <- paste(marker, rnk, sep = "_")
    results[[key]] <- run_marker_rank(
      marker = marker,
      rank_name = rnk,
      pdf_name = sprintf("clr1_%s_%s_%s.pdf", tolower(rnk), marker, SUITE_TAG),
      sel = sel,
      genus_colour_map = genus_colour_map,
      ncol_pdf = ncol_for(marker),
      nrow_pdf = 3L
    )
  }
}

message("Done.")
for (nm in names(results)) {
  message(glue("{nm}: pages={results[[nm]]$n_pages} | {results[[nm]]$pdf}"))
}
