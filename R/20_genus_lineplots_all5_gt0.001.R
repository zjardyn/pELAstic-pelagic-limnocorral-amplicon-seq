# Genus-level figures for ASVs with LOO mean importance > 0.001 in ALL 5 of
#   clr_0.1, clr_0.5, clr_1, rclr, rclr_optspace.
# One PDF per transform (y-axis matches that transform).
# Also writes genus membership CSV per marker.
#
# Run: Rscript R/20_genus_lineplots_all5_gt0.001.R

lines <- readLines("R/16_rclr_family_lineplots.R")
cut <- which(grepl("^# ---- CLR \\(pc=1\\)", lines))[1]
if (is.na(cut) || cut < 2L) {
  stop("Could not find CLR run block in R/16_rclr_family_lineplots.R")
}
eval(parse(text = lines[seq_len(cut - 1L)]), envir = environment())

TRANSFORM_ALL5 <- c("clr_0.1", "clr_0.5", "clr_1", "rclr", "rclr_optspace")
MIN_TRANSFORMS <- length(TRANSFORM_ALL5)

genus_panel_label <- function(genus, family) {
  family_label(genus, family)
}

# Extend transforms to all 5 RF methods (incl. rclr_optspace).
compute_one_transform <- function(X, transform_id) {
  transform_id <- as.character(transform_id)
  if (startsWith(transform_id, "clr_")) {
    pc <- as.numeric(sub("^clr_", "", transform_id))
    Z <- vegan::decostand(X, method = "clr", MARGIN = 1L, pseudocount = pc)
  } else if (identical(transform_id, "rclr")) {
    Z <- vegan::decostand(X, method = "rclr", MARGIN = 1L, impute = FALSE)
    Z[!is.finite(Z)] <- 0
  } else if (identical(transform_id, "rclr_optspace")) {
    Z <- vegan::decostand(X, method = "rclr", MARGIN = 1L, impute = TRUE)
  } else {
    stop("Unsupported transform_id: ", transform_id)
  }
  finish_Z(Z)
}

compute_transform <- function(phy, method) {
  method <- as.character(method)
  if (!method %in% TRANSFORM_ALL5) {
    stop("Unsupported method: ", method)
  }
  X <- otu_samples_x_taxa(phy)
  compute_one_transform(X, method)
}

ylab_for_method <- function(method) {
  method <- as.character(method)
  if (startsWith(method, "clr_")) {
    sprintf("CLR (pseudocount=%s)", sub("^clr_", "", method))
  } else if (identical(method, "rclr")) {
    "rCLR (impute=FALSE)"
  } else if (identical(method, "rclr_optspace")) {
    "rCLR optspace (impute=TRUE)"
  } else {
    method
  }
}

title_tag_for_method <- function(method) {
  method <- as.character(method)
  if (startsWith(method, "clr_")) {
    sprintf("CLR (pc=%s)", sub("^clr_", "", method))
  } else if (identical(method, "rclr")) {
    "rCLR (impute=FALSE)"
  } else if (identical(method, "rclr_optspace")) {
    "rCLR optspace"
  } else {
    method
  }
}

transform_asv_long <- function(phy, asv_ids, method) {
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
      FamilyLabel = genus_panel_label(Genus, Family),
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

write_genus_table <- function(marker, sel_asvs, phy) {
  tax <- as.data.frame(phyloseq::tax_table(phy)) %>%
    tibble::rownames_to_column("asv") %>%
    dplyr::filter(asv %in% sel_asvs) %>%
    dplyr::mutate(
      GenusLabel = genus_panel_label(Genus, Family),
      Family = as.character(Family)
    )

  tab <- tax %>%
    dplyr::group_by(GenusLabel, Family) %>%
    dplyr::summarise(
      n_ASV = dplyr::n(),
      asvs = paste(asv, collapse = ";"),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(n_ASV), GenusLabel)

  out <- file.path(out_dir, sprintf("genus_all5_gt0.001_%s.csv", marker))
  utils::write.csv(tab, out, row.names = FALSE)
  message(glue(
    "{marker}: n_ASV={length(sel_asvs)} | n_Genus={nrow(tab)} | wrote {out}"
  ))
  tab
}

pdf_name_for <- function(method, marker) {
  file.path(
    out_dir,
    sprintf("%s_genus_%s_all5_gt0.001.pdf", method, marker)
  )
}

run_marker_genus <- function(marker, transforms, min_transforms, method,
                             pdf_name, ncol_pdf = 4L, nrow_pdf = 3L) {
  message(glue(
    "=== {marker} | {method} | genus facets | min_transforms={min_transforms} of {length(transforms)} ==="
  ))
  imp_long <- load_importance(marker, transforms)
  sel <- select_asvs(imp_long, min_transforms)
  if (!nrow(sel)) {
    stop(glue("No ASVs passed for {marker}"))
  }
  phy <- load_filtered_phy(marker)

  asv_long <- transform_asv_long(phy, sel$asv, method = method)
  order_tbl <- sample_order_tbl(phy)
  plot_df <- asv_long %>% dplyr::inner_join(order_tbl, by = "sample_id")

  n_asv <- nrow(sel)
  n_gen <- dplyr::n_distinct(plot_df$FamilyLabel)
  title <- glue(
    "{marker} {title_tag_for_method(method)} by Genus | mean_loo_imp > {IMP_THRESH} in {min_transforms}/{length(transforms)} | n_ASV={n_asv} | n_Genus={n_gen}"
  )
  pages <- arrange_pages_by_n(
    plot_df, order_tbl, title = title, method = method,
    ncol = ncol_pdf, nrow = nrow_pdf
  )
  pdf_out <- save_multipage_pdf(pages$plots, pdf_name)
  message(glue("{marker}/{method}: n_ASV={n_asv} | n_Genus={n_gen} | pages={pages$n_pages}"))
  message("  PDF: ", pdf_out)
  invisible(list(n_asv = n_asv, n_genus = n_gen, pdf = pdf_out))
}

run_all_plots_marker <- function(marker, ncol_pdf, nrow_pdf) {
  imp_long <- load_importance(marker, TRANSFORM_ALL5)
  sel <- select_asvs(imp_long, MIN_TRANSFORMS)
  if (!nrow(sel)) {
    stop(glue("No ASVs passed for {marker}"))
  }
  phy <- load_filtered_phy(marker)
  write_genus_table(marker, sel$asv, phy)

  results <- lapply(TRANSFORM_ALL5, function(method) {
    run_marker_genus(
      marker = marker,
      transforms = TRANSFORM_ALL5,
      min_transforms = MIN_TRANSFORMS,
      method = method,
      pdf_name = pdf_name_for(method, marker),
      ncol_pdf = ncol_pdf,
      nrow_pdf = nrow_pdf
    )
  })
  names(results) <- TRANSFORM_ALL5
  results
}

message(glue(
  "Subset: mean_loo_imp > {IMP_THRESH} in {MIN_TRANSFORMS}/{length(TRANSFORM_ALL5)} ",
  "(clr_0.1, clr_0.5, clr_1, rclr, rclr_optspace); facets = Genus"
))

res_16s <- run_all_plots_marker("16S", ncol_pdf = 4L, nrow_pdf = 3L)
res_18s <- run_all_plots_marker("18S", ncol_pdf = 3L, nrow_pdf = 3L)

message("Done.")
message(glue(
  "16S: n_ASV={res_16s[[1]]$n_asv} | n_Genus={res_16s[[1]]$n_genus} | 5 PDFs written"
))
message(glue(
  "18S: n_ASV={res_18s[[1]]$n_asv} | n_Genus={res_18s[[1]]$n_genus} | 5 PDFs written"
))
