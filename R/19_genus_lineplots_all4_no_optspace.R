# Genus-level figures for ASVs with LOO mean importance > 0.001 in ALL 4 of
#   clr_0.1, clr_0.5, clr_1, rclr  (excludes rclr_optspace).
# Facets by Genus (Family fallback if Genus missing); ASV lines within panels.
# Also writes a genus membership CSV per marker.
#
# Run: Rscript R/19_genus_lineplots_all4_no_optspace.R

lines <- readLines("R/16_rclr_family_lineplots.R")
cut <- which(grepl("^# ---- CLR \\(pc=1\\)", lines))[1]
if (is.na(cut) || cut < 2L) {
  stop("Could not find CLR run block in R/16_rclr_family_lineplots.R")
}
eval(parse(text = lines[seq_len(cut - 1L)]), envir = environment())

TRANSFORM_CORE4 <- c("clr_0.1", "clr_0.5", "clr_1", "rclr")
MIN_TRANSFORMS <- length(TRANSFORM_CORE4)

genus_panel_label <- function(genus, family) {
  # Prefer Genus; fall back to Family (same style as family_label but swapped).
  family_label(genus, family)
}

# Override: facet key = Genus (stored as FamilyLabel for reuse of plot helpers)
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

  out <- file.path(
    out_dir,
    sprintf("genus_all4_gt0.001_no_optspace_%s.csv", marker)
  )
  utils::write.csv(tab, out, row.names = FALSE)
  message(glue(
    "{marker}: n_ASV={length(sel_asvs)} | n_Genus={nrow(tab)} | wrote {out}"
  ))
  tab
}

# Patch run_marker title wording to say Genus
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
  write_genus_table(marker, sel$asv, phy)

  asv_long <- transform_asv_long(phy, sel$asv, method = method)
  order_tbl <- sample_order_tbl(phy)
  plot_df <- asv_long %>% dplyr::inner_join(order_tbl, by = "sample_id")

  n_asv <- nrow(sel)
  n_gen <- dplyr::n_distinct(plot_df$FamilyLabel)
  title <- glue(
    "{marker} {title_tag_for_method(method)} by Genus | mean_loo_imp > {IMP_THRESH} in {min_transforms}/{length(transforms)} (no optspace) | n_ASV={n_asv} | n_Genus={n_gen}"
  )
  pages <- arrange_pages_by_n(
    plot_df, order_tbl, title = title, method = method,
    ncol = ncol_pdf, nrow = nrow_pdf
  )
  pdf_path <- file.path(out_dir, pdf_name)
  pdf_out <- save_multipage_pdf(pages$plots, pdf_path)
  message(glue("{marker}/{method}: n_ASV={n_asv} | n_Genus={n_gen} | pages={pages$n_pages}"))
  message("  PDF: ", pdf_out)
  invisible(list(n_asv = n_asv, n_genus = n_gen, pdf = pdf_out))
}

message(glue(
  "Subset: mean_loo_imp > {IMP_THRESH} in {MIN_TRANSFORMS}/{length(TRANSFORM_CORE4)} ",
  "(no optspace); facets = Genus"
))

# Lead figure: mean4 (matches primary transform consensus)
res_16s <- run_marker_genus(
  marker = "16S",
  transforms = TRANSFORM_CORE4,
  min_transforms = MIN_TRANSFORMS,
  method = "clr_sum4",
  pdf_name = "clr_mean4_genus_16S_all4_gt0.001_no_optspace.pdf",
  ncol_pdf = 4L,
  nrow_pdf = 3L
)

res_18s <- run_marker_genus(
  marker = "18S",
  transforms = TRANSFORM_CORE4,
  min_transforms = MIN_TRANSFORMS,
  method = "clr_sum4",
  pdf_name = "clr_mean4_genus_18S_all4_gt0.001_no_optspace.pdf",
  ncol_pdf = 3L,
  nrow_pdf = 3L
)

message("Done.")
message(glue("16S: {res_16s$pdf}"))
message(glue("18S: {res_18s$pdf}"))
