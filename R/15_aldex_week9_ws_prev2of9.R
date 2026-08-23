# ALDEx2 on week-9 wall-strip samples (ASV or Genus).
#
# Load: R/01_load_files2.R (no full-table tax_filter)
# Samples: week-9 WS, all 9 primary (keep WS_H_9)
# Prevalence:
#   ASV (default): ASV >=3 reads in >=2/9
#   Genus: tax_glom to Genus FIRST, then genus >=3 reads in >=2/9
#
# Primary bridge (parallel Spearman correlation arms, shared taxon universe):
#   A) Nominal loading: log1p(plastic_concentration)
#   B) Measured retention: log1p(particles_total_d20)
# Secondary: four-level plastic_level KW (low power at n=9)
#
# Binary (optional): retention high (top 2 by particles_total_d20) vs rest (7)
#   — aldex.ttest (Wilcoxon / Welch) + aldex.effect; suffix _binary_high2
#
# Sensitivity (optional): retention arm without WS_H_9 — sample subset only, no refilter
#
# Server:
#   Rscript R/15_aldex_week9_ws_prev2of9.R
#   ALDEX_MC_SAMPLES=128 ALDEX_DATASET=both Rscript R/15_aldex_week9_ws_prev2of9.R
#   ALDEX_TAX_LEVEL=Genus ALDEX_DATASET=both Rscript R/15_aldex_week9_ws_prev2of9.R
#   ALDEX_RUN_RETENTION_NO_H=0 Rscript R/15_aldex_week9_ws_prev2of9.R  # skip H sensitivity
#   ALDEX_RUN_BINARY_HIGH2=0 Rscript R/15_aldex_week9_ws_prev2of9.R     # skip binary arm

suppressPackageStartupMessages({
  library(phyloseq)
  library(ALDEx2)
  library(dplyr)
  library(tibble)
})

source("R/01_load_files2.R")

mc_samples <- as.integer(Sys.getenv("ALDEX_MC_SAMPLES", unset = "128"))
dataset <- toupper(Sys.getenv("ALDEX_DATASET", unset = "both"))
min_reads <- as.integer(Sys.getenv("ALDEX_MIN_READS", unset = "3"))
min_samples <- as.integer(Sys.getenv("ALDEX_MIN_SAMPLES", unset = "2"))
out_dir <- Sys.getenv("ALDEX_OUT_DIR", unset = "output")
retention_drop <- trimws(Sys.getenv("ALDEX_RETENTION_DROP", unset = "WS_H_9"))
run_retention_no_h <- !identical(Sys.getenv("ALDEX_RUN_RETENTION_NO_H", unset = "1"), "0")
run_binary_high2 <- !identical(Sys.getenv("ALDEX_RUN_BINARY_HIGH2", unset = "1"), "0")
binary_top_n <- as.integer(Sys.getenv("ALDEX_BINARY_TOP_N", unset = "2"))
tax_level_raw <- trimws(Sys.getenv("ALDEX_TAX_LEVEL", unset = "ASV"))
tax_level <- ifelse(
  toupper(tax_level_raw) %in% c("GENUS", "G"),
  "Genus",
  "ASV"
)

stopifnot(dataset %in% c("16S", "18S", "BOTH"))
stopifnot(tax_level %in% c("ASV", "Genus"))
stopifnot(is.finite(binary_top_n), binary_top_n >= 1L)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

unit_label <- if (identical(tax_level, "Genus")) "genus" else "ASV"
tag_suffix <- if (identical(tax_level, "Genus")) "_genus" else ""

message(sprintf(
  "ALDEx2 week-9 WS | mc.samples=%d | dataset=%s | tax=%s | filter=>=%d reads in >=%d/9 %s | retention_no_H=%s | binary_high%d=%s",
  mc_samples, dataset, tax_level, min_reads, min_samples, unit_label,
  if (run_retention_no_h) retention_drop else "off",
  binary_top_n,
  if (run_binary_high2) "on" else "off"
))

subset_week9_ws <- function(phy) {
  phy %>%
    phyloseq::subset_samples(as.character(Date) == "9") %>%
    phyloseq::subset_samples(as.character(Location) == "WS")
}

aggregate_to_genus <- function(phy) {
  phyloseq::tax_glom(phy, taxrank = "Genus", NArm = FALSE)
}

filter_prevalence <- function(phy, min_reads = 3L, min_samples = 2L, taxonomic_unit = "ASV") {
  n_samp <- phyloseq::nsamples(phy)
  otu <- as(phyloseq::otu_table(phy), "matrix")
  if (!phyloseq::taxa_are_rows(phy)) {
    otu <- t(otu)
  }
  keep <- rowSums(otu >= min_reads) >= min_samples
  list(
    counts = otu[keep, , drop = FALSE],
    keep_mask = keep,
    n_taxa_before = nrow(otu),
    n_taxa_after = sum(keep),
    n_samples = n_samp,
    min_reads = min_reads,
    min_samples = min_samples,
    prevalence = min_samples / n_samp,
    taxonomic_unit = taxonomic_unit
  )
}

annotate_taxa <- function(taxon_ids, phy) {
  tax <- as.data.frame(phyloseq::tax_table(phy)) %>%
    tibble::rownames_to_column("asv")
  if (!"Genus" %in% names(tax)) {
    tax$Genus <- NA_character_
  }
  tax %>%
    dplyr::filter(asv %in% taxon_ids) %>%
    dplyr::select(asv, Genus, dplyr::any_of(c("Species", "Family")))
}

prep_counts <- function(phy, label) {
  phy_sub <- subset_week9_ws(phy)
  if (identical(tax_level, "Genus")) {
    phy_work <- aggregate_to_genus(phy_sub)
    filt <- filter_prevalence(
      phy_work,
      min_reads = min_reads,
      min_samples = min_samples,
      taxonomic_unit = "Genus"
    )
    message(sprintf(
      "%s week-9 WS: n=%d | genera %d -> %d (aggregate then >= %d reads in >= %d/9, %.1f%%)",
      label, filt$n_samples, filt$n_taxa_before, filt$n_taxa_after,
      filt$min_reads, filt$min_samples, 100 * filt$prevalence
    ))
  } else {
    phy_work <- phy_sub
    filt <- filter_prevalence(
      phy_work,
      min_reads = min_reads,
      min_samples = min_samples,
      taxonomic_unit = "ASV"
    )
    message(sprintf(
      "%s week-9 WS: n=%d | ASVs %d -> %d (>= %d reads in >= %d/9, %.1f%%)",
      label, filt$n_samples, filt$n_taxa_before, filt$n_taxa_after,
      filt$min_reads, filt$min_samples, 100 * filt$prevalence
    ))
  }

  sd <- as(phyloseq::sample_data(phy_work), "data.frame")
  sd$sample_id <- rownames(sd)

  list(
    counts = filt$counts,
    filter = filt,
    sample_data = sd,
    phy = phy_work
  )
}

run_aldex_arm <- function(counts, sample_ids, conds, cont_var, arm_label) {
  counts <- counts[, sample_ids, drop = FALSE]
  conds <- unname(conds[sample_ids])
  cont_var <- unname(cont_var[sample_ids])

  set.seed(2026L)
  clr <- ALDEx2::aldex.clr(
    counts,
    conds = conds,
    mc.samples = mc_samples,
    denom = "all",
    gamma = NULL,
    verbose = FALSE
  )
  corr <- ALDEx2::aldex.corr(clr, cont.var = cont_var)
  kw <- ALDEx2::aldex.kw(clr)

  list(
    arm = arm_label,
    sample_ids = sample_ids,
    n_samples = length(sample_ids),
    continuous = cont_var,
    conditions = conds,
    clr = clr,
    correlation = corr,
    kw = kw
  )
}

# Binary high (top-N by particles_total_d20) vs rest: Wilcoxon/Welch + effect.
assign_binary_retention <- function(sample_ids, retention_raw, top_n = binary_top_n) {
  retention_raw <- setNames(as.numeric(retention_raw[sample_ids]), sample_ids)
  if (anyNA(retention_raw)) {
    stop("binary retention: NA in particles_total_d20 for: ",
         paste(names(retention_raw)[is.na(retention_raw)], collapse = ", "))
  }
  ord <- order(retention_raw, decreasing = TRUE, method = "radix")
  ranked_ids <- sample_ids[ord]
  ranked_vals <- unname(retention_raw[ord])
  if (top_n >= length(sample_ids)) {
    stop("binary retention: top_n (", top_n, ") must be < n_samples (", length(sample_ids), ")")
  }
  # Tie at the high/rest cut: include all samples tied with the Nth value.
  cut_val <- ranked_vals[top_n]
  high_ids <- ranked_ids[ranked_vals >= cut_val]
  if (length(high_ids) != top_n) {
    warning(sprintf(
      "binary retention: tie at cut (particles_total_d20 >= %.6g) expands high group to n=%d (requested top_n=%d)",
      cut_val, length(high_ids), top_n
    ))
  }
  conds <- setNames(rep("rest", length(sample_ids)), sample_ids)
  conds[high_ids] <- "high"
  list(
    conditions = conds,
    high_ids = high_ids,
    rest_ids = setdiff(sample_ids, high_ids),
    retention_raw = retention_raw,
    cut_value = cut_val,
    top_n_requested = top_n
  )
}

run_aldex_binary_arm <- function(counts, sample_ids, binary_conds, arm_label) {
  counts <- counts[, sample_ids, drop = FALSE]
  conds <- unname(binary_conds[sample_ids])
  n_levels <- length(unique(conds))
  if (n_levels != 2L) {
    stop("binary arm requires exactly 2 condition levels, got ", n_levels)
  }

  set.seed(2026L)
  clr <- ALDEx2::aldex.clr(
    counts,
    conds = conds,
    mc.samples = mc_samples,
    denom = "all",
    gamma = NULL,
    verbose = FALSE
  )
  ttest <- ALDEx2::aldex.ttest(clr, paired.test = FALSE, verbose = FALSE)
  effect <- ALDEx2::aldex.effect(clr, verbose = FALSE)

  list(
    arm = arm_label,
    sample_ids = sample_ids,
    n_samples = length(sample_ids),
    conditions = conds,
    clr = clr,
    ttest = ttest,
    effect = effect
  )
}

save_arm_tables <- function(arm_obj, tax_annot, tag, arm_suffix) {
  corr_cols <- c(
    "spearman.erho", "spearman.ep", "spearman.eBH",
    "pearson.ecor", "pearson.ep", "pearson.eBH"
  )
  corr_tbl <- arm_obj$correlation %>%
    tibble::rownames_to_column("asv") %>%
    dplyr::left_join(tax_annot, by = "asv")

  kw_cols <- intersect(c("kw.ep", "kw.eBH", "glm.ep", "glm.eBH"), names(arm_obj$kw))
  kw_tbl <- arm_obj$kw %>%
    tibble::rownames_to_column("asv") %>%
    dplyr::left_join(tax_annot, by = "asv")

  out_path <- file.path(out_dir, sprintf("aldex_%s_%s.rds", tag, arm_suffix))
  saveRDS(
    list(
      arm = arm_obj$arm,
      sample_ids = arm_obj$sample_ids,
      n_samples = arm_obj$n_samples,
      filter_tag = tag,
      tax_level = tax_level,
      correlation_table = corr_tbl,
      kw_table = kw_tbl,
      raw = list(
        correlation = arm_obj$correlation,
        kw = arm_obj$kw
      )
    ),
    out_path
  )
  message("Wrote ", normalizePath(out_path, winslash = "/", mustWork = FALSE))

  list(
    path = out_path,
    correlation = corr_tbl %>% dplyr::select(asv, Genus, dplyr::any_of(corr_cols)),
    kw = kw_tbl %>% dplyr::select(asv, Genus, dplyr::any_of(kw_cols))
  )
}

save_binary_arm_tables <- function(arm_obj, tax_annot, tag, arm_suffix, binary_meta) {
  ttest_cols <- intersect(
    c("we.ep", "we.eBH", "wi.ep", "wi.eBH"),
    names(arm_obj$ttest)
  )
  effect_prefer <- c(
    "rab.all", "rab.win.high", "rab.win.rest",
    "diff.btw", "diff.win", "effect", "overlap"
  )
  effect_cols <- intersect(effect_prefer, names(arm_obj$effect))
  if (!length(effect_cols)) {
    effect_cols <- names(arm_obj$effect)
  }

  ttest_tbl <- arm_obj$ttest %>%
    tibble::rownames_to_column("asv") %>%
    dplyr::left_join(tax_annot, by = "asv")
  effect_tbl <- arm_obj$effect %>%
    tibble::rownames_to_column("asv") %>%
    dplyr::left_join(tax_annot, by = "asv")
  combined <- dplyr::full_join(
    ttest_tbl,
    effect_tbl %>% dplyr::select(asv, dplyr::all_of(names(arm_obj$effect))),
    by = "asv"
  )

  out_path <- file.path(out_dir, sprintf("aldex_%s_%s.rds", tag, arm_suffix))
  saveRDS(
    list(
      arm = arm_obj$arm,
      sample_ids = arm_obj$sample_ids,
      n_samples = arm_obj$n_samples,
      filter_tag = tag,
      tax_level = tax_level,
      binary = list(
        high_ids = binary_meta$high_ids,
        rest_ids = binary_meta$rest_ids,
        cut_value = binary_meta$cut_value,
        top_n_requested = binary_meta$top_n_requested,
        retention_raw = binary_meta$retention_raw,
        conditions = arm_obj$conditions
      ),
      ttest_table = ttest_tbl,
      effect_table = effect_tbl,
      combined_table = combined,
      raw = list(
        ttest = arm_obj$ttest,
        effect = arm_obj$effect
      )
    ),
    out_path
  )
  message("Wrote ", normalizePath(out_path, winslash = "/", mustWork = FALSE))

  list(
    path = out_path,
    high_ids = binary_meta$high_ids,
    rest_ids = binary_meta$rest_ids,
    ttest = ttest_tbl %>% dplyr::select(asv, Genus, dplyr::any_of(ttest_cols)),
    effect = effect_tbl %>% dplyr::select(asv, Genus, dplyr::any_of(effect_cols))
  )
}

run_marker <- function(phy, label) {
  prep <- prep_counts(phy, label)
  x9f <- prep$counts
  sd <- prep$sample_data
  all_ids <- colnames(x9f)

  tag <- sprintf(
    "%s_9_ws_prev%dof%d%s",
    label, min_samples, length(all_ids), tag_suffix
  )
  tax_annot <- annotate_taxa(rownames(x9f), prep$phy)

  conds <- setNames(as.character(sd[all_ids, "plastic_level"]), all_ids)
  nominal <- setNames(log1p(as.numeric(sd[all_ids, "plastic_concentration"])), all_ids)
  retention <- setNames(log1p(as.numeric(sd[all_ids, "particles_total_d20"])), all_ids)
  retention_raw <- setNames(as.numeric(sd[all_ids, "particles_total_d20"]), all_ids)

  message(sprintf("=== %s loading (nominal) n=%d tax=%s ===", label, length(all_ids), tax_level))
  arm_loading <- run_aldex_arm(
    x9f, all_ids, conds, nominal,
    arm_label = "nominal_loading_log1p_plastic_concentration"
  )
  out_loading <- save_arm_tables(arm_loading, tax_annot, tag, "loading_corr_kw")

  message(sprintf("=== %s retention (measured) n=%d tax=%s ===", label, length(all_ids), tax_level))
  arm_retention <- run_aldex_arm(
    x9f, all_ids, conds, retention,
    arm_label = "measured_retention_log1p_particles_total_d20"
  )
  out_retention <- save_arm_tables(arm_retention, tax_annot, tag, "retention_corr_kw")

  out_sens <- NULL
  if (isTRUE(run_retention_no_h) && nzchar(retention_drop)) {
    keep_ids <- setdiff(all_ids, retention_drop)
    if (length(keep_ids) < 3L) {
      warning(label, ": too few samples after retention sensitivity drop")
    } else {
      message(sprintf(
        "=== %s retention sensitivity (drop %s) n=%d tax=%s ===",
        label, retention_drop, length(keep_ids), tax_level
      ))
      arm_sens <- run_aldex_arm(
        x9f, keep_ids, conds, retention,
        arm_label = paste0("retention_sensitivity_drop_", retention_drop)
      )
      out_sens <- save_arm_tables(
        arm_sens, tax_annot, tag,
        sprintf("retention_corr_kw_no_%s", gsub("[^A-Za-z0-9]+", "_", retention_drop))
      )
    }
  }

  out_binary <- NULL
  if (isTRUE(run_binary_high2)) {
    bin <- assign_binary_retention(all_ids, retention_raw, top_n = binary_top_n)
    message(sprintf(
      "=== %s retention binary high%d vs rest | high=%s (cut=%.4g) | rest n=%d | tax=%s ===",
      label, binary_top_n,
      paste(bin$high_ids, collapse = ","),
      bin$cut_value,
      length(bin$rest_ids),
      tax_level
    ))
    arm_binary <- run_aldex_binary_arm(
      x9f, all_ids, bin$conditions,
      arm_label = sprintf(
        "retention_binary_high%d_vs_rest_particles_total_d20",
        binary_top_n
      )
    )
    out_binary <- save_binary_arm_tables(
      arm_binary, tax_annot, tag,
      sprintf("retention_binary_high%d", binary_top_n),
      binary_meta = bin
    )
  }

  summary_path <- file.path(out_dir, sprintf("aldex_%s_summary.rds", tag))
  saveRDS(
    list(
      marker = label,
      tax_level = tax_level,
      filter = prep$filter,
      tag = tag,
      arms = list(
        loading = out_loading,
        retention = out_retention,
        retention_sensitivity = out_sens,
        retention_binary_high2 = out_binary
      ),
      key_columns = list(
        correlation = c("spearman.erho", "spearman.ep", "spearman.eBH"),
        kw = c("kw.ep", "kw.eBH"),
        binary_ttest = c("wi.ep", "wi.eBH", "we.ep", "we.eBH"),
        binary_effect = c("effect", "diff.btw", "overlap")
      )
    ),
    summary_path
  )
  message("Wrote ", normalizePath(summary_path, winslash = "/", mustWork = FALSE))

  invisible(list(
    loading = out_loading,
    retention = out_retention,
    retention_sensitivity = out_sens,
    retention_binary_high2 = out_binary
  ))
}

if (dataset %in% c("16S", "BOTH")) {
  run_marker(phy_16S, "16S")
}
if (dataset %in% c("18S", "BOTH")) {
  run_marker(phy_18S, "18S")
}

message("Done.")
