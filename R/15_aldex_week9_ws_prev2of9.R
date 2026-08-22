# ALDEx2 on week-9 wall-strip samples (ASV level).
#
# Load: R/01_load_files2.R (no full-table tax_filter)
# Samples: week-9 WS, all 9 primary (keep WS_H_9)
# Prevalence: ASV >=3 reads in >=2/9 (same rule as RF; filter once on all 9)
#
# Primary bridge (parallel Spearman correlation arms, shared ASV universe):
#   A) Nominal loading: log1p(plastic_concentration)
#   B) Measured retention: log1p(particles_total_d20)
# Secondary: four-level plastic_level KW (low power at n=9)
#
# Sensitivity (optional): retention arm without WS_H_9 — sample subset only, no refilter
#
# Server:
#   Rscript R/15_aldex_week9_ws_prev2of9.R
#   ALDEX_MC_SAMPLES=128 ALDEX_DATASET=both Rscript R/15_aldex_week9_ws_prev2of9.R
#   ALDEX_RUN_RETENTION_NO_H=0 Rscript R/15_aldex_week9_ws_prev2of9.R  # skip H sensitivity

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

stopifnot(dataset %in% c("16S", "18S", "BOTH"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message(sprintf(
  "ALDEx2 week-9 WS | mc.samples=%d | dataset=%s | ASV filter=>=%d reads in >=%d/9 | retention_no_H=%s",
  mc_samples, dataset, min_reads, min_samples,
  if (run_retention_no_h) retention_drop else "off"
))

subset_week9_ws <- function(phy) {
  phy %>%
    phyloseq::subset_samples(as.character(Date) == "9") %>%
    phyloseq::subset_samples(as.character(Location) == "WS")
}

filter_prevalence_asv <- function(phy, min_reads = 3L, min_samples = 2L) {
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
    taxonomic_unit = "ASV"
  )
}

annotate_asv <- function(asv_ids, phy) {
  tax <- as.data.frame(phyloseq::tax_table(phy)) %>%
    tibble::rownames_to_column("asv")
  if (!"Genus" %in% names(tax)) {
    tax$Genus <- NA_character_
  }
  tax %>%
    dplyr::filter(asv %in% asv_ids) %>%
    dplyr::select(asv, Genus, dplyr::any_of(c("Species", "Family")))
}

prep_counts <- function(phy, label) {
  phy_sub <- subset_week9_ws(phy)
  filt <- filter_prevalence_asv(phy_sub, min_reads = min_reads, min_samples = min_samples)
  sd <- as(phyloseq::sample_data(phy_sub), "data.frame")
  sd$sample_id <- rownames(sd)

  message(sprintf(
    "%s week-9 WS: n=%d | ASVs %d -> %d (>= %d reads in >= %d/9, %.1f%%)",
    label, filt$n_samples, filt$n_taxa_before, filt$n_taxa_after,
    filt$min_reads, filt$min_samples, 100 * filt$prevalence
  ))

  list(
    counts = filt$counts,
    filter = filt,
    sample_data = sd,
    phy = phy_sub
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

run_marker <- function(phy, label) {
  prep <- prep_counts(phy, label)
  x9f <- prep$counts
  sd <- prep$sample_data
  all_ids <- colnames(x9f)

  tag <- sprintf("%s_9_ws_prev%dof%d", label, min_samples, length(all_ids))
  tax_annot <- annotate_asv(rownames(x9f), prep$phy)

  conds <- setNames(as.character(sd[all_ids, "plastic_level"]), all_ids)
  nominal <- setNames(log1p(as.numeric(sd[all_ids, "plastic_concentration"])), all_ids)
  retention <- setNames(log1p(as.numeric(sd[all_ids, "particles_total_d20"])), all_ids)

  message(sprintf("=== %s loading (nominal) n=%d ===", label, length(all_ids)))
  arm_loading <- run_aldex_arm(
    x9f, all_ids, conds, nominal,
    arm_label = "nominal_loading_log1p_plastic_concentration"
  )
  out_loading <- save_arm_tables(arm_loading, tax_annot, tag, "loading_corr_kw")

  message(sprintf("=== %s retention (measured) n=%d ===", label, length(all_ids)))
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
        "=== %s retention sensitivity (drop %s) n=%d ===",
        label, retention_drop, length(keep_ids)
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

  summary_path <- file.path(out_dir, sprintf("aldex_%s_summary.rds", tag))
  saveRDS(
    list(
      marker = label,
      filter = prep$filter,
      tag = tag,
      arms = list(
        loading = out_loading,
        retention = out_retention,
        retention_sensitivity = out_sens
      ),
      key_columns = list(
        correlation = c("spearman.erho", "spearman.ep", "spearman.eBH"),
        kw = c("kw.ep", "kw.eBH")
      )
    ),
    summary_path
  )
  message("Wrote ", normalizePath(summary_path, winslash = "/", mustWork = FALSE))

  invisible(list(
    loading = out_loading,
    retention = out_retention,
    retention_sensitivity = out_sens
  ))
}

if (dataset %in% c("16S", "BOTH")) {
  run_marker(phy_16S, "16S")
}
if (dataset %in% c("18S", "BOTH")) {
  run_marker(phy_18S, "18S")
}

message("Done.")
