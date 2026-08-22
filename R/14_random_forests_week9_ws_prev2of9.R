# Random forests (week-9 wall strip) — ASV-level transform + bootstrap importance.
#
# Load: R/01_load_files2.R (no full-table tax_filter)
# Samples: Date == 9, Location == WS; keep all 9 (including WS_H_9)
# Prevalence (primary): >=3 reads in >=2 of 9 samples (~22%)
# Response: log1p(particles_total_d20)
# Features: ASVs (not genus-aggregated)
#
# Robustness design:
#   - No outer CV (no 5-fold / LOOCV)
#   - Fixed RF params (mtry = floor(sqrt(p)), min.node.size = 5)
#   - 1000 outer bootstrap resamples of the 9 samples (sample composition)
#   - ranger permutation importance within each bootstrap fit
#   - OOB error recorded as a descriptive diagnostic only
#   - Transform grid measures zero-handling sensitivity (not the main robustness test)
#
# Transform grid:
#   clr_1, clr_0.1  — vegan CLR with QUBS-matched pseudocounts
#   rclr                     — vegan rCLR without matrix completion (impute=FALSE)
#   rclr_optspace            — vegan rCLR with optspace completion (impute=TRUE)
#
# Per transform: mean importance ± SD, topk_freq across bootstraps
# Across transforms: n_transforms_in_topk (not only strict intersection)
#
# Server:
#   RF_N_CL=40 RF_DATASET=both Rscript R/14_random_forests_week9_ws_prev2of9.R

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ranger)
  library(vegan)
  library(glue)
})

source("R/01_load_files2.R")

# --- runtime knobs ---
num_cores <- as.integer(Sys.getenv("RF_N_CL", unset = "40"))
n_boot <- as.integer(Sys.getenv("RF_N_BOOT", unset = Sys.getenv("RF_N_RUNS", unset = "1000")))
num_trees <- as.integer(Sys.getenv("RF_NUM_TREES", unset = "5000"))
dataset <- toupper(Sys.getenv("RF_DATASET", unset = "both")) # 16S | 18S | both
min_reads <- as.integer(Sys.getenv("RF_MIN_READS", unset = "3"))
min_samples <- as.integer(Sys.getenv("RF_MIN_SAMPLES", unset = "2"))
top_k <- as.integer(Sys.getenv("RF_TOP_K", unset = "20"))
out_dir <- Sys.getenv("RF_OUT_DIR", unset = "output")
drop_samples <- strsplit(
  Sys.getenv("RF_DROP_SAMPLES", unset = ""),
  ",",
  fixed = TRUE
)[[1]] |>
  trimws()
drop_samples <- drop_samples[nzchar(drop_samples)]

stopifnot(dataset %in% c("16S", "18S", "BOTH"))
options(ranger.num.threads = num_cores)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

drop_label <- if (length(drop_samples)) paste(drop_samples, collapse = ",") else "none"

message(sprintf(
  "RF week-9 WS | cores=%d | outer_boot=%d | trees=%d | dataset=%s | filter=>=%d reads in >=%d samples | drop=%s | top_k=%d | no outer CV",
  num_cores, n_boot, num_trees, dataset, min_reads, min_samples,
  drop_label, top_k
))

# --- prep ---
subset_week9_ws_rf <- function(phy, drop = drop_samples) {
  phy <- phy %>%
    phyloseq::subset_samples(as.character(Date) == "9") %>%
    phyloseq::subset_samples(as.character(Location) == "WS") %>%
    phyloseq::subset_samples(!is.na(particles_total_d20))
  if (length(drop)) {
    keep <- setdiff(phyloseq::sample_names(phy), drop)
    phy <- phyloseq::prune_samples(keep, phy)
  }
  if (phyloseq::nsamples(phy) < 3L) {
    stop("Too few samples after exclusion: ", phyloseq::nsamples(phy))
  }
  phy
}

filter_prevalence_subset <- function(phy, min_reads = 3L, min_samples = 2L) {
  n_samp <- phyloseq::nsamples(phy)
  if (min_samples > n_samp) {
    stop(sprintf("min_samples=%d > n_samples=%d", min_samples, n_samp))
  }
  otu <- as(phyloseq::otu_table(phy), "matrix")
  if (!phyloseq::taxa_are_rows(phy)) {
    otu <- t(otu)
  }
  keep <- rowSums(otu >= min_reads) >= min_samples
  phy_f <- phyloseq::prune_taxa(keep, phy)
  phy_f <- phyloseq::prune_samples(phyloseq::sample_sums(phy_f) > 0, phy_f)
  list(
    phy = phy_f,
    n_taxa_before = nrow(otu),
    n_taxa_after = sum(keep),
    n_samples = phyloseq::nsamples(phy_f),
    min_reads = min_reads,
    min_samples = min_samples,
    prevalence = min_samples / n_samp
  )
}

counts_samples_by_taxa <- function(phy) {
  X <- as(phyloseq::otu_table(phy), "matrix")
  if (phyloseq::taxa_are_rows(phy)) {
    X <- t(X)
  }
  X <- X[, colSums(X, na.rm = TRUE) > 0, drop = FALSE]
  storage.mode(X) <- "double"
  X
}

# Transforms: return samples x taxa matrix (finite columns only)
transform_counts <- function(X, transform) {
  transform <- as.character(transform)
  if (startsWith(transform, "clr_")) {
    pc <- as.numeric(sub("^clr_", "", transform))
    if (!is.finite(pc) || pc <= 0) {
      stop("Bad CLR pseudocount in transform id: ", transform)
    }
    Z <- vegan::decostand(X, method = "clr", MARGIN = 1L, pseudocount = pc)
  } else if (identical(transform, "rclr")) {
    # Observed-parts geometric mean only; zeros are NA after log.
    # For dense RF input, set those entries to 0 on the rCLR scale
    # (no optspace / DEICODE completion).
    Z <- vegan::decostand(X, method = "rclr", MARGIN = 1L, impute = FALSE)
    Z[!is.finite(Z)] <- 0
  } else if (identical(transform, "rclr_optspace")) {
    # vegan 2.7+ optspace matrix completion (DEICODE-style; default impute=TRUE)
    Z <- vegan::decostand(X, method = "rclr", MARGIN = 1L, impute = TRUE)
  } else {
    stop("Unknown transform: ", transform)
  }
  Z <- as.matrix(Z)
  bad_col <- vapply(seq_len(ncol(Z)), function(j) any(!is.finite(Z[, j])), logical(1))
  if (any(bad_col)) {
    Z <- Z[, !bad_col, drop = FALSE]
  }
  bad_row <- apply(Z, 1, function(v) any(!is.finite(v)))
  if (any(bad_row)) {
    Z <- Z[!bad_row, , drop = FALSE]
  }
  Z
}

drop_near_zero_var <- function(X, freq_cut = 95 / 5, unique_cut = 10) {
  # caret::nearZeroVar logic without pulling in caret/trainControl
  if (ncol(X) < 2L) {
    return(X)
  }
  drop <- vapply(seq_len(ncol(X)), function(j) {
    x <- X[, j]
    if (all(!is.finite(x))) {
      return(TRUE)
    }
    u <- length(unique(x[is.finite(x)]))
    if (u <= unique_cut) {
      tab <- sort(table(x), decreasing = TRUE)
      if (length(tab) == 1L) {
        return(TRUE)
      }
      return((tab[[1]] / tab[[2]]) >= freq_cut)
    }
    FALSE
  }, logical(1))
  if (all(drop)) {
    return(X)
  }
  X[, !drop, drop = FALSE]
}

TRANSFORM_IDS <- c("clr_1", "clr_0.1", "rclr", "rclr_optspace")

rf_importance_one <- function(
    phy,
    transform,
    n_boot = 1000L,
    num_trees = 5000L,
    min_node_size = 5L,
    top_k = 20L,
    verbose = TRUE
) {
  X0 <- counts_samples_by_taxa(phy)
  X <- transform_counts(X0, transform)
  X <- drop_near_zero_var(X)
  if (ncol(X) < 2L) {
    stop("Fewer than 2 taxa remain after transform/NZV for ", transform)
  }

  samp_df <- as.data.frame(phyloseq::sample_data(phy))
  Y <- log1p(as.numeric(samp_df[rownames(X), "particles_total_d20"]))
  names(Y) <- rownames(X)
  if (any(!is.finite(Y))) {
    stop("Non-finite response values for transform ", transform)
  }

  n <- nrow(X)
  p <- ncol(X)
  mtry <- max(1L, floor(sqrt(p)))
  sample_ids <- rownames(X)

  if (isTRUE(verbose)) {
    message(sprintf(
      "  transform=%s | n=%d | p=%d | outer_boot=%d | trees=%d | mtry=%d (fixed)",
      transform, n, p, n_boot, num_trees, mtry
    ))
  }

  # Shared bootstrap index matrix (reproducible across transforms if seed set upstream)
  boot_idx <- matrix(NA_integer_, nrow = n_boot, ncol = n)
  for (b in seq_len(n_boot)) {
    set.seed(100000L + b)
    boot_idx[b, ] <- sample.int(n, size = n, replace = TRUE)
  }

  imps <- vector("list", n_boot)
  oob_diag <- vector("list", n_boot)

  for (b in seq_len(n_boot)) {
    if (isTRUE(verbose) && (b == 1L || b %% 50L == 0L || b == n_boot)) {
      message(sprintf("  bootstrap %d/%d", b, n_boot))
    }
    idx <- boot_idx[b, ]
    Xb <- X[idx, , drop = FALSE]
    Yb <- Y[idx]
    # Unique rownames for ranger when outer bootstrap duplicates samples
    rownames(Xb) <- make.unique(sample_ids[idx], sep = "__b")

    dat <- data.frame(y = Yb, Xb, check.names = FALSE)
    set.seed(200000L + b)
    fit <- ranger::ranger(
      dependent.variable.name = "y",
      data = dat,
      num.trees = num_trees,
      mtry = mtry,
      min.node.size = min_node_size,
      importance = "permutation",
      replace = TRUE,
      sample.fraction = 1,
      num.threads = num_cores,
      seed = 200000L + b,
      verbose = FALSE
    )

    vi <- fit$variable.importance
    imps[[b]] <- tibble::tibble(
      asv = names(vi),
      Overall = as.numeric(vi),
      boot = b
    )
    oob_diag[[b]] <- tibble::tibble(
      boot = b,
      oob_mse = fit$prediction.error,
      oob_r2 = fit$r.squared
    )
  }

  imp_df <- dplyr::bind_rows(imps) %>%
    dplyr::group_by(boot) %>%
    dplyr::mutate(rank = dplyr::dense_rank(dplyr::desc(Overall))) %>%
    dplyr::ungroup()

  oob_df <- dplyr::bind_rows(oob_diag)

  agg <- imp_df %>%
    dplyr::group_by(asv) %>%
    dplyr::summarise(
      mean_importance = mean(Overall, na.rm = TRUE),
      sd_importance = stats::sd(Overall, na.rm = TRUE),
      iqr_importance = stats::IQR(Overall, na.rm = TRUE),
      topk_freq = mean(rank <= top_k, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(mean_importance))

  tax <- as.data.frame(phyloseq::tax_table(phy)) %>%
    tibble::rownames_to_column("asv")
  if (!"Genus" %in% names(tax)) {
    tax$Genus <- NA_character_
  }

  out <- tax %>%
    dplyr::select(asv, Genus) %>%
    dplyr::inner_join(agg, by = "asv") %>%
    dplyr::arrange(dplyr::desc(mean_importance))

  list(
    transform = transform,
    n_samples = n,
    n_taxa = p,
    mtry = mtry,
    top_k = top_k,
    n_boot = n_boot,
    importance = out,
    boot_level = imp_df,
    oob_diagnostic = oob_df,
    oob_summary = tibble::tibble(
      mean_oob_mse = mean(oob_df$oob_mse, na.rm = TRUE),
      sd_oob_mse = stats::sd(oob_df$oob_mse, na.rm = TRUE),
      mean_oob_r2 = mean(oob_df$oob_r2, na.rm = TRUE),
      sd_oob_r2 = stats::sd(oob_df$oob_r2, na.rm = TRUE)
    )
  )
}

# Count how many transforms place each ASV in the top-K (by mean importance)
transform_topk_membership <- function(result_list, top_k = 20L) {
  sets <- lapply(result_list, function(res) {
    head(res$importance$asv, top_k)
  })
  names(sets) <- vapply(result_list, `[[`, character(1), "transform")
  all_asvs <- unique(unlist(sets, use.names = FALSE))
  memb <- tibble::tibble(asv = all_asvs) %>%
    dplyr::mutate(
      n_transforms_in_topk = vapply(
        asv,
        function(a) sum(vapply(sets, function(s) a %in% s, logical(1))),
        integer(1)
      ),
      transforms_in_topk = vapply(
        asv,
        function(a) paste(names(sets)[vapply(sets, function(s) a %in% s, logical(1))], collapse = ";"),
        character(1)
      )
    ) %>%
    dplyr::arrange(dplyr::desc(n_transforms_in_topk), asv)

  list(
    top_k = top_k,
    n_transforms = length(sets),
    n_per_transform = vapply(sets, length, integer(1)),
    sets = sets,
    membership = memb,
    # strict intersection kept as a convenience summary, not the primary metric
    intersection_asvs = Reduce(intersect, sets)
  )
}

run_marker <- function(phy, label) {
  phy_sub <- subset_week9_ws_rf(phy)
  filt <- filter_prevalence_subset(phy_sub, min_reads = min_reads, min_samples = min_samples)
  phy_f <- filt$phy

  message(sprintf(
    "%s RF subset: n_samples=%d | taxa %d -> %d (filter >=%d reads in >=%d samples, prevalence=%.1f%%) | dropped=%s",
    label, filt$n_samples, filt$n_taxa_before, filt$n_taxa_after,
    filt$min_reads, filt$min_samples, 100 * filt$prevalence,
    drop_label
  ))
  sd <- as(phyloseq::sample_data(phy_f), "data.frame")
  print(sd[, c("plastic_level", "particles_total_d20"), drop = FALSE])

  tag <- sprintf(
    "%s_9_ws_n%d_prev%dof%d",
    label, filt$n_samples, min_samples, filt$n_samples
  )

  results <- vector("list", length(TRANSFORM_IDS))
  names(results) <- TRANSFORM_IDS

  for (tr in TRANSFORM_IDS) {
    message(sprintf("=== %s | %s ===", label, tr))
    res <- rf_importance_one(
      phy_f,
      transform = tr,
      n_boot = n_boot,
      num_trees = num_trees,
      top_k = top_k,
      verbose = TRUE
    )
    results[[tr]] <- res
    path_i <- file.path(out_dir, sprintf("rf_%s_%s.rds", tag, tr))
    saveRDS(
      list(
        marker = label,
        filter = filt,
        drop_samples = drop_samples,
        transform = tr,
        design = list(
          outer = "bootstrap_samples",
          n_boot = n_boot,
          num_trees = num_trees,
          mtry = res$mtry,
          outer_cv = "none",
          importance = "ranger_permutation",
          oob = "diagnostic_only"
        ),
        result = res
      ),
      path_i
    )
    message("Wrote ", normalizePath(path_i, winslash = "/", mustWork = FALSE))
    message(sprintf(
      "  OOB diagnostic: mean R2=%.3f (SD=%.3f); mean MSE=%.4f",
      res$oob_summary$mean_oob_r2, res$oob_summary$sd_oob_r2, res$oob_summary$mean_oob_mse
    ))
  }

  sens <- transform_topk_membership(results, top_k = top_k)
  ref <- results[["clr_1"]]$importance
  memb_tbl <- sens$membership %>%
    dplyr::left_join(ref, by = "asv") %>%
    dplyr::arrange(
      dplyr::desc(n_transforms_in_topk),
      dplyr::desc(mean_importance)
    )

  summary_path <- file.path(out_dir, sprintf("rf_%s_transform_topk%d_membership.rds", tag, top_k))
  saveRDS(
    list(
      marker = label,
      filter = filt,
      drop_samples = drop_samples,
      transforms = TRANSFORM_IDS,
      top_k = top_k,
      n_boot = n_boot,
      transform_sensitivity = sens,
      membership_table = memb_tbl,
      per_transform_paths = file.path(out_dir, sprintf("rf_%s_%s.rds", tag, TRANSFORM_IDS))
    ),
    summary_path
  )
  message(sprintf(
    "%s top-%d transform membership: %d ASVs in >=1 transform; strict intersection n=%d",
    label, top_k, nrow(memb_tbl), length(sens$intersection_asvs)
  ))
  message("Wrote ", normalizePath(summary_path, winslash = "/", mustWork = FALSE))
  print(
    memb_tbl %>%
      dplyr::select(asv, Genus, n_transforms_in_topk, mean_importance, sd_importance, topk_freq, transforms_in_topk),
    n = 50
  )

  invisible(list(results = results, transform_sensitivity = sens, membership_table = memb_tbl))
}

if (dataset %in% c("16S", "BOTH")) {
  run_marker(phy_16S, "16S")
}
if (dataset %in% c("18S", "BOTH")) {
  run_marker(phy_18S, "18S")
}

message("Done.")
