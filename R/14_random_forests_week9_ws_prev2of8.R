# Random forests (week-9 wall strip) — ASV-level transform sensitivity grid.
#
# Load: R/01_load_files2.R (no full-table tax_filter)
# Samples: Date == 9, Location == WS; drop bad none (default WS_H_9) → n = 8
# Prevalence (primary): >=3 reads in >=2 of the 8 samples (~25%)
#   Sensitivity: set RF_MIN_SAMPLES=3 for >=3 of 8
# Response: log1p(particles_total_d20)
# Features: ASVs (not genus-aggregated)
#
# Transform grid (1000 runs each by default):
#   clr_1, clr_0.5, clr_0.1  — vegan CLR with QUBS-matched pseudocounts
#   rclr                     — vegan rCLR without matrix completion (impute=FALSE)
#   rclr_optspace            — vegan rCLR with optspace completion (impute=TRUE;
#                              DEICODE-style; vegan 2.7+ default)
#
# Per transform: mean importance ± SD, topk_freq
# After all transforms: intersection of top-K ASV sets
#
# Server:
#   RF_N_CL=40 RF_DATASET=both Rscript R/14_random_forests_week9_ws_prev2of8.R
#   RF_N_CL=40 RF_DATASET=16S RF_MIN_SAMPLES=3 Rscript R/14_random_forests_week9_ws_prev2of8.R

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(caret)
  library(ranger)
  library(vegan)
  library(glue)
})

source("R/01_load_files2.R")

# --- runtime knobs ---
num_cores <- as.integer(Sys.getenv("RF_N_CL", unset = "40"))
n_runs <- as.integer(Sys.getenv("RF_N_RUNS", unset = "1000"))
num_trees <- as.integer(Sys.getenv("RF_NUM_TREES", unset = "10000"))
dataset <- toupper(Sys.getenv("RF_DATASET", unset = "both")) # 16S | 18S | both
min_reads <- as.integer(Sys.getenv("RF_MIN_READS", unset = "3"))
min_samples <- as.integer(Sys.getenv("RF_MIN_SAMPLES", unset = "2")) # primary=2; sens=3
top_k <- as.integer(Sys.getenv("RF_TOP_K", unset = "20"))
out_dir <- Sys.getenv("RF_OUT_DIR", unset = "output")
drop_samples <- strsplit(
  Sys.getenv("RF_DROP_SAMPLES", unset = "WS_H_9"),
  ",",
  fixed = TRUE
)[[1]] |>
  trimws()
drop_samples <- drop_samples[nzchar(drop_samples)]

stopifnot(dataset %in% c("16S", "18S", "BOTH"))
options(ranger.num.threads = num_cores)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list = ls(name = env), pos = env)
}
unregister_dopar()

message(sprintf(
  "RF week-9 WS | cores=%d | runs=%d | trees=%d | dataset=%s | filter=>=%d reads in >=%d samples | drop=%s | top_k=%d",
  num_cores, n_runs, num_trees, dataset, min_reads, min_samples,
  paste(drop_samples, collapse = ","), top_k
))

# --- prep ---
subset_week9_ws_rf <- function(phy, drop = drop_samples) {
  phy <- phy %>%
    phyloseq::subset_samples(as.character(Date) == "9") %>%
    phyloseq::subset_samples(as.character(Location) == "WS") %>%
    phyloseq::subset_samples(!is.na(particles_total_d20))
  keep <- setdiff(phyloseq::sample_names(phy), drop)
  phy <- phyloseq::prune_samples(keep, phy)
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

TRANSFORM_IDS <- c("clr_1", "clr_0.5", "clr_0.1", "rclr", "rclr_optspace")

rf_importance_one <- function(
    phy,
    transform,
    seeds = NULL,
    n_runs = 1000L,
    num_trees = 10000L,
    min_node_size = 5L,
    splitrule = "variance",
    metric = "RMSE",
    top_k = 20L,
    verbose = TRUE
) {
  X0 <- counts_samples_by_taxa(phy)
  X <- transform_counts(X0, transform)

  if (ncol(X) > 1L) {
    nzv_idx <- caret::nearZeroVar(as.data.frame(X))
    if (length(nzv_idx) > 0L) {
      X <- X[, -nzv_idx, drop = FALSE]
    }
  }
  if (ncol(X) < 2L) {
    stop("Fewer than 2 taxa remain after transform/NZV for ", transform)
  }

  samp_df <- as.data.frame(phyloseq::sample_data(phy))
  Y <- log1p(as.numeric(samp_df[rownames(X), "particles_total_d20"]))
  if (any(!is.finite(Y))) {
    stop("Non-finite response values for transform ", transform)
  }

  p <- ncol(X)
  mtry <- max(1L, floor(sqrt(p)))
  tune_grid <- expand.grid(
    mtry = mtry,
    splitrule = splitrule,
    min.node.size = min_node_size
  )
  ctrl <- caret::trainControl(method = "cv", number = 5L, savePredictions = "final")

  if (is.null(seeds)) {
    seeds <- seq_len(n_runs)
  }
  if (isTRUE(verbose)) {
    message(sprintf(
      "  transform=%s | n=%d | p=%d | runs=%d | mtry=%d",
      transform, nrow(X), p, length(seeds), mtry
    ))
  }

  imps <- lapply(seq_along(seeds), function(i) {
    s <- seeds[[i]]
    if (isTRUE(verbose) && (i == 1L || i %% 50L == 0L || i == length(seeds))) {
      message(sprintf("  run %d/%d (seed=%s)", i, length(seeds), as.character(s)))
    }
    set.seed(s)
    m <- caret::train(
      x = X,
      y = Y,
      method = "ranger",
      trControl = ctrl,
      tuneGrid = tune_grid,
      metric = metric,
      importance = "permutation",
      num.trees = num_trees,
      num.threads = num_cores,
      replace = TRUE,
      sample.fraction = 1
    )
    vi <- caret::varImp(m)$importance
    tibble::tibble(
      asv = rownames(vi),
      Overall = vi$Overall,
      run = i,
      seed = s
    )
  })

  imp_df <- dplyr::bind_rows(imps) %>%
    dplyr::group_by(run) %>%
    dplyr::mutate(rank = dplyr::dense_rank(dplyr::desc(Overall))) %>%
    dplyr::ungroup()

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
    n_samples = nrow(X),
    n_taxa = ncol(X),
    top_k = top_k,
    importance = out,
    run_level = imp_df
  )
}

intersect_topk <- function(result_list, top_k = 20L) {
  sets <- lapply(result_list, function(res) {
    head(res$importance$asv, top_k)
  })
  names(sets) <- vapply(result_list, `[[`, character(1), "transform")
  common <- Reduce(intersect, sets)
  list(
    top_k = top_k,
    n_per_transform = vapply(sets, length, integer(1)),
    n_intersection = length(common),
    asvs = common,
    sets = sets
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
    paste(drop_samples, collapse = ",")
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
      n_runs = n_runs,
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
        n_runs = n_runs,
        num_trees = num_trees,
        result = res
      ),
      path_i
    )
    message("Wrote ", normalizePath(path_i, winslash = "/", mustWork = FALSE))
  }

  ix <- intersect_topk(results, top_k = top_k)
  # Annotated intersection table using clr_1 importances when available
  ref <- results[["clr_1"]]$importance
  ix_tbl <- tibble::tibble(asv = ix$asvs) %>%
    dplyr::left_join(ref, by = "asv") %>%
    dplyr::arrange(dplyr::desc(mean_importance))

  summary_path <- file.path(out_dir, sprintf("rf_%s_intersection_topk%d.rds", tag, top_k))
  saveRDS(
    list(
      marker = label,
      filter = filt,
      drop_samples = drop_samples,
      transforms = TRANSFORM_IDS,
      top_k = top_k,
      intersection = ix,
      intersection_table = ix_tbl,
      per_transform_paths = file.path(out_dir, sprintf("rf_%s_%s.rds", tag, TRANSFORM_IDS))
    ),
    summary_path
  )
  message(sprintf(
    "%s intersection of top-%d across %d transforms: n=%d",
    label, top_k, length(TRANSFORM_IDS), ix$n_intersection
  ))
  message("Wrote ", normalizePath(summary_path, winslash = "/", mustWork = FALSE))
  if (nrow(ix_tbl) > 0L) {
    print(ix_tbl %>% dplyr::select(asv, Genus, mean_importance, sd_importance, topk_freq), n = 50)
  }

  invisible(list(results = results, intersection = ix, intersection_table = ix_tbl))
}

if (dataset %in% c("16S", "BOTH")) {
  run_marker(phy_16S, "16S")
}
if (dataset %in% c("18S", "BOTH")) {
  run_marker(phy_18S, "18S")
}

message("Done.")
