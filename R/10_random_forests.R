source("R/01_load_files.R")
source("R/functions.R")

library(caret)
library(ranger)

num_cores <- 40
options(ranger.num.threads = num_cores)
# num_cores <- parallel::detectCores() - 1

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister_dopar()
# if the variable "total particles/m2" is used, remove NAs
phy_16S_9_ws <- phy_16S %>%
    subset_samples(Date == 9) %>%
    subset_samples(Location == "WS") %>%
    subset_samples(!is.na(particles_total_d20)) %>%
    ps_mutate(plastic_retained = if_else(particles_total_d20 > 10,  "retained", "not_retained")) 

phy_18S_9_ws <- phy_18S %>%
    subset_samples(Date == 9) %>%
    subset_samples(Location == "WS") %>%
    subset_samples(!is.na(particles_total_d20)) %>%
    ps_mutate(plastic_retained = if_else(particles_total_d20 > 10,  "retained", "not_retained"))

# For full setup, run with: 
# n_runs: 1000, 
# num_trees:  10,000, 

##
rf_importance <- function(
    phy,
    seeds = NULL,
    n_runs = 5,
    num_trees = 10000,
    # resamples = 100,
    mtry_grid = floor(sqrt(ncol(X))),              # e.g., NULL (auto), or c( sqrt(p)/2, sqrt(p), 2*sqrt(p) )
    min_node_grid = 5,          # p >> n  
    splitrule = "variance",        # regression default
    metric = "RMSE",               # regression metric
    verbose = TRUE                 # print progress and chosen params
) {
    X <- as.data.frame(as.matrix(otu_table(phy)))
    if (taxa_are_rows(phy)) {
        X <- t(X)
    }
    # Align to samples present in sample_data
    samp_df <- as.data.frame(sample_data(phy))
    keep_samps <- intersect(rownames(X), rownames(samp_df))
    X <- X[keep_samps, , drop = FALSE]
    # Drop degenerate samples/taxa before transform (columns only, keep all rows)
    X <- X[, colSums(X, na.rm = TRUE) > 0, drop = FALSE]
    # rclr per sample
    X <- decostand(X, method = "rclr", MARGIN = 1)
    X <- as.matrix(X)
    # Remove columns with any non-finite introduced by transform
    bad_col <- vapply(seq_len(ncol(X)), function(j) any(!is.finite(X[, j])), logical(1))
    if (any(bad_col)) X <- X[, !bad_col, drop = FALSE]
    # Remove near-zero variance taxa to stabilize RF (columns only)
    if (ncol(X) > 1) {
        nzv_idx <- caret::nearZeroVar(as.data.frame(X))
        if (length(nzv_idx) > 0) X <- X[, -nzv_idx, drop = FALSE]
    }
    # Remove non-finite rows introduced by transform to avoid empty bootstrap samples
    bad_row <- apply(X, 1, function(v) any(!is.finite(v)))
    if (any(bad_row)) X <- X[!bad_row, , drop = FALSE]
    # Target (regression on log1p particles) — build from data.frame column to avoid list coercion
    samp_df <- as.data.frame(sample_data(phy))
    Y <- log1p(samp_df[rownames(X), , drop = FALSE]$particles_total_d20)

    # caret bootstrap control (resampling inside each run)
    # ctrl <- trainControl(method = "boot", number = resamples, savePredictions = "final")
    ctrl <- trainControl(method = "cv", number = 5, savePredictions = "final")
    # Set up tuning grid: FIX mtry to sqrt(p)
    p <- ncol(X)
    base <- max(1, floor(sqrt(p)))
    tune_grid <- expand.grid(mtry = base, splitrule = splitrule, min.node.size = min_node_grid)

    if (isTRUE(verbose)) {
        message(sprintf("rf_importance_5: n_runs=%d, num_trees=%d, metric=%s",
                        n_runs, num_trees, metric))
        message(sprintf("mtry_grid=c(%s); min_node_grid=c(%s); splitrule=%s",
                        paste(mtry_grid, collapse=","), paste(min_node_grid, collapse=","), splitrule))
        message(sprintf("p=%d (taxa after preprocessing)", p))
    }

    # establish seeds sequence from n_runs if not provided
    if (is.null(seeds)) seeds <- seq_len(n_runs)
    imps <- lapply(seq_along(seeds), function(i) {
        s <- seeds[i]
        if (isTRUE(verbose)) message(sprintf("run %d/%d (seed=%s)", i, length(seeds), as.character(s)))
        set.seed(s)
        m <- train(
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
            sample.fraction = 1 # bootstrap sampling
        )
        if (isTRUE(verbose)) message(sprintf("  bestTune: mtry=%s, min.node.size=%s",
                                            as.character(m$bestTune$mtry), as.character(m$bestTune$min.node.size)))
        vi <- varImp(m)$importance
        vi$asv <- rownames(vi)
        vi$run <- i
        vi
    })

    imp_df <- bind_rows(imps)
    # Rank within each run (higher Overall = better rank 1)
    imp_df <- imp_df %>% group_by(run) %>% mutate(rank = dplyr::dense_rank(dplyr::desc(Overall))) %>% ungroup()
    top_k <- 20L
    agg <- imp_df %>%
        group_by(asv) %>%
        summarise(
            mean_importance = mean(Overall, na.rm = TRUE),
            sd_importance   = sd(Overall, na.rm = TRUE),
            iqr_importance  = IQR(Overall, na.rm = TRUE),
            topk_freq       = mean(rank <= top_k, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        arrange(desc(mean_importance))

    out <- tax_table(phy)[, "Genus"] %>%
        as.data.frame() %>%
        rownames_to_column("asv") %>%
        inner_join(agg, by = "asv") %>%
        arrange(desc(mean_importance))
    out
}

rf_runs_16S <- rf_importance(phy_16S_9_ws, n_runs = 1000)
rf_runs_16S %>% as_tibble() %>% print(n = 50)
saveRDS(rf_runs_16S, file = "output/rf_runs_1000_16S.rds")
rf_runs_18S <- rf_importance(phy_18S_9_ws, n_runs = 1000)
rf_runs_18S %>% as_tibble() %>% print(n = 50)
saveRDS(rf_runs_18S, file = "output/rf_runs_1000_18S.rds")


