source("R/01_load_files.R")
source("R/functions.R")

library(caret)
library(ranger)

num_cores <- parallel::detectCores() - 1

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

# Lean setup. For full setup, run with: 
# n_runs: 1000, 
# num_trees:  10,000, 
# resamples:  1000, 
# mtry_grid: seq(2, floor(p/3), by = max(1, floor(p/30)))
# min_node_grid: c(3, 5, 10, 20)


rf_importance <- function(
    phy,
    seeds = NULL,
    n_runs = 5,
    num_trees = 1000,
    resamples = 100,
    mtry_grid = NULL,              # e.g., NULL (auto), or c( sqrt(p)/2, sqrt(p), 2*sqrt(p) )
    min_node_grid = c(5),          # e.g., c(3,5,10)
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
    cl <- colnames(X)
    rl <- rownames(X)
    # rclr per sample
    X <- decostand(X, method = "rclr", MARGIN = 1)
    colnames(X) <- cl
    rownames(X) <- rl
    X <- as.matrix(X)
    # Target (regression on log1p particles) — subset via data.frame to avoid S4 zero-dimension issues
    samp_df <- as.data.frame(sample_data(phy))
    Y <- log1p(samp_df[rownames(X), "particles_total_d20", drop = TRUE]) %>% pull(particles_total_d20)

    ctrl <- trainControl(method = "boot", number = resamples, savePredictions = "final")

    # Set up tuning grid
    p <- ncol(X)
    if (is.null(mtry_grid)) {
        base <- max(1, floor(sqrt(p)))
        mtry_grid <- unique(pmax(1, floor(c(base/2, base, base*2))))
        mtry_grid <- mtry_grid[mtry_grid <= p]
    } else {
        mtry_grid <- pmax(1, pmin(p, floor(mtry_grid)))
    }
    tune_grid <- expand.grid(mtry = mtry_grid, splitrule = splitrule, min.node.size = min_node_grid)

    if (isTRUE(verbose)) {
        message(sprintf("rf_importance_5: n_runs=%d, num_trees=%d, resamples=%d, metric=%s",
                        n_runs, num_trees, resamples, metric))
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
            num.threads = num_cores
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

rf_runs_16S <- rf_importance(phy_16S_9_ws)
rf_runs_16S %>% as_tibble() %>% print(n = 50)
saveRDS(rf_runs_16S, file = "out/rf_runs_16S.rds")

rf_runs_18S <- rf_importance(phy_18S_9_ws)
rf_runs_18S %>% as_tibble() %>% print(n = 50)
saveRDS(rf_runs_18S, file = "out/rf_runs_18S.rds")
