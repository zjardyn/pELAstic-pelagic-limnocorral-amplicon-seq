# RF figures with LOO importance as PRIMARY ranking.
#
# Workflow reminder:
#   - Fixed ranger params (ntree=5000, mtry=sqrt(p), min.node.size=5)
#   - Leave-one-out CV: each fold fits on n-1 samples
#   - Primary taxon score = mean permutation importance across LOO folds
#   - Full-n importance kept only as secondary
#   - Consensus: mean_loo_importance > 0.001 across transforms
#
# Writes under figures/rf_loocv_importance/loo_primary/{16S,18S}/
#   - importance rank PDFs (per marker, all transforms)
#   - Jaccard matrices (CSV + PDF) for LOO-imp > 0.001 sets
#   - family/genus line plots (all4 + all5) via R/18–R/21
#
# Run: Rscript R/22_rf_loo_primary_figures.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(glue)
})

out_dir <- Sys.getenv(
  "RF_LOO_PRIMARY_OUT_DIR",
  unset = "figures/rf_loocv_importance/loo_primary"
)
rf_dir <- Sys.getenv("RF_OUT_DIR", unset = "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
for (mk in c("16S", "18S")) {
  dir.create(file.path(out_dir, mk), showWarnings = FALSE, recursive = TRUE)
}
marker_dir <- function(marker) file.path(out_dir, marker)

# Same 8-transform grid as R/14 / R/23 (line plots R/18–21 still use 4/5-set scripts)
TRANSFORM_IDS <- c(
  "clr_1", "clr_0.5", "clr_0.1", "clr_0.01",
  "rclr", "rclr_optspace",
  "gbm", "czm"
)
IMP_THRESH <- 0.001

rf_path <- function(marker, transform) {
  file.path(
    rf_dir,
    sprintf("rf_%s_9_ws_n9_prev2of9_%s.rds", marker, transform)
  )
}

load_imp <- function(marker, transform) {
  x <- readRDS(rf_path(marker, transform))
  imp <- x$result$importance
  stopifnot("mean_loo_importance" %in% names(imp))
  imp %>%
    mutate(
      transform = transform,
      mean_importance = mean_loo_importance
    ) %>%
    arrange(desc(mean_loo_importance))
}

# --- Importance rank plots (LOO primary) ---
make_imp_rank_plot <- function(imp, title, tag) {
  df <- imp %>%
    mutate(
      rank = row_number(),
      ymin = pmax(0, mean_loo_importance - sd_importance),
      ymax = mean_loo_importance + sd_importance
    )
  ggplot(df, aes(x = rank, y = mean_loo_importance)) +
    geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "grey80", alpha = 0.8) +
    geom_line(linewidth = 0.4, colour = "grey30") +
    labs(
      title = title,
      tag = tag,
      x = "Rank",
      y = "Mean LOO permutation importance ± SD"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.tag = element_text(face = "bold", size = 16),
      plot.title = element_text(size = 10)
    )
}

for (marker in c("16S", "18S")) {
  plots <- vector("list", length(TRANSFORM_IDS))
  names(plots) <- TRANSFORM_IDS
  tags <- LETTERS[seq_along(TRANSFORM_IDS)]
  for (i in seq_along(TRANSFORM_IDS)) {
    tr <- TRANSFORM_IDS[[i]]
    imp <- load_imp(marker, tr)
    plots[[i]] <- make_imp_rank_plot(
      imp,
      title = glue("{marker} | {tr} | n={nrow(imp)} ASVs"),
      tag = tags[[i]]
    )
    # per-transform top-50 bars
    top <- imp %>%
      slice_head(n = 50) %>%
      mutate(
        label = ifelse(is.na(Genus) | Genus == "", substr(asv, 1, 8), Genus),
        label = make.unique(as.character(label)),
        ymin = pmax(0, mean_loo_importance - sd_importance),
        ymax = mean_loo_importance + sd_importance
      )
    p_bar <- ggplot(top, aes(x = reorder(label, mean_loo_importance), y = mean_loo_importance)) +
      geom_col(fill = "grey40") +
      geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.25) +
      coord_flip() +
      labs(
        title = glue("{marker} {tr}: top 50 by mean LOO importance"),
        x = NULL,
        y = "Mean LOO permutation importance ± SD"
      ) +
      theme_bw(base_size = 10)
    ggsave(
      file.path(marker_dir(marker), sprintf("rf_imp_%s_%s_top50.pdf", marker, tr)),
      p_bar,
      width = 8,
      height = 10,
      useDingbats = FALSE
    )
  }
  p_all <- wrap_plots(plots, nrow = 1)
  ggsave(
    file.path(marker_dir(marker), sprintf("rf_imp_%s_all_transforms_all_taxa.pdf", marker)),
    p_all,
    width = 28,
    height = 4.5,
    useDingbats = FALSE
  )
  message("Wrote importance ranks for ", marker)
}

# --- Jaccard of LOO-imp > threshold sets ---
jaccard <- function(a, b) {
  a <- unique(a)
  b <- unique(b)
  inter <- length(intersect(a, b))
  uni <- length(union(a, b))
  if (uni == 0) return(NA_real_)
  inter / uni
}

for (marker in c("16S", "18S")) {
  sets <- lapply(TRANSFORM_IDS, function(tr) {
    imp <- load_imp(marker, tr)
    imp$asv[imp$mean_loo_importance > IMP_THRESH]
  })
  names(sets) <- TRANSFORM_IDS
  mat <- outer(
    TRANSFORM_IDS,
    TRANSFORM_IDS,
    Vectorize(function(i, j) jaccard(sets[[i]], sets[[j]]))
  )
  dimnames(mat) <- list(TRANSFORM_IDS, TRANSFORM_IDS)
  csv_path <- file.path(marker_dir(marker), sprintf("jaccard_%s_loo_imp_gt0.001.csv", marker))
  write.csv(mat, csv_path)
  df <- as.data.frame(as.table(mat)) %>%
    rename(transform_a = Var1, transform_b = Var2, jaccard = Freq)
  p_jac <- ggplot(df, aes(transform_a, transform_b, fill = jaccard)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", jaccard)), size = 3) +
    scale_fill_viridis_c(limits = c(0, 1), name = "Jaccard") +
    labs(
      title = glue("{marker}: Jaccard of ASVs with mean LOO importance > {IMP_THRESH}"),
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(
    file.path(marker_dir(marker), sprintf("jaccard_%s_loo_imp_gt0.001.pdf", marker)),
    p_jac,
    width = 11,
    height = 9,
    useDingbats = FALSE
  )
  message(
    marker, " set sizes LOO>0.001: ",
    paste(sprintf("%s=%d", names(sets), lengths(sets)), collapse = ", ")
  )
}

# --- Family / genus line plots into the same folder ---
Sys.setenv(CLR_FAMILY_OUT_DIR = normalizePath(out_dir, mustWork = TRUE))
message("CLR_FAMILY_OUT_DIR=", Sys.getenv("CLR_FAMILY_OUT_DIR"))

plot_scripts <- c(
  "R/18_family_lineplots_all4_no_optspace.R",
  "R/19_genus_lineplots_all4_no_optspace.R",
  "R/20_genus_lineplots_all5_gt0.001.R",
  "R/21_family_lineplots_all5_gt0.001.R"
)
for (sc in plot_scripts) {
  message("=== Running ", sc, " ===")
  # Isolate each script so sourced helpers don't collide
  status <- system2(
    "Rscript",
    args = c(sc),
    env = c(
      sprintf("CLR_FAMILY_OUT_DIR=%s", normalizePath(out_dir, mustWork = TRUE)),
      sprintf("RF_OUT_DIR=%s", normalizePath(rf_dir, mustWork = TRUE)),
      sprintf("PATH=%s", Sys.getenv("PATH"))
    )
  )
  if (!identical(as.integer(status), 0L)) {
    stop("Plot script failed: ", sc, " (status=", status, ")")
  }
}

message("Done. All LOO-primary RF figures in: ", normalizePath(out_dir, mustWork = TRUE))
