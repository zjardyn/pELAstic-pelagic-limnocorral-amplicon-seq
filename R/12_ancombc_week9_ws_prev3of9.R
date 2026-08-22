# ANCOM-BC2 on week-9 wall-strip samples with subset-specific prevalence filter.
#
# Load: R/01_load_files2.R (NO tax_filter — prevalence applied here after subsetting)
# Filter: >= min_reads in >= min_samples of the week-9 WS analysis set
#   default: >=3 reads in >=3 of 9 samples (33% prevalence)
# Design: plastic_level ~ none / low / medium / high (keep both "none" samples)
#
# Runs (per marker):
#   1) global test once
#   2) trend test n_trend times (default 100) with independent RNG seeds
#
# Server usage examples:
#   Rscript R/12_ancombc_week9_ws_prev3of9.R
#   ANCOMBC_N_CL=40 ANCOMBC_DATASET=16S Rscript R/12_ancombc_week9_ws_prev3of9.R
#   ANCOMBC_N_CL=40 ANCOMBC_DATASET=both ANCOMBC_N_TREND=100 Rscript R/12_ancombc_week9_ws_prev3of9.R

suppressPackageStartupMessages({
  library(phyloseq)
  library(ANCOMBC)
  library(dplyr)
})

source("R/01_load_files2.R")

# --- runtime knobs (env overrides for the server) ---
num_cores <- as.integer(Sys.getenv("ANCOMBC_N_CL", unset = "40"))
n_trend <- as.integer(Sys.getenv("ANCOMBC_N_TREND", unset = "100"))
dataset <- toupper(Sys.getenv("ANCOMBC_DATASET", unset = "both")) # 16S | 18S | both
min_reads <- as.integer(Sys.getenv("ANCOMBC_MIN_READS", unset = "3"))
min_samples <- as.integer(Sys.getenv("ANCOMBC_MIN_SAMPLES", unset = "3"))
out_dir <- Sys.getenv("ANCOMBC_OUT_DIR", unset = "output")
tax_level <- Sys.getenv("ANCOMBC_TAX_LEVEL", unset = "Genus")

stopifnot(dataset %in% c("16S", "18S", "BOTH"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message(sprintf(
  "ANCOMBC week-9 WS | cores=%d | trend_runs=%d | dataset=%s | filter=>=%d reads in >=%d samples | tax=%s",
  num_cores, n_trend, dataset, min_reads, min_samples, tax_level
))

# --- helpers ---
subset_week9_ws <- function(phy) {
  phy %>%
    phyloseq::subset_samples(as.character(Date) == "9") %>%
    phyloseq::subset_samples(as.character(Location) == "WS")
}

# Prevalence on the analysis subset only (not the full experiment table).
filter_prevalence_subset <- function(phy, min_reads = 3L, min_samples = 3L) {
  n_samp <- phyloseq::nsamples(phy)
  if (min_samples > n_samp) {
    stop(sprintf(
      "min_samples=%d exceeds n_samples=%d in analysis subset",
      min_samples, n_samp
    ))
  }
  otu <- as(phyloseq::otu_table(phy), "matrix")
  if (!phyloseq::taxa_are_rows(phy)) {
    otu <- t(otu)
  }
  keep <- rowSums(otu >= min_reads) >= min_samples
  phy_f <- phyloseq::prune_taxa(keep, phy)
  # drop empty samples if any (should not happen for week-9 WS)
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

ensure_plastic_level <- function(phy) {
  sd <- as(phyloseq::sample_data(phy), "data.frame")
  if (!"plastic_level" %in% names(sd)) {
    stop("sample_data missing plastic_level")
  }
  sd$plastic_level <- factor(
    as.character(sd$plastic_level),
    levels = c("none", "low", "medium", "high")
  )
  if (any(is.na(sd$plastic_level))) {
    stop("plastic_level has values outside none/low/medium/high")
  }
  phyloseq::sample_data(phy) <- phyloseq::sample_data(sd)
  phy
}

# Shared trend contrasts (none < low < medium < high)
monotonic_increase_matrix <- matrix(
  c(
    1, 0, 0,
    -1, 1, 0,
    0, -1, 1
  ),
  nrow = 3, byrow = TRUE
)
monotonic_decrease_matrix <- matrix(
  c(
    -1, 0, 0,
    1, -1, 0,
    0, 1, -1
  ),
  nrow = 3, byrow = TRUE
)
umbrella_peak_matrix <- matrix(
  c(
    1, 0, 0,
    -1, 1, 0,
    0, 1, -1
  ),
  nrow = 3, byrow = TRUE
)

trend_control <- list(
  contrast = list(
    monotonic_increase_matrix,
    monotonic_decrease_matrix,
    umbrella_peak_matrix
  ),
  node = list(3, 3, 2),
  solver = "ECOS",
  B = 100
)

run_ancombc2 <- function(phy, do_global = FALSE, do_trend = FALSE) {
  ANCOMBC::ancombc2(
    data = phy,
    tax_level = tax_level,
    fix_formula = "plastic_level",
    rand_formula = NULL,
    p_adj_method = "holm",
    pseudo_sens = TRUE,
    prv_cut = 0, # prevalence already applied on the analysis subset
    lib_cut = 0,
    s0_perc = 0.05,
    group = "plastic_level",
    struc_zero = FALSE,
    neg_lb = FALSE,
    alpha = 0.05,
    n_cl = num_cores,
    verbose = TRUE,
    global = isTRUE(do_global),
    pairwise = FALSE,
    dunnet = FALSE,
    trend = isTRUE(do_trend),
    trend_control = if (isTRUE(do_trend)) trend_control else NULL
  )
}

prep_one <- function(phy, label) {
  phy_sub <- subset_week9_ws(phy)
  filt <- filter_prevalence_subset(phy_sub, min_reads = min_reads, min_samples = min_samples)
  phy_f <- ensure_plastic_level(filt$phy)
  message(sprintf(
    "%s week-9 WS: n_samples=%d | taxa %d -> %d (filter >=%d reads in >=%d samples, prevalence=%.1f%%)",
    label, filt$n_samples, filt$n_taxa_before, filt$n_taxa_after,
    filt$min_reads, filt$min_samples, 100 * filt$prevalence
  ))
  print(table(phyloseq::sample_data(phy_f)$plastic_level))
  list(phy = phy_f, filter = filt)
}

run_marker <- function(phy, label) {
  prep <- prep_one(phy, label)
  phy_f <- prep$phy
  tag <- sprintf(
    "%s_9_ws_prev%dof%d",
    label, min_samples, phyloseq::nsamples(phy_f)
  )

  # 1) global once
  message(sprintf("=== %s GLOBAL (single run) ===", label))
  set.seed(1L)
  out_global <- run_ancombc2(phy_f, do_global = TRUE, do_trend = FALSE)
  global_path <- file.path(out_dir, sprintf("ancombc_global_%s.rds", tag))
  saveRDS(
    list(
      marker = label,
      filter = prep$filter,
      seed = 1L,
      result = out_global
    ),
    global_path
  )
  message("Wrote ", normalizePath(global_path, winslash = "/", mustWork = FALSE))

  # 2) trend 100x (independent seeds)
  message(sprintf("=== %s TREND (%d runs) ===", label, n_trend))
  RNGkind("L'Ecuyer-CMRG")
  set.seed(20260822L)
  rng_vec <- sample.int(899999L, n_trend, replace = FALSE) + 100000L

  list_trend <- vector("list", n_trend)
  names(list_trend) <- as.character(rng_vec)

  for (i in seq_len(n_trend)) {
    seed_i <- rng_vec[[i]]
    message(sprintf("%s trend %d/%d  seed=%d", label, i, n_trend, seed_i))
    set.seed(seed_i)
    list_trend[[as.character(seed_i)]] <- run_ancombc2(
      phy_f,
      do_global = FALSE,
      do_trend = TRUE
    )
    # checkpoint every 10 runs in case the job is killed
    if (i %% 10L == 0L || i == n_trend) {
      trend_path <- file.path(out_dir, sprintf("ancombc_trend_stability_%s.rds", tag))
      saveRDS(
        list(
          marker = label,
          filter = prep$filter,
          seeds = rng_vec,
          n_completed = i,
          results = list_trend[seq_len(i)]
        ),
        trend_path
      )
      message("Checkpoint wrote ", normalizePath(trend_path, winslash = "/", mustWork = FALSE))
    }
  }

  invisible(list(global = out_global, trend = list_trend, filter = prep$filter))
}

if (dataset %in% c("16S", "BOTH")) {
  run_marker(phy_16S, "16S")
}
if (dataset %in% c("18S", "BOTH")) {
  run_marker(phy_18S, "18S")
}

message("Done.")
