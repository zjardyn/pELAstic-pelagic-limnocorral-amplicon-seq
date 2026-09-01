# One-off smoke: GBM/CZM/clr_0.01 transforms only (no RF).
suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(vegan)
})
stopifnot(requireNamespace("zCompositions", quietly = TRUE))
message("zCompositions ", as.character(packageVersion("zCompositions")))
source("R/01_load_files2.R")

subset_week9_ws_rf <- function(phy) {
  phy %>%
    phyloseq::subset_samples(as.character(Date) == "9") %>%
    phyloseq::subset_samples(as.character(Location) == "WS") %>%
    phyloseq::subset_samples(!is.na(particles_total_d20))
}
filter_prevalence_subset <- function(phy, min_reads = 3L, min_samples = 2L) {
  otu <- as(phyloseq::otu_table(phy), "matrix")
  if (!phyloseq::taxa_are_rows(phy)) {
    otu <- t(otu)
  }
  keep <- rowSums(otu >= min_reads) >= min_samples
  phy_f <- phyloseq::prune_taxa(keep, phy)
  phyloseq::prune_samples(phyloseq::sample_sums(phy_f) > 0, phy_f)
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
transform_zcomp <- function(X, method) {
  X_imp <- zCompositions::cmultRepl(
    X,
    label = 0,
    method = method,
    output = "p-counts",
    z.warning = 1,
    z.delete = FALSE,
    suppress.print = TRUE
  )
  X_imp <- as.matrix(X_imp)
  storage.mode(X_imp) <- "double"
  vegan::decostand(X_imp, method = "clr", MARGIN = 1L)
}

for (lab in c("16S", "18S")) {
  phy <- if (identical(lab, "16S")) phy_16S else phy_18S
  phy_f <- filter_prevalence_subset(subset_week9_ws_rf(phy))
  X <- counts_samples_by_taxa(phy_f)
  message(sprintf("--- %s n=%d p=%d ---", lab, nrow(X), ncol(X)))
  for (m in c("GBM", "CZM")) {
    Z <- transform_zcomp(X, m)
    message(sprintf(
      "  %s -> dim=%dx%d finite=%s range=[%.3f,%.3f]",
      m, nrow(Z), ncol(Z), all(is.finite(Z)), min(Z), max(Z)
    ))
  }
  Zc <- vegan::decostand(X, method = "clr", MARGIN = 1L, pseudocount = 0.01)
  message(sprintf(
    "  clr_0.01 -> dim=%dx%d finite=%s",
    nrow(Zc), ncol(Zc), all(is.finite(Zc))
  ))
}
message("SMOKE_OK")
