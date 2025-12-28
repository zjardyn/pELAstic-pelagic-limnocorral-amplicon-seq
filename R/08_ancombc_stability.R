source("R/01_load_files.R")
source("R/functions.R")

# ANCOM-BC2 analysis trend tests with stability analysis 
# This should be run on a server 
# num_cores <- parallel::detectCores() - 1
num_cores <-  40 

# if the variable "total particles/m2" is used, remove NAs
cat("Running ANCOM-BC2 analysis on 16S and 18S ws data.\n")

# phy_16S_9_ws_levels <- phy_16S %>%
#     subset_samples(Date == 9) %>%
#     subset_samples(Location == "WS")
# cat("16S data subsetted.\n")

phy_18S_9_ws_levels <- phy_18S %>%
    subset_samples(Date == 9) %>%
    subset_samples(Location == "WS") 
cat("18S data subsetted.\n")

RNGkind("L'Ecuyer-CMRG")
rng_sampler <- function(){
    function() sample(100000:999999, 1)
}
rng_vec <- vapply(seq_len(100), function(i) rng_sampler()(), integer(1))

# Monotonic increasing: testing pattern none < low < medium < high
monotonic_increase_matrix <- matrix(c(
    1, 0, 0,     # First constraint: low-none > 0
    -1, 1, 0,    # Second constraint: medium-none > low-none => medium-low > 0 
    0, -1, 1     # Third constraint: high-none > medium-none => high-medium > 0
), nrow = 3, byrow = TRUE)

# Monotonic decreasing: testing pattern none > low > medium > high
monotonic_decrease_matrix <- matrix(c(
    -1, 0, 0,    # First constraint: low-none < 0
    1, -1, 0,    # Second constraint: medium-none < low-none => medium-low < 0
    0, 1, -1     # Third constraint: high-none < medium-none => high-medium < 0
), nrow = 3, byrow = TRUE)

# Umbrella shape (peak at medium): none < low < medium > high
umbrella_peak_matrix <- matrix(c(
    1, 0, 0,     # First constraint: low-none > 0 (increasing to peak)
    -1, 1, 0,    # Second constraint: medium-none > low-none => medium-low > 0
    0, 1, -1     # Third constraint: medium-none > high-none => medium-high > 0
), nrow = 3, byrow = TRUE)

run_ancombc <- function(rng_seed){
    set.seed(rng_seed)
    output_16S_trend <- ancombc2(
        data = phy_18S_9_ws_levels,
        tax_level = "Genus",
        fix_formula = "plastic_level", 
        rand_formula = NULL,
        p_adj_method = "holm",
        pseudo_sens = TRUE,
        prv_cut = 0,
        lib_cut = 0,
        s0_perc = 0.05,
        group = "plastic_level",
        struc_zero = FALSE,
        neg_lb = FALSE,
        alpha = 0.05,
        n_cl = num_cores,
        verbose = TRUE,
        global = TRUE,    # Global test
        trend = TRUE,     # Trend analysis
        trend_control = list(
            contrast = list(
                monotonic_increase_matrix,  # Increasing trend
                monotonic_decrease_matrix,  # Decreasing trend
                umbrella_peak_matrix        # Umbrella shape (peak at low)
            ),
            node = list(3, 3, 2),  # For increasing/decreasing use num_groups-1, for umbrella use position of peak
            solver = "ECOS",
            B = 100
        )
    )
}

list_output <- vector("list", length(rng_vec))
names(list_output) <- as.character(rng_vec)

for(i in 1:length(rng_vec)){
    cat("RNG seed:", rng_vec[i], "\n")
    cat("Iteration:", i, "\n")
    output <- run_ancombc(rng_vec[i])
    list_output[[as.character(rng_vec[i])]] <- output
}
saveRDS(list_output, "output/ancombc_stability_18S_9_ws_trend.rds")
