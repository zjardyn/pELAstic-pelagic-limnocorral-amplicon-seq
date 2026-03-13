library(phyloseq)
library(magrittr)
suppressPackageStartupMessages(library(microViz))

files <- list.files("data", pattern = "\\.(rds|rdata|rda)$", full.names = TRUE, ignore.case = TRUE)
for(file in files){
    varname <- tools::file_path_sans_ext(basename(file))
    obj <- NULL
    ok <- FALSE
    try({
        obj <- readRDS(file)
        ok <- TRUE
    }, silent = TRUE)
    if (isTRUE(ok)) {
        assign(varname, obj)
    } else {
        env <- new.env()
        lr <- try(load(file, envir = env), silent = TRUE)
        if (!inherits(lr, "try-error")) {
            objs <- ls(env)
            if (length(objs) == 1) {
                assign(varname, get(objs[[1]], envir = env))
            } else if (length(objs) > 1) {
                for (nm in objs) {
                    assign(nm, get(nm, envir = env))
                }
            }
        } else {
            stop(paste("Failed to load:", file))
        }
    }
}
# Clean up
to_remove <- c("file", "files", "varname", "obj", "ok", "env", "lr", "objs")
rm(list = to_remove[to_remove %in% ls()])
# Reuse the robust loader to get phy_16S and phy_18S into the environment

sample_data(phy_16S)

# Helper to extract an ASV table with taxa as rows
get_asv_table <- function(phy) {
  phy <- ps_filter(phy, Location == "WS")
  otu_mat <- as(otu_table(phy), "matrix")
  if (!taxa_are_rows(phy)) {
    otu_mat <- t(otu_mat)
  }
  # merge WS_E_6 and WS_E_6-2 samples (sum counts) if both are present
  if (all(c("WS_E_6", "WS_E_6-2") %in% colnames(otu_mat))) {
    otu_mat[, "WS_E_6"] <- otu_mat[, "WS_E_6"] + otu_mat[, "WS_E_6-2"]
    otu_mat <- otu_mat[, setdiff(colnames(otu_mat), "WS_E_6-2"), drop = FALSE]
  }
  # merge WS_I_3 and WS_I_3-2 samples (sum counts) if both are present
  if (all(c("WS_I_3", "WS_I_3-2") %in% colnames(otu_mat))) {
    otu_mat[, "WS_I_3"] <- otu_mat[, "WS_I_3"] + otu_mat[, "WS_I_3-2"]
    otu_mat <- otu_mat[, setdiff(colnames(otu_mat), "WS_I_3-2"), drop = FALSE]
  }
  # order taxa by total abundance (most abundant first)
  ord <- order(rowSums(otu_mat), decreasing = TRUE)
  otu_mat <- otu_mat[ord, , drop = FALSE]
  asv_ids <- rownames(otu_mat)
  df <- as.data.frame(otu_mat, check.names = FALSE)
  df <- cbind(ASV = asv_ids, df)
  df
}

asv_16S <- get_asv_table(phy_16S)
asv_18S <- get_asv_table(phy_18S)

if (!dir.exists("output")) {
  dir.create("output", recursive = TRUE)
}

write.table(
  asv_16S,
  file = "output/asv_table_16S.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  asv_18S,
  file = "output/asv_table_18S.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)



