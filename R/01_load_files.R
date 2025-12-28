library(phyloseq)
library(magrittr)
suppressPackageStartupMessages(library(microViz))


# Load all files in the data directory
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
# rm(to_remove, phy_dir)

# Filter and fix taxa
# Report taxa counts before filtering
before_16S_taxa <- if (exists("phy_16S")) phyloseq::ntaxa(phy_16S) else NA_integer_
before_18S_taxa <- if (exists("phy_18S")) phyloseq::ntaxa(phy_18S) else NA_integer_

phy_16S <- subset_taxa(phy_16S, Order != "Chloroplast")
phy_16S <- subset_taxa(phy_16S, Family != "Mitochondria")

after_16S_rm_chl_mito <- phyloseq::ntaxa(phy_16S)

phy_16S <- phy_16S %>%
    tax_filter(min_prevalence = 0.2, min_sample_abundance = 3) 

after_16S_tax_filter <- phyloseq::ntaxa(phy_16S)
phy_18S <- phy_18S %>%
    tax_filter(min_prevalence = 0.2, min_sample_abundance = 3)

# Report intermediate for 18S before other filters
after_18S_tax_filter <- phyloseq::ntaxa(phy_18S)

# Just don't do it!
# phy_16S <- phyloseq::rarefy_even_depth(phy_16S)
# phy_18S <- phyloseq::rarefy_even_depth(phy_18S)
phy_18S <- subset_taxa(phy_18S, Genus != "Incertae_sedis Family")

after_18S_rm_incertae <- phyloseq::ntaxa(phy_18S)

phy_18S <- phy_18S %>% tax_fix(unknowns = c("eukaryotic_picoplankton_environmental_sample", "uncultured", "uncultured_Eimeriidae", "uncultured_eukaryote", "uncultured_freshwater_eukaryote","Incertae_Sedis Family", "uncultured_Chlorophyta", "uncultured_fungus"),
                               anon_unique = TRUE)

# after_18S_tax_fix <- phyloseq::ntaxa(phy_18S)
# phy_ws_18S <- phy_ws_18S %>% tax_fix(unknowns = c("eukaryotic_picoplankton_environmental_sample", "uncultured", "uncultured_Eimeriidae", "uncultured_eukaryote", "uncultured_freshwater_eukaryote"))
# phy_ms_18S <- phy_ms_18S %>% tax_fix(unknowns = c("eukaryotic_picoplankton_environmental_sample", "uncultured", "uncultured_Eimeriidae", "uncultured_eukaryote", "uncultured_freshwater_eukaryote"))

phy_16S <- phy_16S %>% tax_fix(unknowns = c("Incertae Sedis", "endosymbionts"))

# Print concise summary
cat(sprintf(
    "16S taxa: before=%s, after_rm_chl_mito=%s, after_tax_filter=%s\n",
    before_16S_taxa, after_16S_rm_chl_mito, after_16S_tax_filter
))
cat(sprintf(
    "18S taxa: before=%s, after_tax_filter=%s, after_rm_incertae=%s\n",
    before_18S_taxa, after_18S_tax_filter, after_18S_rm_incertae 
))
