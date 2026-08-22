source("R/01_load_files.R")
source("R/functions.R")

# phy <- phy_16S
phy <- phy_18S
phy2 <- phy %>%
    ps_mutate(week_3 = if_else(Date == "1", true = 1, false = 0),
                week_6 = if_else(Date == "2", true = 1, false = 0),
                week_9 = if_else(Date == "3", true = 1, false = 0),
                wall_strip = if_else(Location == "WS", true = 1, false = 0)) %>%
    tax_fix() %>%
    tax_transform("clr", rank = "Genus")

otu_table <- as(otu_table(phy2), "matrix")
otu_table <- t(otu_table)
meta_data <- as(sample_data(phy2), "data.frame") %>%
    mutate(log_plastic = log10(plastic_concentration + 1))

env1 <- meta_data$week_9
env2 <- meta_data$wall_strip
env3 <- meta_data$log_plastic

vp <- varpart(otu_table, env1, env2, env3)
vp$part$indfract
vp$part

# Unique fraction of time (week 9 vs weeks 3/6)
rda_X1 <- rda(
    otu_table ~ week_9 + Condition(wall_strip) + Condition(log_plastic),
    data = meta_data
)
anova.cca(rda_X1, permutations = 999)

# Unique fraction of location (wall strip vs microscope slide)
rda_X2 <- rda(
    otu_table ~ wall_strip + Condition(week_9) + Condition(log_plastic),
    data = meta_data
)
anova.cca(rda_X2, permutations = 999)

# Unique fraction of log10(plastic + 1)
rda_X3 <- rda(
    otu_table ~ log_plastic + Condition(week_9) + Condition(wall_strip),
    data = meta_data
)
anova.cca(rda_X3, permutations = 999)
