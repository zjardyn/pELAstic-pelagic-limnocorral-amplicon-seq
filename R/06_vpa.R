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
meta_data <- as(sample_data(phy2), "data.frame")
env1 <- meta_data %>% dplyr::select(week_9) %>% dplyr::pull()
env2 <- meta_data %>% dplyr::select(wall_strip) %>% dplyr::pull()

v3 <- "plastic_concentration"
env3 <- meta_data[[v3]]
env3 <- log10(env3 + 1)

vp <- varpart(otu_table, env1, env2, env3)
vp$part$indfract
vp$part

X1 <- meta_data %>% dplyr::select(week_9)
X2 <- meta_data %>% dplyr::select(wall_strip)
X3 <- meta_data %>% dplyr::select(plastic_concentration) %>%
  mutate(plastic_concentration = log1p(plastic_concentration))

# test time point alone
rda_X1 <- rda(otu_table ~ week_9 + Condition(wall_strip) + Condition(plastic_concentration),
              data = meta_data)
anova.cca(rda_X1, step = 1000)

# test location alone
rda_X2 <- rda(otu_table ~ wall_strip + Condition(week_9) + Condition(plastic_concentration),
              data = meta_data)
anova.cca(rda_X2, step = 1000)

# test log plastic concentration alone
rda_X3 <- rda(otu_table ~ plastic_concentration + Condition(week_9) + Condition(wall_strip),
              data = meta_data)
anova.cca(rda_X3, step = 1000)