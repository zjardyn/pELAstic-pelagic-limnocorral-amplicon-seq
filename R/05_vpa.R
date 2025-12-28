source("R/01_load_files.R")
source("R/functions.R")

 # look at individual fractions
vpa_plot(phy_16S, week_9, wall_strip, label1 = "Date", label2 = "Location", output = "table")

# vpa_plot(phy_16S, week_9, wall_strip, "plastic_conc_log",
#          label1 = "Date", label2 = "Location", label3 = "Plastic Concentration", output = "plot")

vpa_plot(phy_18S, week_9, wall_strip, label1 = "Date", label2 = "Location", output = "table")

# vpa_plot(phy_18S, week_9, wall_strip, "plastic_conc_log",
#          label1 = "Date", label2 = "Location", label3 = "Plastic Concentration", output = "plot")
