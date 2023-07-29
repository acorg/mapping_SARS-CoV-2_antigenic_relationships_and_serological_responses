
# Setup workspace
rm(list = ls())
set.seed(100)
library(tidyverse)
library(Racmacs)
library(patchwork)

# Load functions
source("functions/scales.R")
source("functions/map_longinfo.R")
source("functions/sr_group_labels.R")
source("functions/ag_labels.R")
source("functions/plot_map_titers.R")
source("functions/plot_joint_titerplot.R")

# Read the map
map <- read.acmap("data/maps/map_full.ace")

# Read in serum group gmts
ag_means_individual_effects <- readRDS("data/generated_data/ag_gmt_individual_effects.rds")
sr_group_gmts <- ag_means_individual_effects %>% rename(logtiter = ag_gmt_accounting_for_individual_effect)

# Read in outlier data
outliers <- readRDS("data/generated_data/outliers.rds")

# Drop levels referring to potential second infections
srGroups(map) <- droplevels(
  factor(
    gsub(" ?2nd", "", srGroups(map), fixed = T),
    levels(srGroups(map))
  )
)

# Do the joint titerplot
plot_joint_titerplot(
  map = map,
  y_aes = "logtiter",
  highlighted_sera = list(
    outliers$higher_than_homologous, 
    outliers$manual
  ),
  highlighted_sera_col = c("red", "blue"),
  sr_group_gmts = sr_group_gmts,
  plot_sr_group_gmts = FALSE,
  output_pdf = "figures/som/figS2-outlier_titerplot.pdf",
  output_png = "figures/som/figS2-outlier_titerplot.png"
)

