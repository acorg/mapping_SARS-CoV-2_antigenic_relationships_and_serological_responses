
# Setup workspace
rm(list = ls())
library(tidyverse)
library(Racmacs)
library(patchwork)

# Load functions
source("functions/scales.R")
source("functions/map_longinfo.R")
source("functions/sr_group_labels.R")
source("functions/ag_labels.R")
source("functions/impute_titers.R")
source("functions/plot_map_titers.R")
source("functions/plot_joint_titerplot.R")

# Read the map
map <- read.acmap("data/maps/map_full_no_outliers.ace")

# Read in serum group gmts
ag_means_individual_effects <- readRDS("data/generated_data/ag_gmt_individual_effects.rds")
sr_group_gmts <- ag_means_individual_effects %>% rename(logtiter = ag_gmt_accounting_for_individual_effect)

# Subset to not include mutants
map <- subsetMap(
  map, 
  antigens = agGroups(map) == "wildtype"
)

# Do a joint titerplot
plot_joint_titerplot(
  map = map,
  y_aes = "imputed_logtiter_minus_individual_effect",
  highlighted_sera = c(),
  highlighted_sera_cols = c(),
  sr_group_gmts = sr_group_gmts,
  sr_individual_effects = readRDS("data/generated_data/individual_effects.rds"),
  output_pdf = "figures/som/figS5-titerplot_minus_individual_effects.pdf",
  output_png = "figures/som/figS5-titerplot_minus_individual_effects.png"
)

