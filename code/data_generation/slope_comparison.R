
# Setup workspace
rm(list = ls())
library(tidyverse)
library(Racmacs)
library(patchwork)
library(titertools)
set.seed(100)

# Read the map
map <- read.acmap("data/maps/map_full_no_outliers.ace")

# Remove non wildtype antigens
map <- subsetMap(
  map = map,
  antigens = agGroups(map) == "wildtype"
)

# Function to do a slope comparison
compare_slopes <- function(
    map,
    sr_groups_subset,
    homologous_ag
) {
  
  
  # Set sera subset
  sr_subset <- srGroups(map) %in% sr_groups_subset
  
  compare_crossreactivity_breadth(
    titers = titerTable(map)[,sr_subset],
    sr_groups = srGroups(map)[sr_subset],
    homologous_ag = homologous_ag,
    reference_sr_group = sr_groups_subset[1],
    dilution_stepsize = 0
  )
  
}

# Run model for WT-vaccinees
wt_result <- compare_slopes(
  map = map,
  sr_groups = c(
    "D614G",
    "2x mRNA-1273",
    "3x mRNA-1273 BD01",
    "3x mRNA-1273 BD29",
    "3x mRNA-1273 (6 month)"
  ),
  homologous_ag = "D614G"
)

# Run model for B.1.351 vaccinees
beta_result <- compare_slopes(
  map = map,
  sr_groups = c("B.1.351", "2x mRNA-1273.351"),
  homologous_ag = "B.1.351"
)

# Save results
saveRDS(wt_result, "data/generated_data/slope_comparison_wt.rds")
saveRDS(beta_result, "data/generated_data/slope_comparison_beta.rds")

