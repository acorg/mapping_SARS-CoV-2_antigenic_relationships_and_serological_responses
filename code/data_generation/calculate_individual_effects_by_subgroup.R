
# Setup workspace
rm(list = ls())
library(tidyverse)
library(ggplot2)
library(Racmacs)
library(posterior)

# Read the map
map <- read.acmap("data/maps/map_full_no_outliers.ace")

# Subset serum groups by source
srGroups(map) <- factor(paste(srGroups(map), srExtra(map), sep = "; "))

# Remove non wildtype antigens
map <- subsetMap(
  map = map,
  antigens = agGroups(map) == "wildtype"
)

# Get serum group data
individual_effects <- tibble(
  sr_name = character(0),
  sr_individual_effect = numeric(0)
)

ag_means_individual_effects <- tibble(
  sr_group = factor(NULL, levels(srGroups(map))),
  ag_name = character(0),
  ag_gmt_accounting_for_individual_effect = numeric(0)
)

# Calculate sr group individual effects
for (sr_group in levels(srGroups(map))) {
  
  message(sr_group)
  
  sr_group_sr_names <- srNames(map)[srGroups(map) == sr_group]
  titers <- titerTable(map)[,srGroups(map) == sr_group, drop = F]
  
  # Calculate the fit
  individual_effects_result <- titertools::estimate_sr_effects(
    titers = titers,
    dilution_stepsize = 0,
    sigma = 0.62,
    ci_method = "quap"
  )
  
  # Construct tables
  sr_group_individual_effects <- tibble(
    sr_name = sr_group_sr_names,
    sr_individual_effect = individual_effects_result[sprintf("sr_effects[%s]", seq_len(ncol(titers))), "estimate"]
  )
  
  sr_group_ag_means_individual_effects <- tibble(
    ag_name = agNames(map),
    sr_group = factor(sr_group, levels(srGroups(map))),
    ag_gmt_accounting_for_individual_effect = individual_effects_result[sprintf("ag_means[%s]", seq_len(nrow(titers))), "estimate"]
  )
  
  individual_effects <- bind_rows(
    individual_effects,
    sr_group_individual_effects
  )
  
  ag_means_individual_effects <- bind_rows(
    ag_means_individual_effects,
    sr_group_ag_means_individual_effects
  )
  
}

# Save the individual effects
saveRDS(individual_effects, file = "data/generated_data/individual_effects_by_subgroup.rds")
saveRDS(ag_means_individual_effects, file = "data/generated_data/ag_gmt_individual_effects_by_subgroup.rds")

