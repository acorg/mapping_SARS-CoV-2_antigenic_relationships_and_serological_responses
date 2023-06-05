
# Setup workspace
rm(list = ls())
library(Racmacs)
library(tidyverse)
library(patchwork)
set.seed(300)

source("functions/scales.R")
source("functions/homologous_ags.R")

# Read the map
map <- read.acmap("data/maps/map_full_with_extras_no_outliers.ace")

# Get the titer table
ags <- agNames(map)
all_sr_groups <- srGroups(map)
sr_groups <- unique(srGroups(map))
titer_table <- titerTable(map)

# Setup table for folddrop data
expand.grid(
  sr_group = sr_groups,
  ag1 = ags,
  ag2 = ags
) %>%
  as_tibble() %>%
  mutate(
    ag1 = as.character(ag1),
    ag2 = as.character(ag2),
    folddiff = 0,
    folddiff_upper = 0,
    folddiff_lower = 0
  ) -> folddrop_data

# Calculate and add folddrop data
for (sr_group in sr_groups) {
  message(sr_group)
  for (ag1num in 1:length(ags)) {
    for (ag2num in ag1num:length(ags)) {
      if (ag1num != ag2num) {
        
        folddrop <- titertools::log2diff(
          titers1 = titer_table[ag1num, all_sr_groups == sr_group],
          titers2 = titer_table[ag2num, all_sr_groups == sr_group],
          sigma = 0.87,
          ci_method = "HDI",
          dilution_stepsize = dilutionStepsize(map)
        )
        
        # Add folddrop data
        folddrop_data$folddiff[
          folddrop_data$ag1 == ags[ag1num] &
          folddrop_data$ag2 == ags[ag2num] & 
          folddrop_data$sr_group == sr_group
        ] <- folddrop["mean", "estimate"]
        
        folddrop_data$folddiff_lower[
          folddrop_data$ag1 == ags[ag1num] &
            folddrop_data$ag2 == ags[ag2num] & 
            folddrop_data$sr_group == sr_group
        ] <- folddrop["mean", "lower"]
        
        folddrop_data$folddiff_upper[
          folddrop_data$ag1 == ags[ag1num] &
            folddrop_data$ag2 == ags[ag2num] & 
            folddrop_data$sr_group == sr_group
        ] <- folddrop["mean", "upper"]
        
        # Add reverse folddrop data
        folddrop_data$folddiff[
          folddrop_data$ag2 == ags[ag1num] &
            folddrop_data$ag1 == ags[ag2num] & 
            folddrop_data$sr_group == sr_group
        ] <- -folddrop["mean", "estimate"]
        
        folddrop_data$folddiff_lower[
          folddrop_data$ag2 == ags[ag1num] &
            folddrop_data$ag1 == ags[ag2num] & 
            folddrop_data$sr_group == sr_group
        ] <- -folddrop["mean", "lower"]
        
        folddrop_data$folddiff_upper[
          folddrop_data$ag2 == ags[ag1num] &
            folddrop_data$ag1 == ags[ag2num] & 
            folddrop_data$sr_group == sr_group
        ] <- -folddrop["mean", "upper"]
        
      }
    }
  }
}

# Save folddrop results
saveRDS(folddrop_data, "data/generated_data/all_folddrops.rds")
