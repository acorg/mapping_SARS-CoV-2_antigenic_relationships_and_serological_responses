
# Setup workspace
rm(list = ls())
set.seed(100)
library(tidyverse)
library(Racmacs)
library(patchwork)

# Load functions
source("functions/homologous_ags.R")

# Read the map
map <- read.acmap("data/maps/map_full.ace")

# Exclude certain antigens from higher than homologous comparison
excluded_comp_ags <- c(
  "P.1", 
  "BA.1.1", "BA.2", "BA.2.12.1", "BA.3", "BA.4/BA.5",
  "B.1.617.2+K417N", "B.1.617.2 (AY.1)+K417N", "B.1.617.2 (AY.2)+K417N", "B.1.617.2 (AY.3)+E484Q",
  agNames(map)[agGroups(map) != "wildtype"] # Exclude mutants
)

# Determine those individuals with more than 2-fold higher than homologous titers
sr_homologous_logtiters <- unlist(srHomologousLogTiters(map))
sr_max_logtiters <- unname(apply(logtiterTable(map)[!agNames(map) %in% excluded_comp_ags,], 2, max, na.rm = TRUE))
excluded_sr_higher_than_homologous <- srNames(map)[sr_max_logtiters > sr_homologous_logtiters + 1]

# Set manually judged possible second infections
omicron_second <- c(
  "0177-02V8LC00-001",
  "0177-02VG6C00-001",
  "BA1-35",
  "BA1-38",
  "BA1-41"
)

other_second <- c(
  "0410-0GY9YA00-001",
  "Ser_051",
  "P1-0055",
  "B117-003",
  "204755836",
  "B1617-74"
)

manual_excluded_sera <- c(
  omicron_second,
  other_second
)

# Save the excluded sera
saveRDS(
  object = list(
    higher_than_homologous = excluded_sr_higher_than_homologous,
    manual = manual_excluded_sera
  ),
  file = "data/generated_data/outliers.rds"
)

# Create a map with outlier data indicated
possible_second <- c(excluded_sr_higher_than_homologous, manual_excluded_sera)
map_full_no_outliers <- subsetMap(
  map = map, 
  sera = !srNames(map) %in% possible_second
)

# Also remove outliers for map with extra mutant titrations
map_full_with_extras <- read.acmap("data/maps/map_full_with_extras.ace")
map_full_with_extras_no_outliers <- subsetMap(
  map = map_full_with_extras, 
  sera = !srNames(map_full_with_extras) %in% possible_second
)

save.acmap(map_full_no_outliers, "data/maps/map_full_no_outliers.ace")
save.acmap(map_full_with_extras_no_outliers, "data/maps/map_full_with_extras_no_outliers.ace")

