
# Setup workspace
rm(list = ls())
library(tidyverse)
library(Racmacs)
set.seed(100)

# source("functions/homologous_ags.R")
source("functions/adjust_titers_by_slope.R")
source("functions/set_pt_drawing_order.R")

# Set optimization parameters
num_optimizations <- 500

# Read in the full map
map_full_with_outliers <- read.acmap("data/maps/map_full.ace")
map_full_no_outliers <- read.acmap("data/maps/map_full_no_outliers.ace")

# Subset the map to sera with at least 3 detectable titers, but do include the single BA.2 serum (which has 2)
ndetectable <- colSums(titerTable(map_full_no_outliers) != "*" & !grepl("<", titerTable(map_full_no_outliers)))
sera_subset <- srNames(map_full_no_outliers)[ndetectable >= 3 | srGroups(map_full_no_outliers) == "BA.2"]

map_subset_no_outliers <- subsetMap(
  map_full_no_outliers, 
  sera = sera_subset
)

# Subset the map to sera with at least 3 detectable titers, but do include the single BA.2 serum (which has 2)
ndetectable <- colSums(titerTable(map_full_with_outliers) != "*" & !grepl("<", titerTable(map_full_with_outliers)))
sera_subset <- srNames(map_full_with_outliers)[ndetectable >= 3 | srGroups(map_full_with_outliers) == "BA.2"]

map_subset_with_outliers <- subsetMap(
  map_full_with_outliers, 
  sera = sera_subset
)

# Adjust the titers to account for slope differences
slope_data <- readRDS("data/generated_data/slope_comparison_wt.rds")$slope_factors %>%
  select(sr_group, estimate)

map_subset_no_outliers_slope_adjusted <- adjust_titers_by_slope(
  map = map_subset_no_outliers,
  slopes = slope_data
)

map_subset_with_outliers_slope_adjusted <- adjust_titers_by_slope(
  map = map_subset_with_outliers,
  slopes = slope_data
)

# Make a map with outliers from wildtype variants alone
map_subset_with_outliers_slope_adjusted %>%
  subsetMap(
    antigens = agGroups(map_subset_with_outliers_slope_adjusted) == "wildtype"
  ) %>% optimizeMap(
    number_of_dimensions = 2,
    number_of_optimizations = num_optimizations,
    minimum_column_basis = "none",
    options = list(ignore_disconnected = TRUE)
  ) -> map_subset_with_outliers_slope_adjusted_wt_only

# Make a map from wildtype variants alone
map_subset_no_outliers_slope_adjusted %>%
  subsetMap(
    antigens = agGroups(map_subset_no_outliers_slope_adjusted) == "wildtype"
  ) %>% optimizeMap(
    number_of_dimensions = 2,
    number_of_optimizations = num_optimizations,
    minimum_column_basis = "none",
    options = list(ignore_disconnected = TRUE)
  ) -> map_subset_no_outliers_slope_adjusted_wt_only

map_subset_no_outliers_slope_adjusted %>%
  subsetMap(
    antigens = agGroups(map_subset_no_outliers_slope_adjusted) == "wildtype"
  ) %>% optimizeMap(
    number_of_dimensions = 3,
    number_of_optimizations = num_optimizations,
    minimum_column_basis = "none",
    options = list(ignore_disconnected = TRUE)
  ) -> map_subset_no_outliers_slope_adjusted_wt_only_3d

# Make a map from wildtype variants plus d614g mutants
map_subset_no_outliers_slope_adjusted %>%
  subsetMap(
    antigens = agGroups(map_subset_no_outliers_slope_adjusted) %in% c("wildtype", "d614g_mutant")
  ) %>% optimizeMap(
    number_of_dimensions = 2,
    number_of_optimizations = num_optimizations,
    minimum_column_basis = "none",
    options = list(ignore_disconnected = TRUE)
  ) -> map_subset_no_outliers_slope_adjusted_wt_plus_d614g_mutants

map_subset_no_outliers_slope_adjusted %>%
  subsetMap(
    antigens = agGroups(map_subset_no_outliers_slope_adjusted) %in% c("wildtype", "d614g_mutant")
  ) %>% optimizeMap(
    number_of_dimensions = 3,
    number_of_optimizations = num_optimizations,
    minimum_column_basis = "none",
    options = list(ignore_disconnected = TRUE)
  ) -> map_subset_no_outliers_slope_adjusted_wt_plus_d614g_mutants_3d

# Make a map from wildtype variants plus 417 mutants
map_subset_no_outliers_slope_adjusted %>%
  subsetMap(
    antigens = agGroups(map_subset_no_outliers_slope_adjusted) %in% c("wildtype", "417_mutant")
  ) %>% optimizeMap(
    number_of_dimensions = 2,
    number_of_optimizations = num_optimizations,
    minimum_column_basis = "none",
    options = list(ignore_disconnected = TRUE)
  ) -> map_subset_no_outliers_slope_adjusted_wt_plus_417_mutants

map_subset_no_outliers_slope_adjusted %>%
  subsetMap(
    antigens = agGroups(map_subset_no_outliers_slope_adjusted) %in% c("wildtype", "417_mutant")
  ) %>% optimizeMap(
    number_of_dimensions = 3,
    number_of_optimizations = num_optimizations,
    minimum_column_basis = "none",
    options = list(ignore_disconnected = TRUE)
  ) -> map_subset_no_outliers_slope_adjusted_wt_plus_417_mutants_3d

# Make a map from wildtype variants BA.1+A484K mutant
map_subset_no_outliers_slope_adjusted %>%
  subsetMap(
    antigens = agGroups(map_subset_no_outliers_slope_adjusted) %in% c("wildtype", "BA1_mutant")
  ) %>% optimizeMap(
    number_of_dimensions = 2,
    number_of_optimizations = num_optimizations,
    minimum_column_basis = "none",
    options = list(ignore_disconnected = TRUE)
  ) -> map_subset_no_outliers_slope_adjusted_wt_plus_ba1_mutant

map_subset_no_outliers_slope_adjusted %>%
  subsetMap(
    antigens = agGroups(map_subset_no_outliers_slope_adjusted) %in% c("wildtype", "BA1_mutant")
  ) %>% optimizeMap(
    number_of_dimensions = 3,
    number_of_optimizations = num_optimizations,
    minimum_column_basis = "none",
    options = list(ignore_disconnected = TRUE)
  ) -> map_subset_no_outliers_slope_adjusted_wt_plus_ba1_mutant_3d


# Reorient the maps
mapTransformation(map_subset_no_outliers_slope_adjusted_wt_only) <- rbind(
  c(-0.6931977, -0.7207475),
  c(-0.7207475, 0.6931977)
)

map_subset_with_outliers_slope_adjusted_wt_only <- realignMap(map_subset_with_outliers_slope_adjusted_wt_only, map_subset_no_outliers_slope_adjusted_wt_only)
map_subset_no_outliers_slope_adjusted_wt_plus_d614g_mutants <- realignMap(map_subset_no_outliers_slope_adjusted_wt_plus_d614g_mutants, map_subset_no_outliers_slope_adjusted_wt_only)
map_subset_no_outliers_slope_adjusted_wt_plus_417_mutants   <- realignMap(map_subset_no_outliers_slope_adjusted_wt_plus_417_mutants, map_subset_no_outliers_slope_adjusted_wt_only)
map_subset_no_outliers_slope_adjusted_wt_plus_ba1_mutant    <- realignMap(map_subset_no_outliers_slope_adjusted_wt_plus_ba1_mutant, map_subset_no_outliers_slope_adjusted_wt_only)

map_subset_no_outliers_slope_adjusted_wt_only <- set_pt_drawing_order(map_subset_no_outliers_slope_adjusted_wt_only)
map_subset_with_outliers_slope_adjusted_wt_only <- set_pt_drawing_order(map_subset_with_outliers_slope_adjusted_wt_only)
map_subset_no_outliers_slope_adjusted_wt_plus_d614g_mutants <- set_pt_drawing_order(map_subset_no_outliers_slope_adjusted_wt_plus_d614g_mutants)
map_subset_no_outliers_slope_adjusted_wt_plus_417_mutants <- set_pt_drawing_order(map_subset_no_outliers_slope_adjusted_wt_plus_417_mutants)
map_subset_no_outliers_slope_adjusted_wt_plus_ba1_mutant <- set_pt_drawing_order(map_subset_no_outliers_slope_adjusted_wt_plus_ba1_mutant)

# Save the maps
save.acmap(keepSingleOptimization(map_subset_no_outliers_slope_adjusted_wt_only), "data/maps/map_ndsubset_no_outliers_slope_adjusted.ace")
save.acmap(keepSingleOptimization(map_subset_with_outliers_slope_adjusted_wt_only), "data/maps/map_ndsubset_with_outliers_slope_adjusted.ace")
save.acmap(keepSingleOptimization(map_subset_no_outliers_slope_adjusted_wt_plus_d614g_mutants), "data/maps/map_ndsubset_no_outliers_plus_d614g_mutants_slope_adjusted.ace")
save.acmap(keepSingleOptimization(map_subset_no_outliers_slope_adjusted_wt_plus_417_mutants), "data/maps/map_ndsubset_no_outliers_plus_417_mutants_slope_adjusted.ace")
save.acmap(keepSingleOptimization(map_subset_no_outliers_slope_adjusted_wt_plus_ba1_mutant), "data/maps/map_ndsubset_no_outliers_plus_ba1_mutant_slope_adjusted.ace")

save.acmap(keepSingleOptimization(map_subset_no_outliers_slope_adjusted_wt_only_3d), "data/maps/map_ndsubset_no_outliers_slope_adjusted_3d.ace")
save.acmap(keepSingleOptimization(map_subset_no_outliers_slope_adjusted_wt_plus_d614g_mutants_3d), "data/maps/map_ndsubset_no_outliers_plus_d614g_mutants_slope_adjusted_3d.ace")
save.acmap(keepSingleOptimization(map_subset_no_outliers_slope_adjusted_wt_plus_417_mutants_3d), "data/maps/map_ndsubset_no_outliers_plus_417_mutants_slope_adjusted_3d.ace")
save.acmap(keepSingleOptimization(map_subset_no_outliers_slope_adjusted_wt_plus_ba1_mutant_3d), "data/maps/map_ndsubset_no_outliers_plus_ba1_mutant_slope_adjusted_3d.ace")

