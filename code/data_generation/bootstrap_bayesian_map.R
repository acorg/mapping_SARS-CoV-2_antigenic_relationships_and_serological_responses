
# Clean workspace
rm(list = ls())
set.seed(100)
library(Racmacs)

# Set variables
bootstrap_repeats <- 500
optimizations_per_repeat <- 500

# Load the map
map <- read.acmap('data/maps/map_ndsubset_no_outliers_slope_adjusted.ace')

# Resample bootstrap the map
bootstrap_map_ags <- bootstrapMap(
  map,
  method = "bayesian",
  bootstrap_repeats = bootstrap_repeats,
  bootstrap_ags = TRUE,
  bootstrap_sr = FALSE,
  reoptimize = TRUE,
  optimizations_per_repeat = optimizations_per_repeat
) |> bootstrapBlobs(gridspacing = 0.1)

bootstrap_map_srs <- bootstrapMap(
  map,
  method = "bayesian",
  bootstrap_repeats = bootstrap_repeats,
  bootstrap_ags = FALSE,
  bootstrap_sr = TRUE,
  reoptimize = TRUE,
  optimizations_per_repeat = optimizations_per_repeat
) |> bootstrapBlobs(gridspacing = 0.1)

bootstrap_map_ags_srs <- bootstrapMap(
  map,
  method = "bayesian",
  bootstrap_repeats = bootstrap_repeats,
  bootstrap_ags = TRUE,
  bootstrap_sr = TRUE,
  reoptimize = TRUE,
  optimizations_per_repeat = optimizations_per_repeat
) |> bootstrapBlobs(gridspacing = 0.1)

# Save the results
saveRDS(bootstrap_map_ags, "data/generated_data/map_bootstrap_bayesian_ags.rds")
saveRDS(bootstrap_map_srs, "data/generated_data/map_bootstrap_bayesian_srs.rds")
saveRDS(bootstrap_map_ags_srs, "data/generated_data/map_bootstrap_bayesian_ags_srs.rds")


