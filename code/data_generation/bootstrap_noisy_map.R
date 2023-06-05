
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
  method = "noisy",
  bootstrap_repeats = bootstrap_repeats,
  ag_noise_sd = 0.4,
  titer_noise_sd = 0,
  reoptimize = TRUE,
  optimizations_per_repeat = optimizations_per_repeat
) |> bootstrapBlobs(gridspacing = 0.1)

bootstrap_map_titers <- bootstrapMap(
  map,
  method = "noisy",
  bootstrap_repeats = bootstrap_repeats,
  ag_noise_sd = 0,
  titer_noise_sd = 0.62,
  reoptimize = TRUE,
  optimizations_per_repeat = optimizations_per_repeat
) |> bootstrapBlobs(gridspacing = 0.1)

bootstrap_map_ags_titers <- bootstrapMap(
  map,
  method = "noisy",
  bootstrap_repeats = bootstrap_repeats,
  ag_noise_sd = 0.4,
  titer_noise_sd = 0.62,
  reoptimize = TRUE,
  optimizations_per_repeat = optimizations_per_repeat
) |> bootstrapBlobs(gridspacing = 0.1)

# Save the results
saveRDS(bootstrap_map_ags, "data/generated_data/map_bootstrap_noisy_ags.rds")
saveRDS(bootstrap_map_titers, "data/generated_data/map_bootstrap_noisy_titers.rds")
saveRDS(bootstrap_map_ags_titers, "data/generated_data/map_bootstrap_noisy_ags_titers.rds")


