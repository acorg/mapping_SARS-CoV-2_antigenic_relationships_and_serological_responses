
# Clear generated data directories
unlink("data/generated_data", recursive = TRUE)
unlink("data/maps", recursive = TRUE)
dir.create("data/generated_data")
dir.create("data/maps")

# Process the map data
message("==== Processing map data ====")
source("code/data_generation/process_map_data.R")

# Determining outliers
message("==== Determining outliers ====")
source("code/data_generation/determine_outliers.R")

# Calculating individual effect sizes
message("==== Calculating individual effect sizes ====")
source("code/data_generation/calculate_individual_effects.R")
source("code/data_generation/calculate_individual_effects_by_subgroup.R")

# Calculate slope comparisons
message("==== Calculating slope comparisons ====")
source("code/data_generation/slope_comparison.R")

# Calculate folddrop data
message("==== Calculating fold drop differences ====")
source("code/data_generation/calculate_folddrops.R")
source("code/data_generation/calculate_all_folddrops.R")

# Optimize the maps
message("==== Creating the antigenic maps ====")
source("code/data_generation/making_the_maps.R")

# Impute titers
message("==== Imputing titers ====")
source("code/data_generation/impute_titers.R")

# Bootstrap the map
message("==== Bootstrapping the map ====")
source("code/data_generation/bootstrap_noisy_map.R")
source("code/data_generation/bootstrap_resample_map.R")

# Do cross-validation on the maps
message("==== Cross-validating the map ====")
source("code/data_generation/cross_validate_map.R")

# Dimension test the maps
message("==== Dimension testing the map ====")
source("code/data_generation/dimension_test_map.R")

