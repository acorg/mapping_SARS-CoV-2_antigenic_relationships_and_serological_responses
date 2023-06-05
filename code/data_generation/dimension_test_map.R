
# Clean workspace
rm(list = ls())
set.seed(100)
library(Racmacs)

# Load the map
map <- read.acmap('data/maps/map_ndsubset_no_outliers_slope_adjusted.ace')

# Run the dimensionality testing
dimtest_results <- dimensionTestMap(
  map,
  dimensions_to_test = 1:5,
  test_proportion = 0.1,
  minimum_column_basis = "none",
  fixed_column_bases = colBases(map),
  number_of_optimizations = 500,
  replicates_per_dimension = 200
)

# Save the results
saveRDS(dimtest_results, "data/generated_data/dimtest_results.rds")
