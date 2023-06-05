
# Setup workspace
rm(list = ls())
library(tidyverse)
library(Racmacs)
set.seed(100)

# Set variables
test_proportion <- 0.1
number_of_optimizations <- 500
number_of_replicates <- 500

# Read the map
map <- read.acmap("data/maps/map_ndsubset_no_outliers_slope_adjusted.ace")

# Function to cross-validate map predictions
crossvalidateMap <- function(
  map,
  test_proportion = 0.1,
  number_of_optimizations = 1000,
  number_of_replicates = 1000,
  optimization_number = 1,
  options = list()
) {
  
  # Perform the CV testing
  cv_results <- Racmacs:::runDimensionTestMap(
    map = map,
    dimensions_to_test = mapDimensions(map, optimization_number),
    test_proportion = test_proportion,
    minimum_column_basis = minColBasis(map, optimization_number),
    fixed_column_bases = colBases(map),
    number_of_optimizations = number_of_optimizations,
    replicates_per_dimension = number_of_replicates,
    options = options
  )
  
  # Summarise the results
  do.call(
    bind_rows,
    lapply(seq_along(cv_results$results), \(n) {
      tibble(
        measured_titer = cv_results$titers[cv_results$results[[n]]$test_indices],
        predicted_logtiter = cv_results$results[[n]]$predictions[[1]],
        titer_index = as.vector(cv_results$results[[n]]$test_indices),
        run = n
      )
    })
  ) %>% 
    mutate(
      ag_num = Racmacs:::agNumMatrix(map)[titer_index],
      sr_num = Racmacs:::srNumMatrix(map)[titer_index]
    )
  
}

results <- crossvalidateMap(
  map, 
  number_of_optimizations = number_of_optimizations,
  number_of_replicates = number_of_replicates,
  test_proportion = test_proportion
)

results %>%
  mutate(
    ag_name = agNames(map)[ag_num],
    sr_group = srGroups(map)[sr_num],
    measured_logtiter = Racmacs:::log_titers(measured_titer, 0),
    predicted_titer = as.character(round(2^predicted_logtiter*10, 6)),
    measured_titer_type = Racmacs:::titer_types_int(measured_titer),
    residual = measured_logtiter - predicted_logtiter
  ) -> results

saveRDS(results, "data/generated_data/cv_testing.rds")


