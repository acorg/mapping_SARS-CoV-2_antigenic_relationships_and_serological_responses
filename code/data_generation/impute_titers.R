
rm(list = ls())
library(Racmacs)
library(tidyverse)

source("functions/map_longinfo.R")
source("functions/diagnostics.R")
source("functions/imputation.R")

# Impute titers for 2D map
map <- read.acmap("data/maps/map_ndsubset_no_outliers_slope_adjusted.ace")
residual_data <- residualErrorTable(map)

# Add antigen and sera info
residual_data %>% 
  left_join(long_ag_info(map), by = "ag_num") %>%
  left_join(long_sr_info(map), by = "sr_num") -> residual_data

# Impute data
residual_data %>% 
  filter(
    measured_titer != "*"
  ) %>%
  group_by(
    ag_name,
    sr_group
  ) %>%
  group_modify(
    .f = impute_table_residuals
  ) %>%
  filter(
    residual_error_imputed != Inf
  ) -> residual_data

saveRDS(residual_data, 'data/generated_data/imputed_titers.rds')


# Impute titers for 3D map
map <- read.acmap("data/maps/map_ndsubset_no_outliers_slope_adjusted_3d.ace")
residual_data <- residualErrorTable(map)

# Add antigen and sera info
residual_data %>% 
  left_join(long_ag_info(map), by = "ag_num") %>%
  left_join(long_sr_info(map), by = "sr_num") -> residual_data

# Impute data
residual_data %>% 
  filter(
    measured_titer != "*",
    !is.na(predicted_logtiter)
  ) %>%
  group_by(
    ag_name,
    sr_group
  ) %>%
  group_modify(
    .f = impute_table_residuals
  ) %>%
  filter(
    residual_error_imputed != Inf
  ) -> residual_data

saveRDS(residual_data, 'data/generated_data/imputed_titers_3d.rds')


