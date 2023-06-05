
# Setup workspace
rm(list = ls())
library(Racmacs)
library(tidyverse)
library(patchwork)
set.seed(300)

# Load functions
source("functions/homologous_ags.R")

# Read the map
map <- read.acmap("data/maps/map_full_no_outliers.ace")

# Get the titer table
titer_table <- titerTable(map)

# Setup to store results
all_group_results <- tibble(
  sr_group = character(0),
  ag = character(0),
  diff = numeric(0),
  diff_upper = numeric(0),
  diff_lower = numeric(0),
  homologous = logical(0)
)

# Append results for each group
for (sr_group in levels(srGroups(map))) {
  
  message(sr_group)
  sr <- which(srGroups(map) == sr_group)
  homo_ag_name <- homologous_ags[sr_group]
  homo_ag <- match(homo_ag_name, agNames(map))
  sr_group_results <- tibble(
    sr_group = rep(sr_group, numAntigens(map)),
    ag = character(numAntigens(map)),
    diff = numeric(numAntigens(map)),
    diff_upper = numeric(numAntigens(map)),
    diff_lower = numeric(numAntigens(map)),
    homologous = logical(numAntigens(map)),
    n = numeric(numAntigens(map))
  )
  
  for (ag in seq_len(numAntigens(map))) {
    homologous_titers <- titer_table[homo_ag, sr]
    ag_titers <- titer_table[ag, sr]
    titer_diff_est <- titertools::log2diff(
      titers1 = homologous_titers,
      titers2 = ag_titers,
      sigma = 0.87,
      ci_method = "HDI",
      dilution_stepsize = 0
    )
    
    # Populate results
    sr_group_results$ag[ag] <- agNames(map)[ag]
    if (ag == homo_ag) {
      sr_group_results$diff[ag] <- 0
      sr_group_results$diff_upper[ag] <- 0
      sr_group_results$diff_lower[ag] <- 0
      sr_group_results$homologous[ag] <- TRUE
    } else {
      sr_group_results$diff[ag] <- titer_diff_est["mean", "estimate"]
      sr_group_results$diff_upper[ag] <- titer_diff_est["mean", "upper"]
      sr_group_results$diff_lower[ag] <- titer_diff_est["mean", "lower"]
      sr_group_results$homologous[ag] <- FALSE
    }
    
    sr_group_results$n[ag] <- sum(
      homologous_titers != "*" & ag_titers != "*"
    )
    
  }
  
  # Append the results
  all_group_results <- bind_rows(all_group_results, sr_group_results)
  
}

# Save folddrop results
saveRDS(all_group_results, "data/generated_data/folddrops.rds")


