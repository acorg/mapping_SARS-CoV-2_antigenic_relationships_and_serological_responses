
library(Racmacs)
library(tidyverse)
map <- read.acmap("data/maps/map_full.ace")
map_no_outliers <- read.acmap("data/maps/map_full_no_outliers.ace")
map_lndscps <- read.acmap("data/maps/map_ndsubset_no_outliers_slope_adjusted.ace")
map_outliers <- subsetMap(map, sera = !srNames(map) %in% srNames(map_no_outliers))

# Summary function
summary_tibble <- function(x) {
  x <- summary(x)
  tibble::tibble(
    sr_group = names(x),
    n = x
  )
}

# Get number of second infections
num_excluded_sera <- numSera(map_outliers)
num_excluded_sera_by_group <- summary_tibble(droplevels(srGroups(map_outliers)))

# Merge serum groups
sr_groups <- as.character(srGroups(map))
srGroups(map) <- factor(sr_groups)

# Get additional numbers
num_sera <- numSera(map)
num_included_sera <- numSera(map_no_outliers)
num_excluded_sera <- num_sera - num_included_sera
num_sera_by_group <- summary_tibble(srGroups(map))
num_sera_by_group_lndscps <- summary_tibble(srGroups(map_lndscps))
num_sera_groups <- length(unique(srGroups(map)))
num_wt_variants <- sum(agGroups(map) == "wildtype")

outliers <- readRDS("data/generated_data/outliers.rds")
num_outliers_by_criteria <- length(outliers$higher_than_homologous)
num_outliers_by_manual_selection <- length(outliers$manual)
