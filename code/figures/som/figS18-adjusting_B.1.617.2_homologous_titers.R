

# Setup workspace
rm(list = ls())
library(Racmacs)
library(tidyverse)
library(patchwork)
set.seed(300)

source("functions/scales.R")
source("functions/sr_group_labels.R")

# Read the map
map <- read.acmap("data/maps/map_ndsubset_no_outliers_slope_adjusted.ace")

# Adjust homologous titers down by four fold
titerTable(map)["B.1.617.2", srGroups(map) == "B.1.617.2"] <- Racmacs:::reactivity_adjust_titers(
  titers = titerTable(map)["B.1.617.2", srGroups(map) == "B.1.617.2"],
  adjustment = -2
)

map_lowhomologous <- optimizeMap(
  map = map,
  number_of_dimensions = 2,
  number_of_optimizations = 500,
  minimum_column_basis = "none"
)

map_lowhomologous <- realignMap(
  map_lowhomologous,
  map
)

gpA <- ggplot(procrustesMap(map, map_lowhomologous, antigens = FALSE))
gpB <- ggplot(procrustesMap(map, map_lowhomologous, sera = FALSE))

gp <- gpA + gpB + plot_annotation(tag_levels = "A") & theme(plot.tag.position = c(0, 1), plot.tag = element_text(size = 42, hjust = 0, vjust = 1, margin = margin(t = 0.5, l = 0.5, unit = "cm")))
ggsave(plot = gp, filename = "figures/som/figS18-adjusting_B.1.617.2_homologous.pdf", width = 20, height = 9, units = "in")
ggsave(plot = gp, filename = "figures/som/figS18-adjusting_B.1.617.2_homologous.png", width = 20, height = 9, units = "in")
