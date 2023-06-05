
# Setup workspace
rm(list = ls())
library(tidyverse)
library(Racmacs)

map <- read.acmap("data/maps/map_ndsubset_no_outliers_slope_adjusted.ace")
agSize(map) <- agSize(map)*0.7
srSize(map) <- srSize(map)*0.7

srOutline(map) <- adjustcolor(srOutline(map), alpha.f = 0.4)

gp <- ggplot(map, padding = 0.5) + theme(legend.position = "none")

ggsave(
  filename = "figures/som/figS22-full_map.pdf",
  plot = gp,
  width = 10,
  height = 9,
  units = "in"
)

ggsave(
  filename = "figures/som/figS22-full_map.png",
  plot = gp,
  width = 10,
  height = 9,
  units = "in"
)

