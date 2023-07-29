
rm(list = ls())
library(Racmacs)

map_no_outliers <- read.acmap("data/maps/map_ndsubset_no_outliers_slope_adjusted.ace")
map_with_outliers <- read.acmap('data/maps/map_ndsubset_with_outliers_slope_adjusted.ace')

# Get map limits
maplims <- Racmacs:::mapPlotLims(map_with_outliers, sera = F)

gp <- ggplot(
  procrustesMap(map_with_outliers, map_no_outliers, sera = FALSE),
  xlim = maplims$xlim,
  ylim = maplims$ylim
)

ggsave(
  plot = gp,
  filename = "figures/som/figS19-pc_map_with_outliers.pdf",
  width = 9,
  height = 8,
  units = "in"
)

ggsave(
  plot = gp,
  filename = "figures/som/figS19-pc_map_with_outliers.png",
  width = 9,
  height = 8,
  units = "in"
)
