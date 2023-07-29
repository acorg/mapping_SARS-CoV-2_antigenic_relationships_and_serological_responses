
# Setup workspace
rm(list = ls())
library(tidyverse)
library(Racmacs)
library(ggtext)
source("functions/ag_designations.R")
source("functions/viewer_plotting.R")
source("functions/screenshot.R")

# Load the maps
map2d <- read.acmap("data/maps/map_ndsubset_no_outliers_slope_adjusted.ace")
map3d <- read.acmap("data/maps/map_ndsubset_no_outliers_slope_adjusted_3d.ace")

# Plot the 2D map
agSize(map2d) <- agSize(map2d)*0.7
srSize(map2d) <- srSize(map2d)*0.7
maplims <- Racmacs:::mapPlotLims(map2d, sera = FALSE, padding = 0.5)
srOutline(map2d) <- adjustcolor(srOutline(map2d), alpha.f = 0.4)
maplims$xlim[2] <- maplims$xlim[2]
gp <- ggplot(map2d, xlim = maplims$xlim, ylim = maplims$ylim, fill.alpha = NULL, outline.alpha = NULL)
ggsave("figures/main/fig3a-main_map.pdf", gp, width = 7.65, height = 6, units = "in")

# Plot the map with serum points hidden
map2d_nosera <- map2d
srShown(map2d_nosera) <- FALSE
gp <- ggplot(map2d_nosera, xlim = maplims$xlim, ylim = maplims$ylim, fill.alpha = NULL, outline.alpha = NULL)
ggsave("figures/summary/fig0A-main_map.pdf", gp, width = 7.65, height = 6, units = "in")

# Screenshot the 3d map
unlink("figures/main/fig3b-3d_map.png")
screenshotWebpage(
  url = "figures/main/maps/3d_map.html",
  file = "figures/main/fig3b-3d_map.png",
  vwidth = 538*2,
  vheight = 437*2
)


