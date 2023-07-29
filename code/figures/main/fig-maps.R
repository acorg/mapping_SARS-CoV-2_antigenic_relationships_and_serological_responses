
# Setup workspace
rm(list = ls())
library(tidyverse)
library(Racmacs)
library(ggtext)
source("functions/ag_designations.R")
source("functions/viewer_plotting.R")

# Setup plotting directory
unlink("figures/main/maps", recursive = TRUE)
dir.create("figures/main/maps")
maps_dir <- "figures/main/maps"

# Load the maps
map2d <- read.acmap("data/maps/map_ndsubset_no_outliers_slope_adjusted.ace")
map3d <- read.acmap("data/maps/map_ndsubset_no_outliers_slope_adjusted_3d.ace")

# Plot the 3D map in different angles
map3d_viewer <- plot_map(
  map = procrustes_2d_to_3d(map2d, map3d, flip_z_axis = F, excluded_ags = "BA.4/BA.5"),
  scaling_size = 0.15,
  sr_scaling_size = 0.1,
  alter_opacity = F,
  sera_opacity = 0.2,
  grid.col = NA,
  # rotation = c(-1.2409, 0.0066, 0.1270),
  rotation = c(-1.2408, 0.0255, 0.1205),
  translation = c(-0.0799, 0.0153, 0.0057),
  zoom = 1.4,
  datafile = file.path(maps_dir, "3d_map.rds")
)

map3d_viewer_side <- plot_map(
  map = procrustes_2d_to_3d(map2d, map3d, flip_z_axis = F, excluded_ags = "BA.4/BA.5"),
  scaling_size = 0.15,
  sr_scaling_size = 0.1,
  alter_opacity = F,
  sera_opacity = 0.2,
  grid.col = NA,
  rotation = c(-1.5332, 1.5658, 1.7594),
  translation = c(-0.0091, 0.0009, -0.0802),
  zoom = 1.5
)

map3d_viewer_front <- plot_map(
  map = procrustes_2d_to_3d(map2d, map3d, flip_z_axis = F, excluded_ags = "BA.4/BA.5"),
  scaling_size = 0.15,
  sr_scaling_size = 0.1,
  alter_opacity = F,
  sera_opacity = 0.2,
  grid.col = NA,
  rotation = c(0, 0, 0),
  translation = c(0, 0, 0),
  zoom = 1.5
)

htmlwidgets::saveWidget(
  widget = map3d_viewer,
  file = file.path(maps_dir, "3d_map.html"),
  libdir = ".lib"
)

htmlwidgets::saveWidget(
  widget = map3d_viewer_side,
  file = file.path(maps_dir, "3d_map_side.html"),
  libdir = ".lib"
)

htmlwidgets::saveWidget(
  widget = map3d_viewer_front,
  file = file.path(maps_dir, "3d_map_front.html"),
  libdir = ".lib"
)

