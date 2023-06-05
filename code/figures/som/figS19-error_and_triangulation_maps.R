
# Setup workspace
rm(list = ls())
library(Racmacs)
library(patchwork)

# Read map
map <- read.acmap('data/maps/map_ndsubset_no_outliers_slope_adjusted.ace')

# Calculate triangulation blobs
triangulation <- triangulationBlobs(map, grid_spacing = 0.05, stress_lim = 1)

# Plot maps
maplims <- Racmacs:::mapPlotLims(map, sera = TRUE)
gp_triangulation <- ggplot(triangulation, xlim = maplims$xlim, ylim = maplims$ylim)

agSize(map) <- agSize(map)*0.5
srSize(map) <- srSize(map)*0.5
gp_error <- ggplot(map, xlim = maplims$xlim, ylim = maplims$ylim, show_error_lines = TRUE)

# Save the plot
gp <- gp_error + gp_triangulation + plot_annotation(tag_levels = "A") & theme(plot.tag.position = c(0.05, 0.95), plot.tag = element_text(size = 24))

ggsave(
  filename = "figures/som/figS19-error_and_triangulation_maps.png",
  plot = gp,
  width = 12,
  height = 5.8,
  units = "in"
)

ggsave(
  filename = "figures/som/figS19-error_and_triangulation_maps.pdf",
  plot = gp,
  width = 12,
  height = 5.8,
  units = "in"
)

# Save the map with triangulation blobs data
saveRDS(
  object = triangulation,
  file = "figures/som/figS19b-triangulation_blobs.rds"
)


