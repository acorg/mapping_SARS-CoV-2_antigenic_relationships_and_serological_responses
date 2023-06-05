
# Setup workspace
rm(list = ls())
library(tidyverse)
library(Racmacs)
library(ablandscapes)
library(r3js)
source("functions/fit_lndscp_val.R")
source("functions/variables.R")

# Setup plotting directory
unlink("figures/main/lndscps", recursive = TRUE)
dir.create("figures/main/lndscps")
lndscps_dir <- "figures/main/lndscps"

# Load map
map <- read.acmap("data/maps/map_ndsubset_no_outliers_slope_adjusted.ace")

# Load gmts after accounting for individual variation
ag_means_individual_effects <- readRDS("data/generated_data/ag_gmt_individual_effects.rds")

# Set limits
plot_lims <- Racmacs:::mapPlotLims(map, padding = 0.5, sera = FALSE)
plot_xlim <- plot_lims$xlim
plot_ylim <- plot_lims$ylim
plot_zlim <- c(0, 10)

# Work out which sera are out of bounds
sr_x_coords <- srCoords(map)[,1]
sr_y_coords <- srCoords(map)[,2]
margin <- 0.2
sr_out_of_bound <- sr_x_coords < plot_xlim[1] + margin | 
  sr_x_coords > plot_xlim[2] - margin | 
  sr_y_coords < plot_ylim[1] + margin |
  sr_y_coords > plot_ylim[2] - margin

lndscp_xlim <- range(agCoords(map)[,1])
lndscp_ylim <- range(agCoords(map)[,2])

# Set viewing angles
angle <- list(
  rotation = c(-1.3934, 0.0020, -0.2476),
  translation = c(0.0459, 0.0303, 0.0745),
  zoom = 1.4
)

# Function to plot serum group landscape
plot_sr_group_lndscp <- function(
    sr_group, 
    sr_slope = 1,
    sr_group_color = NULL,
    add_mean_surface = TRUE,
    add_individual = TRUE,
    add_titers = TRUE,
    titer_plane = log2(5),
    return_mean_surface = FALSE,
    return_data3js = FALSE,
    fix_surface = NA,
    filename
  ) {
  
  # Set sera group details
  if (is.null(sr_group_color)) sr_group_color <- unique(srOutline(map)[srGroups(map) == sr_group])
  
  sr_group_logtiters <- logtiterTable(map)[ ,srGroups(map) == sr_group, drop = F]
  sr_group_mean_logtiters <- rowMeans(sr_group_logtiters, na.rm = T)
  sr_group_ag_means_individual_effects <- ag_means_individual_effects[ag_means_individual_effects$sr_group == sr_group,]
  sr_group_mean_logtiters <- sr_group_ag_means_individual_effects$ag_gmt_accounting_for_individual_effect[match(agNames(map), sr_group_ag_means_individual_effects$ag_name)]
  sr_group_mean_logtiters[sr_group_mean_logtiters < plot_zlim[1]] <- plot_zlim[1]
  
  sr_group_coords <- srCoords(map)[srGroups(map) == sr_group,,drop = F]
  sr_group_colbases <- colBases(map)[srGroups(map) == sr_group]
  
  # Set map subset
  map_subset <- map
  srShown(map_subset)[sr_out_of_bound] <- FALSE
  srShown(map_subset)[srGroups(map_subset) != sr_group] <- FALSE
  agSize(map_subset) <- 0.8*agSize(map_subset)
  
  # Plot the base plot
  data3js <- lndscp_r3jsplot(
    fit = list(acmap = map_subset),
    aspect.z = 0.5,
    show.surface = FALSE,
    show.titers = FALSE,
    output = "data3js",
    xlim = plot_xlim,
    ylim = plot_ylim,
    zlim = plot_zlim,
    show.sidegrid = TRUE,
    show.axis = FALSE,
    options = list(
      opacity.basemap.ags = 1,
      cex.basemap.ags = 1,
      cex.basemap.sr = 1.5,
      lwd.grid             = 1.5,
      opacity.basemap.sr = 1,
      sidegrid.lwd = 2,
      sidegrid.at = list(z = 1:floor(plot_zlim[2])),
      sidegrid.col = "grey95"
    )
  )
  
  # Get fitted surface
  grid_x_coords <- seq(from = lndscp_xlim[1], to = lndscp_xlim[2], by = 0.25)
  grid_y_coords <- seq(from = lndscp_ylim[1], to = lndscp_ylim[2], by = 0.25)
  grid_x_matrix <- matrix(grid_x_coords, length(grid_y_coords), length(grid_x_coords), byrow = T)
  grid_y_matrix <- matrix(grid_y_coords, length(grid_y_coords), length(grid_x_coords), byrow = F)
  grid_z_matrix <- matrix(NA, length(grid_y_coords), length(grid_x_coords))
  
  grid_z_matrix[] <- vapply(
    seq_len(length(grid_z_matrix)),
    \(n) {
      fit_lndscp_val(
        x = grid_x_matrix[n], 
        y = grid_y_matrix[n],
        sr_coords = sr_group_coords,
        sr_colbases = sr_group_colbases,
        sr_slope = sr_slope
      )
    }, numeric(1)
  )
  
  # Fix z matrix to a given max point
  if (!is.na(fix_surface)) {
    # max_z <- max(grid_z_matrix)
    max_z <- fit_lndscp_val(
      x = agCoords(map)["D614G", 1], 
      y = agCoords(map)["D614G", 2],
      sr_coords = sr_group_coords,
      sr_colbases = sr_group_colbases,
      sr_slope = sr_slope
    )
    grid_z_matrix[] <- grid_z_matrix + (fix_surface - max_z)
  }
  
  # Add the surface
  if (return_mean_surface) data3js <- r3js::plot3js.new()
  if (add_mean_surface) {
    
    data3js <- r3js::surface3js(
      data3js,
      x = grid_x_matrix,
      y = grid_y_matrix,
      z = grid_z_matrix,
      col = sr_group_color,
      opacity = 0.8,
      toggle = "Mean landscape",
      wireframe = FALSE,
      doubleSide = TRUE
    )
    
    data3js <- r3js::surface3js(
      data3js,
      x = grid_x_matrix,
      y = grid_y_matrix,
      z = grid_z_matrix,
      col = adjustcolor(
        sr_group_color,
        red.f = 0.25,
        green.f = 0.25,
        blue.f = 0.25
      ),
      opacity = 0.8,
      toggle = "Mean landscape",
      wireframe = TRUE
    )
  
  }
  
  if (return_mean_surface) return(data3js)
  
  # Add individual landscapes
  if (add_individual) {
    for (i in seq_along(sr_group_colbases)) {
      
      # Calculate the individual surface
      sr_coord <- sr_group_coords[i,]
      sr_colbase <- sr_group_colbases[i]
      
      grid_dists <- as.matrix(dist(rbind(
        sr_coord,
        cbind(
          as.vector(grid_x_matrix),
          as.vector(grid_y_matrix)
        )
      )))[1, -1]
      
      grid_z_matrix <- matrix(NA, length(grid_y_coords), length(grid_x_coords))
      grid_z_matrix[] <- sr_colbase - grid_dists*sr_slope
      
      # Add the individual surface
      data3js <- r3js::surface3js(
        data3js,
        x = grid_x_matrix,
        y = grid_y_matrix,
        z = grid_z_matrix,
        col = "grey70",
        opacity = 0.2,
        toggle = "Individual landscapes",
        wireframe = FALSE,
        doubleSide = TRUE
      )
      
      data3js <- r3js::surface3js(
        data3js,
        x = grid_x_matrix,
        y = grid_y_matrix,
        z = grid_z_matrix,
        col = "grey70",
        opacity = 0.4,
        toggle = "Individual landscapes",
        wireframe = TRUE
      )
      
    }
  }
  
  # Add the titers
  if (add_titers) {
  
    data3js <- lndscp3d_titers(
      data3js = data3js,
      object = list(
        coords = agCoords(map)[!is.na(sr_group_mean_logtiters),],
        logtiters = sr_group_mean_logtiters[!is.na(sr_group_mean_logtiters)],
        acmap_indices = which(!is.na(sr_group_mean_logtiters)),
        acmap = map
      ),
      zlim = plot_zlim,
      options = list(
        cex.titer = 1.5,
        cex.titer.impulse = 1.2,
        col.impulse = "grey60",
        lwd.impulse = 0.4,
        lwd.titer.outline = 0.5
      )
    )
    
  }
  
  # Add the titer 50 plane
  if (!is.na(titer_plane)) {
    
    x_grid <- c(plot_xlim[1], seq(from = plot_xlim[1], to = plot_xlim[2], by = 0.25), plot_xlim[2])
    y_grid <- c(plot_ylim[1], seq(from = plot_ylim[1], to = plot_ylim[2], by = 0.25), plot_ylim[2])
    z_grid <- matrix(titer_plane, length(x_grid), length(y_grid))
    data3js <- r3js::surface3js(
      data3js,
      x = x_grid,
      y = y_grid,
      z = z_grid,
      col = "grey50",
      opacity = 0.35,
      toggle = "Titer 50"
    )
    
    # data3js <- r3js::surface3js(
    #   data3js,
    #   x = x_grid,
    #   y = y_grid,
    #   z = z_grid,
    #   col = "grey50",
    #   opacity = 1,
    #   toggle = "Titer 50",
    #   wireframe = TRUE
    # )
    
  }
  
  # Draw border
  data3js <- r3js::lines3js(
    data3js,
    x = c(plot_xlim[1], plot_xlim[1], plot_xlim[2], plot_xlim[2], plot_xlim[1]),
    y = c(plot_ylim[1], plot_ylim[2], plot_ylim[2], plot_ylim[1], plot_ylim[1]),
    z = rep(plot_zlim[1], 5),
    lwd = 2,
    col = "grey70"
  )
  
  # Return data3js
  if (return_data3js) return(data3js)
  
  # Save r3js data
  saveRDS(
    list(
      data3js = data3js,
      rotation = angle$rotation,
      translation = angle$translation,
      zoom = angle$zoom,
      title = sr_group
    ),
    file = gsub("html$", "rds", filename)
  )
  
  # Create html widget
  widget <- r3js(
    data3js = data3js,
    rotation = angle$rotation,
    translation = angle$translation,
    zoom = angle$zoom,
    title = sr_group
  )
  
  # Add additional code
  widget <- htmlwidgets::onRender(
    x = widget,
    jsCode = paste0("function(el, x, data){
      if(window.top !== window.self) {
        var divsToHide = document.getElementsByClassName('toggle');
        for(var i = 0; i < divsToHide.length; i++){
            divsToHide[i].style.display = 'none';
        }
      }
    }")
  )
  
  # Save widget
  htmlwidgets::saveWidget(
    widget = widget,
    file = filename,
    selfcontained = FALSE,
    libdir = "lib"
  )
  
}

# Read slope data
slope_data <- readRDS("data/generated_data/slope_comparison_wt.rds")
slope_data <- slope_data$slope_factors

## Plot interactive basemap
data3js <- ablandscapes:::lndscp3d_setup(
  xlim = plot_xlim,
  ylim = plot_ylim,
  zlim = plot_zlim,
  aspect.z = 0.5,
  show.sidegrid = FALSE,
  show.axis = FALSE,
  options = list(
    lwd.grid = 1.5
  )
)

data3js <- ablandscapes:::lndscp3d_map(
  data3js = data3js,
  fit = list(acmap = map),
  xlim = plot_xlim,
  ylim = plot_ylim,
  zlim = plot_zlim,
  show.map.sera = FALSE,
  options = list(
    opacity.basemap.ags = 1,
    cex.basemap.ags = 1,
    cex.basemap.sr = 1.5
  )
)

data3js <- r3js::light3js(
  data3js,
  type = "ambient"
)

save3js(
  data3js = data3js, 
  rotation = c(0, 0, 0), 
  zoom = 0.6677, 
  title = "basemap",
  file = file.path(lndscps_dir, "basemap.html")
)

## WT comparison
sr_group_cols <- wt_slope_comparison_cols

data3js <- plot_sr_group_lndscp(
  sr_group = "D614G",
  sr_group_color = sr_group_cols["D614G"],
  add_individual = FALSE,
  add_titers = FALSE,
  titer_plane = NA,
  return_data3js = TRUE
)

data3js_3x <- plot_sr_group_lndscp(
  sr_group = "3x mRNA-1273",
  sr_group_color = sr_group_cols["3x mRNA-1273"],
  sr_slope = slope_data$estimate[slope_data$sr_group == "3x mRNA-1273"],
  return_mean_surface = TRUE
)

data3js_3xBD01 <- plot_sr_group_lndscp(
  sr_group = "3x mRNA-1273 BD01",
  sr_group_color = sr_group_cols["3x mRNA-1273 BD01"],
  sr_slope = slope_data$estimate[slope_data$sr_group == "3x mRNA-1273 BD01"],
  return_mean_surface = TRUE
)

data3js_2x <- plot_sr_group_lndscp(
  sr_group = "2x mRNA-1273",
  sr_group_color = sr_group_cols["2x mRNA-1273"],
  sr_slope = slope_data$estimate[slope_data$sr_group == "2x mRNA-1273"],
  return_mean_surface = TRUE
)

data3js_3x_6month <- plot_sr_group_lndscp(
  sr_group = "3x mRNA-1273 (6 month)",
  sr_group_color = sr_group_cols["3x mRNA-1273 (6 month)"],
  sr_slope = slope_data$estimate[slope_data$sr_group == "3x mRNA-1273 (6 month)"],
  return_mean_surface = TRUE
)

data3js$plot <- c(
  data3js$plot, 
  data3js_3x$plot, 
  data3js_3xBD01$plot, 
  data3js_2x$plot,
  data3js_3x_6month$plot
)

# Save widget data
saveRDS(
  object = list(
    data3js = data3js,
    rotation = angle$rotation,
    translation = angle$translation,
    zoom = angle$zoom,
    title = "WT slope comparison"
  ),
  file = file.path(lndscps_dir, "wt_comparison.rds")
)

# Create html widget
widget <- r3js(
  data3js = data3js,
  rotation = angle$rotation,
  translation = angle$translation,
  zoom = angle$zoom,
  title = "WT slope comparison"
)

widget <- htmlwidgets::onRender(
  widget,
  jsCode = paste0("function(el, x, data){
  el.style.outline = 'solid 2px #eeeeee';
  }")
)

htmlwidgets::saveWidget(
  widget = widget,
  file = file.path(lndscps_dir, "wt_comparison.html"),
  selfcontained = FALSE,
  libdir = "lib"
)

## WT slope comparison
sr_group_cols <- wt_slope_comparison_cols
surface_height <- 8

data3js <- plot_sr_group_lndscp(
  sr_group = "D614G",
  sr_group_color = sr_group_cols["D614G"],
  add_individual = FALSE,
  add_titers = FALSE,
  return_data3js = TRUE,
  titer_plane = NA,
  fix_surface = surface_height
)

data3js_3x <- plot_sr_group_lndscp(
  sr_group = "3x mRNA-1273",
  sr_group_color = sr_group_cols["3x mRNA-1273"],
  sr_slope = slope_data$estimate[slope_data$sr_group == "3x mRNA-1273"],
  return_mean_surface = TRUE,
  fix_surface = surface_height
)

data3js_3xBD01 <- plot_sr_group_lndscp(
  sr_group = "3x mRNA-1273 BD01",
  sr_group_color = sr_group_cols["3x mRNA-1273 BD01"],
  sr_slope = slope_data$estimate[slope_data$sr_group == "3x mRNA-1273 BD01"],
  return_mean_surface = TRUE,
  fix_surface = surface_height
)

data3js_2x <- plot_sr_group_lndscp(
  sr_group = "2x mRNA-1273",
  sr_group_color = sr_group_cols["2x mRNA-1273"],
  sr_slope = slope_data$estimate[slope_data$sr_group == "2x mRNA-1273"],
  return_mean_surface = TRUE,
  fix_surface = surface_height
)

data3js_3x_6month <- plot_sr_group_lndscp(
  sr_group = "3x mRNA-1273 (6 month)",
  sr_group_color = sr_group_cols["3x mRNA-1273 (6 month)"],
  sr_slope = slope_data$estimate[slope_data$sr_group == "3x mRNA-1273 (6 month)"],
  return_mean_surface = TRUE,
  fix_surface = surface_height
)


data3js$plot <- c(
  data3js$plot,
  data3js_3x$plot,
  data3js_3xBD01$plot,
  data3js_2x$plot,
  data3js_3x_6month$plot
)

# Save widget data
saveRDS(
  object = list(
    data3js = data3js,
    rotation = angle$rotation,
    translation = angle$translation,
    zoom = angle$zoom,
    title = "WT slope comparison fixed height"
  ),
  file = file.path(lndscps_dir, "wt_slope_comparison.rds")
)

# Create html widget
widget <- r3js(
  data3js = data3js,
  rotation = angle$rotation,
  translation = angle$translation,
  zoom = angle$zoom,
  title = "WT slope comparison fixed height"
)

widget <- htmlwidgets::onRender(
  widget,
  jsCode = paste0("function(el, x, data){
  el.style.outline = 'solid 2px #eeeeee';
  }")
)

htmlwidgets::saveWidget(
  widget = widget,
  file = file.path(lndscps_dir, "wt_slope_comparison.html"),
  selfcontained = FALSE,
  libdir = "lib"
)

for (sr_group in unique(srGroups(map))) {
  
  message(sr_group)
  if (sr_group %in% slope_data$sr_group) sr_slope <- slope_data$estimate[slope_data$sr_group == sr_group]
  else                                   sr_slope <- 1
  
  plot_sr_group_lndscp(
    sr_group = sr_group, 
    sr_slope = sr_slope,
    filename = file.path(lndscps_dir, paste0(sr_group, ".html"))
  )
  
}
