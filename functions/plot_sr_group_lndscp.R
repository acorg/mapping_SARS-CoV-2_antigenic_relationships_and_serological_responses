
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
    lndscp_width = "400px",
    lndscp_height = "300px"
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
    
    x_grid <- seq(from = plot_xlim[1], to = plot_xlim[2], by = 0.5)
    y_grid <- seq(from = plot_ylim[1], to = plot_ylim[2], by = 0.5)
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
    
    data3js <- r3js::surface3js(
      data3js,
      x = x_grid,
      y = y_grid,
      z = z_grid,
      col = colorspace::darken("grey50", amount = 0.1),
      opacity = 0.6,
      toggle = "Titer 50",
      wireframe = TRUE
    )
    
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
  
  # Create html widget
  r3js(
    data3js = data3js,
    rotation = angle$rotation,
    translation = angle$translation,
    zoom = angle$zoom,
    title = sr_group,
    width = lndscp_width,
    height = lndscp_height
  )
  
}

