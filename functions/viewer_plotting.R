
# Function to plot the antigenic map
plot_map <- function(
    map, 
    scaling_size = 0.6, 
    ag_scaling_size = scaling_size, 
    sr_scaling_size = scaling_size, 
    alter_opacity = TRUE, 
    sera_opacity = 0.6,
    grid.col = "#cfcfcf",
    width = NULL,
    height = NULL,
    datafile = NULL,
    ...
) {
  
  # Scale antigens and sera
  agSize(map) <- agSize(map)*ag_scaling_size
  srSize(map) <- srSize(map)*sr_scaling_size
  
  # Alter point opacity
  srOutline(map) <- adjustcolor(srOutline(map), alpha.f = sera_opacity)
  if (alter_opacity) {
    agFill(map) <- adjustcolor(agFill(map), alpha.f = 0.9)
    agOutline(map) <- adjustcolor(agOutline(map), alpha.f = 0.9)
  }
  
  # View the map
  plotlims <- Racmacs:::mapPlotLims(map, sera = FALSE)
  
  # Save map data to a file
  if (!is.null(datafile)) {
    saveRDS(
      object = list(
        map = map, 
        height = height,
        width = width,
        options = RacViewer.options(
          point.opacity = "inherit",
          xlim = plotlims$xlim,
          ylim = plotlims$ylim,
          grid.col = grid.col,
          ...
        )
      ),
      file = datafile
    )
  }
  
  # Return the viewer
  RacViewer(
    map = map, 
    height = height,
    width = width,
    options = RacViewer.options(
      point.opacity = "inherit",
      xlim = plotlims$xlim,
      ylim = plotlims$ylim,
      grid.col = grid.col,
      ...
    )
  )
  
}

highlight_ags <- function(map, ags) {
  
  agFill(map)[!agNames(map) %in% ags] <- adjustcolor(agFill(map)[!agNames(map) %in% ags], alpha.f = 0.3)
  agOutline(map)[!agNames(map) %in% ags] <- adjustcolor(agOutline(map)[!agNames(map) %in% ags], alpha.f = 0.3)
  map
  
}

output_map_divs <- function(output1, output2) {
  
  htmltools::div(
    htmltools::div(
      htmltools::h4("2 dimensions"),
      output1,
      style = "margin-right:20px;"
    ),
    htmltools::div(
      htmltools::h4("3 dimensions"),
      output2
    ),
    style = "display:flex;"
  )
  
}

linkags <- function(map, linked_ag_pairs) {
  
  for (linked_ag_pair in linked_ag_pairs) {
    
    agcoords1 <- agCoords(map)[agNames(map) == linked_ag_pair[1],]
    agcoords2 <- agCoords(map)[agNames(map) == linked_ag_pair[2],]
    while (length(agcoords1) < 3) agcoords1 <- c(agcoords1, 0)
    while (length(agcoords2) < 3) agcoords2 <- c(agcoords2, 0)
    
    map <- r3js::lines3js(
      data3js = map,
      x = c(agcoords1[1], agcoords2[1]),
      y = c(agcoords1[2], agcoords2[2]),
      z = c(agcoords1[3], agcoords2[3]),
      lwd = 8,
      geometry = FALSE
    )
    
  }
  
  # Return the map
  map
  
}

add_procrustes_grid <- function(
    map3d, 
    map2d, 
    excluded_ags = c(),
    grid_col = "grey80",
    flip_z_axis = FALSE
) {
  
  srCoords(map2d)[] <- NA
  agCoords(map2d)[agNames(map2d) %in% excluded_ags, ] <- NA
  map3d <- realignMap(map3d, map2d)
  if (flip_z_axis) {
    agCoords(map3d)[,3] <- -agCoords(map3d)[,3]
    srCoords(map3d)[,3] <- -srCoords(map3d)[,3]
  }
  agCoords(map3d) <- agCoords(map3d) # Bake transformation
  plotlims <- Racmacs:::mapPlotLims(map3d, sera = FALSE)
  for (x in seq(from = plotlims$xlim[1], to = plotlims$xlim[2])) {
    map3d <- r3js::lines3js(
      data3js = map3d,
      x = c(x, x),
      y = plotlims$ylim,
      z = c(0, 0),
      col = grid_col
    )
  }
  
  for (y in seq(from = plotlims$ylim[1], to = plotlims$ylim[2])) {
    map3d <- r3js::lines3js(
      data3js = map3d,
      x = plotlims$xlim,
      y = c(y, y),
      z = c(0, 0),
      col = grid_col
    )
  }
  
  map3d
  
}

add_2d_points <- function(
    map3d, 
    map2d, 
    point_scaling = 0.1,
    point_opacity = 1
  ) {
  
  r3js::points3js(
    data3js = map3d,
    x = agCoords(map2d)[,1],
    y = agCoords(map2d)[,2],
    z = rep(0, numAntigens(map2d)),
    shape = "circle filled",
    col = agOutline(map2d),
    lwd = 0.3,
    opacity = point_opacity,
    fill = agFill(map2d),
    size = agSize(map2d)*point_scaling
  )
  
}

link_2d_points <- function(
    map3d, 
    map2d, 
    col = "black"
  ) {
  
  for (i in seq_len(numAntigens(map3d))) {
    map3d <- r3js::lines3js(
      data3js = map3d,
      x = c(agCoords(map2d)[i, 1], agCoords(map3d)[i, 1]),
      y = c(agCoords(map2d)[i, 2], agCoords(map3d)[i, 2]),
      z = c(0, agCoords(map3d)[i, 3]),
      lwd = 2,
      col = col
    )
  }
  map3d
  
}

procrustes_2d_to_3d <- function(
    map2d, 
    map3d,
    link_col = "black",
    point_opacity_2d = 1,
    excluded_ags = c(),
    flip_z_axis = FALSE
) {
  
  map3d <- add_procrustes_grid(map3d, map2d, excluded_ags = excluded_ags, flip_z_axis = flip_z_axis)
  map3d <- add_2d_points(map3d, map2d, point_scaling = 0.2, point_opacity = point_opacity_2d)
  map3d <- link_2d_points(map3d, map2d, col = link_col)
  map3d
  
}

