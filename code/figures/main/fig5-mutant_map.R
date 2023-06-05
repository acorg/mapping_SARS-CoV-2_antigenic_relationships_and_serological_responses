
# Setup workspace
rm(list = ls())
library(tidyverse)
library(Racmacs)
library(patchwork)
source("functions/viewer_plotting.R")
source("functions/screenshot.R")
set.seed(100)

# Function to add an arrow
add_arrow <- function(from, to, from_offset, to_offset) {
  
  unit_vector <- c(to[1] - from[1], to[2] - from[2])
  unit_vector <- unit_vector / sqrt(sum(unit_vector^2))
  list(
    annotate(
      geom = "segment",
      x = from[1] + unit_vector[1]*from_offset,
      y = from[2] + unit_vector[2]*from_offset,
      xend = to[1] - unit_vector[1]*(to_offset + 0.1),
      yend = to[2] - unit_vector[2]*(to_offset + 0.1),
      linewidth = 1.5
    ),
    annotate(
      geom = "segment",
      x = from[1] + unit_vector[1]*from_offset,
      y = from[2] + unit_vector[2]*from_offset,
      xend = to[1] - unit_vector[1]*to_offset,
      yend = to[2] - unit_vector[2]*to_offset,
      linewidth = 0,
      arrow = arrow(
        type = "closed",
        length = unit(0.04, "npc"),
        angle = 20
      )
    )
  )
  
}

# Function for plotting
doplot <- function(
    map, 
    highlighted_ags, 
    ag_label_adjustments,
    arrow_pairs = list()
  ) {

  mutants <- agNames(map) %in% highlighted_ags
  
  # Do not show the L452R variant since it was only titrated against WT sera
  agShown(map)[agNames(map) == "D614G+L452R"] <- FALSE
  agSize(map) <- agSize(map)*0.7
  srSize(map) <- srSize(map)*0.7
  
  # Adjust colors
  agFill(map)[!mutants] <- adjustcolor(agFill(map)[!mutants], alpha.f = 0.4)
  agOutline(map)[!mutants] <- adjustcolor(agOutline(map)[!mutants], alpha.f = 0.4)
  srOutline(map) <- adjustcolor(srOutline(map), alpha.f = 0.4)
  
  # Plot the map
  lims <- Racmacs:::mapPlotLims(map, sera = FALSE, padding = 0.5)
  
  drawing_order <- ptDrawingOrder(map)
  for (ag in highlighted_ags) {
    agnum <- match(ag, agNames(map))
    drawing_order <- c(drawing_order[drawing_order != agnum], agnum)
  }
  ptDrawingOrder(map) <- drawing_order
  agOutlineWidth(map)[mutants] <- 2
  
  # Do base plot
  par(mar = rep(1, 4))
  gp <- ggplot(
    map, 
    fill.alpha = NULL, 
    outline.alpha = NULL,
    xlim = lims$xlim - 0.3,
    ylim = lims$ylim
  )
  
  # Label the mutants
  label_adjustments <- agCoords(map)
  label_adjustments[] <- 0
  
  for (i in seq_along(highlighted_ags)) {
    label_adjustments[agNames(map) == highlighted_ags[i]] <- ag_label_adjustments[[i]]
  }
  
  text_coords <- agCoords(map)[mutants,] + label_adjustments[mutants,]
  gp <- gp + annotate(
    "text",
    x = text_coords[,1],
    y = text_coords[,2],
    hjust = 0,
    label = agNames(map)[mutants],
    fontface = "bold",
    size = 3.5
  )
  
  # Annotate other variants
  variants <- c("P.1", "B.1.351", "B.1.617.2", "BA.1", "B.1.1.7", "D614G")
  xnudge <- rep(0, length(variants))
  ynudge <- rep(0.4, length(variants))
  ynudge[3] <- -0.4
  ynudge[6] <- 0
  xnudge[6] <- -0.675
  
  for (i in seq_along(variants)) {
    if (!variants[i] %in% highlighted_ags) {
      gp <- gp + annotate(
        "text",
        x = agCoords(map)[variants[i], 1] + xnudge[i],
        y = agCoords(map)[variants[i], 2] + ynudge[i],
        label = variants[i],
        fontface = "bold",
        size = 4,
        color = adjustcolor(agFill(map)[agNames(map) == variants[i]], alpha.f = 1.2)
      )
    }
  }
  
  # Annotate arrows
  for (arrow_pair in arrow_pairs) {
    
    gp <- gp + add_arrow(
      from = agCoords(map)[arrow_pair[1],],
      to = agCoords(map)[arrow_pair[2],],
      from_offset = agSize(map)[arrow_pair[1]]*0.022,
      to_offset = agSize(map)[arrow_pair[2]]*0.026
    )
    
  }
  
  # Return the plot
  gp
  
}

gpA <- doplot(
  map = read.acmap("data/maps/map_ndsubset_no_outliers_plus_d614g_mutants_slope_adjusted.ace"),
  highlighted_ags = c(
    "D614G+N501Y",
    "D614G",
    "D614G+E484K",
    "D614G+L452R+E484Q",
    "D614G+E484Q"
  ),
  ag_label_adjustments = list(
    c(0.3, 0),
    c(-1.1, 0),
    c(0.25, 0),
    c(0.2, -0.2),
    c(0.2, -0.2)
  ),
  arrow_pairs = list(
    c("D614G", "D614G+E484K"),
    c("D614G", "D614G+E484Q"),
    c("D614G", "D614G+L452R+E484Q")
  )
)

gpB <- doplot(
  map = read.acmap("data/maps/map_ndsubset_no_outliers_plus_417_mutants_slope_adjusted.ace"),
  highlighted_ags = c(
    "P.1",
    "B.1.429",
    "B.1.351",
    "B.1.617.1",
    "P.1+T417K",
    "B.1.429+K417N",
    "B.1.351+N417K",
    "B.1.617.1+K417N"
  ),
  ag_label_adjustments = list(
    c(-0.7, 0),
    c(-1.1, 0),
    c(0.4, 0),
    c(-1.3, 0),
    c(-1.35, 0),
    c(-1.85, 0),
    c(0.3, 0),
    c(-2, 0)
  ),
  arrow_pairs = list(
    c("P.1", "P.1+T417K"),
    c("B.1.429", "B.1.429+K417N"),
    c("B.1.351", "B.1.351+N417K"),
    c("B.1.617.1", "B.1.617.1+K417N")
  )
)

gpC <- doplot(
  map = read.acmap("data/maps/map_ndsubset_no_outliers_plus_ba1_mutant_slope_adjusted.ace"),
  highlighted_ags = c(
    "BA.1",
    "BA.1+A484K"
  ),
  ag_label_adjustments = list(
    c(-0.8, 0),
    c(-1, 0.35)
  ),
  arrow_pairs = list(
    c("BA.1", "BA.1+A484K")
  )
)

# Plot the 3d mutant map
main_map <- read.acmap("data/maps/map_ndsubset_no_outliers_slope_adjusted.ace")
map2d_ba1e484K <- read.acmap("data/maps/map_ndsubset_no_outliers_plus_ba1_mutant_slope_adjusted.ace")
map3d_ba1e484K <- read.acmap("data/maps/map_ndsubset_no_outliers_plus_ba1_mutant_slope_adjusted_3d.ace")
map3d_ba1e484K <- highlight_ags(map3d_ba1e484K, c("BA.1", "BA.1+A484K"))
srOutlineWidth(map3d_ba1e484K) <- srOutlineWidth(map3d_ba1e484K)*0.5

# Align maps
map3dpc2d_ba1e484K <- procrustes_2d_to_3d(
  map2d = map2d_ba1e484K, 
  map3d = map3d_ba1e484K, 
  link_col = "grey80",
  point_opacity_2d = 0.4,
  excluded_ags = c("BA.4/BA.5"),
  flip_z_axis = TRUE
)

# Link the mutant pairs
map3dpc2d_ba1e484K <- linkags(
  map = map3dpc2d_ba1e484K, 
  linked_ag_pairs = list(c("BA.1", "BA.1+A484K"))
)

# Save map data for panel A
read.acmap("data/maps/map_ndsubset_no_outliers_plus_d614g_mutants_slope_adjusted.ace") |> 
  linkags(
    linked_ag_pairs = list(
      c("D614G", "D614G+N501Y"),
      c("D614G", "D614G+L452R"),
      c("D614G", "D614G+E484Q"),
      c("D614G", "D614G+E484K"),
      c("D614G", "D614G+L452R+E484Q")
    )
  ) |> 
  saveRDS("figures/main/fig5a-mutant_maps.rds")

# Save map data for panel B
read.acmap("data/maps/map_ndsubset_no_outliers_plus_417_mutants_slope_adjusted.ace") |> 
  linkags(
    linked_ag_pairs = list(
      c("P.1", "P.1+T417K"),
      c("B.1.429", "B.1.429+K417N"),
      c("B.1.351", "B.1.351+N417K"),
      c("B.1.617.1", "B.1.617.1+K417N")
    )
  ) |> 
  saveRDS("figures/main/fig5b-mutant_maps.rds")

# Create a web viewer
map3dpc2d_ba1e484K_viewer <- plot_map(
  map3dpc2d_ba1e484K, 
  scaling_size = 0.2, 
  alter_opacity = F, 
  sera_opacity = 0.2,
  grid.col = NA,
  rotation = c(-1.1883, 0.0143, -0.1414),
  translation = c(-0.0037, 0.0270, 0.0064),
  zoom = 1.2,
  datafile = "figures/main/fig5c-mutant_maps.rds"
)

# Save the web viewer
htmlwidgets::saveWidget(
  widget = map3dpc2d_ba1e484K_viewer,
  file = "figures/main/fig5c-mutant_maps.html",
  libdir = ".lib"
)

# Screenshot an image
screenshotWebpage(
  url = "figures/main/fig5c-mutant_maps.html",
  file = "figures/main/fig5c-mutant_maps.png",
  vwidth = 1200,
  vheight = 1000,
  cliprect = c(100, 50, 1000, 920)
)

# Get an image of the 3d map
gpC_3d <- cowplot::ggdraw() + 
  cowplot::draw_image("figures/main/fig5c-mutant_maps.png") + 
  cowplot::draw_label("C", 0.05, 0.94, size = 28) +
  cowplot::draw_label("BA.1", 0.755, 0.58, hjust = 0, vjust = 0, size = 10, fontface = "bold") +
  cowplot::draw_label("BA.1+A484K", 0.705, 0.375, hjust = 0, vjust = 0, size = 10, fontface = "bold") +
  cowplot::draw_line(
    x = c(0.7, 0.89),
    y = c(0.59, 0.395),
    arrow = arrow(angle = 20, length = unit(0.5, "cm"), type = "closed"),
    linewidth = 0
  ) +
  theme(
    panel.border = element_rect(colour = "grey70", linewidth = 0.5),
    plot.margin = margin(0.4, 0.4, 0.4, 0.4, unit = "cm")
  )

# Save the main plot
gp <- gpA + gpB + gpC_3d +
  plot_layout(widths = c(1, 1, 0.95)) +
  plot_annotation(tag_levels = list(c("A", "B", ""))) & 
  theme(plot.tag = element_text(size = 28), plot.tag.position = c(0.05, 0.94))

ggsave("figures/main/fig5-mutant_maps.pdf", gp, width = 17, height = 5.2, units = "in")
ggsave("figures/main/fig5-mutant_maps.png", gp, width = 17, height = 5.2, units = "in")

# Save the SOM figure plot
ggsave("figures/som/figS29-BA1_A484K_2d.pdf", gpC, width = 6, height = 4.8, units = "in")
ggsave("figures/som/figS29-BA1_A484K_2d.png", gpC, width = 6, height = 4.8, units = "in")

