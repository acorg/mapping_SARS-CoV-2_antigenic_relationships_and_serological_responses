
# Setup workspace
rm(list = ls())
library(tidyverse)
library(Racmacs)
library(patchwork)
set.seed(100)

# Load the maps
map <- read.acmap("data/maps/map_ndsubset_no_outliers_slope_adjusted.ace")
map_no_vaccine_sera <- subsetMap(
  map,
  sera = !srGroups(map) %in% c(
    "2x mRNA-1273",
    "3x mRNA-1273 BD01",
    "3x mRNA-1273 BD29",
    "3x mRNA-1273 (6 month)",
    "2x mRNA-1273.351"
  )
)

# Reoptimize the map without vaccine sera
map_no_vaccine_sera <- optimizeMap(
  map_no_vaccine_sera,
  number_of_dimensions = 2,
  number_of_optimizations = 1000,
  minimum_column_basis = "none"
) |> 
  realignMap(map)

# Calculate map limits
lims <- Racmacs:::mapPlotLims(map, sera = FALSE)

full_comparison <- ggplot(
  procrustesMap(
    map = map, 
    comparison_map = map_no_vaccine_sera,
    sera = FALSE
  ),
  xlim = lims$xlim,
  ylim = lims$ylim
)

subset_comparison <- ggplot(
  procrustesMap(
    map = map, 
    comparison_map = map_no_vaccine_sera, 
    antigens = agNames(map)[!agNames(map) %in% c(
      "BA.1.1",
      "BA.2.12.1",
      "BA.3",
      "BA.4/BA.5"
    )],
    sera = FALSE
  ),
  xlim = lims$xlim,
  ylim = lims$ylim
)

# Combine the plots
gp <- full_comparison + subset_comparison + plot_annotation(tag_levels = "A") &
  theme(
    plot.tag.position = c(0, 1), 
    plot.tag = element_text(
      hjust = 0,
      vjust = 1,
      size = 36,
      margin = margin(
        t = 0.3,
        l = 0.4,
        unit = "cm"
      )
    )
  )

# Save the plot
ggsave(
  plot = gp,
  filename = "figures/som/figS17-map_made_from_infection_only_sera.pdf",
  width = 11*1.4,
  height = 5*1.4,
  units = "in"
)

ggsave(
  plot = gp,
  filename = "figures/som/figS17-map_made_from_infection_only_sera.png",
  width = 11*1.4,
  height = 5*1.4,
  units = "in"
)


