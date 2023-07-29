
# Setup workspace
rm(list = ls())
library(tidyverse)
library(patchwork)
library(Racmacs)
set.seed(100)

# Load map
map <- read.acmap('data/maps/map_ndsubset_no_outliers_slope_adjusted.ace')

# Style map
agSize(map) <- agSize(map) * 0.5
srSize(map) <- srSize(map) * 0.5

# Keep column bases fixed
colbases <- colBases(map)

# Set plot limits
maplims <- Racmacs:::mapPlotLims(map, sera = FALSE)

maps <- lapply(
  agNames(map), function(variant) {
    
    # Remove sera from serum group
    message(variant)
    map %>% removeAntigens(
      agNames(map) == variant
    ) %>%
      optimizeMap(2, 500, "none", fixed_column_bases = colbases) %>%
      realignMap(map) %>%
      procrustesMap(map, sera = FALSE) -> map_subset
    
    ggplot(
      map_subset,
      xlim = maplims$xlim,
      ylim = maplims$ylim,
      plot_stress = TRUE
    ) + 
      labs(
        title = paste("-", variant)
      )
  }
)

# Remove all B.1.617.2 variants
map %>% removeAntigens(
  grepl("B.1.617.2", agNames(map), fixed = T)
) %>%
  optimizeMap(2, 500, "none", fixed_column_bases = colbases) %>%
  realignMap(map) %>%
  procrustesMap(map, sera = FALSE) -> map_subset

gpA <- ggplot(
  map_subset,
  xlim = maplims$xlim,
  ylim = maplims$ylim,
  plot_stress = TRUE
) + 
  labs(
    title = "- all B.1.617.2 variants"
  )


gp <- wrap_plots(c(maps, list(gpA)), ncol = 4) + plot_annotation(tag_levels = "A") & theme(plot.tag.position = c(0.1, 0.8), legend.position = "none")
ggsave("figures/som/figS15-excluding_variants.pdf", gp, width = 8.5*1.3, height = 12.5*1.3, units = "in")
ggsave("figures/som/figS15-excluding_variants.png", gp, width = 8.5*1.3, height = 12.5*1.3, units = "in")

