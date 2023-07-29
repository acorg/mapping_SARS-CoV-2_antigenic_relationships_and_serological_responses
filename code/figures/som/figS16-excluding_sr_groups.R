
# Setup workspace
rm(list = ls())
library(tidyverse)
library(patchwork)
library(Racmacs)
source("functions/sr_group_labels.R")

# Load map
map <- read.acmap("data/maps/map_ndsubset_no_outliers_slope_adjusted.ace")

# Style map
agSize(map) <- agSize(map) * 0.5
srSize(map) <- srSize(map) * 0.5

# Set serum group order
srGroups(map) <- factor(
  srGroups(map),
  c(
    "2x mRNA-1273",
    "3x mRNA-1273 BD01",
    "3x mRNA-1273 BD29",
    "3x mRNA-1273 (6 month)",
    "2x mRNA-1273.351",
    "D614G",
    "B.1.1.7",
    "P.1",
    "B.1.351",
    "B.1.526+E484K",
    "B.1.617.2",
    "B.1.637",
    "C.37",
    "BA.1",
    "BA.2"
  )
)

plot_remove_srgroups <- function(map, sr_groups, plottitle) {
  
  map %>% removeSera(
    srGroups(map) %in% sr_groups
  ) %>%
    optimizeMap(2, 500, "none") %>%
    realignMap(map) %>%
    procrustesMap(map, sera = FALSE) -> map_subset
  
  ggplot(
    map_subset,
    xlim = maplims$xlim,
    ylim = maplims$ylim,
    plot_stress = TRUE
  ) + 
    labs(
      title = plottitle
    )
  
}

# Set plot limits
maplims <- Racmacs:::mapPlotLims(map, sera = FALSE)

maps <- lapply(
  levels(srGroups(map)), function(sr_group) {
    
    # Remove sera from serum group
    message(sr_group)
    plot_remove_srgroups(map, sr_group, sr_group_labels_multiline[sr_group])
    
  }
)

gp <- wrap_plots(maps, ncol = 4) + plot_annotation(tag_levels = "A") & theme(plot.tag.position = c(0.08, 0.8), legend.position = "none")
ggsave("figures/som/figS16-excluding_sr_groups.pdf", gp, width = 12, height = 12, units = "in")
ggsave("figures/som/figS16-excluding_sr_groups.png", gp, width = 12, height = 12, units = "in")

