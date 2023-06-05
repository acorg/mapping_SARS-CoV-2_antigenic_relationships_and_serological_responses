
# Setup workspace
rm(list = ls())
library(Racmacs)
library(tidyverse)
library(patchwork)
set.seed(300)

source("functions/scales.R")
source("functions/ag_labels.R")
source("functions/sr_group_labels.R")
source("functions/map_longinfo.R")

# Read GMT data
ag_gmts <- readRDS("data/generated_data/ag_gmt_individual_effects_by_subgroup.rds") |> separate(sr_group, c("sr_group", "source"), "; ")

# Read the map
map <- read.acmap("data/maps/map_full_no_outliers.ace")

# Get individual titer data
long_map_info(map) |> 
  rename(source = sr_extra) -> long_titer_data

# Get n per source
long_titer_data |> 
  distinct(sr_name, source, sr_group) |> 
  group_by(sr_group, source) |> 
  count() -> n_by_source

# Determine serum groups with multiple sources
ag_gmts |> 
  distinct(sr_group, source) |> 
  group_by(sr_group) |> 
  count() |> 
  filter(n > 1) |> 
  pluck("sr_group") -> included_sr_groups

# Subset by serum group
plots <- lapply(included_sr_groups, function(sr_group_subset) {

  ag_gmts |> 
    filter(sr_group %in% sr_group_subset) -> plotdata
  
  # Determine antigens to include and antigen order
  plotdata |> 
    group_by(sr_group, ag_name) |> 
    summarise(
      ag_gmt_accounting_for_individual_effect = mean(ag_gmt_accounting_for_individual_effect),
      .groups = "drop"
    ) |> 
    filter(
      !is.na(ag_gmt_accounting_for_individual_effect)
    ) |> 
    arrange(
      -ag_gmt_accounting_for_individual_effect
    ) |> 
    pluck("ag_name") -> ag_order
  
  plotdata |> 
    mutate(
      ag_gmt_accounting_for_individual_effect = pmax(0, ag_gmt_accounting_for_individual_effect)
    ) -> plotdata
  
  long_titer_data |> 
    filter(
      sr_group == sr_group_subset
    ) -> plotdata_individual
  
  # Determine n by source
  source_n <- n_by_source |> filter(sr_group == sr_group_subset)
  source_labels <- sprintf(
    "%s (n=%s)", 
    source_n$source, 
    source_n$n
  )
  names(source_labels) <- source_n$source
  source_labels <- gsub("Moderna", "mRNA-1273", source_labels)
  
  plotdata |> 
    ggplot(
      aes(
        x = ag_name,
        y = ag_gmt_accounting_for_individual_effect,
        group = source,
        color = source
      )
    ) + 
    geom_line(
      data = plotdata_individual,
      mapping = aes(group = sr_name, y = logtiter),
      alpha = 0.2
    ) +
    geom_line(
      linewidth = 1,
      alpha = 0.8
    ) +
    geom_point(
      shape = "circle",
      alpha = 0.8
    ) +
    scale_x_discrete(
      limits = ag_order,
      labels = ag_who_labels
    ) + 
    scale_y_titer(
      ymin = 0
    ) +
    scale_color_brewer(
      type = "qual",
      palette = "Dark2",
      labels = source_labels
    ) +
    coord_cartesian(
      ylim = c(-0.4, 14)
    ) +
    annotate(
      geom = "tile",
      fill = agFill(map)[ag_order],
      x = ag_order,
      y = -0.75,
      color = NA,
      width = 1,
      height = 0.5
    ) +
    labs(
      title = sr_group_labels[sr_group_subset],
      x = "",
      y = "Titer",
      color = "Source"
    ) + 
    titerplot_theme() + 
    theme(
      plot.title = element_text(size = 11),
      legend.key = element_blank(),
      axis.text.x = ggtext::element_markdown()
    )
  
})

# Wrap and save plots
gp <- wrap_plots(plots[c(3, 4, 5, 1, 2)], nrow = 3, widths = c(1, 0.3), byrow = F) + 
  plot_annotation(tag_levels = "A") &
  theme(plot.tag.position = c(0, 1), plot.tag = element_text(hjust = 0, vjust = 0.5, size = 24))

ggsave(
  plot = gp,
  filename = "figures/som/figS2-comparing_sample_sources.pdf",
  width = 10,
  height = 12,
  units = "in"
)

ggsave(
  plot = gp,
  filename = "figures/som/figS2-comparing_sample_sources.png",
  width = 10,
  height = 12,
  units = "in"
)
