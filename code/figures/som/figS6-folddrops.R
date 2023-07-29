
# Setup workspace
rm(list = ls())
library(Racmacs)
library(tidyverse)
library(patchwork)
set.seed(300)

source("functions/scales.R")
source("functions/annotate_across_facets.R")
source("functions/sr_group_labels.R")
source("functions/ag_labels.R")

# Read the map
map <- read.acmap("data/maps/map_full_no_outliers.ace")

# Subset to only wildtypes
map <- subsetMap(map, antigens = agGroups(map) == "wildtype")

# Get the folddrop results
all_group_results <- readRDS("data/generated_data/folddrops.rds")

# Cycle through serum groups
plots <- list()

# Function to plot foldchange
plot_folddiff <- function(
    all_group_results, 
    sr_groups, 
    ag_names,
    ag_order = NULL,
    show.legend = FALSE, 
    coloring,
    add_pointrange = TRUE,
    xlim = NULL,
    ylim = c(-10, 2),
    ncol = 2
  ) {
  
  # Subset data
  all_group_results %>%
    filter(
      sr_group %in% sr_groups,
      ag %in% ag_names
    ) -> sr_group_results
  
  # Define order
  if (is.null(ag_order)) {
    sr_group_gmts <- readRDS("data/generated_data/ag_gmt_individual_effects.rds")
    ag_order <- sr_group_gmts %>% filter(sr_group == "D614G") %>% arrange(-ag_gmt_accounting_for_individual_effect) %>% select(ag_name) %>% unlist() %>% unname()
    ag_order <- ag_order[ag_order %in% agNames(map)[agGroups(map) == "wildtype"]]
    omicron_order <- c("BA.1", "BA.1.1", "BA.2", "BA.2.12.1", "BA.3", "BA.4/BA.5")
    ag_order <- ag_order[!ag_order %in% omicron_order]
    ag_order <- c(ag_order, omicron_order)
  }
  
  # Refactor serum groups
  sr_group_results$sr_group <- factor(
    sr_group_results$sr_group,
    sr_groups
  )
  
  sr_group_results %>%
    ggplot(
      aes(
        x = ag,
        y = diff,
        color = sr_group
      )
    ) -> gp
  
  if (add_pointrange) {
    
    gp <- gp + geom_pointrange(
      mapping = aes(
        ymin = diff_lower,
        ymax = diff_upper
      ),
      show.legend = show.legend,
      position = position_dodge2(width = 0.75)
    )
    
  }
  
  ag_fill_colors <- agFill(map)
  names(ag_fill_colors) <- agNames(map)
  
  gp <- gp + scale_x_discrete(
    limits = ag_order,
    labels = ag_who_labels
  ) + 
    scale_y_continuous(
      breaks = floor(ylim[2]):ceiling(ylim[1]),
      labels = foldchange
    ) +
    scale_color_manual(
      values = srGroupOutline(map)
    ) +
    geom_text(
      data = distinct(sr_group_results, sr_group),
      mapping = aes(
        x = 21,
        y = 1.75,
        label = sr_group_labels[as.character(sr_group)]
      ),
      show.legend = FALSE,
      color = "black",
      # fontface = "bold",
      hjust = 1,
      size = 3
    ) +
    geom_text(
      data = filter(sr_group_results, homologous),
      mapping = aes(
        x = length(ag_order) + 0.2,
        y = 0.7,
        label = paste("n =", n)
      ),
      size = 3,
      hjust = 1,
      color = "black"
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey40"
    ) +
    annotate_across_facets(
      geom = "tile",
      x = distinct(sr_group_results, sr_group, ag, .keep_all = TRUE)$ag,
      y = ylim[1] - 0.5,
      facet_var = "sr_group",
      facet_value = distinct(sr_group_results, sr_group, ag, .keep_all = TRUE)$sr_group,
      fill = ag_fill_colors[distinct(sr_group_results, sr_group, ag, .keep_all = TRUE)$ag]
    ) +
    coord_cartesian(
      xlim = xlim,
      ylim = ylim
    ) +
    labs(
      x = "",
      y = "Estimated fold drop"
    ) + 
    facet_wrap(
      vars(sr_group),
      ncol = ncol
    ) +
    titerplot_theme() +
    theme(
      axis.text.x = ggtext::element_markdown(
        angle = 45,
        hjust = 1
      ),
      panel.grid.minor = element_blank(),
      strip.text = element_blank()
    ) -> gp
  
  gp
  
}

gpA <- plot_folddiff(
  all_group_results = all_group_results, 
  sr_groups = c(
    "D614G",
    "B.1.1.7",
    "B.1.351",
    "P.1",
    "B.1.617.2",
    "B.1.526+E484K",
    "B.1.637",
    "C.37",
    "BA.1",
    "BA.2"
  ),
  ag_names = agNames(map),
  show.legend = FALSE, 
  coloring = "ag",
  ncol = 2
)

gpB <- plot_folddiff(
  all_group_results = all_group_results, 
  sr_groups = c("2x mRNA-1273", "3x mRNA-1273 BD01", "3x mRNA-1273 BD29", "3x mRNA-1273 (6 month)", "2x mRNA-1273.351"),
  ag_names = agNames(map),
  show.legend = FALSE, 
  coloring = "ag",
  ncol = 1
)

gpA <- gpA + theme(plot.margin = margin(l = 0.5, unit = "cm"), plot.tag.position = c(-0.01, 0.98))
gpB <- gpB + theme(plot.margin = margin(l = 1, unit = "cm"), plot.tag.position = c(-0.02, 0.98))

gp <- gpA + gpB +
  plot_layout(
    widths = c(2, 1)
  ) + plot_annotation(
    tag_levels = list(c("A", "B "))
  ) & 
  theme(
    plot.tag = element_text(size = 28)
  )

# Save the plot
ggsave("figures/som/figS6-folddrops.pdf", gp, width = 14, height = 12, units = "in")
ggsave("figures/som/figS6-folddrops.png", gp, width = 14, height = 12, units = "in")

