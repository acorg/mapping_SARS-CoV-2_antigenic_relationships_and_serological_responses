
# Setup workspace
rm(list = ls())
library(tidyverse)
library(Racmacs)
library(ablandscapes)
library(r3js)
library(patchwork)

# Load functions
source("functions/scales.R")
source("functions/fit_lndscp_val.R")
source("functions/sr_group_labels.R")
source("functions/ag_labels.R")

# Load map
map <- read.acmap("data/maps/map_ndsubset_no_outliers_slope_adjusted.ace")

# Read slope data
slope_data <- readRDS("data/generated_data/slope_comparison_wt.rds")
slope_data <- slope_data$slope_factors

# Load gmts after accounting for individual variation
ag_means_individual_effects <- readRDS("data/generated_data/ag_gmt_individual_effects.rds")
sr_group_gmts <- ag_means_individual_effects %>% rename(logtiter = ag_gmt_accounting_for_individual_effect)

# Work out antigen order
ag_order <- sr_group_gmts %>% filter(sr_group == "D614G") %>% arrange(-logtiter) %>% select(ag_name) %>% unlist() %>% unname()
ag_order <- ag_order[ag_order %in% agNames(map)[agGroups(map) == "wildtype"]]
omicron_order <- c("BA.1", "BA.1.1", "BA.2", "BA.2.12.1", "BA.3", "BA.4/BA.5")
ag_order <- ag_order[!ag_order %in% omicron_order]
ag_order <- c(ag_order, omicron_order)

# Setup to store all group fits
all_group_fits <- tibble(
  ag_name = character(0),
  sr_group = factor(NULL, levels = levels(srGroups(map))),
  sr_group_mean_logtiters = numeric(0),
  fitted_ag_gmts = numeric(0)
)

for (sr_group in levels(srGroups(map))) {

  # Set sera details
  sr_group_ag_means_individual_effects <- ag_means_individual_effects[ag_means_individual_effects$sr_group == sr_group,]
  sr_group_mean_logtiters <- sr_group_ag_means_individual_effects$ag_gmt_accounting_for_individual_effect[match(agNames(map), sr_group_ag_means_individual_effects$ag_name)]
  sr_group_coords <- srCoords(map)[srGroups(map) == sr_group,]
  sr_group_colbases <- colBases(map)[srGroups(map) == sr_group]
  
  sr_slope <- 1
  if (sr_group %in% slope_data$sr_group) {
    sr_slope <- slope_data$estimate[slope_data$sr_group == sr_group]
  }
  
  fitted_ag_gmts <- vapply(
    seq_len(numAntigens(map)),
    \(n) {
      fit_lndscp_val(
        x = agCoords(map)[n,1], 
        y = agCoords(map)[n,2],
        sr_coords = sr_group_coords,
        sr_colbases = sr_group_colbases,
        sr_slope = sr_slope
      )
    }, numeric(1)
  )
  
  sr_group_fits <- tibble(
    ag_name = agNames(map),
    sr_group = factor(sr_group, levels = levels(srGroups(map))),
    sr_group_mean_logtiters = sr_group_mean_logtiters,
    fitted_ag_gmts = fitted_ag_gmts
  )
  
  all_group_fits <- bind_rows(all_group_fits, sr_group_fits)

}

all_group_fits$sr_group_mean_logtiters[all_group_fits$sr_group_mean_logtiters < 1] <- 0.5
all_group_fits$fitted_ag_gmts[all_group_fits$fitted_ag_gmts < 1] <- 0.5

plot_group_fits <- function(sr_groups_subset, ncol) {
  
  plotdata <- all_group_fits %>% filter(sr_group %in% sr_groups_subset) 
  plotdata %>%
    ggplot() + 
    geom_line(
      data = filter(plotdata, !is.na(sr_group_mean_logtiters)),
      aes(
        x = ag_name,
        y = sr_group_mean_logtiters,
        group = sr_group
      ),
      color = "blue",
      linetype = "22"
    ) + 
    geom_line(
      aes(
        x = ag_name,
        y = sr_group_mean_logtiters,
        group = sr_group
      ),
      color = "blue"
    ) + 
    geom_line(
      aes(
        x = ag_name,
        y = sr_group_mean_logtiters,
        group = sr_group
      ),
      color = "blue"
    ) + 
    geom_point(
      aes(
        x = ag_name,
        y = sr_group_mean_logtiters,
        shape = "measured",
        color = "measured"
      )
    ) + 
    geom_point(
      aes(
        x = ag_name,
        y = fitted_ag_gmts,
        shape = "predicted",
        color = "predicted"
      ),
    ) + 
    scale_x_discrete(
      limits = ag_order,
      labels = ag_who_labels
    ) +
    scale_y_titer(ymin = 1) +
    scale_color_manual(
      values = c(
        measured = "blue",
        predicted = "red"
      )
    ) + 
    scale_shape_manual(
      values = c(
        measured = "circle",
        predicted = "circle open"
      )
    ) + 
    facet_wrap(
      vars(sr_group),
      labeller = as_labeller(sr_group_labels),
      ncol = ncol
    ) + 
    coord_cartesian(
      ylim = shrinkrange(c(-0.5, 12.5), 0.05)
    ) +
    labs(
      shape = "GMT",
      colour = "GMT",
      x = ""
    ) +
    titerplot_theme() +
    theme(
      axis.text.x = ggtext::element_markdown()
    ) -> gp
  
  # Annotate region below detectable
  gp <- gp +
    annotate(
      "tile",
      x = agNames(map),
      y = 0,
      height = 2,
      fill = "grey50",
      color = NA,
      alpha = 0.3
    )
  
  # Annotate colors for each antigen
  for (n in seq_len(numAntigens(map))) {
    gp <- gp +
      annotate(
        "tile",
        x = agNames(map)[n],
        y = -0.75,
        height = 1,
        fill = agFill(map)[n],
        color = NA
      ) 
  }
  
  # Return the plot
  gp
  
}

infection_sera <- c(
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
)

vaccine_sera <- c(
  "2x mRNA-1273",
  "3x mRNA-1273 BD01",
  "3x mRNA-1273 BD29",
  "3x mRNA-1273 (6 month)",
  "2x mRNA-1273.351"
)

gpA <- plot_group_fits(infection_sera, ncol = 2)
gpB <- plot_group_fits(vaccine_sera, ncol = 1)

gp <- gpA + gpB +
  plot_layout(
    widths = c(2, 1),
    guides = "collect"
  ) + plot_annotation(
    tag_levels = list(c("A", "B "))
  ) & 
  theme(
    plot.tag = element_text(size = 24)
  )

# Save the plot
ggsave("figures/som/figS24-ablandscapes_fitted_vs_measured.pdf", gp, width = 14, height = 12, units = "in")
ggsave("figures/som/figS24-ablandscapes_fitted_vs_measured.png", gp, width = 14, height = 12, units = "in")

