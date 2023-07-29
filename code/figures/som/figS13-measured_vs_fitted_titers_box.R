
rm(list = ls())
library(Racmacs)
library(tidyverse)

source("functions/scales.R")
source("functions/sr_group_labels.R")

map <- read.acmap("data/maps/map_ndsubset_no_outliers_slope_adjusted.ace")
residual_data <- readRDS('data/generated_data/imputed_titers.rds')

# Plot the residual error with imputed values for < cases
## Boxplots split by antigens and sera
residual_data %>%
  mutate(
    ag_name = factor(ag_name, levels = c(
      "D614G", "B.1.1.7", "B.1.429", "B.1.617.2 (AY.1)+K417N", "B.1.617.2", "B.1.617.2+K417N", 
      "B.1.617.2 (AY.2)+K417N", "C.37", "P.1", "B.1.1.7+E484K", "B.1.526+E484K", "B.1.617.1", "B.1.351", 
      "B.1.617.2 (AY.3)+E484Q", "B.1.621", "BA.1", "BA.1.1", "BA.2", "BA.2.12.1", "BA.3", "BA.4/BA.5"))
  ) %>%
  ggplot(
    aes(
      x = ag_name,
      y = residual_error_imputed,
      fill = ag_name
    )
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_boxplot(
    lwd = 0.25,
    outlier.shape = NA,
    alpha = 0.8,
    show.legend = FALSE
  ) +
  ggbeeswarm::geom_beeswarm(
    aes(
      color = titer_type
    ),
    size = 0.5,
    cex = 0.5,
    alpha = 0.8,
    show.legend = FALSE
  ) +
  facet_wrap(
    vars(sr_group),
    labeller = as_labeller(sr_group_labels)
  ) +
  scale_color_manual(
    values = c(
      "=" = "grey20",
      "<" = "lightblue",
      ">" = "pink"
    )
  ) +
  scale_fill_manual(
    values = agFillScale(map)
  ) +
  scale_x_discrete(
    limits = agNames(map)
  ) +
  scale_y_continuous(
    breaks = -5:5
  ) +
  coord_cartesian(
    ylim = c(-5, 5)
  ) +
  titerplot_theme() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(
    x = "",
    y = "Measured log titer - Fitted log titer"
  ) -> gp

ggsave('figures/som/figS13-measured_vs_fitted_titers_box.png', width = 11, height = 10)
ggsave('figures/som/figS13-measured_vs_fitted_titers_box.pdf', width = 11, height = 10)

