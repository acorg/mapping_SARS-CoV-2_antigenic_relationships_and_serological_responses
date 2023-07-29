
rm(list = ls())
library(Racmacs)
library(tidyverse)

source("functions/diagnostics.R")
source("functions/imputation.R")
source("functions/scales.R")
source("functions/sr_group_labels.R")

map <- read.acmap('data/maps/map_ndsubset_no_outliers_slope_adjusted.ace')

cv_data <- readRDS("data/generated_data/cv_testing.rds")

dilution_stepsize <- 0

cv_data$measured_logtiter_lower <- cv_data$measured_logtiter - dilution_stepsize / 2
cv_data$measured_logtiter_upper <- cv_data$measured_logtiter + dilution_stepsize / 2

cv_data$residual_lower = cv_data$measured_logtiter_lower - cv_data$predicted_logtiter
cv_data$residual_upper = cv_data$measured_logtiter_upper - cv_data$predicted_logtiter

# Impute data
cv_data %>%
  filter(
    measured_titer != "*"
  ) %>%
  group_by(
    ag_name,
    sr_group
  ) %>%
  group_modify(
    .f = impute_table_residuals
  ) -> cv_data

# Plot the residual error with imputed values for < cases
## Boxplots split by antigens and sera
cv_data %>%
  mutate(
    sr_group = factor(sr_group, levels = c(
      "2x mRNA-1273", "3x mRNA-1273 BD01", "3x mRNA-1273 BD29", "3x mRNA-1273 (6 month)", 
      "2x mRNA-1273.351", "D614G", "B.1.1.7", "P.1", "B.1.351", "B.1.526+E484K", 
      "B.1.617.2", "B.1.637", "C.37", "BA.1", "BA.2")),
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
  facet_wrap(
    vars(sr_group),
    labeller = as_labeller(sr_group_labels)
  ) +
  scale_fill_manual(
    values = agFillScale(map)
  ) +
  scale_color_manual(
    values = c(
      "=" = "grey20",
      "<" = "lightblue",
      ">" = "pink"
    )
  ) +
  scale_y_continuous(
    breaks = -10:10
  ) +
  coord_cartesian(
    ylim = c(-10, 10)
  ) +
  titerplot_theme() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(
    x = "",
    y = "Measured log titer - Predicted log titer"
  ) -> gp

ggsave(
  plot = gp,
  filename = 'figures/som/figS23-measured_vs_predicted_titers_box.png', 
  width = 11, 
  height = 14
)

ggsave(
  plot = gp,
  filename = 'figures/som/figS23-measured_vs_predicted_titers_box.pdf', 
  width = 11, 
  height = 14
)
