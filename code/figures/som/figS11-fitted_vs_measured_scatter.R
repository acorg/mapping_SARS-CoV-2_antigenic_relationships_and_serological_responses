
rm(list = ls())
library(Racmacs)
library(tidyverse)

source("functions/map_longinfo.R")
source("functions/diagnostics.R")
source("functions/scales.R")

# Load map
map <- read.acmap("data/maps/map_ndsubset_no_outliers_slope_adjusted.ace")

# Get residual error table
residual_data <- residualErrorTable(map)

# Add antigen and sera info
residual_data %>% 
  left_join(long_ag_info(map), by = "ag_num") %>%
  left_join(long_sr_info(map), by = "sr_num") -> residual_data

# Plot the fitted against the measured titers as a scatter plot
residual_data %>%
  ggplot(
    aes(
      x = pmax(1, predicted_logtiter),
      y = pmax(1, measured_logtiter_upper)
    )
  ) +
  geom_point(
    alpha = 0.4
  ) +
  theme_bw() +
  coord_cartesian(
    xlim = c(1, 15),
    ylim = c(1, 15)
  ) +
  xlab('Fitted titers') + 
  ylab('Measured titers') + 
  scale_x_continuous(
    breaks = 1:15,
    labels = c("<20", 2^(2:15)*10)
  ) +
  scale_y_continuous(
    breaks = 1:15,
    labels = c("<20", 2^(2:15)*10)
  ) +
  geom_abline(intercept = 0, lty=2, col = 'grey') + 
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) 

ggsave("figures/som/figS11-fitted_vs_measured_scatter.png", width=5.3, height=5)
ggsave("figures/som/figS11-fitted_vs_measured_scatter.pdf", width=5.3, height=5)
