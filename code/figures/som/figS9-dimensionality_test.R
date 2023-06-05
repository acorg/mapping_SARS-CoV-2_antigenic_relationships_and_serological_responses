
# Setup workspace
rm(list = ls())
library(Racmacs)
library(patchwork)
library(tidyverse)
source("functions/viewer_plotting.R")
source("functions/screenshot.R")

# Load the dimension test results
dimtest_results <- readRDS("data/generated_data/dimtest_results.rds")

# Calculate confidence intervals
dimtest_results |> 
  mutate(
    mean_rmse_detectable_upper = mean_rmse_detectable + sqrt(var_rmse_detectable / replicates)*qt(0.975, df = replicates - 1),
    mean_rmse_detectable_lower = mean_rmse_detectable + sqrt(var_rmse_detectable / replicates)*qt(0.025, df = replicates - 1)
  ) -> dimtest_results

# Plot the results
dimtest_results |> 
  ggplot(
    aes(
      x = dimensions,
      y = mean_rmse_detectable,
      ymin = mean_rmse_detectable_lower,
      ymax = mean_rmse_detectable_upper
    )
  ) + 
  geom_line() +
  geom_linerange() +
  geom_point() +
  theme_bw() +
  labs(
    x = 'Dimensions',
    y = 'Mean RMSE of detectable titers'
  ) -> gpA

# Plot a 3D map
map3D <- read.acmap('data/maps/map_ndsubset_no_outliers_slope_adjusted_3d.ace')


# Screenshot the maps
screenshotWebpage(
  url = "figures/main/maps/3d_map_front.html",
  file = "figures/som/figS9b-dimensionality_test.png",
  vwidth = 500,
  vheight = 520,
  cliprect = c(10, 20, 480, 470)
)

screenshotWebpage(
  url = "figures/main/maps/3d_map_side.html",
  file = "figures/som/figS9c-dimensionality_test.png",
  vwidth = 500,
  vheight = 520,
  cliprect = c(10, 20, 480, 470)
)


# Combine the plots
gpB <- cowplot::ggdraw() + 
  cowplot::draw_image("figures/som/figS9b-dimensionality_test.png") +
  theme(
    panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.5),
    plot.margin = margin(0.4, 0.4, 0.4, 0.4, unit = "cm"),
    plot.tag.position = c(0, 1)
  )

gpC <- cowplot::ggdraw() + 
  cowplot::draw_image("figures/som/figS9c-dimensionality_test.png") +
  theme(
    panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.5),
    plot.margin = margin(0.4, 0.4, 0.4, 0.4, unit = "cm"),
    plot.tag.position = c(0, 1)
  )

gpA <- gpA + theme(
  plot.tag.position = c(-0.15, 1), 
  plot.margin = margin(0.4, 0.4, 0.4, 1, unit = "cm")
)
gp <- gpA + gpB + gpC + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 16, face = "bold", hjust = 0, vjust = 1, margin = margin(0.2, 0.2, 0.2, 0.2, unit = "cm")))


# Save the combined plot
ggsave(
  filename = "figures/som/figS9-dimensionality_test.pdf",
  plot = gp,
  width = 8,
  height = 3,
  units = "in"
)

ggsave(
  filename = "figures/som/figS9-dimensionality_test.png",
  plot = gp,
  width = 8,
  height = 3,
  units = "in"
)
