
# Setup workspace
rm(list = ls())
library(Racmacs)
library(tidyverse)
library(patchwork)
set.seed(300)

source("functions/scales.R")
source("functions/ag_labels.R")
source("functions/sr_group_labels.R")
source("functions/variables.R")
source("functions/homologous_ags.R")
source("functions/screenshot.R")

# Screenshot the slope comparison
file.remove("figures/summary/fig0B-lndscps_comparison.png")
screenshotWebpage(
  url = "figures/main/lndscps/wt_comparison_no_D614G.html",
  file = "figures/summary/fig0B-lndscps_comparison.png",
  vwidth = 1600,
  vheight = 880,
  cliprect = c(200, 0, 1200, 880)
)

sr_group_cols <- wt_slope_comparison_cols

gp <- readRDS("code/figures/summary/fig2B.rds")
gp <- gp + 
  theme(
    legend.position = "right",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "cm")
  ) + 
  # geom_vline(
  #   xintercept = 6,
  #   linetype = "solid",
  #   colour = "grey70"
  # ) + 
  scale_x_continuous(
    breaks = c(1, 6, 12)
  ) +
  scale_color_manual(
    values = sr_group_cols,
    labels = sr_group_labels_no_sera[c(
      "2x mRNA-1273",
      "3x mRNA-1273 BD01",
      "3x mRNA-1273 BD29",
      "3x mRNA-1273 (6 month)"
    )],
    limits = c(
      "2x mRNA-1273",
      "3x mRNA-1273 BD01",
      "3x mRNA-1273 BD29",
      "3x mRNA-1273 (6 month)"
    )
  ) + 
  theme(
    panel.grid.major.y = element_line(color = "#f6f6f6"),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 7.5),
    axis.text.y = element_text(size = 7.5),
    legend.text = element_text(size = 7.5),
    legend.title = element_text(size = 8),
    legend.key.height = unit(0.42, "cm")
  )

ggsave("figures/summary/fig0B-ntd_and_immunodominance.pdf", gp, width = 4, height = 3, units = "in")
ggsave("figures/summary/fig0B-ntd_and_immunodominance.png", gp, width = 4, height = 3, units = "in")
