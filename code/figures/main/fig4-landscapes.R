
# Setup workspace
rm(list = ls())
library(tidyverse)
source("functions/sr_group_labels.R")
lndscps_dir <- "figures/main/lndscps"

# Create a composite figure
landscape_panel <- function(
    panel_label,
    serum_group, 
    img_filepath 
) {
  
  cowplot::ggdraw() + 
    cowplot::draw_image(img_filepath) + 
    cowplot::draw_label(sr_group_labels_no_sera[serum_group], 0.03, 0.04, hjust = 0, vjust = 0, size = 12, fontface = "bold") +
    cowplot::draw_label(panel_label, 0.03, 0.96, hjust = 0, vjust = 1, size = 18, fontface = "bold") +
    theme(
      plot.margin = margin(0.15, 0.15, 0.15, 0.15, unit = "cm"),
      panel.border = element_rect(colour = "grey80", linewidth = 0.5)
    )
  
}

composite_plot <- cowplot::plot_grid(
  landscape_panel("A", "D614G", file.path(lndscps_dir, "D614G.png")),
  landscape_panel("B", "B.1.1.7", file.path(lndscps_dir, "B.1.1.7.png")),
  landscape_panel("K", "2x mRNA-1273", file.path(lndscps_dir, "2x mRNA-1273.png")),
  landscape_panel("C", "B.1.351", file.path(lndscps_dir, "B.1.351.png")),
  landscape_panel("D", "P.1", file.path(lndscps_dir, "P.1.png")),
  landscape_panel("L", "3x mRNA-1273 BD01", file.path(lndscps_dir, "3x mRNA-1273 BD01.png")),
  landscape_panel("E", "B.1.617.2", file.path(lndscps_dir, "B.1.617.2.png")),
  landscape_panel("F", "B.1.526+E484K", file.path(lndscps_dir, "B.1.526+E484K.png")),
  landscape_panel("M", "3x mRNA-1273 BD29", file.path(lndscps_dir, "3x mRNA-1273 BD29.png")),
  landscape_panel("G", "B.1.637", file.path(lndscps_dir, "B.1.637.png")),
  landscape_panel("H", "C.37", file.path(lndscps_dir, "C.37.png")),
  landscape_panel("N", "3x mRNA-1273 (6 month)", file.path(lndscps_dir, "3x mRNA-1273 (6 month).png")),
  landscape_panel("I", "BA.1", file.path(lndscps_dir, "BA.1.png")),
  landscape_panel("J", "BA.2", file.path(lndscps_dir, "BA.2.png")),
  landscape_panel("O", "2x mRNA-1273.351", file.path(lndscps_dir, "2x mRNA-1273.351.png")),
  ncol = 3
)

ggsave(
  filename = "figures/main/fig4_lndscps.pdf",
  plot = composite_plot,
  width = 10,
  height = 12,
  units = "in"
)

ggsave(
  filename = "figures/main/fig4_lndscps.png",
  plot = composite_plot,
  width = 10,
  height = 12,
  units = "in"
)
