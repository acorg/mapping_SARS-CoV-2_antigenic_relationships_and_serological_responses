
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

pt_scaling <- 1.5

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
    ylab,
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
        ymin = diff_lower,
        ymax = diff_upper,
        color = sr_group
      )
    ) -> gp
  
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
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey40"
    ) +
    coord_cartesian(
      xlim = xlim,
      ylim = ylim
    ) +
    labs(
      x = "",
      y = ylab
    ) +
    titerplot_theme() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      panel.grid.minor = element_blank(),
      strip.text = element_blank()
    ) -> gp
  
  gp
  
}

# Plot vaccine vs non-vaccine plots
vac1_plot <- plot_folddiff(
  all_group_results = all_group_results, 
  sr_groups = c(
    "D614G",
    "2x mRNA-1273", 
    "3x mRNA-1273 BD01", 
    "3x mRNA-1273", 
    "3x mRNA-1273 (6 month)"
  ),
  ag_names = agNames(map),
  show.legend = TRUE,
  coloring = "sr_group",
  ylim = c(-8.8, 0.5),
  ylab = "Estimated fold drop vs D614G"
)

# Get estimates of slope difference
wt_slope_est_data <- readRDS("data/generated_data/slope_comparison_wt.rds")
wt_slope_est_data$ag_folddrops %>%
  rename(
    diff = estimate
  ) %>%
  mutate(
    diff_lower = NA_real_,
    diff_upper = NA_real_
  ) -> folddrop_data

folddrop_data <- bind_rows(
  folddrop_data,
  tibble(ag = "D614G", diff = 0)
)

# Set antigen and serum order
ag_order <- folddrop_data$ag[order(-folddrop_data$diff)]

slope_data <- wt_slope_est_data$slope_factors

size_closed <- 0.3
size_open <- 3

# Function to add slope-based fold-drop for subjects
add_slope_folddrop <- function(
    sr_group_subset,
    nudge
) {
  
  if (sr_group_subset %in% c("D614G", "B.1.351")) slope_est <- 1
  else slope_est <- slope_data$estimate[slope_data$sr_group == sr_group_subset]
  
  list(
    geom_line(
      mapping = aes(group = sr_group_subset),
      data = mutate(
        folddrop_data, 
        diff = diff*slope_est,
        sr_group = sr_group_subset
      ),
      linetype = "11",
      position = position_nudge(x = nudge),
      na.rm = TRUE
    ),
    geom_point(
      mapping = aes(group = sr_group_subset),
      data = mutate(
        folddrop_data, 
        diff = diff*slope_est,
        sr_group = sr_group_subset
      ),
      position = position_nudge(x = nudge),
      size = 0.9*pt_scaling,
      alpha = 0.6,
      na.rm = TRUE
    ),
    geom_point(
      mapping = aes(group = sr_group_subset),
      data = mutate(
        folddrop_data, 
        diff = diff*slope_est,
        sr_group = sr_group_subset
      ),
      position = position_nudge(x = nudge),
      color = "white",
      size = 0.3*pt_scaling,
      na.rm = TRUE
    )
  )
  
}

add_measured_folddrop <- function(
    sr_group_subset,
    nudge
) {
  
  list(
    geom_linerange(
      data = filter(all_group_results, sr_group == sr_group_subset),
      size = 0.25*pt_scaling,
      position = position_nudge(x = nudge),
      alpha = 0.5,
      na.rm = TRUE
    ),
    geom_point(
      data = filter(all_group_results, sr_group == sr_group_subset),
      size = 1*pt_scaling,
      position = position_nudge(x = nudge),
      na.rm = TRUE
    )
  )
  
}

nudges <- seq(from = -0.4, to = 0.4, length.out = 5)

sr_group_cols <- wt_slope_comparison_cols

# Add number per group to group legend labels
all_group_results |> 
  filter(homologous) -> group_n_tibble
group_n <- group_n_tibble$n
names(group_n) <- group_n_tibble$sr_group

sr_group_labels_no_sera_with_n <- sprintf(
  "%s (n = %s)",
  sr_group_labels_no_sera,
  group_n[names(sr_group_labels_no_sera)]
)
names(sr_group_labels_no_sera_with_n) <- names(sr_group_labels_no_sera)

vac1_plot + 
  add_slope_folddrop(
    sr_group_subset = "D614G",
    nudge = nudges[1]
  ) +
  add_slope_folddrop(
    sr_group_subset = "2x mRNA-1273",
    nudge = nudges[2]
  ) +
  add_slope_folddrop(
    sr_group_subset = "3x mRNA-1273 BD01",
    nudge = nudges[3]
  ) +
  add_slope_folddrop(
    sr_group_subset = "3x mRNA-1273 BD29",
    nudge = nudges[4]
  ) +
  add_slope_folddrop(
    sr_group_subset = "3x mRNA-1273 (6 month)",
    nudge = nudges[5]
  ) +
  add_measured_folddrop(
    sr_group_subset = "D614G",
    nudge = nudges[1]
  ) +
  add_measured_folddrop(
    sr_group_subset = "2x mRNA-1273",
    nudge = nudges[2]
  ) +
  add_measured_folddrop(
    sr_group_subset = "3x mRNA-1273 BD01",
    nudge = nudges[3]
  ) +
  add_measured_folddrop(
    sr_group_subset = "3x mRNA-1273 BD29",
    nudge = nudges[4]
  ) +
  add_measured_folddrop(
    sr_group_subset = "3x mRNA-1273 (6 month)",
    nudge = nudges[5]
  ) +
  scale_x_discrete(
    limits = ag_order,
    labels = ag_who_labels
  ) + 
  scale_color_manual(
    values = sr_group_cols,
    labels = sr_group_labels_no_sera_with_n[c(
      "D614G",
      "2x mRNA-1273",
      "3x mRNA-1273 BD01",
      "3x mRNA-1273 BD29",
      "3x mRNA-1273 (6 month)"
    )],
    limits = c(
      "D614G",
      "2x mRNA-1273",
      "3x mRNA-1273 BD01",
      "3x mRNA-1273 BD29",
      "3x mRNA-1273 (6 month)"
    )
  ) +
  annotate(
    "tile",
    color = agFill(map)[match(ag_order, agNames(map))],
    fill = agFill(map)[match(ag_order, agNames(map))],
    x = ag_order,
    y = -9.25,
    height = 0.5
  ) +
  labs (
    color = "Serum group"
  ) +
  theme(
    axis.title = element_text(size = 9),
    # legend.position = c(0.22, 0.25) - 0.03,
    legend.position = "none",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.background = element_rect(
      # colour = "grey20",
      colour = NA,
      linewidth = 0.2
    ),
    axis.text.x = ggtext::element_markdown(size = 8)
  ) -> gpA

# Plot breadth estimate over time
slope_data %>%
  mutate(
    num_months = c(1, 6, 7, 12)
  ) %>%
  ggplot(
    aes(
      x = num_months,
      y = estimate,
      ymin = lower,
      ymax = upper,
      color = sr_group
    )
  ) + 
  # geom_vline(
  #   xintercept = 6,
  #   linetype = "dashed",
  #   colour = "grey40"
  # ) +
  geom_line(
    color = "grey50",
    size = 0.4,
    linetype = "42"
  ) +
  geom_linerange(
    size = 0.4
  ) +
  geom_point() + 
  scale_x_continuous(
    breaks = c(1, 6, 7, 12)
  ) + 
  scale_y_continuous(
    breaks = seq(from = 0, to = 1, by = 0.1)
  ) + 
  scale_color_manual(
    values = sr_group_cols,
    labels = sr_group_labels_no_sera_with_n[c(
      "D614G",
      "2x mRNA-1273",
      "3x mRNA-1273 BD01",
      "3x mRNA-1273 BD29",
      "3x mRNA-1273 (6 month)"
    )],
    limits = c(
      "D614G",
      "2x mRNA-1273",
      "3x mRNA-1273 BD01",
      "3x mRNA-1273 BD29",
      "3x mRNA-1273 (6 month)"
    )
  ) +
  labs(
    x = "~Months since 2nd\nvaccine dose",
    y = "Fold-drop magnitude\nrelative to D614G sera",
    color = "Serum group"
  ) + 
  coord_cartesian(
    xlim = c(0, 13),
    ylim = c(0, 1),
    expand = FALSE
  ) + 
  theme(
    panel.border = element_rect(
      colour = "grey60",
      fill = NA
    ),
    # axis.text.x = element_text(
    #   angle = 90,
    #   vjust = 0.5,
    #   hjust = 1,
    #   size = 8
    # ),
    panel.background = element_blank(),
    axis.title = element_text(size = 9, lineheight = 1.1),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(
      # linetype = "dotted",
      colour = "grey95"
    ),
    legend.position = c(0.3, -0.5),
    plot.margin = margin(0.25, 0.25, 4.5, 0.5, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(0.2, "cm"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.key = element_blank(),
    # axis.text.x.bottom = element_blank(),
    axis.title.x.bottom = element_text(margin = margin(t = 0.25, unit = "cm"))
  ) -> gpB

# Show full folddrop data
homo_logtiter_table <- matrix(
  unlist(srHomologousLogTiters(map)),
  nrow = numAntigens(map),
  ncol = numSera(map),
  byrow = TRUE
)

# Screenshot slope comparison landscapes
file.remove("figures/main/lndscps/wt_comparison.png")
screenshotWebpage(
  url = "figures/main/lndscps/wt_comparison.html",
  file = "figures/main/lndscps/wt_comparison.png",
  vwidth = 1600,
  vheight = 880,
  cliprect = c(200, 0, 1200, 880)
)

file.remove("figures/main/lndscps/wt_slope_comparison.png")
screenshotWebpage(
  url = "figures/main/lndscps/wt_slope_comparison.html",
  file = "figures/main/lndscps/wt_slope_comparison.png",
  vwidth = 1600,
  vheight = 880,
  cliprect = c(200, 0, 1200, 880)
)

# Add in extra panels
gpC <- cowplot::ggdraw() + cowplot::draw_image("figures/main/lndscps/wt_comparison.png")
gpD <- cowplot::ggdraw() + cowplot::draw_image("figures/main/lndscps/wt_slope_comparison.png")

gpC <- gpC + theme(panel.border = element_rect(colour = "grey80", linewidth = 0.5, fill = NA), plot.margin = margin(0, 0.25, 0.5, 0.25, "cm"))
gpD <- gpD + theme(panel.border = element_rect(colour = "grey80", linewidth = 0.5, fill = NA), plot.margin = margin(0, 0.25, 0.5, 0.25, "cm"))

# Combine plot in panels
saveRDS(
  object = gpA,
  file = "figures/main/fig2a.rds"
)

saveRDS(
  object = gpB,
  file = "figures/main/fig2b.rds"
)

gp <- cowplot::plot_grid(
  cowplot::plot_grid(
    gpA, gpB, 
    ncol = 2,
    labels = c("A", "B"),
    label_x = c(0.02, 0.04),
    label_y = 0.98,
    rel_widths = c(1, 0.4)
  ),
  cowplot::plot_grid(
    gpC, gpD, 
    ncol = 2,
    labels = c("C", "D"),
    label_x = 0.04,
    label_y = 0.98
  ),
  nrow = 2,
  rel_heights = c(1, 0.6),
  label_x = 0.02,
  label_y = 0.99
)

# Save the plot
ggsave(
  plot = gp, 
  filename = "figures/main/fig2-wt_folddrop_comparison.pdf",
  width = 8, 
  height = 8, 
  units = "in"
)

ggsave(
  plot = gp, 
  filename = "figures/main/fig2-wt_folddrop_comparison.png",
  width = 8, 
  height = 8, 
  units = "in"
)

# Plotting the B.1.351 results
vac2_plot <- plot_folddiff(
  all_group_results = all_group_results, 
  sr_groups = c("2x mRNA-1273.351", "B.1.351"),
  ag_names = agNames(map),
  show.legend = TRUE,
  coloring = "sr_group",
  ylim = c(-8.8, 1.5),
  ylab = "Estimated fold drop vs B.1.351"
)

# Get estimates of slope difference
beta_slope_est_data <- readRDS("data/generated_data/slope_comparison_beta.rds")
beta_slope_est_data$ag_folddrops %>%
  rename(
    diff = estimate
  ) %>%
  mutate(
    diff_lower = NA_real_,
    diff_upper = NA_real_
  ) -> folddrop_data

folddrop_data <- bind_rows(
  folddrop_data,
  tibble(ag = "B.1.351", diff = 0)
)

# Set antigen order
ag_order <- folddrop_data$ag[order(-folddrop_data$diff)]

slope_data <- beta_slope_est_data$slope_factors

vac2_plot + 
  add_slope_folddrop(
    sr_group_subset = "B.1.351",
    nudge = -0.15
  ) +
  add_slope_folddrop(
    sr_group_subset = "2x mRNA-1273.351",
    nudge = 0.15
  ) +
  add_measured_folddrop(
    sr_group_subset = "B.1.351",
    nudge = -0.15
  ) +
  add_measured_folddrop(
    sr_group_subset = "2x mRNA-1273.351",
    nudge = 0.15
  ) +
  scale_x_discrete(
    limits = ag_order,
    labels = ag_who_labels
  ) + 
  scale_color_manual(
    values = c(
      "B.1.351" = "#e41a1c",
      "2x mRNA-1273.351" = "#0000FF"
    ),
    labels = sr_group_labels_no_sera[c(
      "B.1.351",
      "2x mRNA-1273.351"
    )]
  ) +
  annotate(
    "tile",
    color = agFill(map)[match(ag_order, agNames(map))],
    fill = agFill(map)[match(ag_order, agNames(map))],
    x = ag_order,
    y = -9.25,
    height = 0.5
  ) +
  labs (
    color = "Serum group"
  ) +
  theme(
    axis.text.x = ggtext::element_markdown(),
    legend.position = c(0.2, 0.21),
    # legend.position = "none",
    legend.background = element_rect(
      colour = "grey20",
      size = 0.2
    )
  ) -> gp

ggsave(
  plot = gp, 
  filename = "figures/som/figS7-beta_foldrop_comparison.pdf", 
  width = 7, 
  height = 4.4, 
  units = "in"
)

ggsave(
  plot = gp, 
  filename = "figures/som/figS7-beta_foldrop_comparison.png", 
  width = 7, 
  height = 4.4, 
  units = "in"
)

# Plot confidence intervals
wt_slope_est_data <- readRDS(
  "data/generated_data/slope_comparison_wt.rds"
)$slope_factors

wt_slope_est_data %>%
  ggplot(
    aes(
      x = factor(sr_group, sr_group),
      y = estimate,
      ymin = lower,
      ymax = upper
    )
  ) + 
  geom_pointrange() + 
  titerplot_theme() + 
  labs(
    x = "Serum group",
    y = "Slope"
  )
