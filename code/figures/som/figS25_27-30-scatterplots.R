
# Setup workspace
rm(list = ls())
library(tidyverse)
library(Racmacs)
library(titertools)
library(patchwork)

# Load functions
source("functions/scales.R")
source("functions/map_longinfo.R")
source("functions/sr_group_labels.R")

# Read map data
map <- read.acmap("data/maps/map_full_with_extras_no_outliers.ace")

# Get long info
plotdata <- long_map_info(map)

# Convert log2 to fold change
foldchange <- \(x) {
  
  if (is.na(x)) return("")
  xabs <- abs(x)
  foldchange <- 2^xabs
  if (x < 0) foldchange <- -foldchange
  round(foldchange, 1)
  
}

# Function to plot scatterplot
scatterplot_ags <- function(ag1, ag2, sr_group_subset) {

  group_col_data <- distinct(plotdata, sr_group, sr_group_cols)
  sr_group_colors <- group_col_data$sr_group_cols
  names(sr_group_colors) <- group_col_data$sr_group
  
  xlim <- c(0, 10)
  ylim <- xlim
  
  titers1 <- plotdata %>% filter(ag_name == ag1) %>% select(sr_group, sr_num, titer, logtiter) %>% rename(titer1 = titer, logtiter1 = logtiter)
  titers2 <- plotdata %>% filter(ag_name == ag2) %>% select(sr_group, sr_num, titer, logtiter) %>% rename(titer2 = titer, logtiter2 = logtiter)
  
  left_join(
      titers1, 
      titers2, 
      by = c("sr_group", "sr_num")
    ) %>%
    filter(
      sr_group %in% sr_group_subset
    ) -> titerdata
  
  titerdata$sr_group <- factor(
    titerdata$sr_group,
    levels = sr_group_subset
  )
  
  folddiffdata <- readRDS("data/generated_data/all_folddrops.rds") %>%
    filter(
      sr_group %in% sr_group_subset,
      ag1 == !!ag1,
      ag2 == !!ag2
    )
  
  titerdata %>% 
    filter(titer1 != "*" & titer2 != "*") %>% 
    ggplot() + 
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed"
    ) +
    geom_point(
      aes(
        x = logtiter1,
        y = logtiter2,
        color = sr_group
      )
    ) + 
    facet_grid(
      cols = vars(sr_group),
      labeller = as_labeller(sr_group_labels_multiline_with_sera)
    ) +
    scale_x_titer(
      axisname = ag1
    ) + 
    scale_y_titer(
      axisname = ag2
    ) +
    coord_cartesian(
      xlim = xlim,
      ylim = ylim
    ) +
    scale_color_manual(
      values = sr_group_colors
    ) +
    scale_fill_manual(
      values = sr_group_colors
    ) +
    titerplot_theme() + 
    theme(
      legend.position = "none",
      plot.tag.position = c(-0.05, 0.98),
      plot.margin = margin(l = 0.75, b = 0.5, r = 0.2, unit = "cm"),
      strip.text = element_text(size = 8),
      axis.title = element_text(size = 10)
    ) -> gp
  
  label_foldchange <- function(diff, diff_upper, diff_lower) {
    if (is.na(diff)) return("")
    sprintf(
      "%s [%s,%s]",
      foldchange(diff),
      foldchange(diff_upper),
      foldchange(diff_lower)
    )
  }
  
  for (sr_group_name in folddiffdata$sr_group) {
    sr_group_diff <- filter(folddiffdata, sr_group == sr_group_name)
    gp <- gp + geom_ribbon(
      data = tibble(
        sr_group = factor(sr_group_name, levels(titerdata$sr_group)),
        x = extendrange(xlim),
        ymin = extendrange(xlim) + sr_group_diff$folddiff_lower,
        ymax = extendrange(xlim) + sr_group_diff$folddiff_upper
      ),
      aes(
        x = x,
        ymin = ymin,
        ymax = ymax,
        fill = sr_group
      ),
      alpha = 0.1
    ) + 
    geom_line(
      data = tibble(
        sr_group = factor(sr_group_name, levels(titerdata$sr_group)),
        x = extendrange(xlim),
        y = extendrange(xlim) + sr_group_diff$folddiff
      ),
      aes(
        x = x,
        y = y,
        color = sr_group
      ),
      linetype = "dotted"
    ) + 
    geom_text(
      data = sr_group_diff,
      aes(
        x = xlim[1],
        y = ylim[2] - 0.5,
        label = label_foldchange(
          folddiff, folddiff_upper, folddiff_lower
        )
      ),
      size = 3,
      hjust = 0
    )
  }
  
  # Return the plot
  gp
  
}

# Helper functions for plotting
remove_srlabel <- function(x) x + theme(strip.text = element_blank())
space_srlabel <- function(x) x + theme(strip.text = element_text(size = 11, margin = margin(b = 0.5, unit = "cm")))
remove_xaxis <- function(x) x + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) + labs(x = NULL)
  
# Effect of BA.4/BA.5
scatterplotsBA45 <- scatterplot_ags("BA.1", "BA.4/BA.5", c("2x mRNA-1273", "3x mRNA-1273 BD29", "2x mRNA-1273.351", "B.1.617.2", "BA.1"))
ggsave("figures/som/figS25-BA45_scatters.pdf", plot = scatterplotsBA45, width = 10.5, height = 2.8, units = "in")
ggsave("figures/som/figS25-BA45_scatters.png", plot = scatterplotsBA45, width = 10.5, height = 2.8, units = "in")

# Effect of D614G mutation plots
mutants_d614g_A <- scatterplot_ags("D614G", "D614G+E484K", c("D614G", "B.1.351")) + theme(plot.tag.position = c(-0.05, 0.82))
mutants_d614g_B <- scatterplot_ags("D614G", "D614G+E484Q", c("D614G", "B.1.351"))
mutants_d614g_C <- scatterplot_ags("D614G", "D614G+L452R", c("D614G", "B.1.351"))
mutants_d614g_D <- scatterplot_ags("D614G", "D614G+L452R+E484Q", c("D614G", "B.1.351"))
mutants_d614g <- remove_xaxis(space_srlabel(mutants_d614g_A)) / remove_xaxis(remove_srlabel(mutants_d614g_B)) / remove_xaxis(remove_srlabel(mutants_d614g_C)) / remove_srlabel(mutants_d614g_D)
mutants_d614g <- mutants_d614g + patchwork::plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 16))
ggsave("figures/som/figS27-D614G_mutant_scatters.pdf", plot = mutants_d614g, width = 4.6, height = 8.5, units = "in")
ggsave("figures/som/figS27-D614G_mutant_scatters.png", plot = mutants_d614g, width = 4.6, height = 8.5, units = "in")

# Effect of 417 mutation plots
mutants417_A <- scatterplot_ags("P.1", "P.1+T417K", c("D614G", "B.1.351")) + theme(plot.tag.position = c(-0.05, 0.86))
mutants417_B <- scatterplot_ags("B.1.351", "B.1.351+N417K", c("D614G", "B.1.351"))
mutants417_C <- scatterplot_ags("B.1.617.1", "B.1.617.1+K417N", c("D614G", "B.1.351"))
mutants417_D <- scatterplot_ags("B.1.429", "B.1.429+K417N", c("D614G", "B.1.351"))
mutants417 <- space_srlabel(mutants417_A) / remove_srlabel(mutants417_B) / remove_srlabel(mutants417_C) / remove_srlabel(mutants417_D)
mutants417 <- mutants417 + patchwork::plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 16))
ggsave("figures/som/figS28-417_mutant_scatters.pdf", plot = mutants417, width = 4.6, height = 10, units = "in")
ggsave("figures/som/figS28-417_mutant_scatters.png", plot = mutants417, width = 4.6, height = 10, units = "in")

# Effect of 501 mutation plots
mutants501 <- scatterplot_ags("D614G", "D614G+N501Y", c("2x mRNA-1273", "D614G", "B.1.1.7", "P.1", "B.1.351"))
ggsave("figures/som/figS29-501_mutant_scatters.pdf", plot = mutants501, width = 9.5, height = 2.8, units = "in")
ggsave("figures/som/figS29-501_mutant_scatters.png", plot = mutants501, width = 9.5, height = 2.8, units = "in")

# Effect of BA.1+A484K mutation plots
mutantsBA1 <- scatterplot_ags("BA.1", "BA.1+A484K", c("2x mRNA-1273", "2x mRNA-1273.351", "B.1.617.2", "BA.1"))
ggsave("figures/som/figS30-BA1_A484K_scatter.pdf", plot = mutantsBA1, width = 8.5, height = 2.8, units = "in")
ggsave("figures/som/figS30-BA1_A484K_scatter.png", plot = mutantsBA1, width = 8.5, height = 2.8, units = "in")

