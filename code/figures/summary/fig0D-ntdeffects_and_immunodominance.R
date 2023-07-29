
# Setup workspace
rm(list = ls())
set.seed(100)
library(tidyverse)
library(Racmacs)
library(titertools)
library(bayestestR)
library(patchwork)

# Load functions
source("functions/map_longinfo.R")
source("functions/scales.R")
source("functions/homologous_ags.R")
source("functions/sr_group_labels.R")

# Read substitution effect size data
titercomp <- readRDS("data/generated_data/substitution_effect_sizes.rds")

titercomp |> 
  filter(
    sr_group %in% c(
      "2x mRNA-1273",
      "D614G",
      "B.1.617.2",
      "B.1.1.7",
      "P.1",
      "B.1.351"
    ),
    position_diff %in% c(
      # "RBD=", 
      "E484K", 
      "E484Q", 
      "N501Y",
      "K417N", 
      "L452R"
    )
  ) |> 
  mutate(
    position_diff = factor(
      position_diff,
      levels = c(
        # "RBD=", 
        "E484K", 
        "E484Q", 
        "N501Y",
        "K417N", 
        "L452R" 
      )
    )
  ) |> 
  mutate(
    sr_group = factor(
      sr_group,
      rev(c(
        "2x mRNA-1273",
        "D614G",
        "B.1.617.2",
        "B.1.1.7",
        "P.1",
        "B.1.351"
      ))
    )
  ) -> titercomp_subset

titercomp_plotdata <- titercomp_subset |> 
  filter(
    !grepl("Overall", pairname)
  ) |> 
  group_by(
    position_diff, sr_group
  ) |> 
  group_modify(
    ~{
      
      # Calculate the overall folddiff
      folddiff_result <- titertools::log2diff(
        titers1 = unlist(.x$titers1),
        titers2 = unlist(.x$titers2),
        dilution_stepsize = 0,
        # sigma = 0.87,
        ci_method = "HDI"
      )
      
      # Return the results as a tibble
      tibble(
        folddiff = folddiff_result["mean", "estimate"],
        folddiff_upper = folddiff_result["mean", "upper"],
        folddiff_lower = folddiff_result["mean", "lower"],
        sr_group_aas = unique(.x$sr_group_aas)
      )
      
    }
  )

# titercomp_subset |> 
#   filter(
#     !grepl("Overall", pairname)
#   ) |> 
#   group_by(
#     position_diff, sr_group
#   ) |> 
#   summarise(
#     folddiff_min = min(folddiff),
#     folddiff_max = max(folddiff)
#   ) -> titercomp_range
# 
# titercomp_subset |> 
#   filter(
#     grepl("Overall", pairname)
#   ) -> titercomp_overall
# 
# titercomp_overall |> 
#   left_join(
#     titercomp_range,
#     by = c("position_diff", "sr_group")
#   ) -> titercomp_plotdata

titercomp_plotdata |> 
  ggplot(
    aes(
      x = folddiff,
      xmin = folddiff_lower,
      xmax = folddiff_upper,
      y = sr_group,
      color = sr_group_aas
    )
  ) + 
  geom_linerange() + 
  geom_point() + 
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  facet_grid(
    cols = vars(position_diff)
  ) + scale_color_manual(
    values = c(
      "A" = "#1b9e77",
      "E" = "#d95f02",
      "K" = "#7570b3",
      "L" = "#e7298a",
      "N" = "#66a61e",
      "R" = "#e6ab02",
      "T" = "#a6761d",
      "Y" = "#666666"
    )
  ) + 
  coord_cartesian(
    xlim = c(-3.5, 3.5)
  ) +
  scale_x_continuous(
    breaks = -3:3,
    labels = log2foldchange
  ) +
  scale_y_discrete(
    labels = c(
      "2x mRNA-1273" = "2x mRNA-1273",
      "D614G" = "B.1 (D614G)",
      "B.1.617.2" = "Delta (B.1.617.2)",
      "B.1.1.7" = "Alpha (B.1.1.7)",
      "P.1" = "Gamma (P.1)",
      "B.1.351" = "Beta (B.1.351)"
    )
  ) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(fill = NA, colour = "grey60"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.major.x = element_line(colour = "#f6f6f6"),
    strip.background = element_blank(),
    legend.key.size = unit(0.4, "cm"),
    axis.title.x = element_text(margin = margin(t = 0.3, unit = "cm"), size = 10),
    axis.title.y = element_text(margin = margin(r = 0.3, unit = "cm"), size = 10),
    legend.key = element_rect(fill = "white"),
    legend.title = element_text(size = 9, margin = margin(t = 0.3, unit = "cm"))
  ) + 
  labs(
    x = "Mean fold-difference effect of RBD substitution",
    y = "Serum group",
    color = "Amino acid in\neliciting variant"
  ) -> gp

ggsave("figures/summary/fig0D-ntd_and_immunodominance.pdf", gp, width = 9, height = 2, units = "in")
ggsave("figures/summary/fig0D-ntd_and_immunodominance.png", gp, width = 9, height = 2, units = "in")
