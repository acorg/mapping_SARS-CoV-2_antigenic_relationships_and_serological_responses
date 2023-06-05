
# Setup workspace
rm(list = ls())
library(tidyverse)
source("functions/scales.R")
source("functions/sr_group_labels.R")

# Read in foldchange data
folddrops <- readRDS("data/generated_data/all_folddrops.rds")
comparison_folddrops <- readxl::read_excel("data/raw_data/comparison_to_other_studies.xlsx", range = "B3:J43")

# Set dataset assay types
assay_types <- c(
  "Bekliz et al. (2022)" = "live virus",
  "Dejnirattisai et al. (2022)" = "live virus",
  "Rössler et al. (2022)" = "live virus",
  "van der Straten et al. (2022)" = "lentiviral",
  "Evans et al. (2022)" = "lentiviral",
  "Zhou et al. (2022)" = "lentiviral",
  "Wang et al. (2022)" = "lentiviral"
)

TMPRSS2 <- c(
  "Bekliz et al. (2022)" = F,
  "Dejnirattisai et al. (2022)" = F,
  "Rössler et al. (2022)" = T,
  "van der Straten et al. (2022)" = F,
  "Evans et al. (2022)" = F,
  "Zhou et al. (2022)" = F,
  "Wang et al. (2022)" = T
)

comparison_folddrops |> 
  pivot_longer(
    cols = 3:9,
    names_to = "dataset",
    values_to = "folddiff"
  ) |> 
  separate_wider_delim(
    Variants,
    delim = " vs ",
    names = c("ag1", "ag2")
  ) |> 
  rename(
    sr_group = "Serum group"
  ) |> 
  mutate(
    assay_type = assay_types[dataset]
  ) -> comparison_folddrops_long

# Refactor serum groups
comparison_folddrops_long$sr_group <- fct_recode(
  comparison_folddrops_long$sr_group,
  "3x mRNA-1273 BD29" = "4 weeks post 3x mRNA-1273",
  "2x mRNA-1273" = "4 weeks post 2x mRNA-1273"
)

# Filter folddrop data
folddrops |> 
  filter(
    paste(ag1, ag2, sr_group) %in% paste(comparison_folddrops_long$ag1, comparison_folddrops_long$ag2, comparison_folddrops_long$sr_group),
    !is.na(folddiff)
  ) -> folddrops_subset

comparison_folddrops_long |> 
  filter(
    paste(ag1, ag2, sr_group) %in% paste(folddrops_subset$ag1, folddrops_subset$ag2, folddrops_subset$sr_group)
  ) -> comparison_folddrops_long

# Set x axis order
comparison_folddrops_long <- comparison_folddrops_long |> mutate(yorder = paste(ag2, sr_group))
folddrops_subset <- folddrops_subset |> mutate(yorder = paste(ag2, sr_group))
yaxis_order <- folddrops_subset |> arrange(folddiff) |> pluck("yorder")
comparison_folddrops_long <- comparison_folddrops_long |> mutate(yorder = factor(yorder, yaxis_order))
folddrops_subset <- folddrops_subset |> mutate(yorder = factor(yorder, yaxis_order))

dataset_order <- c(
  "Bekliz et al. (2022)",
  "Dejnirattisai et al. (2022)*",
  "Rössler et al. (2022)",
  "van der Straten et al. (2022)",
  "Evans et al. (2022)",
  "Zhou et al. (2022)",
  "Wang et al. (2022)",
  "Our data"
)

comparison_folddrops_long |> 
  mutate(
    sr_group = factor(
      sr_group,
      c(
        "2x mRNA-1273",
        "3x mRNA-1273 BD29",
        "B.1.1.7",
        "B.1.351",
        "P.1",
        "B.1.617.2",
        "BA.1"
      )
    ),
    dataset = replace(dataset, dataset == "Dejnirattisai et al. (2022)", "Dejnirattisai et al. (2022)*")
  ) |> 
  ggplot(
    aes(
      x = folddiff,
      y = yorder,
      color = factor(dataset, dataset_order),
      shape = factor(dataset, dataset_order),
      size = factor(dataset, dataset_order)
    )
  ) + 
  geom_point(
    shape = "circle open",
    size = 2,
    stroke = 1
  ) + 
  geom_errorbar(
    data = mutate(folddrops_subset, dataset = "Our data"),
    mapping = aes(
      x = folddiff,
      xmin = folddiff_lower,
      xmax = folddiff_upper
    ),
    color = "black",
    alpha = 0.8,
    width = 0.4,
    linewidth = 0.5,
    show.legend = FALSE
  ) +
  geom_point(
    data = mutate(folddrops_subset, dataset = "Our data"),
    mapping = aes(
      x = folddiff
    )
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  coord_cartesian(
    xlim = range(c(
      comparison_folddrops_long$folddiff,
      folddrops_subset$folddiff
    ), na.rm = T)
  ) +
  facet_grid(
    rows = vars(sr_group),
    cols = vars(assay_type),
    scales = "free_y",
    space = "free_y",
    margins = "assay_type",
    labeller = labeller(.rows = as_labeller(sr_group_labels_multiline_with_sera))
  ) +
  scale_x_continuous(
    labels = log2foldchange,
    breaks = \(x) ceiling(min(x)):floor(max(x))
  ) +
  scale_y_discrete(
    labels = \(x) {
      gsub("(^[^ ]*).*$", "\\1", x)
    }
  ) +
  scale_color_manual(
    values = c(
      "Bekliz et al. (2022)" = "#e41a1c",
      "Dejnirattisai et al. (2022)*" = "#377eb8",
      "Rössler et al. (2022)" = "#4daf4a",
      "van der Straten et al. (2022)" = "#984ea3",
      "Evans et al. (2022)" = "#ff7f00",
      "Zhou et al. (2022)" = "#a65628",
      "Wang et al. (2022)" = "#f781bf",
      "Our data" = "#000000"
    ),
    limits = dataset_order
  ) +
  scale_shape_manual(
    values = c(
      "Bekliz et al. (2022)" = "circle open",
      "Dejnirattisai et al. (2022)*" = "circle open",
      "Rössler et al. (2022)" = "circle open",
      "van der Straten et al. (2022)" = "circle open",
      "Evans et al. (2022)" = "circle open",
      "Zhou et al. (2022)" = "circle open",
      "Wang et al. (2022)" = "circle open",
      "Our data" = "circle"
    ),
    limits = dataset_order
  ) +
  scale_size_manual(
    values = c(
      "Bekliz et al. (2022)" = 2,
      "Dejnirattisai et al. (2022)*" = 2,
      "Rössler et al. (2022)" = 2,
      "van der Straten et al. (2022)" = 2,
      "Evans et al. (2022)" = 2,
      "Zhou et al. (2022)" = 2,
      "Wang et al. (2022)" = 2,
      "Our data" = 1.5
    ),
    limits = dataset_order
  ) +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "grey60", fill = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey95"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.text.y = element_text(angle = 0),
    strip.background = element_rect(fill = "#eeeeee"),
    axis.title.x = element_text(margin = margin(t = 0.2, unit = "cm")),
    legend.position = "top"
  ) + 
  labs(
    x = "Fold difference relative to D614G (*or 614D where D614G not available)",
    y = "Variant",
    color = "Dataset",
    shape = "Dataset",
    size = "Dataset"
  ) -> gp

ggsave(
  plot = gp,
  filename = "figures/som/figS6-forestplot_comparison.pdf",
  width = 10,
  height = 8.5,
  units = "in"
)

ggsave(
  plot = gp,
  filename = "figures/som/figS6-forestplot_comparison.png",
  width = 10,
  height = 8.5,
  units = "in"
)
