
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

# Read map data
map <- read.acmap("data/maps/map_full_no_outliers.ace")

# Get homologous sequences for each serum group
sr_group_seqs <- tibble(
  sr_group = names(homologous_ags),
  sr_sequences = unname(as.list(as.data.frame(t(agSequences(map)[match(homologous_ags, agNames(map)),]))))
) %>% 
  filter(
    sr_group %in% srGroups(map)
  ) %>%
  mutate(
    sr_group = factor(sr_group, levels(srGroups(map)))
  )

# Convert to long data form
map_longdata <- long_map_info(map)

# Fetch sequence matrix
ag_seq_matrix <- agSequences(map)
rownames(ag_seq_matrix) <- agNames(map)
colnames(ag_seq_matrix) <- seq_len(ncol(ag_seq_matrix))

# Get pairs that differ by only a given position in the RBD
rbd_sites <- 331:531
ntd_sites <- 14:330

# Set base sequence comparison
base_seq <- ag_seq_matrix["D614G",]

  
# Work out pairs of comparison strains
comp <- list()
for (i in 1:(nrow(ag_seq_matrix) - 1)) {
  for (j in (i+1):nrow(ag_seq_matrix)) {
    
    variant1 <- rownames(ag_seq_matrix)[i]
    variant2 <- rownames(ag_seq_matrix)[j]
    
    seq1 <- ag_seq_matrix[variant1, rbd_sites]
    seq2 <- ag_seq_matrix[variant2, rbd_sites]
    
    mismatches <- which(seq1 != seq2)
    mismatch_positions <- names(mismatches)
    
    ntd_diffs <- sum(ag_seq_matrix[variant1, ntd_sites] != ag_seq_matrix[variant2, ntd_sites]) > 0
    
    if (length(mismatch_positions) == 0) {
      
      comp$null <- rbind(comp$null, c(variant1, variant2, ntd_diffs))
      
    } else if (length(mismatch_positions) == 1) {
      
      aa1 <- seq1[mismatch_positions]
      aa2 <- seq2[mismatch_positions]
      
      variants <- c(variant1, variant2)
      variant_aas <- c(aa1, aa2)
      
      # Set ordering of aas
      base_aa <- base_seq[mismatch_positions]
      if (aa1 == base_aa || aa2 == base_aa) {
        pair_order <- rev(order(c(
          aa1 == base_aa, 
          aa2 == base_aa
        )))
      } else {
        pair_order <- c(1, 2)
      }
      
      label <- paste0(
        variant_aas[pair_order][1],
        mismatch_positions,
        variant_aas[pair_order][2]
      )
      
      comp[[label]] <- rbind(
        comp[[label]],
        c(variants[pair_order], ntd_diffs)
      )
    }
    
  }
}

comporder <- c("null", sort(names(comp)[names(comp) != "null"]))
comp <- comp[comporder]

# Assemble table of data
titercompdata <- do.call(
  bind_rows,
  lapply(seq_along(comp), \(x) {
    tibble(
      position_diff = names(comp)[x],
      ag1 = comp[[x]][,1],
      ag2 = comp[[x]][,2],
      ntd_diffs = as.logical(comp[[x]][,3]),
      pairname = paste(ag2, ag1, sep = " - ")
    )
  })
)

# Bind in folddrop data
folddrop_data <- readRDS("data/generated_data/all_folddrops.rds")
titercompdata %>%
  left_join(
    folddrop_data,
    by = c("ag1", "ag2"),
    multiple = "all"
  ) %>%
  filter(
    sr_group %in% srGroups(map)
  ) -> titercompdata

# Add in raw titer data
titer_table <- titerTable(map)
sr_groups <- as.character(srGroups(map))
titercompdata %>%
  mutate(
    titers1 = lapply(seq_along(ag1), \(i) titer_table[ag1[i], sr_groups == sr_group[i]]),
    titers2 = lapply(seq_along(ag2), \(i) titer_table[ag2[i], sr_groups == sr_group[i]])
  ) -> titercompdata

# Summarise overall fold differences by substitution
titercompdata %>%
  group_by(
    position_diff,
    sr_group
  ) %>%
  summarise(
    .groups = "drop",
    foldchange_stat = list(
      titertools::log2diff(
        titers1 = unlist(titers1),
        titers2 = unlist(titers2),
        dilution_stepsize = 0,
        ci_method = "quap"
      )
    )
  ) %>%
  mutate(
    folddiff = vapply(foldchange_stat, \(x) x["mean", "estimate"], numeric(1)),
    folddiff_upper = vapply(foldchange_stat, \(x) x["mean", "upper"], numeric(1)),
    folddiff_lower = vapply(foldchange_stat, \(x) x["mean", "lower"], numeric(1))
  ) %>%
  mutate(
    pairname = "Overall",
    ntd_diffs = TRUE
  ) -> titercompdata_overall

# Bind both
titercompdata <- bind_rows(
  titercompdata,
  titercompdata_overall
)

# Add additional info and ordering
titercompdata %>%
  mutate(
    position_diff_positions = as.numeric(stringr::str_extract(position_diff, "[0-9]{3}"))
  ) -> titercompdata


# Rename null to RBD equal
titercompdata$position_diff <- forcats::fct_recode(
  titercompdata$position_diff,
  `RBD=` = "null"
)

# Add additional info on overall
overall_subset <- titercompdata$pairname == "Overall" & titercompdata$position_diff != "RBD="
titercompdata$pairname[overall_subset] <- sprintf(
  "Overall (%s%s - %s%s)",
  titercompdata$position_diff_positions[overall_subset],
  substr(titercompdata$position_diff[overall_subset], 1, 1),
  titercompdata$position_diff_positions[overall_subset],
  substr(
    titercompdata$position_diff[overall_subset], 
    nchar(as.character(titercompdata$position_diff[overall_subset])), 
    nchar(as.character(titercompdata$position_diff[overall_subset]))
  )
)

# Order pairnames
titercompdata %>%
  mutate(
    pairname = factor(
      pairname, 
      rev(unique(titercompdata$pairname[order(titercompdata$ntd_diffs, match(titercompdata$ag1, agNames(map)))]))
    )
  ) -> titercompdata

# Order position differences
titercompdata %>%
  mutate(
    position_diff = factor(
      position_diff, 
      unique(titercompdata$position_diff[order(titercompdata$position_diff_positions)])
    )
  ) -> titercompdata

titercompdata %>%
  filter(!is.na(folddiff)) -> titercompdata

# Merge sequence data
titercompdata %>%
  left_join(
    sr_group_seqs,
    by = "sr_group"
  ) -> titercompdata

titercompdata %>%
  mutate(
    sr_group_aas = vapply(
      seq_along(position_diff_positions), function(i) {
        sr_sequences[[i]][position_diff_positions[i]]
      }, character(1)
    )
  ) -> titercompdata

# Set serum group labels
sr_group_labels <- c(
  "D614G" = "D614G\nsera",
  "B.1.1.7" = "B.1.1.7\nsera",
  "P.1" = "P.1\nsera",
  "B.1.351" = "B.1.351\nsera",
  "B.1.526+E484K" = "B.1.526+E484K\nsera",
  "B.1.617.2" = "B.1.617.2\nsera",
  "B.1.637" = "B.1.637\nsera",
  "C.37" = "C.37\nsera",
  "BA.1" = "BA.1\nsera",
  "BA.2" = "BA.2\nsera",
  "2x mRNA-1273" = "4 weeks post\n2x mRNA-1273",
  "3x mRNA-1273 BD01" = ">3 months post\n2x mRNA-1273",
  "3x mRNA-1273 BD29" = "4 weeks post\n3x mRNA-1273",
  "3x mRNA-1273 (6 month)" = ">3 months post\n3x mRNA-1273",
  "2x mRNA-1273.351" = "4 weeks post\n2x mRNA-1273.351"
)

# Function for doing the actual plotting
dosubplot <- function(x, add_color_scale = TRUE, plot_widths = c(0.35, 1)) {
  
  x %>%
    mutate(
      ntd_equal = !ntd_diffs
    ) -> x
  
  x %>% 
    filter(
      grepl("Overall", pairname)
    ) -> x_overall
  
  x %>%
    filter(
      !grepl("Overall", pairname)
    ) -> x
  
  x %>%
    ggplot(
      aes(
        x = folddiff,
        y = pairname,
        xmin = folddiff_lower,
        xmax = folddiff_upper,
        shape = ntd_equal,
        color = sr_group_aas
      )
    ) +
    scale_x_continuous(
      breaks = function(x) ceiling(min(x)):floor(max(x)),
      labels = log2foldchange
    ) + 
    scale_y_discrete() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey20") +
    geom_vline(
      data = x_overall,
      aes(xintercept = folddiff, color = sr_group_aas),
      alpha = 1,
      linetype = "22",
      size = 0.5,
      show.legend = FALSE
    ) +
    geom_linerange(
      aes(alpha = ntd_equal),
      size = 1,
      show.legend = FALSE
    ) + 
    geom_point(
      aes(alpha = ntd_equal),
      size = 2,
      stroke = 0.8
    ) + 
    geom_point(
      data = filter(x, !ntd_equal),
      size = 1,
      color = "white",
      shape = "circle",
      show.legend = FALSE
    ) + 
    facet_grid(
      rows = vars(position_diff),
      cols = vars(sr_group),
      labeller = labeller(.cols = sr_group_labels),
      scales = "free_y",
      space = "free_y"
    ) + 
    coord_cartesian(
      xlim = c(-5, 5)
    ) +
    scale_shape_manual(
      values = c(
        "TRUE" = "circle",
        "FALSE" = "circle open"
      )
    ) +
    scale_alpha_manual(
      guide = guide_none(),
      values = c(
        "TRUE" = 1,
        "FALSE" = 0.4
      )
    ) +
    labs(
      x = "Fold difference (Variant B - Variant A)",
      y = "",
      color = "Amino acid",
      alpha = "NTD equal",
      shape = "NTD equal"
    ) +
    titerplot_theme() +
    theme(
      strip.text.y = element_text(
        size = 10,
        angle = 0
      ),
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5
      ),
      legend.background = element_blank(),
      legend.justification = c(0, 0),
      legend.key = element_blank(),
      legend.key.width = unit(.3, "cm"),
      legend.title = element_text(
        size = 10
      ), 
      legend.text = element_text(
        size = 8
      ),
      legend.box.background = element_rect(
        fill = NA,
        colour = "grey60"
      )
    ) -> gp
  
  if (add_color_scale) {
    
    gp <- gp + scale_color_manual(
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
    )
    
  } else {
    
    gp <- gp + scale_color_manual(
      values = c(
        "NA" = "black"
      )
    )
    
  }
  
  # Add line separating overall
  linedata <- x %>% group_by(position_diff) %>% summarise(
    n_total = length(unique(pairname)),
    n_nontd = length(unique(pairname[!ntd_diffs]))
  )
  gp + geom_hline(
    aes(yintercept = n_nontd + 0.5),
    data = linedata,
    linetype = "longdash", 
    color = "grey90"
  )
  
  # Add separate plot indicating variant A and variant B
  x %>%
    distinct(
      position_diff,
      ag1,
      ag2,
      pairname
    ) %>% 
    pivot_longer(
      cols = c("ag1", "ag2"),
      names_to = "ag_AB",
      values_to = "ag_name"
    ) %>%
    mutate(
      ag_AB = factor(
        ifelse(ag_AB == "ag1", "Variant A", "Variant B"),
        c("Variant A", "Variant B")
      )
    ) -> pairdata
  
  pairdata %>%
    ggplot(
      aes(
        y = pairname,
        label = ag_name
      )
    ) + 
    geom_text(
      aes(
        x = 1
      ),
      hjust = 0,
      size = 2.8
    ) +
    coord_cartesian(
      xlim = c(1, 1.05)
    ) +
    facet_grid(
      rows = vars(position_diff),
      cols = vars(ag_AB),
      scales = "free_y",
      space = "free_y"
    ) + 
    theme(
      strip.background = element_blank(),
      strip.text = element_text(
        size = 10
      ),
      strip.text.y = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.background = element_blank(),
      plot.margin = margin(l=0.25, r = 0, unit = "cm"),
      panel.spacing.x = unit(0, "cm")
    ) -> gp_pairnames
  
  gp_pairnames + 
    (gp + theme(axis.text.y = element_blank()) + labs(y = NULL)) + 
    plot_layout(widths = plot_widths)
  
}

doplot <- function(x, include_figB = TRUE, plotheights = c(1, 1), plot_widths = c(0.35, 1)) {
  
  figA <- dosubplot(filter(x, position_diff != "RBD="), plot_widths = plot_widths)
  figA <- figA & theme(plot.tag.position = c(0, 0.995))
  if (!include_figB) return(figA)
  
  figB <- dosubplot(filter(x, position_diff == "RBD=") %>% mutate(sr_group_aas = "NA"), add_color_scale = FALSE, plot_widths = plot_widths)
  figB <- figB & theme(strip.text.y = element_blank(), legend.position = "none", plot.tag.position = c(0, 0.975))
  fig <- figA / figB + plot_layout(
    heights = plotheights,
    guides = 'keep'
  ) + plot_annotation(tag_levels = list(c("A", "", "B", "")))
  fig <- fig & theme(
    plot.tag = element_text(size = 18)
  )
  fig
  
}

titercompdata %>%
  mutate(
    sr_group = factor(
      sr_group,
      names(sr_group_labels_no_sera)
    )
  ) -> titercompdata

titercompdata %>%
  mutate(
    position_diff = factor(
      position_diff,
      levels = c(
        "RBD=", 
        "E484K", 
        "E484Q", 
        "N501Y",
        "K417N", 
        "L452R" 
      )
    )
  ) %>%
  filter(
    sr_group %in% c(
      "2x mRNA-1273",
      "D614G",
      "B.1.617.2",
      "B.1.1.7",
      "P.1",
      "B.1.351"
    ),
    !is.na(position_diff)
  ) %>%
  mutate(
    sr_group = factor(
      sr_group,
      c(
        "2x mRNA-1273",
        "D614G",
        "B.1.617.2",
        "B.1.1.7",
        "P.1",
        "B.1.351"
      )
    )
  ) -> titercompdata_subset

# Save the effect size estimates
titercompdata_subset %>%
  select(
    position_diff,
    sr_group,
    pairname,
    folddiff,
    folddiff_lower,
    folddiff_upper
  ) %>%
  rename(
    diff = folddiff,
    diff_lower = folddiff_lower,
    diff_upper = folddiff_upper
  ) %>%
  mutate(
    foldchange = log2foldchange(diff),
    foldchange_upper = log2foldchange(diff_upper),
    foldchange_lower = log2foldchange(diff_lower)
  ) %>%
  arrange(
    !grepl("Overall", pairname),
    position_diff
  ) %>%
  readr::write_csv(
    "data/generated_data/single_substitution_effect_sizes.csv"
  )

cropped_values <- titercompdata_subset |> filter(sr_group == "B.1.617.2", ag1 == "B.1.617.2", ag2 == "B.1.617.2 (AY.3)+E484Q")
log2foldchange(cropped_values$folddiff)

gp <- doplot(titercompdata_subset, plotheights = c(1, 0.26))
ggsave("figures/main/fig6-ntd_and_immunodominance.pdf", gp, width = 13.2, height = 9.5, units = "in")
ggsave("figures/main/fig6-ntd_and_immunodominance.png", gp, width = 13.2, height = 9.5, units = "in")

gp_full <- doplot(titercompdata, plot_widths = c(0.15, 1), include_figB = FALSE)
ggsave("figures/som/figS30-ntd_and_immunodominance_full.pdf", gp_full, width = 24, height = 9.5, units = "in")
ggsave("figures/som/figS30-ntd_and_immunodominance_full.png", gp_full, width = 24, height = 9.5, units = "in")
