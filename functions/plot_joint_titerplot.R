
# Make a note of excluded sera
plot_joint_titerplot <- function(
    map,
    y_aes,
    highlighted_sera,
    highlighted_sera_cols,
    plottype = "lineplot",
    ag_means_individual_effects = NULL,
    sr_individual_effects = NULL,
    sr_group_gmts = NULL,
    plot_sr_group_gmts = TRUE,
    output_pdf,
    output_png
) {
  
  # Work out antigen order
  ag_order <- sr_group_gmts %>% filter(sr_group == "D614G") %>% arrange(-logtiter) %>% select(ag_name) %>% unlist() %>% unname()
  ag_order <- ag_order[ag_order %in% agNames(map)[agGroups(map) == "wildtype"]]
  omicron_order <- c(
    "BA.1", 
    "BA.1.1", 
    "BA.2", 
    "BA.2.12.1", 
    "BA.3", 
    "BA.4/BA.5"
  )
  
  ag_order <- ag_order[!ag_order %in% omicron_order]
  ag_order <- c(ag_order, omicron_order)
  
  # Subset to not include mutants
  map <- subsetMap(
    map, 
    antigens = agGroups(map) == "wildtype"
  )
  
  # Subset to not include likely second infections
  map <- subsetMap(
    map, 
    sera = !grepl("?2nd", srGroups(map), fixed = T)
  )
  
  # Plot panel A
  map %>% subsetMap(
    sera = !srGroups(map) %in% c(
      "2x mRNA-1273",
      "3x mRNA-1273 BD01",
      "3x mRNA-1273 BD29",
      "3x mRNA-1273 (6 month)",
      "2x mRNA-1273.351"
    )
  ) -> map_subsetA
  
  # Reorder serum groups
  srGroups(map_subsetA) <- factor(
    srGroups(map_subsetA),
    c(
      "D614G",
      "B.1.1.7",
      "B.1.351",
      "P.1",
      "B.1.617.2",
      "B.1.526+E484K",
      "B.1.637",
      "C.37",
      "BA.1",
      "BA.2"
    )
  )
  
  # Plot the titers
  gpA <- plot_map_titers(
    map_subsetA,
    ag_order = ag_order,
    highlighted_sera = highlighted_sera,
    highlighted_sera_cols = highlighted_sera_cols,
    sr_group_gmts = filter(sr_group_gmts, sr_group %in% srGroups(map_subsetA)),
    sr_individual_effects = sr_individual_effects,
    plottype = plottype,
    plot_sr_group_gmts = plot_sr_group_gmts,
    y_aes = y_aes
  )
  
  # Plot panel B
  map %>% subsetMap(
    sera = srGroups(map) %in% c(
      "2x mRNA-1273",
      "3x mRNA-1273 BD01",
      "3x mRNA-1273 BD29",
      "3x mRNA-1273 (6 month)",
      "2x mRNA-1273.351"
    )
  ) -> map_subsetB
  
  # Reorder serum groups
  srGroups(map_subsetB) <- factor(
    srGroups(map_subsetB),
    c(
      "2x mRNA-1273",
      "3x mRNA-1273 BD01",
      "3x mRNA-1273 BD29",
      "3x mRNA-1273 (6 month)",
      "2x mRNA-1273.351"
    )
  )
  
  gpB <- plot_map_titers(
    map_subsetB,
    ag_order = ag_order,
    ncol = 1,
    highlighted_sera = highlighted_sera,
    highlighted_sera_cols = highlighted_sera_cols,
    sr_group_gmts = filter(sr_group_gmts, sr_group %in% srGroups(map_subsetB)),
    sr_individual_effects = sr_individual_effects,
    plottype = plottype,
    plot_sr_group_gmts = plot_sr_group_gmts,
    y_aes = y_aes
  )
  
  gpA <- gpA + theme(plot.margin = margin(l = 0.5, unit = "cm"), plot.tag.position = c(-0.01, 0.98))
  gpB <- gpB + theme(plot.margin = margin(l = 1, unit = "cm"), plot.tag.position = c(-0.02, 0.98))
  
  design <- "
    12
  "
  gp <- gpA + gpB +
    plot_layout(
      design = design,
      widths = c(2, 1),
      heights = c(1, 0.085),
    ) + plot_annotation(
      tag_levels = list(c("A", "B "))
    ) & 
    theme(
      plot.tag = element_text(size = 28)
    )
  
  # Save the plot
  ggsave(output_pdf, gp, width = 14, height = 12, units = "in")
  ggsave(output_png, gp, width = 14, height = 12, units = "in")
  
}
