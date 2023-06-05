
# Function to plot the map titers
plot_map_titers <- function(
    map, 
    ag_order = NULL, 
    ncol = 2, 
    annotate_sr_groups = TRUE,
    sr_individual_effects,
    highlighted_sera,
    highlighted_sera_cols,
    plottype = "titerplot",
    sr_group_gmts = NULL,
    plot_sr_group_gmts = TRUE,
    y_aes
  ) {
  
  # Set the antigen x axis order
  if (is.null(ag_order)) {
    ag_order <- sr_group_gmts %>% filter(sr_group == "D614G") %>% arrange(-logtiter) %>% select(ag_name) %>% unlist() %>% unname()
    omicron_order <- c("BA.1", "BA.1.1", "BA.2", "BA.2.12.1", "BA.3", "BA.4/BA.5")
    ag_order <- ag_order[!ag_order %in% omicron_order]
    ag_order <- c(ag_order, omicron_order)
    ag_order <- ag_order[ag_order %in% agNames(map)]
  }
  
  # Get long info
  plotdata <- long_map_info(map)
  
  # Add in estimated sera individual effects
  if (!is.null(sr_individual_effects)) {
    
    plotdata %>% 
      left_join(
        sr_individual_effects,
        by = "sr_name"
      ) %>%
      mutate(
        adjusted_titer = Racmacs:::reactivity_adjust_titers(
          titers = titer,
          adjustment = -sr_individual_effect
        )
      ) -> plotdata
    
    # Impute non-detectable titers
    plotdata %>%
      group_by(
        sr_group, ag_name
      ) %>%
      group_modify(
        ~ .x %>%
          mutate(
            imputed_titer = titer %>%
              Racmacs:::reactivity_adjust_titers(-sr_individual_effect) %>%
              impute_titers() %>%
              Racmacs:::reactivity_adjust_titers(sr_individual_effect)
          )
      ) %>%
      mutate(
        imputed_logtiter = Racmacs:::log_titers(imputed_titer, unique(dilution_stepsize))
      ) %>%
      ungroup() -> plotdata
    
    plotdata %>%
      mutate(
        imputed_logtiter_minus_individual_effect = imputed_logtiter - sr_individual_effect
      ) %>%
      mutate(
        logtiter = pmax(0.5, logtiter),
        imputed_logtiter_minus_individual_effect = pmax(0.5, imputed_logtiter_minus_individual_effect)
      ) -> plotdata
    
  }
  
  # Clamp titers to <20 and indicate they are <20
  plotdata$lessthan20 <- FALSE
  plotdata$lessthan20[plotdata$logtiter < 1] <- TRUE
  plotdata$logtiter[plotdata$logtiter < 1] <- 0.5
  
  if (!is.null(sr_group_gmts)) {
    sr_group_gmts$lessthan20 <- FALSE
    sr_group_gmts$lessthan20[sr_group_gmts$logtiter < 1] <- TRUE
    sr_group_gmts$logtiter[sr_group_gmts$logtiter < 1] <- 0.5
  }
  
  sr_group_line_colors <- srGroupOutline(map, 0.2)
  sr_group_line_colors["mRNA-1273"] <- "#808080"
  sr_group_line_colors["B.1.351"] <- "#e6b800"

  
  # Function to add trace for subjects
  add_trace <- function(
    data, 
    linecolor, 
    pointcolor = linecolor,
    alpha,
    size
    ) {
    
    list(
      geom_line(
        data = filter(data, !is.na(logtiter)),
        aes(
          group = sr_name,
          y = .data[[y_aes]]
        ),
        alpha = alpha*0.7,
        linewidth = size,
        color = linecolor,
        linetype = "11",
        na.rm = TRUE
      ),
      geom_line(
        data = data,
        aes(
          group = sr_name,
          y = .data[[y_aes]]
        ),
        alpha = 1,
        linewidth = size,
        color = "white",
        na.rm = TRUE
      ),
      geom_line(
        data = data,
        aes(
          group = sr_name,
          y = .data[[y_aes]]
        ),
        alpha = alpha,
        linewidth = size,
        color = linecolor,
        na.rm = TRUE
      ),
      geom_point(
        data = filter(data, !is.na(logtiter)),
        aes(
          y = .data[[y_aes]]
        ),
        alpha = alpha,
        color = pointcolor
      )
    )
    
  }
  
  # Do the plot
  plotdata %>%
    ggplot(
      aes(
        x = ag_name,
        y = .data[[y_aes]],
        color = sr_group,
        fill = sr_group
      )
    ) -> gp
  
  if (plottype == "lineplot") {
    gp <- gp + add_trace(
      data = filter(plotdata, !sr_name %in% unlist(highlighted_sera)),
      linecolor = "grey70",
      alpha = 0.3,
      size = 0.6
    )
  }
  
  if (plottype == "boxplot") {
    gp <- gp + 
      ggbeeswarm::geom_quasirandom(
        color = "black",
        size = 1,
        alpha = 0.5
      ) + geom_boxplot(
        outlier.shape = NA,
        alpha = 0.6
      )
  }
  
  gp <-  gp +
    scale_color_manual(
      values = sr_group_line_colors
    ) +
    scale_fill_manual(
      values = srGroupOutline(map)
    ) +
    scale_y_titer(
      ymin = 1
    ) + 
    scale_x_discrete(
      limits = ag_order,
      labels = ag_who_labels
    ) +
    facet_wrap(
      vars(sr_group),
      ncol = ncol
    ) + 
    labs(
      x = ""
    ) +
    coord_cartesian(
      ylim = c(-1, 14.5),
      expand = FALSE
    ) +
    titerplot_theme() + 
    theme(
      legend.position = "none",
      strip.text = element_blank(),
      axis.text.x = ggtext::element_markdown()
    ) -> gp
  
  # Highlight any additional sera
  for (i in seq_along(highlighted_sera)) {
    
    gp <- gp + add_trace(
      data = filter(plotdata, sr_name %in% highlighted_sera[[i]]),
      linecolor = highlighted_sera_cols[i],
      alpha = 0.4,
      size = 0.6
    )
    
  }
  
  # Add sera group gmts
  if (plot_sr_group_gmts) {
    
    # Filter averages to only cases where individual titers are shown
    plotdata %>% 
      filter(!is.na(logtiter)) %>% 
      distinct(sr_name, ag_name, sr_group) %>% 
      group_by(sr_group, ag_name) %>% 
      count() -> sr_group_ag_counts
    
    sr_group_gmt_subset <- filter(
      sr_group_gmts, 
      ag_name %in% unique(filter(plotdata, !is.na(logtiter))$ag_name)
    )
    
    sr_group_gmt_subset %>% 
      left_join(sr_group_ag_counts, by = c("sr_group", "ag_name")) %>%
      mutate(
        logtiter = replace(logtiter, n == 0, NA)
      ) -> sr_group_gmt_subset
    
    # Add the averages
    gp <- gp + 
      geom_line(
        data = filter(sr_group_gmt_subset, !is.na(logtiter)),
        aes(group = sr_group, y = logtiter),
        size = 0.8,
        alpha = 0.4,
        linetype = "21",
        na.rm = TRUE
      ) +
      geom_line(
        data = sr_group_gmt_subset,
        aes(group = sr_group, y = logtiter),
        size = 0.8,
        na.rm = TRUE
      ) +
      geom_point(
        data = filter(sr_group_gmt_subset, !is.na(logtiter)),
        aes(y = logtiter),
        size = 2.5
      ) +
      geom_point(
        data = filter(sr_group_gmt_subset, !is.na(logtiter)),
        aes(y = logtiter),
        size = 1.5,
        fill = "white",
        shape = 21
      )
    
  }
  
  # Add n = annotations
  if (annotate_sr_groups) {
    
    gp <- gp +
      geom_text(
        data = distinct(plotdata, sr_group),
        aes(
          x = length(ag_order) + 0.2,
          y = 13.5,
          label = sr_group_labels[as.character(sr_group)]
        ),
        size = 3,
        hjust = 1,
        # fontface = "bold",
        color = "black"
      ) +
      geom_text(
        data = plotdata %>% 
          filter(!sr_name %in% unlist(highlighted_sera)) %>% 
          distinct(sr_name, sr_group) %>% 
          group_by(sr_group) %>% 
          count(),
        aes(
          x = length(ag_order) + 0.2,
          y = 12.25,
          label = paste("n =", n)
        ),
        size = 3,
        hjust = 1,
        color = "black"
      )
    
    for (i in seq_along(highlighted_sera)) {
      
      yval <- 12.25 - 1.1*i
      gp <- gp + geom_text(
        data = plotdata %>% 
          filter(sr_name %in% unlist(highlighted_sera[[i]])) %>% 
          distinct(sr_name, sr_group) %>% 
          group_by(sr_group) %>% 
          count() %>% 
          mutate(logtiter = yval),
        aes(
          x = length(ag_order) + 0.2,
          label = paste("n =", n)
        ),
        size = 3,
        hjust = 1,
        color = highlighted_sera_cols[i]
      )
      
    }
    
  }
  
  # Annotate region below detectable
  gp <- gp +
    annotate(
      "tile",
      x = agNames(map),
      y = 0,
      height = 2,
      fill = "grey50",
      color = NA,
      alpha = 0.3
    )
  
  # Annotate colors for each antigen
  for (n in seq_len(numAntigens(map))) {
    gp <- gp +
      annotate(
        "tile",
        x = agNames(map)[n],
        y = -0.75,
        height = 1,
        fill = agFill(map)[n],
        color = NA
      )
  }
  
  # Return the plot
  gp
  
}

