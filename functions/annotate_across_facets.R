
annotate_across_facets <- function(
    geom,
    x,
    y,
    facet_var,
    facet_value,
    ...
    ) {

  datatibble <- tibble(
    x = x,
    y = y
  )
  
  datatibble[[facet_var]] <- facet_value
  
  layer(
    geom = geom,
    params = list(
      na.rm = FALSE, 
      ...
    ),
    stat = StatIdentity,
    position = PositionIdentity,
    data = datatibble,
    mapping = aes_all(c("x", "y")),
    inherit.aes = FALSE,
    show.legend = FALSE
  )
  
}

