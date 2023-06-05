
adjust_titers_by_slope <- function(map, slopes) {
  
  # Adjust slopes of raw titers
  adjustedtitertable <- titerTable(map)
  logtitertable   <- logtiterTable(map)
  titertypestable <- Racmacs:::titertypesTable(map)
  
  # Adjust sera by slopes
  for (sr_group in slopes$sr_group) {
    message("Adjusting slope for ", sr_group)
    for (srnum in which(srGroups(map) == sr_group)) {
      
      titertypes <- titertypestable[,srnum]
      logtiters  <- logtitertable[,srnum]
      colbase    <- max(logtiters, na.rm = T)
      logdrops   <- logtiters - colbase
      
      sr_group_slope <- slopes$estimate[slopes$sr_group == sr_group]
      adjusted_logdrops  <- logdrops / sr_group_slope
      adjusted_logtiters <- colbase + adjusted_logdrops
      adjusted_titers    <- Racmacs:::make_titers((2^adjusted_logtiters)*10, titertypes)
      
      adjustedtitertable[,srnum] <- adjusted_titers
      
    }
  }
  
  # Apply the adjusted titers
  titerTable(map) <- adjustedtitertable
  map
  
}
