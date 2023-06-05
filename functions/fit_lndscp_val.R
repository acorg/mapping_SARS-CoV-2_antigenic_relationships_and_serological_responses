
fit_lndscp_val <- function(x, y, sr_coords, sr_colbases, sr_slope) {
  
  sr_distances <- as.matrix(dist(rbind(c(x, y), sr_coords)))[1, -1]
  mean(sr_colbases - sr_distances*sr_slope)
  
}
