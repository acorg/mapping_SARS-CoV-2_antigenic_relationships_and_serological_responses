
# Function to impute titers
impute_titers <- function(titers) {
  titertools:::impute_gmt_titers(
    result = titertools::gmt(titers, ci_method = "quap", dilution_stepsize = 0), 
    titers = titers
  )
}
