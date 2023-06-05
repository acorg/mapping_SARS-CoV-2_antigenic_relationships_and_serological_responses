
#' Function to calculate the residual error
#' @export
residualErrorTable <- function(map) {

  # Work out predicted and measured log titers
  predicted_logtiter <- Racmacs:::colBaseMatrix(map) - mapDistances(map)
  measured_logtiter  <- logtiterTable(map)
  predicted_titer    <- as.character(round(2^predicted_logtiter *10, 2))
  titer_types        <- Racmacs:::titertypesTable(map)
  measured_logtiter_lower <- measured_logtiter - dilutionStepsize(map) / 2
  measured_logtiter_upper <- measured_logtiter + dilutionStepsize(map) / 2
  measured_logtiter_lower[titer_types == 2] <- -Inf
  measured_logtiter_lower[titer_types == 3] <- Inf

  titer_category <- factor(
    titer_types,
    levels = -1:3,
    labels = c(".", "*", "=", "<", ">")
  )

  tibble(
    ag_num                  = factor(as.vector(Racmacs:::agNumMatrix(map))),
    sr_num                  = factor(as.vector(Racmacs:::srNumMatrix(map))),
    titer_type              = titer_category,
    predicted_titer         = as.vector(predicted_titer),
    measured_titer          = as.vector(titerTable(map)),
    predicted_logtiter      = as.vector(predicted_logtiter),
    measured_logtiter_lower = as.vector(measured_logtiter_lower),
    measured_logtiter_upper = as.vector(measured_logtiter_upper),
    residual_lower          = measured_logtiter_lower - predicted_logtiter,
    residual_upper          = measured_logtiter_upper - predicted_logtiter
  )

}

impute_table_residuals <- function(residual_data, group) {

  residual_data %>%
    mutate(
      residual_error_imputed = fit_imputed_censored_normal(
        lower_lims = residual_lower,
        upper_lims = residual_upper,
        mean = mean(residual_upper, na.rm = T),
        sd = sd(residual_upper, na.rm = T)
      )
    )

}
