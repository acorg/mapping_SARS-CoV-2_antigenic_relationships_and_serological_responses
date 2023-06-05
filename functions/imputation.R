
#' Function to fit a censored normal distribution
fit_censored_normal <- function(
  lower_lims,
  upper_lims,
  mean,
  sd
) {

  # Set NA result
  result <- list(
    estimate = c(mean = NA, sd = NA),
    sd = c(mean = NA, sd = NA),
    vcov = diag(2)
  )

  try({

    capture.output({
      fit <- fitdistrplus::fitdistcens(
        censdata = data.frame(
          left = lower_lims,
          right = upper_lims
        ),
        start = list(
          mean = mean,
          sd = sd
        ),
        distr = "norm"
      )
    })

    result$estimate <- fit$estimate
    result$sd       <- fit$sd
    result$vcov     <- fit$vcov

  }, silent = TRUE)

  # Return the result
  result

}

#' Function to impute censored data
#' @export
impute_censored_normal <- function(
  lower_lims,
  upper_lims,
  mu,
  sigma,
  covariance = diag(1, 2)
) {

  # Deal with NA mus
  if (is.na(mu)) {
    imputed_values <- lower_lims
    imputed_values[imputed_values == -Inf] <- upper_lims
    return(imputed_values)
  }

  # Infer censored data
  censored_data <- upper_lims != lower_lims & !is.na(upper_lims)

  # Simply return the same data if no censored data was found
  if (sum(censored_data) == 0) return(upper_lims)

  # Impute mean and sd variables
  imputed_pars <- mvtnorm::rmvnorm(sum(censored_data), mean = c(mu, sigma), sigma = covariance)

  # Impute the censored values
  imputed_values <- rep(NA, length(upper_lims))
  imputed_values[!censored_data] <- upper_lims[!censored_data]
  imputed_values[censored_data] <- truncnorm::rtruncnorm(
    n = sum(censored_data),
    a = lower_lims[censored_data],
    b = upper_lims[censored_data],
    mean = imputed_pars[,1],
    sd = imputed_pars[,2]
  )

  # Return the imputed values
  imputed_values

}

#' Function to impute a whole censored normal distribution including mean and sd
#' @export
fit_imputed_censored_normal <- function(
  lower_lims,
  upper_lims,
  mean,
  sd
) {

  # Deal with missing data cases
  if (sum(!is.na(lower_lims)) == 0) {
    return(lower_lims)
  }

  # Deal with cases where lower_lims are all -Inf
  if (length(unique(lower_lims)) == 1 & lower_lims[1] == -Inf) {
    return(Inf)
  }

  # Fit the distribution
  result <- fit_censored_normal(
    lower_lims = lower_lims,
    upper_lims = upper_lims,
    mean,
    sd
  )

  # Impute the censored data
  imputed_censored_data <- impute_censored_normal(
    lower_lims = lower_lims,
    upper_lims = upper_lims,
    mu = result$estimate["mean"],
    sigma = result$estimate["sd"],
    covariance = result$vcov
  )

  imputed_censored_data

}

