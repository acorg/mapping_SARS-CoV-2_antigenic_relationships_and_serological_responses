
rm(list = ls())
library(Racmacs)
library(tidyverse)

source("functions/map_longinfo.R")
source("functions/diagnostics.R")
source("functions/imputation.R")

# Read map and imputed residual data
map <- read.acmap("data/maps/map_ndsubset_no_outliers_slope_adjusted.ace")
residual_data <- readRDS('data/generated_data/imputed_titers.rds')

# Calculate summary parameters for the fit
residual_data %>% 
  ungroup() %>%
  filter(
    measured_titer != "*"
  ) %>%
  summarise(
    n = sum(measured_titer != "*"),
    fit = list(fit_censored_normal(
      lower_lims = residual_lower,
      upper_lims = residual_upper,
      mean = mean(residual_upper, na.rm = T),
      sd = sd(residual_upper, na.rm = T)
    )),
    .groups = "keep"
  ) %>% 
  mutate(
    mean = vapply(fit, \(x) x$estimate["mean"], numeric(1)),
    sd   = vapply(fit, \(x) x$estimate["sd"], numeric(1))
  ) -> residual_summary

residual_data %>% 
  ungroup() %>%
  filter(
    measured_titer != "*",
    grepl("<", measured_titer)
  ) %>%
  summarise(
    n = sum(measured_titer != "*"),
    mean = mean(residual_error_imputed),
    sd = sd(residual_error_imputed),
    .groups = "keep"
  ) -> residual_summary_nd

residual_data %>% 
  ungroup() %>%
  filter(
    measured_titer != "*",
    !grepl("<", measured_titer)
  ) %>%
  summarise(
    n = sum(measured_titer != "*"),
    mean = mean(measured_logtiter_upper - predicted_logtiter),
    sd = sd(measured_logtiter_upper - predicted_logtiter),
    .groups = "keep"
  ) -> residual_summary_detectable

# Plot a histogram with the overall distribution of the residuals
residual_data %>% 
  mutate(
    titer_type = forcats::fct_recode(
      titer_type,
      "measurable" = "=",
      "non-detectable" = "<"
    )
  ) %>%
  ggplot() + 
  geom_histogram(
    aes(
      x = residual_error_imputed,
      fill = titer_type
    ),
    binwidth = 0.2,
    color = "grey50",
    alpha = 0.8
  ) +
  geom_vline(
    xintercept = residual_summary$mean,
    lty = 2,
    lwd = 0.6
  ) +
  geom_function(
    n = 1000,
    fun = \(x) {
      dnorm(x, sd = 0.87)*(sum(!is.na(residual_data$residual_error_imputed))*0.2)
    },
    color = "grey20",
    linetype = "11"
  ) +
  coord_cartesian(
    xlim = c(-8, 8)
  ) +
  theme_bw() +
  labs(
    x = "Measured log titer - Fitted log titer",
    y = "Count",
    fill = "Titer type",
    subtitle = sprintf(
      "µ = %s, σ = %s",
      round(residual_summary$mean, 2),
      round(residual_summary$sd, 2)
    )
  ) -> gp

ggsave("figures/som/figS12-measured_vs_fitted_titers_histogram.png", width = 7, height = 4)
ggsave("figures/som/figS12-measured_vs_fitted_titers_histogram.pdf", width = 7, height = 4)


