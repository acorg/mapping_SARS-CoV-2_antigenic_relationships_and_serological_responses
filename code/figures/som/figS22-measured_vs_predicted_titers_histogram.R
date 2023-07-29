
rm(list = ls())
library(Racmacs)
library(tidyverse)

source("functions/diagnostics.R")
source("functions/imputation.R")

map <- read.acmap('data/maps/map_ndsubset_no_outliers_slope_adjusted.ace')
cv_data <- readRDS("data/generated_data/cv_testing.rds")

dilution_stepsize <- 0

cv_data$measured_logtiter_lower <- cv_data$measured_logtiter - dilution_stepsize / 2
cv_data$measured_logtiter_upper <- cv_data$measured_logtiter + dilution_stepsize / 2

cv_data$residual_lower = cv_data$measured_logtiter_lower - cv_data$predicted_logtiter
cv_data$residual_upper = cv_data$measured_logtiter_upper - cv_data$predicted_logtiter

cv_data$sr_group <- paste(cv_data$sr_group, "sera", sep=" ")

# Impute data
cv_data %>%
  filter(
    measured_titer != "*"
  ) %>%
  group_by(
    ag_name,
    sr_group
  ) %>%
  group_modify(
    .f = impute_table_residuals
  ) -> cv_data


# Calculate summary parameters for the fit
cv_data %>%
  ungroup() %>%
  filter(
    measured_titer != "*",
    is.finite(residual_upper)
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
  ) -> cv_summary

fit_0 <- fit_censored_normal(
  lower_lims = cv_data$residual_lower[is.finite(cv_data$residual_lower)],
  upper_lims = cv_data$residual_upper[is.finite(cv_data$residual_upper)],
  mean = 0,
  sd = sd(cv_data$residual_error_imputed[is.finite(cv_data$residual_error_imputed)], na.rm = T)
)

sprintf('Standard deviation with mean of 0: %s', fit_0$estimate["sd"])
sprintf('Standard deviation for error in both titers: %s', sqrt(fit_0$estimate["sd"] ^ 2 / 2))


# Plot a histogram with the overall distribution of the residuals
cv_data %>%
  ggplot() +
  geom_histogram(
    aes(
      x = residual_error_imputed,
      fill = factor(measured_titer_type, levels = c(1, 2), labels = c("=", "<"))
    ),
    binwidth = 0.2,
    color = "grey50",
    alpha = 0.8
  ) +
  geom_vline(
    xintercept = cv_summary$mean,
    lty = 2,
    lwd = 0.6
  ) +
  geom_function(
    n = 1000,
    fun = \(x) {
      dnorm(x, sd = 0.87)*(sum(!is.na(cv_data$residual_error_imputed))*0.2)
    },
    color = "grey20",
    linetype = "11"
  ) +
  coord_cartesian(
    xlim = c(-10, 10)
  ) +
  scale_fill_manual(
    values = c(
      "<" = "lightblue",
      "=" = "grey80"
    )
  ) +
  theme_bw() +
  labs(
    x = "Measured log titer - Predicted log titer",
    y = "Count",
    fill = "Titer type",
    subtitle = sprintf(
      "µ = %s, σ = %s",
      round(cv_summary$mean, 2),
      round(cv_summary$sd, 2)
    )
  ) -> gp

ggsave("figures/som/figS22-measured_vs_predicted_titers_histogram.png", width=5, height=4)
ggsave("figures/som/figS22-measured_vs_predicted_titers_histogram.pdf", width=5, height=4)

