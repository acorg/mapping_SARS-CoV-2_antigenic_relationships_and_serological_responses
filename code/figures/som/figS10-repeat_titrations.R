
# Setup workspace
rm(list = ls())
library(tidyverse)
set.seed(10)

repeat_data <- readxl::read_excel("data/repeat_titration_data.xlsx")

repeat_data %>%
  distinct() %>%
  group_by(
    Pub.ID, Virus
  ) %>% 
  group_split() -> repeat_groups

repeat_groups <- repeat_groups[sapply(repeat_groups, nrow) > 1]

# Calculate variation between pairwise repeats
pairwise_variation <- c()
for (repeat_group in repeat_groups) {
  for (row1 in 1:(nrow(repeat_group) - 1)) {
    for (row2 in (row1 + 1):nrow(repeat_group)) {
      titer1 <- repeat_group$ID50[row1]
      titer2 <- repeat_group$ID50[row2]
      logtiter1 <- log2(as.numeric(titer1) / 10)
      logtiter2 <- log2(as.numeric(titer2) / 10)
      pairwise_variation <- rbind(
        pairwise_variation, 
        c(logtiter1, logtiter2)
      )
    }
  }
}

repeat_diff <- pairwise_variation[,1] - pairwise_variation[,2]
repeat_diff <- repeat_diff[!is.na(repeat_diff)]
repeat_diff_mean_stat <- Hmisc::smean.cl.normal(repeat_diff)
repeat_diff_mean_stat <- round(repeat_diff_mean_stat, 3)
sd_mean0 <- sqrt(sum(repeat_diff^2) / (length(repeat_diff) - 1))


# Plot a histogram of differences
doplot <- function() {
  layout(matrix(1:2, nrow = 1))
  hist(
    repeat_diff, 
    breaks = seq(from = -3.75, to = 3.75, by = 0.5),
    xlab = "Repeat difference",
    main = ""
  )
  
  shap_result <- shapiro.test((repeat_diff - mean(repeat_diff, na.rm = T)) / sd(repeat_diff, na.rm = TRUE))
  qqplot(
    qnorm(ppoints(length(repeat_diff)), mean = mean(repeat_diff), sd = sd(repeat_diff)),
    repeat_diff,
    xlab = sprintf("Theoretical quantile (µ=%s, σ=%s)", round(mean(repeat_diff), 2), round(sd(repeat_diff), 2)),
    ylab = "Sample quantile"
  )
  abline(0, 1, lty = 2)
}


# Create plots
cairo_pdf("figures/som/figS10-repeat_variation.pdf", 9.8, 4.8)
doplot()
dev.off()

