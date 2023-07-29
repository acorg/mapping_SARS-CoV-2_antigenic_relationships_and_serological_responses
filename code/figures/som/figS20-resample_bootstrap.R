
rm(list = ls())
library(Racmacs)
library(patchwork)

# bootmap_resample_ags_srs <- readRDS("data/generated_data/map_bootstrap_resample_ags_srs.rds")
# bootmap_resample_srs <- readRDS("data/generated_data/map_bootstrap_resample_srs.rds")
# bootmap_resample_ags <- readRDS("data/generated_data/map_bootstrap_resample_ags.rds")
bootmap_resample_ags_srs <- readRDS("data/generated_data/map_bootstrap_bayesian_ags_srs.rds")
bootmap_resample_srs <- readRDS("data/generated_data/map_bootstrap_bayesian_srs.rds")
bootmap_resample_ags <- readRDS("data/generated_data/map_bootstrap_bayesian_ags.rds")

gpA <- ggplot(bootmap_resample_ags_srs)
gpB <- ggplot(bootmap_resample_ags)
gpC <- ggplot(bootmap_resample_srs)

gp <- gpA + gpB + gpC + plot_annotation(tag_levels = "A") & theme(plot.tag.position = c(0.06, 0.94), plot.tag = element_text(size = 24))
ggsave(plot = gp, filename = "figures/som/figS20-resample_bootstrap_map.pdf", width = 12, height = 4, units = "in")
ggsave(plot = gp, filename = "figures/som/figS20-resample_bootstrap_map.png", width = 12, height = 4, units = "in")
