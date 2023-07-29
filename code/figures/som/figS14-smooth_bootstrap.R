
rm(list = ls())
library(Racmacs)
library(patchwork)

bootmap_titer_ag_noise <- readRDS("data/generated_data/map_bootstrap_noisy_ags_titers.rds")
bootmap_titer_noise <- readRDS("data/generated_data/map_bootstrap_noisy_titers.rds")
bootmap_ag_noise <- readRDS("data/generated_data/map_bootstrap_noisy_ags.rds")

gpA <- ggplot(bootmap_titer_ag_noise)
gpB <- ggplot(bootmap_titer_noise)
gpC <- ggplot(bootmap_ag_noise)

gp <- gpA + gpB + gpC + plot_annotation(tag_levels = "A") & theme(plot.tag.position = c(0.06, 0.94), plot.tag = element_text(size = 24))
ggsave(plot = gp, filename = "figures/som/figS14-noisy_bootstrap_map.pdf", width = 12, height = 4, units = "in")
ggsave(plot = gp, filename = "figures/som/figS14-noisy_bootstrap_map.png", width = 12, height = 4, units = "in")
