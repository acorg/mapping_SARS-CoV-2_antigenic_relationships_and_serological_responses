
# Clear figure directories
unlink("figures", recursive = TRUE)
dir.create("figures")
dir.create("figures/main")
dir.create("figures/som")
dir.create("figures/summary")

# ==== MAIN FIGURES =======

# Generate all interactive landscapes
message("==== Rendering LANDSCAPES ========")
source("code/figures/main/fig-landscapes.R")

# Generate landscapes screenshots
message("==== Rendering LANDSCAPES SCREENSHOTS ========")
source("code/figures/main/fig-landscapes_screenshot.R")

# Generate interactive maps
message("==== Rendering MAPS ========")
source("code/figures/main/fig-maps.R")

# Figure 1
message("==== Rendering FIG. 1 ========")
source("code/figures/main/fig1-titerplots.R")

# Figure 2
message("==== Rendering FIG. 2 & FIG. S6 ========")
source("code/figures/main/fig2-folddrop_comp.R")

# Figure 3
message("==== Rendering FIG. 3 ========")
source("code/figures/main/fig3-map.R")

# Figure 4
message("==== Rendering FIG. 4 ========")
source("code/figures/main/fig4-landscapes.R")

# Figure 5
message("==== Rendering FIG. 5 & FIG. S27 ========")
source("code/figures/main/fig5-mutant_map.R")

# Figure 6
message("==== Rendering FIG. 6 & FIG. S28 ========")
source("code/figures/main/fig6-ntdeffects_and_immunodominance.R")


# ==== SOM FIGURES =======

# Figure S2
message("==== Rendering FIG. S2 ========")
source("code/figures/som/figS2-outlier_titerplot.R")

# Figure S3
message("==== Rendering FIG. S3 ========")
source("code/figures/som/figS3-titer_patterns_by_source.R")

# Figure S4
message("==== Rendering FIG. S4 ========")
source("code/figures/som/figS4-titer_boxplots.R")

# Figure S5
message("==== Rendering FIG. S5 ========")
source("code/figures/som/figS5-titerplot_no_individual_effects.R")

# Figure S6
message("==== Rendering FIG. S6 ========")
source("code/figures/som/figS6-folddrops.R")

# Figure S7
message("==== Rendering FIG. S7 ========")
source("code/figures/som/figS7-forestplot_comparison.R")

# Figure S8
message("==== FIG. S8 already rendered with FIG. 2 ========")

# Figure S9
message("==== Rendering FIG. S9 ========")
source("code/figures/som/figS9-dimensionality_test.R")

# Figure S10
message("==== Rendering FIG. S10 ========")
source("code/figures/som/figS10-repeat_titrations.R")

# Figure S11
message("==== Rendering FIG. S11 ========")
source("code/figures/som/figS11-fitted_vs_measured_scatter.R")

# Figure S12
message("==== Rendering FIG. S12 ========")
source("code/figures/som/figS12-measured_vs_fitted_titers_histogram.R")

# Figure S13
message("==== Rendering FIG. S13 ========")
source("code/figures/som/figS13-measured_vs_fitted_titers_box.R")

# Figure S14
message("==== Rendering FIG. S14 ========")
source("code/figures/som/figS14-smooth_bootstrap.R")

# Figure S15
message("==== Rendering FIG. S15 ========")
source("code/figures/som/figS15-excluding_variants.R")

# Figure S16
message("==== Rendering FIG. S16 ========")
source("code/figures/som/figS16-excluding_sr_groups.R")

# Figure S17
message("==== Rendering FIG. S17 ========")
source("code/figures/som/figS17-map_made_from_infection_only_sera.R")

# Figure S18
message("==== Rendering FIG. S18 ========")
source("code/figures/som/figS18-adjusting_B.1.617.2_homologous_titers.R")

# Figure S19
message("==== Rendering FIG. S19 ========")
source("code/figures/som/figS19-map_outlier_removal.R")

# Figure S20
message("==== Rendering FIG. S20 ========")
source("code/figures/som/figS20-resample_bootstrap.R")

# Figure S21
message("==== Rendering FIG. S21 ========")
source("code/figures/som/figS21-error_and_triangulation_maps.R")

# Figure S22
message("==== Rendering FIG. S22 ========")
source("code/figures/som/figS22-measured_vs_predicted_titers_histogram.R")

# Figure S23
message("==== Rendering FIG. S23 ========")
source("code/figures/som/figS23-measured_vs_predicted_titers_box.R")

# Figure S24
message("==== Rendering FIG. S24 ========")
source("code/figures/som/figS24-full_map.R")

# Figure S25, S27, S28, S29, S30
message("==== Rendering FIG. S25, S27, S28, S29, S30 ========")
source("code/figures/som/figS25_27-30-scatterplots.R")

# Figure S26
message("==== Rendering FIG. S26 ========")
source("code/figures/som/figS26-ablandscapes_fitted_vs_measured.R")

# Figure S31
message("==== FIG. S31 already rendered with FIG. 5 ========")

# Figure S32 & S33
message("==== FIG. S32 & S33 already rendered with FIG. 6 ========")

# Generate print summary figures
message("==== Rendering FIG. 0 ========")
source("code/figures/summary/fig0B-folddrop_comp.R")
source("code/figures/summary/fig0D-ntdeffects_and_immunodominance.R")

