#!/usr/bin/env Rscript

# Figure 3: Spatial Confounding in δ²Hwax-δ²Hprecip Relationship
# Shows how spatial structure affects the apparent slope

library(posterior)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(viridis)
library(scales)

source("posterior_helpers.R")

cat("\n=== Figure 3: Spatial Confounding in δ²Hwax-δ²Hprecip Relationship ===\n\n")

# ========================================
# LOAD MODELS AND DATA
# ========================================

cat("Loading models and data via posterior_helpers...\n")

draws_baseline    <- load_draws("baseline")
draws_sp          <- load_draws("baseline_sp")
sediment_data     <- load_sediment()
stan_data_sp      <- load_stan_data("baseline_sp")
d2H_sd            <- stan_data_sp$scaling_params$d2H_sd
cat(sprintf("  baseline: %d vars  baseline_sp: %d vars\n",
            length(variables(draws_baseline)), length(variables(draws_sp))))
cat(sprintf("  d2H_sd (standardization factor): %.3f ‰\n", d2H_sd))
cat("  sediment data: ", nrow(sediment_data), " observations\n")

# ========================================
# EXTRACT MODEL PARAMETERS
# ========================================

cat("\nExtracting model parameters...\n")

# Non-spatial model
beta_oipc_baseline <- as.numeric(
  as_draws_matrix(subset_draws(draws_baseline, variable = "beta_oipc")))
slope_mean_baseline <- mean(beta_oipc_baseline)
slope_ci_baseline   <- quantile(beta_oipc_baseline, c(0.025, 0.975))
cat(sprintf("  Non-spatial model slope: %.3f (95%% CI: %.3f-%.3f)\n",
            slope_mean_baseline, slope_ci_baseline[1], slope_ci_baseline[2]))

# Spatial model
beta_oipc_sp <- as.numeric(
  as_draws_matrix(subset_draws(draws_sp, variable = "beta_oipc")))
slope_mean_sp <- mean(beta_oipc_sp)
slope_ci_sp   <- quantile(beta_oipc_sp, c(0.025, 0.975))
cat(sprintf("  Spatial model slope: %.3f (95%% CI: %.3f-%.3f)\n",
            slope_mean_sp, slope_ci_sp[1], slope_ci_sp[2]))

# ========================================
# EXTRACT SPATIAL INTERCEPT EFFECTS (back-transformed to ‰)
# ========================================
#
# Stan defines `alpha_spatial[i]` on the *standardized* response scale
# (4d_leaf_wax_spatial_model.stan:316, 331-335). The quantity already
# includes the global intercept beta_0 added to a mean-zero GP residual.
# The GP residual contribution at observation i, in standardized units, is:
#
#     alpha_spatial[i] - beta_0
#
# To plot or subtract this in original per-mille units, multiply by the
# d2H_wax standardization factor saved in stan_data$scaling_params$d2H_sd
# (same factor Stan uses at line 488 to build `sigma_intercept_spatial`
# in original units). The prior version of this script multiplied
# `alpha_spatial` (still including beta_0) by `sigma_intercept_spatial`
# — a units error that mixed a standardized-scale quantity with an
# original-scale SD. See manuscript/FIGURES.md for the correction note.

cat("\nExtracting spatial intercept effects (back-transformed to ‰)...\n")

n_obs <- nrow(sediment_data)

# Posterior means of alpha_spatial[i] and beta_0 on the standardized scale.
# Pull both as draws_matrices (iters × chains rows flattened, obs as cols),
# take posterior means per column.
alpha_matrix <- as_draws_matrix(subset_draws(draws_sp, variable = "alpha_spatial"))
beta_0_draws <- as.numeric(
  as_draws_matrix(subset_draws(draws_sp, variable = "beta_0")))

beta_0_mean       <- mean(beta_0_draws)
alpha_spatial_std <- colMeans(alpha_matrix)
cat(sprintf("  beta_0 (standardized): %.4f\n", beta_0_mean))
stopifnot(length(alpha_spatial_std) == n_obs)

# GP residual contribution, in original ‰ units:
#   (alpha_spatial - beta_0) × d2H_sd
alpha_spatial <- (alpha_spatial_std - beta_0_mean) * d2H_sd

cat(sprintf("  Spatial intercept contribution range: %.1f to %.1f ‰\n",
            min(alpha_spatial), max(alpha_spatial)))
cat(sprintf("  Spatial intercept contribution mean:  %.3f ‰\n", mean(alpha_spatial)))
cat(sprintf("  Spatial intercept contribution SD:    %.3f ‰\n", sd(alpha_spatial)))

# ========================================
# ASSIGN CONTINENTS
# ========================================

cat("\nAssigning continents based on coordinates...\n")

assign_continent <- function(lon, lat) {
  # Simple continent assignment based on coordinates
  ifelse(lon >= -180 & lon <= -30 & lat >= -60 & lat <= 15, "South America",
  ifelse(lon >= -180 & lon <= -30 & lat > 15 & lat <= 90, "North America",
  ifelse(lon >= -30 & lon <= 50 & lat >= -40 & lat <= 40, "Africa",
  ifelse(lon >= -30 & lon <= 50 & lat > 40 & lat <= 90, "Europe",
  ifelse(lon >= 50 & lon <= 180 & lat >= -60 & lat <= 90, "Asia",
  ifelse(lon >= 100 & lon <= 180 & lat >= -60 & lat <= -10, "Australia",
         "Other"))))))
}

sediment_data$continent <- assign_continent(sediment_data$longitude, sediment_data$latitude)

# Continent colors
continent_colors <- c(
  "North America" = "#1f77b4",
  "South America" = "#ff7f0e",
  "Africa" = "#2ca02c",
  "Europe" = "#d62728",
  "Asia" = "#9467bd",
  "Australia" = "#8c564b",
  "Other" = "#7f7f7f"
)

# ========================================
# PREPARE DATA FOR PLOTTING
# ========================================

cat("\nPreparing data for plotting...\n")

# Create data frame for plotting
# IMPORTANT: Apply spatial adjustment correctly
plot_data <- data.frame(
  d2H_precip = sediment_data$d2H_precip,
  d2H_wax = sediment_data$d2H_wax,
  d2H_wax_adjusted = sediment_data$d2H_wax - alpha_spatial,  # Remove spatial effects
  longitude = sediment_data$longitude,
  latitude = sediment_data$latitude,
  continent = sediment_data$continent
)

# Remove rows with NA in d2H_precip
n_before <- nrow(plot_data)
plot_data <- plot_data[!is.na(plot_data$d2H_precip), ]
n_after <- nrow(plot_data)
cat(sprintf("  Removed %d rows with NA in d2H_precip, %d observations remaining\n",
            n_before - n_after, n_after))

# Check the adjustment worked
cat("\nVerifying spatial adjustment:\n")
cat(sprintf("  Original d2H_wax range: %.1f to %.1f ‰\n",
            min(plot_data$d2H_wax), max(plot_data$d2H_wax)))
cat(sprintf("  Adjusted d2H_wax range: %.1f to %.1f ‰\n",
            min(plot_data$d2H_wax_adjusted), max(plot_data$d2H_wax_adjusted)))
cat(sprintf("  Mean difference: %.3f ‰\n",
            mean(plot_data$d2H_wax - plot_data$d2H_wax_adjusted)))

# ========================================
# PANEL A: NON-SPATIAL MODEL
# ========================================

cat("\nCreating Panel A (non-spatial model)...\n")

panel_a <- ggplot(plot_data, aes(x = d2H_precip, y = d2H_wax)) +

  # Points colored by continent
  geom_point(aes(color = continent),
             size = 2,
             alpha = 0.6) +

  # Regression line with confidence interval
  geom_smooth(method = "lm",
              formula = y ~ x,
              color = "black",
              linewidth = 1,
              alpha = 0.2,
              se = TRUE) +

  # Color scale
  scale_color_manual(values = continent_colors,
                     name = "Continent") +

  # Labels (no title or subtitle as requested)
  labs(
    x = expression(delta^2*H[precip]~"(‰, VSMOW)"),
    y = expression(delta^2*H[wax]~"(‰, VSMOW)")
  ) +

  # Panel label
  annotate("text",
           x = -180,
           y = -80,
           label = "A",
           hjust = 0,
           vjust = 1,
           size = 6,
           fontface = "bold") +

  # Theme
  theme_classic(base_size = 11) +
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  ) +

  # Axis limits
  coord_cartesian(xlim = c(-190, 10), ylim = c(-290, -70))

# ========================================
# PANEL B: SPATIAL MODEL
# ========================================

cat("\nCreating Panel B (spatial model)...\n")

panel_b <- ggplot(plot_data, aes(x = d2H_precip, y = d2H_wax_adjusted)) +

  # Points colored by continent
  geom_point(aes(color = continent),
             size = 2,
             alpha = 0.6) +

  # Regression line with confidence interval
  geom_smooth(method = "lm",
              formula = y ~ x,
              color = "black",
              linewidth = 1,
              alpha = 0.2,
              se = TRUE) +

  # Color scale
  scale_color_manual(values = continent_colors,
                     name = "Continent") +

  # Labels (no title or subtitle as requested)
  labs(
    x = expression(delta^2*H[precip]~"(‰, VSMOW)"),
    y = expression(delta^2*H[wax]~"– spatial effect (‰, VSMOW)")
  ) +

  # Panel label
  annotate("text",
           x = -180,
           y = -80,
           label = "B",
           hjust = 0,
           vjust = 1,
           size = 6,
           fontface = "bold") +

  # Theme
  theme_classic(base_size = 11) +
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  ) +

  # Axis limits (same as Panel A for comparison)
  coord_cartesian(xlim = c(-190, 10), ylim = c(-290, -70))

# ========================================
# COMBINE PANELS
# ========================================

cat("\nCombining panels...\n")

# Extract legend from panel B
legend <- get_legend(panel_b)

# Remove legend from panel B for combining
panel_b_no_legend <- panel_b + theme(legend.position = "none")

# Combine panels side by side
combined_plot <- plot_grid(
  panel_a,
  panel_b_no_legend,
  ncol = 2,
  align = "hv",
  rel_widths = c(1, 1),
  labels = NULL  # We already have A and B labels in the panels
)

# Add legend at the bottom
final_plot <- plot_grid(
  combined_plot,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.15)
)

# ========================================
# SAVE FIGURE
# ========================================

cat("\nSaving figure...\n")

# Save as PDF (publication quality)
output_pdf <- "../figures/main_figs/figure_03_spatial_confounding.pdf"
ggsave(output_pdf, final_plot, width = 12, height = 6, dpi = 300)
cat(sprintf("  Saved PDF: %s\n", output_pdf))

# Save as PNG (for quick viewing)
output_png <- "../figures/main_figs/figure_03_spatial_confounding.png"
ggsave(output_png, final_plot, width = 12, height = 6, dpi = 300)
cat(sprintf("  Saved PNG: %s\n", output_png))

# ========================================
# SUMMARY STATISTICS
# ========================================

cat("\n=== Summary Statistics ===\n")

# Calculate reduction in slope
slope_reduction <- (1 - slope_mean_sp/slope_mean_baseline) * 100
cat(sprintf("\nSlope reduction: %.1f%%\n", slope_reduction))

# Spatial effect statistics
cat(sprintf("\nSpatial intercept effects:\n"))
cat(sprintf("  Mean: %.3f ‰\n", mean(alpha_spatial)))
cat(sprintf("  SD: %.3f ‰\n", sd(alpha_spatial)))
cat(sprintf("  Range: %.3f to %.3f ‰\n", min(alpha_spatial), max(alpha_spatial)))

# Compare R-squared values
lm_raw <- lm(d2H_wax ~ d2H_precip, data = plot_data)
lm_adjusted <- lm(d2H_wax_adjusted ~ d2H_precip, data = plot_data)

cat(sprintf("\nR-squared values:\n"))
cat(sprintf("  Non-spatial model: R² = %.3f\n", summary(lm_raw)$r.squared))
cat(sprintf("  Spatial model: R² = %.3f\n", summary(lm_adjusted)$r.squared))

# Check variance reduction
var_raw <- var(residuals(lm_raw))
var_adjusted <- var(residuals(lm_adjusted))
cat(sprintf("\nResidual variance:\n"))
cat(sprintf("  Non-spatial model: %.1f\n", var_raw))
cat(sprintf("  Spatial model: %.1f\n", var_adjusted))
cat(sprintf("  Variance reduction: %.1f%%\n", (1 - var_adjusted/var_raw) * 100))

# Visual check - scatter by continent
cat("\nMean d2H_wax by continent (for visual verification):\n")
continent_means_raw <- aggregate(d2H_wax ~ continent, data = plot_data, mean)
continent_means_adj <- aggregate(d2H_wax_adjusted ~ continent, data = plot_data, mean)
comparison <- merge(continent_means_raw, continent_means_adj, by = "continent")
names(comparison) <- c("Continent", "Raw", "Adjusted")
comparison$Difference <- comparison$Raw - comparison$Adjusted
print(comparison)

cat("\n=== Figure 3 generation complete ===\n")