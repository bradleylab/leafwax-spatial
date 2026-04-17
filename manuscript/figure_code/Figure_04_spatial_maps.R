#!/usr/bin/env Rscript

library(posterior)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(tidyverse)
library(viridis)
library(scales)
library(gridExtra)
library(cowplot)

source("posterior_helpers.R")

cat("\n=== Figure 4: Spatial Patterns in Intercepts and Slopes ===\n\n")

# ========================================
# LOAD MODEL AND DATA
# ========================================

cat("Loading full_sp draws and stan_data via helpers...\n")
draws <- load_draws("full_sp")
stan_data <- load_stan_data("full_sp")

n_knots <- stan_data$n_pp_knots
cat(sprintf("  Number of knots: %d  Vars in draws: %d\n",
            n_knots, length(variables(draws))))

# Convenience helper: scalar posterior mean
.mean_of <- function(var) mean(as.numeric(
  as_draws_matrix(subset_draws(draws, variable = var))))

# ========================================
# EXTRACT PARAMETERS
# ========================================

# INTERCEPTS
ls_intercept_km         <- .mean_of("ls_intercept_km")
sigma_intercept_spatial <- .mean_of("sigma_intercept_spatial")
z_intercept_knots       <- colMeans(
  as_draws_matrix(subset_draws(draws, variable = "z_intercept_spatial")))
stopifnot(length(z_intercept_knots) == n_knots)

cat(sprintf("\nIntercept parameters:\n"))
cat(sprintf("  Length scale: %.1f km\n", ls_intercept_km))
cat(sprintf("  Spatial SD: %.1f ‰\n", sigma_intercept_spatial))

# SLOPES
ls_slope_km         <- .mean_of("ls_slope_km")
sigma_slope_spatial <- .mean_of("sigma_slope_spatial")
beta_oipc           <- .mean_of("beta_oipc")
z_slope_knots       <- colMeans(
  as_draws_matrix(subset_draws(draws, variable = "z_slope_spatial")))
stopifnot(length(z_slope_knots) == n_knots)

cat(sprintf("\nSlope parameters:\n"))
cat(sprintf("  Length scale: %.1f km\n", ls_slope_km))
cat(sprintf("  Global slope: %.3f\n", beta_oipc))
cat(sprintf("  Spatial SD: %.3f\n", sigma_slope_spatial))

# ========================================
# GET KNOT COORDINATES
# ========================================

cat("\nProcessing knot coordinates...\n")
knot_coords_std <- stan_data$knot_coords  # Standardized

# Get scaling parameters
lon_sd <- stan_data$coord_scaling[1]
lat_sd <- stan_data$coord_scaling[2]
lon_mean <- mean(stan_data$longitude)
lat_mean <- mean(stan_data$latitude)

# Transform knots to original scale
knot_coords <- matrix(NA, nrow=n_knots, ncol=2)
knot_coords[,1] <- knot_coords_std[,1] * lon_sd + lon_mean  # lon
knot_coords[,2] <- knot_coords_std[,2] * lat_sd + lat_mean  # lat
colnames(knot_coords) <- c("lon", "lat")

# ========================================
# MATÉRN 3/2 KERNEL
# ========================================

matern32_kernel <- function(x1, x2, alpha, rho) {
  d <- sqrt(sum((x1 - x2)^2))
  sqrt_3_d_rho <- sqrt(3) * d / rho

  if (d < 1e-10) {
    return(alpha^2)
  } else {
    return(alpha^2 * (1 + sqrt_3_d_rho) * exp(-sqrt_3_d_rho))
  }
}

# ========================================
# CREATE PREDICTION GRID
# ========================================

cat("Creating 2x2 degree prediction grid...\n")
lon_seq <- seq(-179, 179, by=2)
lat_seq <- seq(-89, 89, by=2)
pred_grid <- expand.grid(lon=lon_seq, lat=lat_seq)
n_pred <- nrow(pred_grid)

cat(sprintf("  Grid dimensions: %d x %d = %d points\n",
            length(lon_seq), length(lat_seq), n_pred))

# ========================================
# PANEL A: SPATIAL INTERCEPTS
# ========================================

cat("\nInterpolating spatial intercepts...\n")

# Convert length scale for intercepts
mean_lat_rad <- lat_mean * pi/180
km_per_deg_lon <- 111 * cos(mean_lat_rad)
km_per_deg_lat <- 111
ls_lon_std <- ls_intercept_km / (km_per_deg_lon * lon_sd)
ls_lat_std <- ls_intercept_km / (km_per_deg_lat * lat_sd)
ls_std_intercept <- sqrt(ls_lon_std^2 + ls_lat_std^2) / sqrt(2)

# Build kernel matrix for intercepts
K_knots_intercept <- matrix(NA, nrow=n_knots, ncol=n_knots)
for (i in 1:n_knots) {
  for (j in 1:n_knots) {
    K_knots_intercept[i, j] <- matern32_kernel(
      knot_coords_std[i,],
      knot_coords_std[j,],
      alpha = 1.0,
      rho = ls_std_intercept
    )
  }
}
K_knots_intercept <- K_knots_intercept + diag(1e-9, n_knots)
K_inv_intercept <- solve(K_knots_intercept)

# Interpolate intercepts
pred_grid$spatial_intercept <- NA
for (i in 1:n_pred) {
  pred_coords_std <- c((pred_grid$lon[i] - lon_mean) / lon_sd,
                       (pred_grid$lat[i] - lat_mean) / lat_sd)

  k_vec <- numeric(n_knots)
  for (j in 1:n_knots) {
    k_vec[j] <- matern32_kernel(
      pred_coords_std,
      knot_coords_std[j,],
      alpha = 1.0,
      rho = ls_std_intercept
    )
  }

  pred_z <- as.numeric(t(k_vec) %*% K_inv_intercept %*% z_intercept_knots)
  pred_grid$spatial_intercept[i] <- pred_z * sigma_intercept_spatial
}

cat(sprintf("  Intercept range: %.1f to %.1f ‰\n",
            min(pred_grid$spatial_intercept, na.rm=TRUE),
            max(pred_grid$spatial_intercept, na.rm=TRUE)))

# ========================================
# PANEL B: SPATIAL SLOPES
# ========================================

cat("\nInterpolating spatial slopes...\n")

# Convert length scale for slopes
ls_lon_std <- ls_slope_km / (km_per_deg_lon * lon_sd)
ls_lat_std <- ls_slope_km / (km_per_deg_lat * lat_sd)
ls_std_slope <- sqrt(ls_lon_std^2 + ls_lat_std^2) / sqrt(2)

# Build kernel matrix for slopes
K_knots_slope <- matrix(NA, nrow=n_knots, ncol=n_knots)
for (i in 1:n_knots) {
  for (j in 1:n_knots) {
    K_knots_slope[i, j] <- matern32_kernel(
      knot_coords_std[i,],
      knot_coords_std[j,],
      alpha = 1.0,
      rho = ls_std_slope
    )
  }
}
K_knots_slope <- K_knots_slope + diag(1e-9, n_knots)
K_inv_slope <- solve(K_knots_slope)

# Interpolate slopes
pred_grid$spatial_slope <- NA
for (i in 1:n_pred) {
  pred_coords_std <- c((pred_grid$lon[i] - lon_mean) / lon_sd,
                       (pred_grid$lat[i] - lat_mean) / lat_sd)

  k_vec <- numeric(n_knots)
  for (j in 1:n_knots) {
    k_vec[j] <- matern32_kernel(
      pred_coords_std,
      knot_coords_std[j,],
      alpha = 1.0,
      rho = ls_std_slope
    )
  }

  pred_z <- as.numeric(t(k_vec) %*% K_inv_slope %*% z_slope_knots)
  pred_grid$spatial_slope[i] <- beta_oipc + pred_z * sigma_slope_spatial
}

cat(sprintf("  Slope range: %.3f to %.3f\n",
            min(pred_grid$spatial_slope, na.rm=TRUE),
            max(pred_grid$spatial_slope, na.rm=TRUE)))

# ========================================
# CREATE PLOTS
# ========================================

cat("\nCreating combined figure...\n")

# Get world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# Sample locations
sample_coords <- data.frame(
  lon = stan_data$longitude,
  lat = stan_data$latitude
)

# PANEL A: Intercepts
p_intercept <- ggplot() +
  # Raster layer
  geom_raster(data = pred_grid,
              aes(x = lon, y = lat, fill = spatial_intercept),
              interpolate = TRUE) +

  # World boundaries
  geom_sf(data = world,
          fill = NA,
          color = "gray30",
          linewidth = 0.2,
          alpha = 0.8) +

  # Sample locations
  geom_point(data = sample_coords,
             aes(x = lon, y = lat),
             shape = 21,
             size = 0.5,
             fill = NA,
             color = "black",
             alpha = 0.3,
             stroke = 0.2) +

  # Knot locations
  geom_point(data = as.data.frame(knot_coords),
             aes(x = lon, y = lat),
             shape = 3,
             size = 1,
             color = "black",
             alpha = 0.6,
             stroke = 0.4) +

  # Color scale
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-70, 50),
    oob = scales::squish,
    name = "Spatial Effect (‰)",
    breaks = seq(-60, 60, by = 30),
    guide = guide_colorbar(barwidth = 10, barheight = 0.3)
  ) +

  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +

  theme_minimal(base_size = 10) +
  theme(
    panel.background = element_rect(fill = "#E8F4F8", color = NA),
    panel.grid = element_line(color = "gray90", linewidth = 0.15),
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9),
    plot.margin = margin(5, 5, 5, 5)
  ) +

  labs(
    title = "A. Spatial Intercept Variations",
    subtitle = sprintf("Range: %.1f to %.1f ‰ | SD: %.1f ‰",
                      min(pred_grid$spatial_intercept, na.rm=TRUE),
                      max(pred_grid$spatial_intercept, na.rm=TRUE),
                      sd(pred_grid$spatial_intercept, na.rm=TRUE))
  )

# PANEL B: Slopes
# Determine appropriate color scale for slopes
slope_min <- min(pred_grid$spatial_slope, na.rm=TRUE)
slope_max <- max(pred_grid$spatial_slope, na.rm=TRUE)
slope_range <- slope_max - slope_min
color_min_slope <- floor((slope_min - 0.05 * slope_range) * 20) / 20
color_max_slope <- ceiling((slope_max + 0.05 * slope_range) * 20) / 20

p_slope <- ggplot() +
  # Raster layer
  geom_raster(data = pred_grid,
              aes(x = lon, y = lat, fill = spatial_slope),
              interpolate = TRUE) +

  # World boundaries
  geom_sf(data = world,
          fill = NA,
          color = "gray30",
          linewidth = 0.2,
          alpha = 0.8) +

  # Sample locations
  geom_point(data = sample_coords,
             aes(x = lon, y = lat),
             shape = 21,
             size = 0.5,
             fill = NA,
             color = "black",
             alpha = 0.3,
             stroke = 0.2) +

  # Knot locations
  geom_point(data = as.data.frame(knot_coords),
             aes(x = lon, y = lat),
             shape = 3,
             size = 1,
             color = "black",
             alpha = 0.6,
             stroke = 0.4) +

  # Color scale
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = beta_oipc,
    limits = c(color_min_slope, color_max_slope),
    oob = scales::squish,
    name = "Slope Coefficient",
    breaks = seq(0.44, 0.46, by = 0.01),
    labels = sprintf("%.2f", seq(0.44, 0.46, by = 0.01)),
    guide = guide_colorbar(barwidth = 10, barheight = 0.3)
  ) +

  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +

  theme_minimal(base_size = 10) +
  theme(
    panel.background = element_rect(fill = "#E8F4F8", color = NA),
    panel.grid = element_line(color = "gray90", linewidth = 0.15),
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9),
    plot.margin = margin(5, 5, 5, 5)
  ) +

  labs(
    title = "B. Spatial Slope Variations",
    subtitle = sprintf("Range: %.3f to %.3f | Global: %.3f",
                      min(pred_grid$spatial_slope, na.rm=TRUE),
                      max(pred_grid$spatial_slope, na.rm=TRUE),
                      beta_oipc)
  )

# ========================================
# COMBINE AND SAVE
# ========================================

# Combine panels
combined_plot <- plot_grid(
  p_intercept, p_slope,
  ncol = 1,
  nrow = 2,
  labels = NULL,
  rel_heights = c(1, 1),
  align = "v"
)

# Add overall title
title_grob <- ggdraw() +
  draw_label(
    "Figure 3: Spatial Patterns in Model Parameters (Full_sp Model)",
    fontface = "bold",
    size = 12,
    x = 0.5,
    hjust = 0.5
  )

subtitle_grob <- ggdraw() +
  draw_label(
    sprintf("Matérn 3/2 GP with %d knots | Length scale: %.0f km | 313:1 ratio (intercept:slope variation)",
            n_knots, ls_intercept_km),
    size = 10,
    x = 0.5,
    hjust = 0.5,
    color = "gray30"
  )

# Final plot with title
final_plot <- plot_grid(
  title_grob,
  subtitle_grob,
  combined_plot,
  ncol = 1,
  rel_heights = c(0.05, 0.03, 1)
)

# Save as PDF
output_pdf <- "Figure_03.pdf"
cat(sprintf("\nSaving figure to %s...\n", output_pdf))
ggsave(output_pdf, final_plot, width = 10, height = 12, dpi = 300)

# Save as PNG for quick viewing
output_png <- "Figure_03.png"
ggsave(output_png, final_plot, width = 10, height = 12, dpi = 150)
cat(sprintf("Also saved as %s\n", output_png))

# ========================================
# REGIONAL SUMMARY
# ========================================

cat("\n=== Regional Summary ===\n")

regions <- list(
  "Western Americas" = list(lon=c(-130, -70), lat=c(-20, 20)),
  "Europe" = list(lon=c(-10, 30), lat=c(35, 60)),
  "Africa" = list(lon=c(-20, 50), lat=c(-35, 35)),
  "Asia" = list(lon=c(60, 140), lat=c(0, 50)),
  "Australia" = list(lon=c(110, 155), lat=c(-40, -10))
)

cat("\nIntercepts by region:\n")
for (region_name in names(regions)) {
  region <- regions[[region_name]]
  mask <- pred_grid$lon >= region$lon[1] & pred_grid$lon <= region$lon[2] &
          pred_grid$lat >= region$lat[1] & pred_grid$lat <= region$lat[2]

  if (sum(mask) > 0) {
    regional_intercepts <- pred_grid$spatial_intercept[mask]
    cat(sprintf("  %-20s: Mean = %6.1f ‰, Range = %6.1f to %6.1f ‰\n",
                region_name,
                mean(regional_intercepts, na.rm=TRUE),
                min(regional_intercepts, na.rm=TRUE),
                max(regional_intercepts, na.rm=TRUE)))
  }
}

cat("\nSlopes by region:\n")
for (region_name in names(regions)) {
  region <- regions[[region_name]]
  mask <- pred_grid$lon >= region$lon[1] & pred_grid$lon <= region$lon[2] &
          pred_grid$lat >= region$lat[1] & pred_grid$lat <= region$lat[2]

  if (sum(mask) > 0) {
    regional_slopes <- pred_grid$spatial_slope[mask]
    cat(sprintf("  %-20s: Mean = %.4f, Range = %.4f to %.4f\n",
                region_name,
                mean(regional_slopes, na.rm=TRUE),
                min(regional_slopes, na.rm=TRUE),
                max(regional_slopes, na.rm=TRUE)))
  }
}

cat("\n=== Figure 3 generation complete ===\n")