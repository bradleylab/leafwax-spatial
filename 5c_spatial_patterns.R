#!/usr/bin/env Rscript
#───────────────────────────────────────────────────────────────────────────────
# 5c_spatial_patterns.R
#
# Analyze and visualize spatial patterns in model outputs
# Maps spatial random effects, identifies spatial clusters, and tests residuals
# Creates diagnostic plots for spatial variation in parameters
#
# Input: results/model_fits/*.rds, results/3_sediment_ready_for_modeling.rds
# Output: results/spatial_patterns/ (maps, variograms, diagnostic plots)
#───────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(cmdstanr)
library(posterior)
library(viridis)
library(sf)
library(rnaturalearth)
library(geoR)
library(fields)
library(ape)
library(cowplot)
library(patchwork)

cat("\nSPATIAL PATTERN DIAGNOSTICS\n")
cat("===========================\n")

# Setup
output_dir <- "model_analysis/spatial_pattern_diagnostics"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Get all fitted models
model_dirs <- list.dirs("model_output", recursive = FALSE)
model_names <- basename(model_dirs)

cat("Found", length(model_names), "fitted models\n\n")

# Initialize storage for results
all_results <- list()
all_spatial_summaries <- list()
all_regional_summaries <- list()
all_knot_data <- list()
all_predictions <- list()
all_variograms <- list()

# Process each model
for (model_name in model_names) {
  
  cat("\n\nAnalyzing:", model_name, "\n")
  cat(strrep("-", 60), "\n")
  
  tryCatch({
    
    # Skip non-spatial models
    if (!grepl("_sp$", model_name)) {
      cat("  Skipping non-spatial model\n")
      next
    }
    
    # Load model and data
    fit_file <- file.path("model_output", model_name, "fit.rds")
    stan_data_file <- file.path("prepared_data", 
                                paste0("stan_data_", model_name, ".rds"))
    
    if (!file.exists(fit_file) || !file.exists(stan_data_file)) {
      cat("  Skipping - files not found\n")
      next
    }
    
    fit <- readRDS(fit_file)
    stan_data <- readRDS(stan_data_file)
    
    # Get all variables
    all_vars <- fit$metadata()$variables
    unique_vars <- unique(gsub("\\[.*\\]", "", all_vars))
    
    cat("  Model has", length(all_vars), "variables\n")
    
    # Identify spatial variables
    spatial_vars <- c(
      grep("spatial", unique_vars, value = TRUE),
      grep("sigma_", unique_vars, value = TRUE),
      grep("ls_", unique_vars, value = TRUE),
      grep("lambda", unique_vars, value = TRUE)
    )
    
    cat("  Spatial-related variables found:", paste(head(spatial_vars, 10), collapse = ", "), "\n")
    
    # Extract dimensions
    N <- stan_data$N
    n_knots <- stan_data$n_pp_knots
    
    # Check if this is truly a spatial model
    has_spatial <- length(spatial_vars) > 2 && n_knots > 0
    
    if (!has_spatial) {
      cat("  No spatial parameters found - skipping detailed analysis\n")
      next
    }
    
    # Identify key variables - CORRECTED NAMES
    alpha_vars <- grep("^z_intercept_spatial", all_vars, value = TRUE)
    slope_vars <- grep("^z_slope_spatial", all_vars, value = TRUE)
    
    if (length(alpha_vars) > 0) {
      alpha_var <- alpha_vars[1]
      cat("  Found intercept variable:", alpha_var, "\n")
    } else {
      cat("  No spatial intercept variable found\n")
    }
    
    if (length(slope_vars) > 0) {
      slope_var <- slope_vars[1]
      cat("  Found slope variable:", slope_var, "\n")
    } else {
      cat("  No spatial slope variable found\n")
    }
    
    # 1. SPATIAL SLOPE ANALYSIS
    cat("\n1. CHECKING FOR SPATIAL PARAMETERS\n")
    
    has_spatial_slopes <- length(slope_vars) > 0
    
    if (has_spatial_slopes) {
      # Extract raw spatial slopes at knots
      slope_summary <- fit$summary(variables = slope_vars)
      z_slope_means <- slope_summary$mean[1:n_knots]
      
      # Get sigma_slope - CORRECTED NAME
      sigma_slope_summary <- fit$summary(variables = "sigma_slope_spatial")
      sigma_slope <- sigma_slope_summary$mean
      
      # Get beta_oipc (global slope)
      beta_oipc_summary <- fit$summary(variables = "beta_oipc")
      beta_oipc <- beta_oipc_summary$mean
      
      # Calculate actual slopes at knots
      knot_slopes <- beta_oipc + z_slope_means * sigma_slope
      
      cat("  Found spatial slope summaries\n")
      cat("  Slope range: [", round(min(knot_slopes), 3), ",", 
          round(max(knot_slopes), 3), "]\n")
      cat("  Slope SD:", round(sd(knot_slopes), 3), "\n")
    } else {
      cat("  No spatial slopes found\n")
      # Use global slope for all knots
      beta_oipc_summary <- fit$summary(variables = "beta_oipc")
      beta_oipc <- beta_oipc_summary$mean
      knot_slopes <- rep(beta_oipc, n_knots)
    }
    
    # 1. SPATIAL INTERCEPT ANALYSIS
    cat("\n1. SPATIAL INTERCEPT ANALYSIS\n")
    
    has_spatial_intercepts <- length(alpha_vars) > 0
    
    if (has_spatial_intercepts) {
      # Extract raw spatial intercepts at knots
      intercept_summary <- fit$summary(variables = alpha_vars)
      z_intercept_means <- intercept_summary$mean[1:n_knots]
      
      # Get sigma_intercept - CORRECTED NAME
      sigma_intercept_summary <- fit$summary(variables = "sigma_intercept_spatial")
      sigma_intercept <- sigma_intercept_summary$mean
      
      # Get beta_0 (global intercept)
      beta_0_summary <- fit$summary(variables = "beta_0")
      beta_0 <- beta_0_summary$mean
      
      # Get scaling parameters from stan_data
      d2h_mean <- stan_data$scaling_params$d2H_mean
      d2h_sd <- stan_data$scaling_params$d2H_sd
      oipc_mean <- stan_data$scaling_params$oipc_mean
      
      # Calculate actual intercepts in standardized space
      intercept_std <- beta_0 + z_intercept_means * (sigma_intercept / d2h_sd)
      
      # Transform to original scale
      intercept_at_mean <- intercept_std * d2h_sd + d2h_mean
      
      # Adjust to intercept at OIPC = 0
      intercept_at_zero <- intercept_at_mean + knot_slopes * (0 - oipc_mean)
      
      cat("  Global intercept at OIPC=0:", round(mean(intercept_at_zero), 1), "‰\n")
    } else {
      cat("  No spatial intercepts found\n")
      cat("  Skipping detailed intercept analysis\n")
      
      # Use global intercept
      beta_0_summary <- fit$summary(variables = "beta_0")
      beta_0 <- beta_0_summary$mean
      d2h_mean <- stan_data$scaling_params$d2H_mean
      d2h_sd <- stan_data$scaling_params$d2H_sd
      oipc_mean <- stan_data$scaling_params$oipc_mean
      
      intercept_at_mean <- beta_0 * d2h_sd + d2h_mean
      intercept_at_zero <- rep(intercept_at_mean + beta_oipc * (0 - oipc_mean), n_knots)
    }
    
    # Extract length scales
    ls_vars <- grep("ls_.*_km", unique_vars, value = TRUE)
    if (length(ls_vars) > 0) {
      ls_summary <- fit$summary(variables = ls_vars[1])
      ls_intercept <- ls_slope <- ls_summary$mean
    } else {
      ls_intercept <- ls_slope <- 0
    }
    
    cat("  Length scales: intercept =", round(ls_intercept, 0), 
        "km, slope =", round(ls_slope, 0), "km\n")
    
    # Sigma values
    if (exists("sigma_intercept")) {
      cat("  Spatial SDs: intercept =", round(sigma_intercept, 1), 
          "‰, slope =", round(sigma_slope, 3), "\n")
    }
    
    # 2. MAPPING SPATIAL PATTERNS
    cat("\n2. MAPPING SPATIAL PATTERNS\n")
    
    # Get knot coordinates
    knot_coords <- stan_data$knot_coords
    
    # Transform back to original scale
    lon_mean <- mean(stan_data$longitude)
    lon_sd <- sd(stan_data$longitude)
    lat_mean <- mean(stan_data$latitude)
    lat_sd <- sd(stan_data$latitude)
    
    knot_lon <- knot_coords[,1] * lon_sd + lon_mean
    knot_lat <- knot_coords[,2] * lat_sd + lat_mean
    
    # Create knot data frame
    knot_df <- data.frame(
      model = model_name,
      knot_id = 1:n_knots,
      lon = knot_lon,
      lat = knot_lat,
      slope = knot_slopes,
      intercept = intercept_at_zero,
      density = stan_data$knot_data_density
    )
    
    all_knot_data[[model_name]] <- knot_df
    
    # Get world map
    world <- ne_countries(scale = "medium", returnclass = "sf")
    
    # Create slope map
    p_slope <- ggplot() +
      geom_sf(data = world, fill = "gray90", color = "white") +
      geom_point(data = knot_df, 
                 aes(x = lon, y = lat, color = slope, size = sqrt(density + 1)),
                 alpha = 0.8) +
      scale_color_viridis(name = "Slope", 
                         limits = range(knot_slopes) + c(-0.05, 0.05) * diff(range(knot_slopes))) +
      scale_size_continuous(range = c(2, 8), guide = "none") +
      coord_sf(xlim = c(-180, 180), ylim = c(-60, 80)) +
      labs(title = paste(model_name, "- Spatial variation in slopes"),
           subtitle = sprintf("Range: %.3f - %.3f, SD: %.3f", 
                             min(knot_slopes), max(knot_slopes), sd(knot_slopes))) +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    # Create intercept map
    p_intercept <- ggplot() +
      geom_sf(data = world, fill = "gray90", color = "white") +
      geom_point(data = knot_df, 
                 aes(x = lon, y = lat, color = intercept, size = sqrt(density + 1)),
                 alpha = 0.8) +
      scale_color_gradient2(name = "Intercept (‰)", 
                           low = "blue", mid = "white", high = "red",
                           midpoint = mean(intercept_at_zero)) +
      scale_size_continuous(range = c(2, 8), guide = "none") +
      coord_sf(xlim = c(-180, 180), ylim = c(-60, 80)) +
      labs(title = paste(model_name, "- Spatial variation in intercepts"),
           subtitle = sprintf("Range: %.1f - %.1f ‰", 
                             min(intercept_at_zero), max(intercept_at_zero))) +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    # Save maps
    dir.create(file.path(output_dir, "maps"), showWarnings = FALSE)
    ggsave(file.path(output_dir, "maps", paste0(model_name, "_slope_map.png")),
           p_slope, width = 12, height = 6, dpi = 300)
    ggsave(file.path(output_dir, "maps", paste0(model_name, "_intercept_map.png")),
           p_intercept, width = 12, height = 6, dpi = 300)
    
    # 3. REGIONAL ANALYSIS
    cat("\n3. REGIONAL ANALYSIS\n")
    
    # Convert to sf object
    knot_sf <- st_as_sf(knot_df, coords = c("lon", "lat"), crs = 4326)
    
    # Join with continents
    continents <- world %>% 
      select(continent) %>%
      filter(!is.na(continent))
    
    knot_continent <- st_join(knot_sf, continents, join = st_nearest_feature)
    
    # Regional summary
    regional_summary <- knot_continent %>%
      st_drop_geometry() %>%
      group_by(continent) %>%
      summarise(
        n_knots = n(),
        mean_slope = mean(slope),
        sd_slope = sd(slope),
        min_slope = min(slope),
        max_slope = max(slope),
        mean_intercept = mean(intercept),
        sd_intercept = sd(intercept),
        .groups = "drop"
      ) %>%
      arrange(desc(n_knots))
    
    print(regional_summary)
    all_regional_summaries[[model_name]] <- regional_summary
    
    # Save regional summary
    write.csv(regional_summary, 
              file.path(output_dir, paste0(model_name, "_regional_summary.csv")),
              row.names = FALSE)
    
    # 4. SPATIAL PREDICTIONS AT DATA LOCATIONS
    cat("\n4. SPATIAL PREDICTIONS AT DATA LOCATIONS\n")
    
    # Extract per-observation spatial effects. Stan exports `alpha_spatial[N]`
    # (intercept at each location) and `beta_oipc_spatial[N]` (slope at each
    # location). The GP *residual* contributions are
    #   intercept residual = alpha_spatial - beta_0
    #   slope residual     = beta_oipc_spatial - beta_oipc
    # (4d_leaf_wax_spatial_model.stan:316,342,504-505). The legacy names
    # `gp_intercept` / `gp_slope` never existed in this model.
    has_alpha <- "alpha_spatial" %in% unique_vars
    has_slope <- "beta_oipc_spatial" %in% unique_vars

    if (has_alpha || has_slope) {
      if (has_alpha) {
        alpha_mean <- fit$summary(variables = "alpha_spatial")$mean
        beta_0_mean <- fit$summary(variables = "beta_0")$mean
        gp_intercept_obs <- alpha_mean[1:N] - beta_0_mean
      } else {
        gp_intercept_obs <- rep(0, N)
      }

      if (has_slope) {
        slope_mean <- fit$summary(variables = "beta_oipc_spatial")$mean
        beta_oipc_mean <- fit$summary(variables = "beta_oipc")$mean
        gp_slope_obs <- slope_mean[1:N] - beta_oipc_mean
      } else {
        gp_slope_obs <- rep(0, N)
      }

      # Observation-level data frame (effects are residuals from the global
      # mean, so they plot as deviations).
      obs_df <- data.frame(
        lon = stan_data$longitude * lon_sd + lon_mean,
        lat = stan_data$latitude * lat_sd + lat_mean,
        gp_intercept = gp_intercept_obs,
        gp_slope = gp_slope_obs,
        observed = stan_data$d2H_wax * d2h_sd + d2h_mean
      )
      
      # Map of GP effects at observations
      p_gp_obs <- ggplot() +
        geom_sf(data = world, fill = "gray90", color = "white") +
        geom_point(data = obs_df, 
                   aes(x = lon, y = lat, color = gp_slope),
                   size = 1.5, alpha = 0.6) +
        scale_color_viridis(name = "GP slope effect") +
        coord_sf(xlim = c(-180, 180), ylim = c(-60, 80)) +
        labs(title = paste(model_name, "- GP effects at observations")) +
        theme_minimal()
      
      ggsave(file.path(output_dir, paste0(model_name, "_gp_observations.png")),
             p_gp_obs, width = 12, height = 6, dpi = 300)
    }
    
    # 5. VARIOGRAM ANALYSIS
    cat("\n5. VARIOGRAM ANALYSIS\n")
    
    if (n_knots > 20) {
      # Calculate distances between knots
      knot_dists <- as.matrix(dist(knot_coords))
      
      # Convert to km (approximate)
      knot_dists_km <- knot_dists * 111 * mean(cos(knot_lat * pi/180))
      
      # Calculate slope differences
      slope_diffs <- outer(knot_slopes, knot_slopes, "-")
      slope_semivar <- 0.5 * slope_diffs^2
      
      # Bin by distance
      dist_breaks <- seq(0, max(knot_dists_km), length.out = 20)
      dist_mids <- (dist_breaks[-1] + dist_breaks[-length(dist_breaks)]) / 2
      
      # Calculate empirical variogram
      variogram_data <- data.frame(
        distance = as.vector(knot_dists_km[upper.tri(knot_dists_km)]),
        semivariance = as.vector(slope_semivar[upper.tri(slope_semivar)])
      ) %>%
        filter(distance > 0) %>%
        mutate(dist_bin = cut(distance, breaks = dist_breaks)) %>%
        group_by(dist_bin) %>%
        summarise(
          dist = mean(distance),
          gamma = mean(semivariance),
          n = n(),
          .groups = "drop"
        ) %>%
        filter(!is.na(dist_bin))
      
      # Plot variogram
      p_variogram <- ggplot(variogram_data, aes(x = dist, y = gamma)) +
        geom_point(aes(size = n), alpha = 0.7) +
        geom_smooth(method = "loess", se = FALSE, color = "red") +
        geom_hline(yintercept = var(knot_slopes), linetype = "dashed", alpha = 0.5) +
        geom_vline(xintercept = ls_slope, linetype = "dashed", color = "blue", alpha = 0.5) +
        scale_size_continuous(range = c(2, 6)) +
        labs(title = paste(model_name, "- Empirical variogram of slopes"),
             x = "Distance (km)", y = "Semivariance",
             subtitle = paste("Estimated range:", round(ls_slope, 0), "km")) +
        theme_minimal()
      
      ggsave(file.path(output_dir, paste0(model_name, "_variogram.png")),
             p_variogram, width = 8, height = 6, dpi = 300)
      
      all_variograms[[model_name]] <- variogram_data
    }
    
    # 6. PREDICTION GRID
    cat("\n6. CREATING PREDICTION GRID\n")
    
    # Create a regular grid for predictions
    lon_range <- range(stan_data$longitude * lon_sd + lon_mean)
    lat_range <- range(stan_data$latitude * lat_sd + lat_mean)
    
    # Extend slightly beyond data
    lon_buffer <- diff(lon_range) * 0.1
    lat_buffer <- diff(lat_range) * 0.1
    
    lon_grid <- seq(lon_range[1] - lon_buffer, lon_range[2] + lon_buffer, by = 5)
    lat_grid <- seq(lat_range[1] - lat_buffer, lat_range[2] + lat_buffer, by = 5)
    
    pred_grid <- expand.grid(lon = lon_grid, lat = lat_grid)
    
    # Remove ocean points (roughly)
    pred_grid_sf <- st_as_sf(pred_grid, coords = c("lon", "lat"), crs = 4326)
    land <- world %>% summarise(geometry = st_union(geometry))
    pred_grid_land <- st_intersection(pred_grid_sf, land)
    pred_coords <- st_coordinates(pred_grid_land)
    
    if (nrow(pred_coords) > 50) {
      cat("  Calculating predictions at", nrow(pred_coords), "grid points\n")
      
      # Standardize grid coordinates
      pred_lon_std <- (pred_coords[,1] - lon_mean) / lon_sd
      pred_lat_std <- (pred_coords[,2] - lat_mean) / lat_sd
      
      # Calculate distances from grid points to knots
      grid_knot_dists <- fields::rdist(
        cbind(pred_lon_std, pred_lat_std),
        knot_coords
      )
      
      # GP kernel function (squared exponential)
      ls_std <- ls_intercept / (111 * mean(cos(lat_mean * pi/180)) * mean(c(lon_sd, lat_sd)))
      K_grid_knot <- exp(-0.5 * (grid_knot_dists / ls_std)^2)
      
      # Normalize kernel
      K_grid_knot <- K_grid_knot / rowSums(K_grid_knot)
      
      # Predict slopes at grid points
      grid_slopes <- as.vector(K_grid_knot %*% matrix(knot_slopes, ncol = 1))
      grid_intercepts <- as.vector(K_grid_knot %*% matrix(intercept_at_zero, ncol = 1))
      
      # Create prediction data frame
      grid_pred_df <- data.frame(
        lon = pred_coords[,1],
        lat = pred_coords[,2],
        slope_pred = grid_slopes,
        intercept_pred = grid_intercepts
      )
      
      all_predictions[[model_name]] <- grid_pred_df
      
      # Create prediction map
      p_pred <- ggplot() +
        geom_sf(data = world, fill = "gray90", color = "white") +
        geom_point(data = grid_pred_df,
                   aes(x = lon, y = lat, color = slope_pred),
                   size = 0.8, alpha = 0.8) +
        scale_color_viridis(name = "Predicted\nSlope") +
        coord_sf(xlim = lon_range + c(-lon_buffer, lon_buffer),
                 ylim = lat_range + c(-lat_buffer, lat_buffer)) +
        labs(title = paste(model_name, "- Spatial predictions")) +
        theme_minimal()
      
      ggsave(file.path(output_dir, paste0(model_name, "_predictions.png")),
             p_pred, width = 10, height = 8, dpi = 300)
    }
    
    # 7. PREDICTION UNCERTAINTY
    cat("\n7. ANALYZING PREDICTION UNCERTAINTY\n")
    
    # Calculate posterior SD of predictions if available
    if ("mu" %in% unique_vars) {
      # Extract mu draws and calculate SD
      mu_draws <- fit$draws("mu", format = "matrix")
      mu_sd <- apply(mu_draws, 2, sd)
      
      # Create uncertainty map
      uncertainty_df <- data.frame(
        lon = stan_data$longitude * lon_sd + lon_mean,
        lat = stan_data$latitude * lat_sd + lat_mean,
        uncertainty = mu_sd * d2h_sd  # Convert to original scale
      )
      
      p_uncertainty <- ggplot() +
        geom_sf(data = world, fill = "gray90", color = "white") +
        geom_point(data = uncertainty_df, 
                   aes(x = lon, y = lat, color = uncertainty),
                   size = 1, alpha = 0.7) +
        scale_color_viridis(name = "Prediction\nSD (‰)", option = "plasma") +
        coord_sf(xlim = c(-180, 180), ylim = c(-60, 80)) +
        labs(title = paste(model_name, "- Prediction uncertainty")) +
        theme_minimal()
      
      ggsave(file.path(output_dir, paste0(model_name, "_uncertainty.png")),
             p_uncertainty, width = 12, height = 6, dpi = 300)
    }
    
    # 8. DENSITY VS SPATIAL EFFECTS
    cat("\n8. DENSITY VS SPATIAL EFFECTS\n")
    
    # Plot relationship between data density and spatial effects
    p_density <- ggplot(knot_df, aes(x = density)) +
      geom_point(aes(y = slope - mean(slope)), color = "darkgreen", 
                 alpha = 0.6, size = 3) +
      geom_smooth(aes(y = slope - mean(slope)), method = "loess", 
                  se = FALSE, color = "darkgreen") +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
      labs(title = paste(model_name, "- Data density vs spatial effects"),
           x = "Number of observations within radius",
           y = "Slope deviation from mean") +
      theme_minimal()
    
    ggsave(file.path(output_dir, paste0(model_name, "_density_effects.png")),
           p_density, width = 8, height = 6, dpi = 300)
    
    # Store all results for this model
    all_results[[model_name]] <- list(
      summary = data.frame(
        model = model_name,
        n_knots = n_knots,
        slope_mean = mean(knot_slopes),
        slope_sd = sd(knot_slopes),
        slope_min = min(knot_slopes),
        slope_max = max(knot_slopes),
        intercept_mean = mean(intercept_at_zero),
        intercept_sd = sd(intercept_at_zero),
        intercept_min = min(intercept_at_zero),
        intercept_max = max(intercept_at_zero),
        length_scale = ls_intercept,
        sigma_slope = if(exists("sigma_slope")) sigma_slope else NA,
        sigma_intercept = if(exists("sigma_intercept")) sigma_intercept else NA
      ),
      knot_data = knot_df,
      regional_summary = regional_summary
    )
    
  }, error = function(e) {
    cat("  Error analyzing", model_name, ":", e$message, "\n")
  })
}

# ======================================================================
# MULTI-MODEL COMPARISON
# ======================================================================

cat("\n\n ======================================================================\n")
cat("SPATIAL PATTERN SUMMARY\n")
cat("======================================================================\n\n")

# Extract summaries for spatial models only
spatial_models <- names(all_results)
summary_df <- bind_rows(lapply(all_results, function(x) x$summary))
all_knots <- bind_rows(all_knot_data)

if (length(spatial_models) > 0) {
  
  cat("Summary of spatial patterns:\n")
  print(summary_df)
  
  # Save overall summaries
  write.csv(summary_df, file.path(output_dir, "spatial_model_comparison.csv"), 
            row.names = FALSE)
  write.csv(all_knots, file.path(output_dir, "all_knots_combined.csv"),
            row.names = FALSE)
  
  # Create comparison plots
  
  # 1. Compare spatial variation across models
  p_compare_sd <- ggplot(summary_df, aes(x = reorder(model, -slope_sd))) +
    geom_col(aes(y = slope_sd), fill = "steelblue", alpha = 0.7) +
    coord_flip() +
    labs(title = "Spatial variation in slopes across models",
         x = "Model", y = "SD of slopes") +
    theme_minimal()
  
  # 2. Compare range of slopes
  p_compare_range <- ggplot(summary_df, aes(x = reorder(model, slope_mean))) +
    geom_pointrange(aes(y = slope_mean, ymin = slope_min, ymax = slope_max)) +
    coord_flip() +
    geom_hline(yintercept = 0.8, linetype = "dashed", alpha = 0.5) +
    labs(title = "Range of slopes across models",
         x = "Model", y = "Slope") +
    theme_minimal()
  
  # 3. Length scale comparison
  p_length_scale <- ggplot(summary_df %>% filter(!is.na(length_scale) & length_scale > 0), 
                          aes(x = reorder(model, -length_scale))) +
    geom_col(aes(y = length_scale), fill = "darkgreen", alpha = 0.7) +
    coord_flip() +
    labs(title = "Spatial correlation length scales",
         x = "Model", y = "Length scale (km)") +
    theme_minimal()
  
  # 4. Sigma comparison
  p_sigma <- ggplot(summary_df %>% filter(!is.na(sigma_slope)), 
                   aes(x = reorder(model, -sigma_slope))) +
    geom_col(aes(y = sigma_slope), fill = "darkred", alpha = 0.7) +
    coord_flip() +
    labs(title = "Spatial standard deviations (slopes)",
         x = "Model", y = "Sigma slope") +
    theme_minimal()
  
  # Combine plots
  p_combined <- (p_compare_sd | p_compare_range) / (p_length_scale | p_sigma)
  
  ggsave(file.path(output_dir, "model_comparison.png"),
         p_combined, width = 14, height = 10, dpi = 300)
  
  # Create best model visualization
  best_model <- "full_interact_sp"  # Based on your LOO results
  if (best_model %in% names(all_knot_data)) {
    
    best_knots <- all_knot_data[[best_model]]
    
    # Combined map showing both slope and intercept
    p_best <- ggplot() +
      geom_sf(data = world, fill = "gray90", color = "white") +
      geom_point(data = best_knots,
                 aes(x = lon, y = lat, color = slope, size = sqrt(density + 1)),
                 alpha = 0.8) +
      scale_color_viridis(name = "Slope") +
      scale_size_continuous(range = c(2, 8), name = "Data density",
                           breaks = sqrt(c(1, 10, 50, 100) + 1),
                           labels = c("0", "9", "49", "99")) +
      coord_sf(xlim = c(-180, 180), ylim = c(-60, 80)) +
      labs(title = paste("Best model:", best_model),
           subtitle = "Spatial variation in leaf wax response to precipitation isotopes") +
      theme_minimal() +
      theme(legend.position = "bottom",
            legend.box = "vertical")
    
    ggsave(file.path(output_dir, "best_model_spatial_pattern.png"),
           p_best, width = 12, height = 8, dpi = 300)
    
    # Regional comparison for best model
    if (best_model %in% names(all_regional_summaries)) {
      regional_best <- all_regional_summaries[[best_model]]
      
      p_regional <- ggplot(regional_best, aes(x = reorder(continent, mean_slope))) +
        geom_pointrange(aes(y = mean_slope, 
                           ymin = mean_slope - sd_slope,
                           ymax = mean_slope + sd_slope),
                       size = 0.8, fatten = 3) +
        geom_text(aes(y = mean_slope, label = paste0("n=", n_knots)),
                 hjust = -0.5, size = 3) +
        coord_flip() +
        labs(title = paste("Regional variation in", best_model),
             x = "Continent", y = "Mean slope ± SD") +
        theme_minimal()
      
      ggsave(file.path(output_dir, "best_model_regional.png"),
             p_regional, width = 8, height = 6, dpi = 300)
    }
  }
  
} else {
  cat("No spatial models found\n")
}

# Create diagnostic summary
issues <- c()
if (any(summary_df$slope_min < 0.1, na.rm = TRUE)) {
  issues <- c(issues, "Some models have unrealistically low slopes (<0.1)")
}
if (any(summary_df$slope_max > 1.5, na.rm = TRUE)) {
  issues <- c(issues, "Some models have unrealistically high slopes (>1.5)")
}
if (any(summary_df$length_scale < 100, na.rm = TRUE)) {
  issues <- c(issues, "Some models have very short length scales (<100 km)")
}

cat("\n\nPOTENTIAL ISSUES:\n")
cat("----------------------------------------\n")

if (length(issues) > 0) {
  for (issue in issues) {
    cat("-", issue, "\n")
  }
} else {
  cat("No major issues detected\n")
}

cat("\n\nDiagnostics complete! Results saved to:", output_dir)