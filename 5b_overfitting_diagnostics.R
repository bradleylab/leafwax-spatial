#!/usr/bin/env Rscript
#───────────────────────────────────────────────────────────────────────────────
# 5b_overfitting_diagnostics.R
#
# Comprehensive overfitting diagnostics for spatial models
# Tests model complexity, effective degrees of freedom, and generalization
# Performs cross-validation and out-of-sample prediction tests
#
# Input: results/model_fits/*.rds, results/3_sediment_ready_for_modeling.rds
# Output: results/overfitting_diagnostics/ (diagnostic plots and metrics)
#───────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(posterior)
source("scripts/posterior_helpers.R")
source("scripts/table_helpers.R")
library(loo)
library(sf)
library(viridis)

# Create output directory
output_dir <- "model_analysis/overfitting_diagnostics"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("\nOVERFITTING DIAGNOSTICS\n")
cat("=======================\n")

#───────────────────────────────────────────────────────────────────────────────
# 1. SPATIAL CROSS-VALIDATION
#───────────────────────────────────────────────────────────────────────────────

spatial_cv_analysis <- function(model_name, prepared_data_dir = NULL) {
  cat("\n1. SPATIAL CROSS-VALIDATION for", model_name, "\n")
  cat(strrep("-", 50), "\n")

  draws     <- load_draws(model_name)
  stan_data <- load_stan_data(model_name)

  y_rep <- as_draws_matrix(subset_draws(draws, variable = "d2H_rep"))
  y_pred_mean <- colMeans(y_rep)
  y_obs <- stan_data$d2H_wax
  
  # Convert back to original scale
  d2H_sd_original <- stan_data$scaling_params$d2H_sd
  y_obs_original <- y_obs * d2H_sd_original
  y_pred_mean_original <- y_pred_mean * d2H_sd_original
  residuals_original <- y_obs_original - y_pred_mean_original
  
  # Get coordinates
  coords <- data.frame(
    lon = stan_data$longitude,
    lat = stan_data$latitude,
    obs = y_obs_original,
    pred = y_pred_mean_original,
    residual = residuals_original
  )
  
  # Define spatial regions for cross-validation - FIX: use coords$ prefix
  coords$region <- case_when(
    coords$lon < -30 ~ "Americas",
    coords$lon >= -30 & coords$lon < 60 & coords$lat > 35 ~ "Europe",
    coords$lon >= -30 & coords$lon < 60 & coords$lat <= 35 ~ "Africa", 
    coords$lon >= 60 & coords$lon < 140 & coords$lat > -10 ~ "Asia",
    coords$lon >= 140 | coords$lat < -10 ~ "Oceania",
    TRUE ~ "Other"  # Add catch-all
  )
  
  # Calculate performance by region
  regional_performance <- coords %>%
    group_by(region) %>%
    summarise(
      n = n(),
      rmse = sqrt(mean(residual^2)),
      mae = mean(abs(residual)),
      bias = mean(residual),
      r_squared = cor(obs, pred)^2,
      .groups = "drop"
    )
  
  cat("\nRegional Performance:\n")
  print(regional_performance)
  
  # Test: Performance in sparse vs dense regions
  # Fix the density quartile calculation
  if (stan_data$include_gp == 1) {
    # For each observation, find nearest knot density
    knot_coords <- stan_data$knot_coords
    obs_coords_std <- cbind((coords$lon - mean(coords$lon))/sd(coords$lon),
                            (coords$lat - mean(coords$lat))/sd(coords$lat))
    
    coords$nearest_knot_density <- NA
    for (i in 1:nrow(coords)) {
      dists <- sqrt((obs_coords_std[i,1] - knot_coords[,1])^2 + 
                    (obs_coords_std[i,2] - knot_coords[,2])^2)
      nearest_knot <- which.min(dists)
      coords$nearest_knot_density[i] <- stan_data$knot_data_density[nearest_knot]
    }
    
    coords$density_group <- cut(coords$nearest_knot_density,
                                breaks = c(-Inf, 10, 50, Inf),
                                labels = c("Sparse (<10)", "Medium (10-50)", "Dense (>50)"))
    
    density_performance <- coords %>%
      group_by(density_group) %>%
      summarise(
        n = n(),
        rmse = sqrt(mean(residual^2)),
        mae = mean(abs(residual)),
        .groups = "drop"
      )
    
    cat("\nPerformance by Data Density:\n")
    print(density_performance)
  } else {
    density_performance <- NULL
  }
  
  # Create diagnostic plot with original scale residuals
  p_spatial <- ggplot(coords, aes(x = lon, y = lat)) +
    geom_point(aes(color = abs(residual), shape = region), size = 3, alpha = 0.8) +
    scale_color_viridis(name = "|Residual| (‰)", limits = c(0, 30)) +  # Updated scale
    theme_minimal() +
    labs(title = paste("Spatial Distribution of Residuals:", model_name),
         subtitle = "Color shows absolute residual in ‰, shape shows region")
  
  ggsave(file.path(output_dir, paste0("spatial_residuals_regions_", model_name, ".png")),
         p_spatial, width = 10, height = 6, dpi = 300)
  
  return(list(regional = regional_performance, density = density_performance))
}

#───────────────────────────────────────────────────────────────────────────────
# 2. COMPLEXITY ANALYSIS
#───────────────────────────────────────────────────────────────────────────────

complexity_analysis <- function(model_name, prepared_data_dir = NULL) {
  cat("\n\n2. MODEL COMPLEXITY ANALYSIS for", model_name, "\n")
  cat(strrep("-", 50), "\n")

  loo_result <- tryCatch(load_loo(model_name),
                         error = function(e) { cat("  ", conditionMessage(e), "\n"); NULL })
  if (is.null(loo_result)) return(NULL)
  
  # Extract diagnostics
  p_loo <- loo_result$estimates["p_loo", "Estimate"]
  n_obs <- length(loo_result$diagnostics$pareto_k)
  
  cat("  Number of observations:", n_obs, "\n")
  cat("  Effective number of parameters (p_loo):", round(p_loo, 1), "\n")
  cat("  Ratio p_loo/n:", round(p_loo/n_obs, 3), "\n")
  
  # Check Pareto k diagnostics
  k_values <- loo_result$diagnostics$pareto_k
  k_summary <- c(
    sum(k_values > 0.5),
    sum(k_values > 0.7),
    sum(k_values > 1.0)
  )
  
  cat("\n  Pareto k diagnostics:\n")
  cat("    k > 0.5:", k_summary[1], "(", round(100*k_summary[1]/n_obs, 1), "%)\n")
  cat("    k > 0.7:", k_summary[2], "(", round(100*k_summary[2]/n_obs, 1), "%)\n")
  cat("    k > 1.0:", k_summary[3], "(", round(100*k_summary[3]/n_obs, 1), "%)\n")
  
  # Identify problematic observations
  if (k_summary[2] > 0) {
    cat("\n  Problematic observations (k > 0.7):\n")
    prob_idx <- which(k_values > 0.7)
    stan_data <- load_stan_data(model_name)
    
    prob_obs <- data.frame(
      idx = prob_idx,
      lon = stan_data$longitude[prob_idx],
      lat = stan_data$latitude[prob_idx],
      k = k_values[prob_idx]
    )
    print(head(prob_obs, 10))
  }
  
  return(list(p_loo = p_loo, n_obs = n_obs, k_summary = k_summary))
}

#───────────────────────────────────────────────────────────────────────────────
# 3. PREDICTION UNCERTAINTY ANALYSIS
#───────────────────────────────────────────────────────────────────────────────

uncertainty_analysis <- function(model_name, prepared_data_dir = NULL) {
  cat("\n\n3. PREDICTION UNCERTAINTY ANALYSIS for", model_name, "\n")
  cat(strrep("-", 50), "\n")

  draws     <- load_draws(model_name)
  stan_data <- load_stan_data(model_name)

  y_rep <- as_draws_matrix(subset_draws(draws, variable = "d2H_rep"))
  
  # Calculate prediction intervals
  pred_mean <- colMeans(y_rep)
  pred_sd <- apply(y_rep, 2, sd)
  pred_lower <- apply(y_rep, 2, quantile, 0.025)
  pred_upper <- apply(y_rep, 2, quantile, 0.975)
  interval_width <- pred_upper - pred_lower
  
  # Back-transform to original scale
  scaling_params <- stan_data$scaling_params
  interval_width_permil <- interval_width * scaling_params$d2H_sd
  pred_sd_permil <- pred_sd * scaling_params$d2H_sd
  
  cat("  Prediction interval widths (‰):\n")
  cat("    Mean:", round(mean(interval_width_permil), 1), "\n")
  cat("    SD:", round(sd(interval_width_permil), 1), "\n")
  cat("    Range:", round(range(interval_width_permil), 1), "\n")
  
  # Check if uncertainty varies with data density
  coords <- data.frame(
    lon = stan_data$longitude,
    lat = stan_data$latitude,
    pred_sd = pred_sd_permil,
    interval_width = interval_width_permil
  )
  
  # For each observation, find distance to nearest other observation
  dist_matrix <- as.matrix(dist(coords[,c("lon", "lat")]))
  diag(dist_matrix) <- NA
  coords$nearest_neighbor_dist <- apply(dist_matrix, 1, min, na.rm = TRUE)
  
  # Plot uncertainty vs isolation
  p_uncertainty <- ggplot(coords, aes(x = nearest_neighbor_dist, y = pred_sd)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "loess", se = TRUE) +
    labs(x = "Distance to Nearest Observation (degrees)",
         y = "Prediction SD (‰)",
         title = paste("Prediction Uncertainty vs Data Sparsity:", model_name)) +
    theme_minimal()
  
  ggsave(file.path(output_dir, paste0("uncertainty_vs_sparsity_", model_name, ".png")),
         p_uncertainty, width = 8, height = 6, dpi = 300)
  
  # Check for suspiciously narrow intervals
  narrow_threshold <- 5  # ‰
  n_narrow <- sum(interval_width_permil < narrow_threshold)
  
  if (n_narrow > 0) {
    cat("\n  WARNING:", n_narrow, "observations have 95% intervals < ", narrow_threshold, "‰\n")
    cat("  This might indicate overconfidence/overfitting\n")
  }
  
  return(coords)
}

#───────────────────────────────────────────────────────────────────────────────
# 4. SPATIAL SCALE ANALYSIS
#───────────────────────────────────────────────────────────────────────────────

spatial_scale_analysis <- function(model_name, prepared_data_dir = NULL) {
  cat("\n\n4. SPATIAL SCALE ANALYSIS for", model_name, "\n")
  cat(strrep("-", 50), "\n")

  config <- load_config(model_name)
  if (!config$include_gp) {
    cat("  No GP component - skipping\n")
    return(NULL)
  }

  draws        <- load_draws(model_name)
  vars_present <- variables(draws)

  if ("ls_intercept_km" %in% vars_present) {
    ls_int   <- as.numeric(as_draws_matrix(subset_draws(draws, variable = "ls_intercept_km")))
    ls_slope <- as.numeric(as_draws_matrix(subset_draws(draws, variable = "ls_slope_km")))
    
    cat("  Length scale for intercept (km):\n")
    cat("    Mean:", round(mean(ls_int)), "\n")
    cat("    95% CI:", round(quantile(ls_int, c(0.025, 0.975))), "\n")
    
    cat("\n  Length scale for slope (km):\n") 
    cat("    Mean:", round(mean(ls_slope)), "\n")
    cat("    95% CI:", round(quantile(ls_slope, c(0.025, 0.975))), "\n")
    
    # Compare to data spacing
    stan_data <- load_stan_data(model_name)
    coords <- cbind(stan_data$longitude, stan_data$latitude)
    dist_matrix <- fields::rdist.earth(coords, coords, miles = FALSE)
    diag(dist_matrix) <- NA
    
    mean_nn_dist <- mean(apply(dist_matrix, 1, min, na.rm = TRUE))
    cat("\n  Mean nearest-neighbor distance:", round(mean_nn_dist), "km\n")
    cat("  Length scale / NN distance ratio:", round(mean(ls_int)/mean_nn_dist, 1), "\n")
    
    if (mean(ls_int) < 2 * mean_nn_dist) {
      cat("\n  WARNING: Length scale might be too short relative to data spacing\n")
      cat("  This could indicate overfitting to individual observations\n")
    }
  }
  
  # Extract scale weights
  if ("scale_weights" %in% vars_present) {
    scale_weights <- as_draws_matrix(subset_draws(draws, variable = "scale_weights"))
    scale_weight_means <- colMeans(scale_weights)
    
    cat("\n  Distance scale weights:\n")
    distances <- stan_data$distance_scales
    for (i in 1:length(distances)) {
      cat("    ", distances[i], "km:", round(scale_weight_means[i], 3), "\n")
    }
    
    # Effective scale
    effective_scale <- sum(distances * scale_weight_means)
    cat("\n  Effective integration scale:", round(effective_scale), "km\n")
  }
}

#───────────────────────────────────────────────────────────────────────────────
# 5. COMPARE REDUCED MODELS
#───────────────────────────────────────────────────────────────────────────────

compare_model_sizes <- function() {
  cat("\n\n5. MODEL SIZE COMPARISON\n")
  cat(strrep("-", 50), "\n")

  # Discover all models with loo.rds + config rds in the April mirror.
  model_dirs <- list.dirs(APRIL_RUN, full.names = FALSE, recursive = FALSE)
  model_dirs <- model_dirs[!grepl("^_", model_dirs)]
  models <- model_dirs[sapply(model_dirs, function(m)
    file.exists(file.path(APRIL_RUN, m, "loo.rds")))]

  if (length(models) == 0) {
    cat("  No models with LOO results found\n")
    return(NULL)
  }

  size_comparison <- map_df(models, function(m) {
    loo_result <- tryCatch(load_loo(m), error = function(e) NULL)
    config     <- tryCatch(load_config(m), error = function(e) NULL)
    if (is.null(loo_result) || is.null(config)) return(NULL)
    data.frame(
      model   = m,
      n_knots = config$n_pp_knots,
      p_loo   = loo_result$estimates["p_loo", "Estimate"],
      looic   = loo_result$estimates["looic", "Estimate"],
      stringsAsFactors = FALSE
    )
  })
  
  if (nrow(size_comparison) == 0) {
    cat("  No complete model results found\n")
    return(NULL)
  }
  
  cat("\nModel complexity comparison:\n")
  print(size_comparison)
  
  # Plot p_loo vs performance if we have data
  if (nrow(size_comparison) > 0) {
    p_size <- ggplot(size_comparison, aes(x = p_loo, y = looic, label = model)) +
      geom_point(size = 4) +
      geom_text(hjust = -0.1, vjust = -0.5) +
      labs(x = "Effective Parameters (p_loo)",
           y = "LOOIC (lower is better)",
           title = "Model Complexity vs Performance") +
      theme_minimal()
    
    ggsave(file.path(output_dir, "complexity_vs_performance.png"),
           p_size, width = 8, height = 6, dpi = 300)
  }
  
  return(size_comparison)
}

#───────────────────────────────────────────────────────────────────────────────
# MAIN EXECUTION
#───────────────────────────────────────────────────────────────────────────────

# Automatically discover all models with widened posterior_draws.rds in the
# April mirror. Skip the _prepared_data subdir.
model_dirs <- list.dirs(APRIL_RUN, full.names = FALSE, recursive = FALSE)
model_dirs <- model_dirs[!grepl("^_", model_dirs)]
fitted_models <- character()
for (model_name in model_dirs) {
  if (file.exists(file.path(APRIL_RUN, model_name, "posterior_draws.rds"))) {
    fitted_models <- c(fitted_models, model_name)
  }
}
if (length(fitted_models) == 0) {
  stop("No fitted models found at ", APRIL_RUN)
}

cat("Found", length(fitted_models), "fitted models:", paste(fitted_models, collapse = ", "), "\n")

# prepared_data_dir retained as a legacy function argument; helpers
# (load_stan_data, load_config) resolve paths internally against APRIL_RUN.
prepared_data_dir <- NULL

all_results <- list()

for (model in fitted_models) {
  cat("\n", strrep("=", 70), "\n")
  cat("ANALYZING MODEL:", model, "\n")
  cat(strrep("=", 70), "\n")

  results <- list()

  tryCatch({
    results$spatial_cv <- spatial_cv_analysis(model, prepared_data_dir)
  }, error = function(e) cat("  Spatial CV error:", e$message, "\n"))

  tryCatch({
    results$complexity <- complexity_analysis(model, prepared_data_dir)
  }, error = function(e) cat("  Complexity error:", e$message, "\n"))

  tryCatch({
    results$uncertainty <- uncertainty_analysis(model, prepared_data_dir)
  }, error = function(e) cat("  Uncertainty error:", e$message, "\n"))

  tryCatch({
    results$spatial_scale <- spatial_scale_analysis(model, prepared_data_dir)
  }, error = function(e) cat("  Spatial scale error:", e$message, "\n"))
  
  all_results[[model]] <- results
}

# Overall comparison
size_comp <- compare_model_sizes()

# Summary report
cat("\n\n", strrep("=", 70), "\n")
cat("OVERFITTING DIAGNOSTIC SUMMARY\n")
cat(strrep("=", 70), "\n")

cat("\nKey indicators to check:\n")
cat("1. Regional performance should be similar across regions\n")
cat("2. Performance shouldn't degrade in data-sparse areas\n")
cat("3. p_loo/n ratio should be < 0.1-0.2\n")
cat("4. Few observations should have Pareto k > 0.7\n")
cat("5. Prediction intervals should widen in data-sparse regions\n")
cat("6. Length scales should be >> nearest neighbor distance\n")

# Save full results
saveRDS(all_results, file.path(output_dir, "overfitting_diagnostics.rds"))
write.csv(size_comp, file.path(output_dir, "model_size_comparison.csv"), row.names = FALSE)

#───────────────────────────────────────────────────────────────────────────────
# TABLE S2 PRODUCER — manuscript/tables/Table_S2_regional_performance_body.tex
#───────────────────────────────────────────────────────────────────────────────
#
# Regional-performance rows (Overall + 5 regions) per model. Numbers come
# straight out of spatial_cv_analysis()'s regional_performance df plus an
# "Overall" row computed from the same residuals across all sites.
#
# The body fragment preserves the prior table's row coloring (lightgray
# for the Overall row) and bolded Model / Overall / N columns, so the
# wrapper .tex can ship as a pure preamble-only document.
#
# Companion CSV (table_S2_regional_performance.csv) feeds the numeric
# audit.

.latex_model_name_s2 <- function(name) {
  if (grepl("_sp$", name)) {
    base <- sub("_sp$", "", name)
    base <- gsub("_", "\\_", base, fixed = TRUE)
    paste0(base, "\\_sp")
  } else {
    gsub("_", "\\_", name, fixed = TRUE)
  }
}

# Preserve config.yaml model order: present spatial and non-spatial pairs
# side by side, with paleo-simplicity models at the end.
s2_model_order <- c(
  "baseline", "baseline_sp",
  "baseline_env", "baseline_env_sp",
  "baseline_veg", "baseline_veg_sp",
  "full", "full_sp",
  "full_interact", "full_interact_sp",
  "elevation_only_sp", "elevation_c4_sp",
  "c4_only_sp", "elevation_c4_interact_sp"
)
s2_model_order <- s2_model_order[s2_model_order %in% names(all_results)]

.compute_overall <- function(model_name) {
  # Load the same residuals the regional breakdown used so Overall is
  # rigorously consistent with the regional rows.
  draws     <- load_draws(model_name)
  stan_data <- load_stan_data(model_name)
  y_rep <- as_draws_matrix(subset_draws(draws, variable = "d2H_rep"))
  y_pred_mean <- colMeans(y_rep)
  y_obs <- stan_data$d2H_wax
  d2H_sd <- stan_data$scaling_params$d2H_sd
  y_obs_o  <- y_obs       * d2H_sd
  y_pred_o <- y_pred_mean * d2H_sd
  data.frame(
    region = "Overall",
    n         = length(y_obs),
    rmse      = sqrt(mean((y_obs_o - y_pred_o)^2)),
    r_squared = cor(y_obs_o, y_pred_o)^2,
    stringsAsFactors = FALSE
  )
}

s2_rows_csv <- list()
s2_lines    <- character()

for (m in s2_model_order) {
  reg <- all_results[[m]]$spatial_cv$regional
  if (is.null(reg)) next
  overall <- .compute_overall(m)

  m_escaped <- .latex_model_name_s2(m)
  # Overall row: lightgray background, bold model/region/N/RMSE/R².
  s2_lines <- c(s2_lines,
    sprintf(
      "\\rowcolor{lightgray} \\textbf{%s} & \\textbf{Overall} & \\textbf{%d} & \\textbf{%.1f} & \\textbf{%.3f} \\\\",
      m_escaped, as.integer(overall$n), overall$rmse, overall$r_squared
    )
  )
  s2_rows_csv[[length(s2_rows_csv) + 1]] <- data.frame(
    model = m, region = "Overall",
    n = as.integer(overall$n), rmse = overall$rmse,
    r_squared = overall$r_squared, stringsAsFactors = FALSE
  )

  # Regional rows: empty model column (visually continuing the model).
  region_order <- c("Africa", "Americas", "Asia", "Europe", "Oceania", "Other")
  reg_sorted <- reg[match(intersect(region_order, reg$region), reg$region), ]
  for (i in seq_len(nrow(reg_sorted))) {
    r <- reg_sorted[i, ]
    s2_lines <- c(s2_lines,
      sprintf(" & %s & %d & %.1f & %.3f \\\\",
              r$region, as.integer(r$n), r$rmse, r$r_squared)
    )
    s2_rows_csv[[length(s2_rows_csv) + 1]] <- data.frame(
      model = m, region = r$region,
      n = as.integer(r$n), rmse = r$rmse,
      r_squared = r$r_squared, stringsAsFactors = FALSE
    )
  }
}

s2_csv <- bind_rows(s2_rows_csv)
write.csv(s2_csv,
          "manuscript/tables/Table_S2_regional_performance.csv",
          row.names = FALSE)

s2_path <- "manuscript/tables/Table_S2_regional_performance_body.tex"
writeLines(c(
  "% Auto-generated by 5b_overfitting_diagnostics.R.",
  "% Do not hand-edit — regenerate via `make tables`.",
  "",
  s2_lines
), s2_path)
cat("wrote", s2_path, "\n")

cat("\n\nDiagnostics complete! Results saved to:", output_dir, "\n")