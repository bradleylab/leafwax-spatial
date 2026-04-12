# Extract vegetation coefficients from all models
library(tidyverse)
library(cmdstanr)

# Create output directory
output_dir <- "model_analysis/tables"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("VEGETATION COEFFICIENT EXTRACTION\n")
cat("=================================\n\n")

# Models to check
models_to_check <- c("baseline_veg", "baseline_veg_sp", "full", "full_sp",
                     "full_interact", "full_interact_sp")

# Initialize results list
all_coefficients <- list()

# Extract coefficients for each model
for (model_name in models_to_check) {
  cat("\nProcessing:", model_name, "\n")
  cat(strrep("-", 40), "\n")

  # Check if model output exists
  fit_file <- paste0("model_output/", model_name, "/fit.rds")
  config_file <- paste0("prepared_data/config_", model_name, ".rds")

  if (!file.exists(fit_file)) {
    cat("  Model fit not found - skipping\n")
    next
  }
  if (!file.exists(config_file)) {
    cat("  Config not found - skipping\n")
    next
  }

  # Load model and config
  fit <- readRDS(fit_file)
  config <- readRDS(config_file)

  cat("  Model configuration:\n")
  cat("    Include C4:", config$include_c4, "\n")
  cat("    Include PFT:", config$include_pft, "\n")
  cat("    Include interactions:", config$include_veg_interactions, "\n")

  # Extract vegetation coefficients
  veg_params <- character()

  # C4 coefficient (all vegetation models have this)
  if (config$include_c4) {
    veg_params <- c(veg_params, "beta_c4")
  }

  # PFT coefficients
  if (config$include_pft) {
    veg_params <- c(veg_params, "beta_tree", "beta_shrub", "beta_grass")
  }

  # Extract draws for each parameter
  for (param in veg_params) {
    if (param %in% fit$metadata()$variables) {
      draws <- as.vector(fit$draws(param, format = "matrix"))

      # Calculate statistics
      param_stats <- data.frame(
        model = model_name,
        parameter = param,
        mean = mean(draws),
        median = median(draws),
        sd = sd(draws),
        lower_95 = quantile(draws, 0.025),
        upper_95 = quantile(draws, 0.975),
        lower_50 = quantile(draws, 0.25),
        upper_50 = quantile(draws, 0.75),
        n_eff = NA,  # Will add if available
        rhat = NA,   # Will add if available
        stringsAsFactors = FALSE
      )

      # Check if significantly different from zero
      param_stats$significant <- (param_stats$lower_95 > 0) | (param_stats$upper_95 < 0)

      # Get convergence diagnostics if available
      summary_df <- fit$summary(param)
      if (nrow(summary_df) > 0) {
        param_stats$n_eff <- summary_df$ess_bulk[1]
        param_stats$rhat <- summary_df$rhat[1]
      }

      all_coefficients[[paste(model_name, param, sep = "_")]] <- param_stats

      # Print summary
      cat("\n  ", param, ":\n")
      cat("    Mean (SD):", sprintf("%.3f (%.3f)", param_stats$mean, param_stats$sd), "\n")
      cat("    Median:", sprintf("%.3f", param_stats$median), "\n")
      cat("    95% CI: [", sprintf("%.3f, %.3f", param_stats$lower_95, param_stats$upper_95), "]\n")
      cat("    Significant:", param_stats$significant, "\n")
      if (!is.na(param_stats$rhat)) {
        cat("    Rhat:", sprintf("%.3f", param_stats$rhat), "\n")
      }
    } else {
      cat("\n  ", param, ": Not found in model\n")
    }
  }
}

# Combine all results
vegetation_table <- bind_rows(all_coefficients)

# Save main table
write.csv(vegetation_table, file.path(output_dir, "vegetation_coefficients.csv"),
          row.names = FALSE)
cat("\n\nVegetation coefficients saved to:",
    file.path(output_dir, "vegetation_coefficients.csv"), "\n")

# Create summary comparison
cat("\n\n")
cat("VEGETATION EFFECT SUMMARY\n")
cat("=========================\n\n")

# 1. Significant effects
cat("Significant Vegetation Effects (95% CI excludes zero):\n")
cat(strrep("-", 60), "\n")
sig_effects <- vegetation_table %>%
  filter(significant) %>%
  arrange(model, parameter)

if (nrow(sig_effects) > 0) {
  for (i in 1:nrow(sig_effects)) {
    cat(sprintf("%-20s %-12s: %.3f [%.3f, %.3f]\n",
                sig_effects$model[i], sig_effects$parameter[i],
                sig_effects$mean[i], sig_effects$lower_95[i], sig_effects$upper_95[i]))
  }
} else {
  cat("No significant vegetation effects found\n")
}

# 2. Compare PFT effects across models
cat("\n\nPFT Effect Comparison (mean ± SD):\n")
cat(strrep("-", 60), "\n")
pft_params <- c("beta_tree", "beta_shrub", "beta_grass")
for (param in pft_params) {
  param_data <- vegetation_table %>% filter(parameter == param)
  if (nrow(param_data) > 0) {
    cat("\n", param, ":\n")
    for (i in 1:nrow(param_data)) {
      sig_marker <- if(param_data$significant[i]) "*" else " "
      cat(sprintf("  %-20s: %7.3f ± %5.3f%s\n",
                  param_data$model[i], param_data$mean[i], param_data$sd[i], sig_marker))
    }
  }
}
cat("\n* = significant at 95% level\n")

# 3. C4 effect comparison
cat("\n\nC4 Effect Comparison:\n")
cat(strrep("-", 60), "\n")
c4_data <- vegetation_table %>% filter(parameter == "beta_c4")
if (nrow(c4_data) > 0) {
  for (i in 1:nrow(c4_data)) {
    sig_marker <- if(c4_data$significant[i]) "*" else " "
    cat(sprintf("%-20s: %7.3f [%7.3f, %7.3f]%s\n",
                c4_data$model[i], c4_data$mean[i],
                c4_data$lower_95[i], c4_data$upper_95[i], sig_marker))
  }
}

# 4. Relative importance within models
cat("\n\nRelative Vegetation Importance (absolute mean effects):\n")
cat(strrep("-", 60), "\n")
for (model in unique(vegetation_table$model)) {
  model_data <- vegetation_table %>%
    filter(model == !!model) %>%
    mutate(abs_mean = abs(mean)) %>%
    arrange(desc(abs_mean))

  if (nrow(model_data) > 0) {
    cat("\n", model, ":\n")
    cat("  Strongest:", model_data$parameter[1],
        sprintf("(|β| = %.3f)", model_data$abs_mean[1]), "\n")
    if (nrow(model_data) > 1) {
      for (i in 2:min(nrow(model_data), 4)) {
        cat("           ", model_data$parameter[i],
            sprintf("(|β| = %.3f)", model_data$abs_mean[i]), "\n")
      }
    }
  }
}

# 5. Extract variance components for spatial models
cat("\n\nVariance Components for Spatial Vegetation Models:\n")
cat(strrep("-", 60), "\n")

spatial_veg_models <- c("baseline_veg_sp", "full_sp", "full_interact_sp")
variance_results <- list()

for (model_name in spatial_veg_models) {
  fit_file <- paste0("model_output/", model_name, "/fit.rds")

  if (file.exists(fit_file)) {
    fit <- readRDS(fit_file)

    # Extract GP variance components
    if ("sigma_gp_intercept" %in% fit$metadata()$variables) {
      sigma_int <- mean(as.vector(fit$draws("sigma_gp_intercept", format = "matrix")))
      sigma_slope <- mean(as.vector(fit$draws("sigma_gp_slope", format = "matrix")))
      sigma_obs <- mean(as.vector(fit$draws("sigma", format = "matrix")))

      # Calculate variance proportions
      total_var <- sigma_int^2 + sigma_slope^2 + sigma_obs^2

      cat("\n", model_name, ":\n")
      cat(sprintf("  Intercept GP variance: %.1f%% (σ = %.2f‰)\n",
                  100 * sigma_int^2 / total_var, sigma_int))
      cat(sprintf("  Slope GP variance:     %.1f%% (σ = %.2f)\n",
                  100 * sigma_slope^2 / total_var, sigma_slope))
      cat(sprintf("  Observation variance:  %.1f%% (σ = %.2f‰)\n",
                  100 * sigma_obs^2 / total_var, sigma_obs))
      cat(sprintf("  Total GP variance:     %.1f%%\n",
                  100 * (sigma_int^2 + sigma_slope^2) / total_var))

      variance_results[[model_name]] <- list(
        sigma_intercept = sigma_int,
        sigma_slope = sigma_slope,
        sigma_obs = sigma_obs,
        prop_gp = (sigma_int^2 + sigma_slope^2) / total_var
      )
    }
  }
}

# Compare to baseline_sp if available
baseline_sp_file <- "model_output/baseline_sp/fit.rds"
if (file.exists(baseline_sp_file)) {
  fit_baseline <- readRDS(baseline_sp_file)
  if ("sigma_gp_intercept" %in% fit_baseline$metadata()$variables) {
    sigma_int_base <- mean(as.vector(fit_baseline$draws("sigma_gp_intercept", format = "matrix")))
    sigma_slope_base <- mean(as.vector(fit_baseline$draws("sigma_gp_slope", format = "matrix")))
    sigma_obs_base <- mean(as.vector(fit_baseline$draws("sigma", format = "matrix")))
    total_var_base <- sigma_int_base^2 + sigma_slope_base^2 + sigma_obs_base^2
    gp_prop_base <- (sigma_int_base^2 + sigma_slope_base^2) / total_var_base

    cat("\n\nComparison to baseline_sp (no vegetation):\n")
    cat(sprintf("  Baseline GP variance: %.1f%%\n", 100 * gp_prop_base))

    for (model_name in names(variance_results)) {
      reduction <- 100 * (gp_prop_base - variance_results[[model_name]]$prop_gp) / gp_prop_base
      cat(sprintf("  %s reduces GP variance by %.1f%%\n", model_name, reduction))
    }
  }
}

# Save variance results
if (length(variance_results) > 0) {
  variance_df <- map_df(names(variance_results), function(m) {
    data.frame(
      model = m,
      sigma_intercept = variance_results[[m]]$sigma_intercept,
      sigma_slope = variance_results[[m]]$sigma_slope,
      sigma_obs = variance_results[[m]]$sigma_obs,
      prop_gp = variance_results[[m]]$prop_gp
    )
  })
  write.csv(variance_df, file.path(output_dir, "vegetation_variance_components.csv"),
            row.names = FALSE)
  cat("\n\nVariance components saved to:",
      file.path(output_dir, "vegetation_variance_components.csv"), "\n")
}

cat("\n\nAnalysis complete!\n")