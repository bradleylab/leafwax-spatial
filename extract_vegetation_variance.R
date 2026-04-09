# Extract spatial variance components for vegetation models
library(cmdstanr)
library(tidyverse)

cat("SPATIAL VARIANCE COMPONENTS FOR VEGETATION MODELS\n")
cat("==================================================\n\n")

# Models to check
spatial_veg_models <- c("baseline_veg_sp", "full_sp", "full_interact_sp")

# Store results
variance_results <- list()

for (model_name in spatial_veg_models) {
  cat("\nProcessing:", model_name, "\n")
  cat(strrep("-", 40), "\n")

  fit_file <- paste0("model_output/", model_name, "/fit.rds")

  if (file.exists(fit_file)) {
    fit <- readRDS(fit_file)

    # Extract spatial variance components (correct parameter names)
    if ("sigma_intercept_spatial" %in% fit$metadata()$variables) {
      sigma_int <- mean(as.vector(fit$draws("sigma_intercept_spatial", format = "matrix")))
      sigma_slope <- mean(as.vector(fit$draws("sigma_slope_spatial", format = "matrix")))
      sigma_obs <- mean(as.vector(fit$draws("sigma", format = "matrix")))

      # Calculate variance proportions
      total_var <- sigma_int^2 + sigma_slope^2 + sigma_obs^2

      cat("  Intercept GP variance: ", sprintf("%.1f%% (σ = %.2f‰)",
                                               100 * sigma_int^2 / total_var, sigma_int), "\n")
      cat("  Slope GP variance:     ", sprintf("%.1f%% (σ = %.2f)",
                                               100 * sigma_slope^2 / total_var, sigma_slope), "\n")
      cat("  Observation variance:  ", sprintf("%.1f%% (σ = %.2f‰)",
                                               100 * sigma_obs^2 / total_var, sigma_obs), "\n")
      cat("  Total GP variance:     ", sprintf("%.1f%%",
                                               100 * (sigma_int^2 + sigma_slope^2) / total_var), "\n")

      variance_results[[model_name]] <- data.frame(
        model = model_name,
        sigma_intercept = sigma_int,
        sigma_slope = sigma_slope,
        sigma_obs = sigma_obs,
        var_intercept = sigma_int^2,
        var_slope = sigma_slope^2,
        var_obs = sigma_obs^2,
        prop_intercept = sigma_int^2 / total_var,
        prop_slope = sigma_slope^2 / total_var,
        prop_obs = sigma_obs^2 / total_var,
        prop_gp_total = (sigma_int^2 + sigma_slope^2) / total_var
      )
    } else {
      cat("  Spatial variance parameters not found\n")
    }
  }
}

# Compare to baseline_sp if available
cat("\n\nComparison to baseline_sp (no vegetation):\n")
cat(strrep("-", 40), "\n")

baseline_sp_file <- "model_output/baseline_sp/fit.rds"
if (file.exists(baseline_sp_file)) {
  fit_baseline <- readRDS(baseline_sp_file)

  if ("sigma_intercept_spatial" %in% fit_baseline$metadata()$variables) {
    sigma_int_base <- mean(as.vector(fit_baseline$draws("sigma_intercept_spatial", format = "matrix")))
    sigma_slope_base <- mean(as.vector(fit_baseline$draws("sigma_slope_spatial", format = "matrix")))
    sigma_obs_base <- mean(as.vector(fit_baseline$draws("sigma", format = "matrix")))
    total_var_base <- sigma_int_base^2 + sigma_slope_base^2 + sigma_obs_base^2
    gp_prop_base <- (sigma_int_base^2 + sigma_slope_base^2) / total_var_base

    cat(sprintf("  Baseline_sp GP variance: %.1f%% (σ_int=%.2f‰, σ_slope=%.2f)\n",
                100 * gp_prop_base, sigma_int_base, sigma_slope_base))

    for (model_name in names(variance_results)) {
      reduction <- 100 * (gp_prop_base - variance_results[[model_name]]$prop_gp_total) / gp_prop_base
      cat(sprintf("  %s reduces GP variance by %.1f%%\n", model_name, reduction))
    }
  }
}

# Save results
if (length(variance_results) > 0) {
  variance_df <- bind_rows(variance_results)
  write.csv(variance_df, "model_analysis/tables/vegetation_variance_components.csv",
            row.names = FALSE)
  cat("\n\nVariance components saved to: model_analysis/tables/vegetation_variance_components.csv\n")

  # Print summary table
  cat("\nSummary Table:\n")
  cat(strrep("-", 80), "\n")
  cat(sprintf("%-20s %8s %8s %8s %8s %8s\n",
              "Model", "σ_int", "σ_slope", "σ_obs", "GP%", "Obs%"))
  cat(strrep("-", 80), "\n")
  for (i in 1:nrow(variance_df)) {
    cat(sprintf("%-20s %8.2f %8.2f %8.2f %8.1f %8.1f\n",
                variance_df$model[i],
                variance_df$sigma_intercept[i],
                variance_df$sigma_slope[i],
                variance_df$sigma_obs[i],
                100 * variance_df$prop_gp_total[i],
                100 * variance_df$prop_obs[i]))
  }
}

cat("\nAnalysis complete!\n")