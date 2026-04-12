# Comprehensive analysis of full models for Section 3.4
library(tidyverse)
library(cmdstanr)

# Create output directory
output_dir <- "model_analysis/tables"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("FULL MODEL ANALYSIS FOR SECTION 3.4\n")
cat("====================================\n\n")

#───────────────────────────────────────────────────────────────────────────────
# 1. MODEL PERFORMANCE COMPARISON
#───────────────────────────────────────────────────────────────────────────────

cat("1. MODEL PERFORMANCE METRICS\n")
cat(strrep("-", 60), "\n\n")

# Models to compare
full_models <- c("full", "full_sp", "full_interact", "full_interact_sp")
comparison_models <- c("baseline", "baseline_sp", "baseline_env_sp", "baseline_veg_sp")
all_models <- c(full_models, comparison_models)

performance_metrics <- list()

for (model_name in all_models) {
  # Check if model exists
  loo_file <- paste0("model_output/", model_name, "/loo.rds")
  fit_file <- paste0("model_output/", model_name, "/fit.rds")

  if (file.exists(loo_file) && file.exists(fit_file)) {
    # Load LOO
    loo_result <- readRDS(loo_file)

    # Load fit for R² and RMSE calculation
    fit <- readRDS(fit_file)
    stan_data <- readRDS(paste0("prepared_data/stan_data_", model_name, ".rds"))

    # Get posterior predictions
    y_rep <- fit$draws("d2H_rep", format = "matrix")
    y_pred_mean <- colMeans(y_rep)
    y_obs <- stan_data$d2H_wax

    # Calculate R² and RMSE in standardized units
    ss_res <- sum((y_obs - y_pred_mean)^2)
    ss_tot <- sum((y_obs - mean(y_obs))^2)
    r_squared <- 1 - (ss_res / ss_tot)
    rmse <- sqrt(mean((y_obs - y_pred_mean)^2))

    # Convert RMSE to original scale
    rmse_permil <- rmse * stan_data$scaling_params$d2H_sd

    performance_metrics[[model_name]] <- data.frame(
      model = model_name,
      looic = loo_result$estimates["looic", "Estimate"],
      elpd = loo_result$estimates["elpd_loo", "Estimate"],
      p_loo = loo_result$estimates["p_loo", "Estimate"],
      r_squared = r_squared,
      rmse_std = rmse,
      rmse_permil = rmse_permil,
      n = length(y_obs),
      stringsAsFactors = FALSE
    )

    cat(sprintf("%-20s: LOOIC = %8.1f, R² = %.3f, RMSE = %.1f‰\n",
                model_name, performance_metrics[[model_name]]$looic,
                r_squared, rmse_permil))
  } else {
    cat(sprintf("%-20s: Not found\n", model_name))
  }
}

performance_df <- bind_rows(performance_metrics)

# Sort by LOOIC
performance_df <- performance_df %>% arrange(looic)

cat("\nRanked by LOOIC (lower is better):\n")
print(performance_df %>% select(model, looic, r_squared, rmse_permil))

#───────────────────────────────────────────────────────────────────────────────
# 2. EXTRACT COEFFICIENTS FOR FULL MODELS
#───────────────────────────────────────────────────────────────────────────────

cat("\n\n2. COEFFICIENT EXTRACTION FOR FULL MODELS\n")
cat(strrep("-", 60), "\n")

coefficient_list <- list()

for (model_name in full_models) {
  cat("\n", model_name, ":\n")

  fit_file <- paste0("model_output/", model_name, "/fit.rds")
  config_file <- paste0("prepared_data/config_", model_name, ".rds")

  if (!file.exists(fit_file)) next

  fit <- readRDS(fit_file)
  config <- readRDS(config_file)

  # Parameters to extract
  params_to_extract <- c("beta_oipc", "beta_c4", "beta_tree", "beta_shrub",
                        "beta_grass", "beta_precip_amount")

  # Add interaction terms if applicable
  if (grepl("interact", model_name)) {
    params_to_extract <- c(params_to_extract, "beta_oipc_c4", "beta_precip_c4",
                          "beta_oipc_tree", "beta_oipc_shrub", "beta_oipc_grass",
                          "beta_precip_tree", "beta_precip_shrub", "beta_precip_grass")
  }

  for (param in params_to_extract) {
    if (param %in% fit$metadata()$variables) {
      draws <- as.vector(fit$draws(param, format = "matrix"))

      coef_stats <- data.frame(
        model = model_name,
        parameter = param,
        mean = mean(draws),
        median = median(draws),
        sd = sd(draws),
        lower_95 = quantile(draws, 0.025),
        upper_95 = quantile(draws, 0.975),
        lower_50 = quantile(draws, 0.25),
        upper_50 = quantile(draws, 0.75),
        significant = (quantile(draws, 0.025) > 0) | (quantile(draws, 0.975) < 0),
        stringsAsFactors = FALSE
      )

      coefficient_list[[paste(model_name, param, sep = "_")]] <- coef_stats

      sig_marker <- if(coef_stats$significant) "*" else " "
      cat(sprintf("  %-20s: %7.3f [%7.3f, %7.3f]%s\n",
                  param, coef_stats$mean, coef_stats$lower_95,
                  coef_stats$upper_95, sig_marker))
    }
  }

  # Extract elevation effects (B-spline coefficients)
  if ("beta_elevation" %in% fit$metadata()$variables) {
    elev_draws <- fit$draws("beta_elevation", format = "matrix")
    # Get mean effect across knots
    elev_mean <- mean(elev_draws)
    elev_sd <- sd(elev_draws)

    cat(sprintf("  Elevation (avg)     : %7.3f (SD = %.3f)\n", elev_mean, elev_sd))
  }
}

coefficients_df <- bind_rows(coefficient_list)

#───────────────────────────────────────────────────────────────────────────────
# 3. VARIANCE DECOMPOSITION FOR SPATIAL MODELS
#───────────────────────────────────────────────────────────────────────────────

cat("\n\n3. VARIANCE DECOMPOSITION ANALYSIS\n")
cat(strrep("-", 60), "\n")

spatial_models <- c("full_sp", "full_interact_sp", "baseline_sp",
                   "baseline_env_sp", "baseline_veg_sp")

variance_components <- list()

for (model_name in spatial_models) {
  fit_file <- paste0("model_output/", model_name, "/fit.rds")

  if (file.exists(fit_file)) {
    fit <- readRDS(fit_file)

    # Extract variance components
    if ("sigma_intercept_spatial" %in% fit$metadata()$variables) {
      sigma_int <- mean(as.vector(fit$draws("sigma_intercept_spatial", format = "matrix")))
      sigma_slope <- mean(as.vector(fit$draws("sigma_slope_spatial", format = "matrix")))
      sigma_obs <- mean(as.vector(fit$draws("sigma", format = "matrix")))

      # Calculate proportions
      total_var <- sigma_int^2 + sigma_slope^2 + sigma_obs^2

      variance_components[[model_name]] <- data.frame(
        model = model_name,
        sigma_intercept = sigma_int,
        sigma_slope = sigma_slope,
        sigma_obs = sigma_obs,
        var_intercept_pct = 100 * sigma_int^2 / total_var,
        var_slope_pct = 100 * sigma_slope^2 / total_var,
        var_obs_pct = 100 * sigma_obs^2 / total_var,
        var_spatial_total_pct = 100 * (sigma_int^2 + sigma_slope^2) / total_var,
        stringsAsFactors = FALSE
      )

      cat(sprintf("\n%s:\n", model_name))
      cat(sprintf("  Spatial variance (total): %.1f%%\n",
                  variance_components[[model_name]]$var_spatial_total_pct))
      cat(sprintf("    - Intercept: %.1f%% (σ = %.2f‰)\n",
                  variance_components[[model_name]]$var_intercept_pct, sigma_int))
      cat(sprintf("    - Slope:     %.1f%% (σ = %.2f)\n",
                  variance_components[[model_name]]$var_slope_pct, sigma_slope))
      cat(sprintf("  Observation variance:     %.1f%% (σ = %.2f‰)\n",
                  variance_components[[model_name]]$var_obs_pct, sigma_obs))
    }
  }
}

variance_df <- bind_rows(variance_components)

# Calculate variance reduction
if ("baseline_sp" %in% names(variance_components)) {
  baseline_spatial_var <- variance_components[["baseline_sp"]]$var_spatial_total_pct

  cat("\nVariance reduction from baseline_sp:\n")
  for (model in c("full_sp", "full_interact_sp")) {
    if (model %in% names(variance_components)) {
      model_spatial_var <- variance_components[[model]]$var_spatial_total_pct
      reduction <- baseline_spatial_var - model_spatial_var
      cat(sprintf("  %s: %.1f%% reduction\n", model, reduction))
    }
  }
}

#───────────────────────────────────────────────────────────────────────────────
# 4. PARAMETER CORRELATION ANALYSIS
#───────────────────────────────────────────────────────────────────────────────

cat("\n\n4. PREDICTOR CORRELATION ANALYSIS\n")
cat(strrep("-", 60), "\n")

# Load a representative dataset
stan_data <- readRDS("prepared_data/stan_data_full.rds")

# Create predictor matrix
predictors <- data.frame(
  oipc = stan_data$oipc,
  c4 = stan_data$c4_mean,
  elevation = stan_data$elevation_z
)

# Add precipitation if available
if (!is.null(stan_data$precip_amount)) {
  predictors$precip = stan_data$precip_amount
}

# Add PFT if available
if (!is.null(stan_data$pft_tree)) {
  predictors$tree = stan_data$pft_tree
  predictors$shrub = stan_data$pft_shrub
  predictors$grass = stan_data$pft_grass
}

# Calculate correlation matrix
cor_matrix <- cor(predictors, use = "pairwise.complete.obs")

cat("\nCorrelation matrix:\n")
print(round(cor_matrix, 3))

# Find high correlations
high_cor <- which(abs(cor_matrix) > 0.6 & cor_matrix != 1, arr.ind = TRUE)
if (nrow(high_cor) > 0) {
  cat("\nHigh correlations (|r| > 0.6):\n")
  for (i in 1:nrow(high_cor)) {
    if (high_cor[i, 1] < high_cor[i, 2]) {  # Avoid duplicates
      cat(sprintf("  %s vs %s: r = %.3f\n",
                  rownames(cor_matrix)[high_cor[i, 1]],
                  colnames(cor_matrix)[high_cor[i, 2]],
                  cor_matrix[high_cor[i, 1], high_cor[i, 2]]))
    }
  }
} else {
  cat("\nNo high correlations (|r| > 0.6) found between predictors\n")
}

#───────────────────────────────────────────────────────────────────────────────
# 5. COEFFICIENT COMPARISON ACROSS MODEL COMPLEXITY
#───────────────────────────────────────────────────────────────────────────────

cat("\n\n5. COEFFICIENT CHANGES WITH MODEL COMPLEXITY\n")
cat(strrep("-", 60), "\n")

# Compare beta_oipc across models
cat("\nbeta_oipc (OIPC slope) changes:\n")
oipc_coefs <- coefficients_df %>%
  filter(parameter == "beta_oipc") %>%
  select(model, mean, lower_95, upper_95)

if (nrow(oipc_coefs) > 0) {
  print(oipc_coefs)
}

# Compare vegetation effects
cat("\nVegetation effect changes:\n")
veg_params <- c("beta_c4", "beta_tree", "beta_shrub", "beta_grass")

for (param in veg_params) {
  param_data <- coefficients_df %>% filter(parameter == param)
  if (nrow(param_data) > 0) {
    cat(sprintf("\n%s:\n", param))
    for (i in 1:nrow(param_data)) {
      sig_marker <- if(param_data$significant[i]) "*" else " "
      cat(sprintf("  %-20s: %7.3f [%7.3f, %7.3f]%s\n",
                  param_data$model[i], param_data$mean[i],
                  param_data$lower_95[i], param_data$upper_95[i], sig_marker))
    }
  }
}

#───────────────────────────────────────────────────────────────────────────────
# 6. RELATIVE IMPORTANCE
#───────────────────────────────────────────────────────────────────────────────

cat("\n\n6. RELATIVE IMPORTANCE OF PREDICTORS\n")
cat(strrep("-", 60), "\n")

for (model_name in c("full", "full_sp")) {
  model_coefs <- coefficients_df %>%
    filter(model == model_name, !grepl("interact|x", parameter))

  if (nrow(model_coefs) > 0) {
    # Sort by absolute effect size
    model_coefs <- model_coefs %>%
      mutate(abs_mean = abs(mean)) %>%
      arrange(desc(abs_mean))

    cat(sprintf("\n%s - Ranked by absolute effect size:\n", model_name))
    for (i in 1:min(nrow(model_coefs), 6)) {
      sig_marker <- if(model_coefs$significant[i]) "*" else " "
      cat(sprintf("  %2d. %-20s: |β| = %.3f%s\n",
                  i, model_coefs$parameter[i], model_coefs$abs_mean[i], sig_marker))
    }
  }
}

# Check which predictors lose significance
cat("\n\nPredictors that lose significance in full models:\n")
baseline_sig <- coefficients_df %>%
  filter(model %in% c("baseline", "baseline_sp"), significant) %>%
  pull(parameter) %>% unique()

full_nonsig <- coefficients_df %>%
  filter(model %in% c("full", "full_sp"), !significant) %>%
  pull(parameter) %>% unique()

lost_sig <- intersect(baseline_sig, full_nonsig)
if (length(lost_sig) > 0) {
  cat("  ", paste(lost_sig, collapse = ", "), "\n")
} else {
  cat("  None\n")
}

#───────────────────────────────────────────────────────────────────────────────
# 7. SAVE OUTPUTS
#───────────────────────────────────────────────────────────────────────────────

# Save coefficient comparison table
write.csv(coefficients_df,
          file.path(output_dir, "full_model_coefficients.csv"),
          row.names = FALSE)

# Save variance decomposition
write.csv(variance_df,
          file.path(output_dir, "full_model_variance.csv"),
          row.names = FALSE)

# Save performance metrics
write.csv(performance_df,
          file.path(output_dir, "full_model_performance.csv"),
          row.names = FALSE)

# Save correlation matrix
write.csv(cor_matrix,
          file.path(output_dir, "predictor_correlations.csv"))

cat("\n\nAll outputs saved to:", output_dir, "\n")

#───────────────────────────────────────────────────────────────────────────────
# SUMMARY REPORT
#───────────────────────────────────────────────────────────────────────────────

cat("\n", strrep("=", 60), "\n")
cat("SUMMARY FOR SECTION 3.4\n")
cat(strrep("=", 60), "\n")

# Best performing model
best_model <- performance_df %>% slice(1)
cat("\nBest performing model:", best_model$model, "\n")
cat(sprintf("  LOOIC: %.1f, R²: %.3f, RMSE: %.1f‰\n",
            best_model$looic, best_model$r_squared, best_model$rmse_permil))

# Key findings
cat("\nKey findings:\n")
cat("1. Significant predictors in full_sp:\n")
full_sp_sig <- coefficients_df %>%
  filter(model == "full_sp", significant, !grepl("interact", parameter))
if (nrow(full_sp_sig) > 0) {
  for (i in 1:nrow(full_sp_sig)) {
    cat(sprintf("   - %s: %.3f [%.3f, %.3f]\n",
                full_sp_sig$parameter[i], full_sp_sig$mean[i],
                full_sp_sig$lower_95[i], full_sp_sig$upper_95[i]))
  }
}

cat("\n2. Variance explained by spatial component in full_sp:\n")
if ("full_sp" %in% names(variance_components)) {
  cat(sprintf("   %.1f%% (mostly intercept variation)\n",
              variance_components[["full_sp"]]$var_spatial_total_pct))
}

cat("\n3. Interaction effects:\n")
interact_params <- coefficients_df %>%
  filter(grepl("interact", model), grepl("x|_c4|_tree|_shrub|_grass", parameter),
         significant)
if (nrow(interact_params) > 0) {
  cat("   Significant interactions found:\n")
  for (i in 1:nrow(interact_params)) {
    cat(sprintf("   - %s in %s\n", interact_params$parameter[i], interact_params$model[i]))
  }
} else {
  cat("   No significant interaction effects\n")
}

cat("\nAnalysis complete!\n")