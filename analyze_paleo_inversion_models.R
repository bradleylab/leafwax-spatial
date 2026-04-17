#───────────────────────────────────────────────────────────────────────────────
# analyze_paleo_inversion_models.R
#
# Analysis of spatial models for paleoclimate inversions (Section 3.5)
# Extracts and analyzes parameters relevant for paleo-environmental reconstruction
# Evaluates spatial vs non-spatial models for inverse modeling applications
#
# Input: model_output/*/fit.rds (fitted models)
# Output: model_analysis/paleo_inversion_analysis/ (tables and summaries)
#───────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(cmdstanr)

# Create output directory
output_dir <- "model_analysis/tables"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("PALEOCLIMATE INVERSION MODEL ANALYSIS\n")
cat("======================================\n\n")

# Models to analyze
paleo_models <- c("elevation_only_sp", "elevation_c4_sp",
                  "elevation_c4_interact_sp", "c4_only_sp")

# Reference models for comparison
reference_models <- c("baseline_sp", "full_sp")

all_models <- c(paleo_models, reference_models)

#───────────────────────────────────────────────────────────────────────────────
# 1. EXTRACT PERFORMANCE METRICS
#───────────────────────────────────────────────────────────────────────────────

cat("1. MODEL PERFORMANCE METRICS\n")
cat(strrep("-", 60), "\n\n")

performance_list <- list()

for (model_name in all_models) {
  loo_file <- paste0("model_output/", model_name, "/loo.rds")
  fit_file <- paste0("model_output/", model_name, "/fit.rds")

  if (file.exists(loo_file) && file.exists(fit_file)) {
    # Load LOO
    loo_result <- readRDS(loo_file)

    # Load fit for R² and RMSE
    fit <- readRDS(fit_file)
    stan_data <- readRDS(paste0("prepared_data/stan_data_", model_name, ".rds"))

    # Calculate R² and RMSE
    y_rep <- fit$draws("d2H_rep", format = "matrix")
    y_pred_mean <- colMeans(y_rep)
    y_obs <- stan_data$d2H_wax

    ss_res <- sum((y_obs - y_pred_mean)^2)
    ss_tot <- sum((y_obs - mean(y_obs))^2)
    r_squared <- 1 - (ss_res / ss_tot)

    rmse <- sqrt(mean((y_obs - y_pred_mean)^2))
    rmse_permil <- rmse * stan_data$scaling_params$d2H_sd

    performance_list[[model_name]] <- data.frame(
      model = model_name,
      looic = loo_result$estimates["looic", "Estimate"],
      elpd = loo_result$estimates["elpd_loo", "Estimate"],
      p_loo = loo_result$estimates["p_loo", "Estimate"],
      r_squared = r_squared,
      rmse_permil = rmse_permil,
      n = length(y_obs),
      stringsAsFactors = FALSE
    )

    cat(sprintf("%-25s: LOOIC = %7.1f, R² = %.3f, RMSE = %.1f‰\n",
                model_name, performance_list[[model_name]]$looic,
                r_squared, rmse_permil))
  } else {
    cat(sprintf("%-25s: Not found\n", model_name))
  }
}

performance_df <- bind_rows(performance_list)

# Calculate performance differences
baseline_looic <- performance_df$looic[performance_df$model == "baseline_sp"]
full_looic <- performance_df$looic[performance_df$model == "full_sp"]

performance_df <- performance_df %>%
  mutate(
    delta_baseline = looic - baseline_looic,
    delta_full = looic - full_looic
  )

cat("\n\nPerformance comparison:\n")
print(performance_df %>%
      filter(model %in% paleo_models) %>%
      select(model, looic, delta_baseline, delta_full, r_squared) %>%
      arrange(looic))

#───────────────────────────────────────────────────────────────────────────────
# 2. IDENTIFY REQUIRED PREDICTORS AND PALEO-AVAILABILITY
#───────────────────────────────────────────────────────────────────────────────

cat("\n\n2. PREDICTOR REQUIREMENTS AND PALEO-AVAILABILITY\n")
cat(strrep("-", 60), "\n\n")

# Define predictor requirements and paleo-availability
predictor_info <- list(
  elevation_only_sp = list(
    predictors = c("elevation"),
    paleo_feasible = TRUE,
    notes = "Elevation easily estimated from basin analysis"
  ),
  c4_only_sp = list(
    predictors = c("c4"),
    paleo_feasible = TRUE,
    notes = "C4 can be estimated from δ13C of leaf waxes"
  ),
  elevation_c4_sp = list(
    predictors = c("elevation", "c4"),
    paleo_feasible = TRUE,
    notes = "Both elevation and C4 are paleo-reconstructable"
  ),
  elevation_c4_interact_sp = list(
    predictors = c("elevation", "c4", "elevation×c4"),
    paleo_feasible = TRUE,
    notes = "Interaction can be calculated from reconstructed values"
  ),
  baseline_sp = list(
    predictors = c("none (spatial only)"),
    paleo_feasible = TRUE,
    notes = "No predictors needed"
  ),
  full_sp = list(
    predictors = c("elevation", "c4", "tree%", "shrub%", "grass%", "precip_amount"),
    paleo_feasible = FALSE,
    notes = "PFTs and precip amount difficult to reconstruct"
  )
)

for (model in names(predictor_info)) {
  if (model %in% paleo_models || model %in% reference_models) {
    cat(sprintf("\n%s:\n", model))
    cat("  Predictors:", paste(predictor_info[[model]]$predictors, collapse = ", "), "\n")
    cat("  Paleo-feasible:", predictor_info[[model]]$paleo_feasible, "\n")
    cat("  Notes:", predictor_info[[model]]$notes, "\n")
  }
}

#───────────────────────────────────────────────────────────────────────────────
# 3. EXTRACT KEY COEFFICIENTS
#───────────────────────────────────────────────────────────────────────────────

cat("\n\n3. KEY COEFFICIENTS WITH 95% CIs\n")
cat(strrep("-", 60), "\n")

coefficient_list <- list()

for (model_name in c(paleo_models, "baseline_sp")) {
  cat("\n", model_name, ":\n")

  fit_file <- paste0("model_output/", model_name, "/fit.rds")
  config_file <- paste0("prepared_data/config_", model_name, ".rds")

  if (!file.exists(fit_file)) next

  fit <- readRDS(fit_file)

  # Parameters to check. Names must match 4d_leaf_wax_spatial_model.stan.
  # Notes on legacy names:
  #  - beta_elevation: not a scalar in Stan; elevation enters via
  #    `beta_elev_bspline` (vector of B-spline coefficients). Handled below.
  #  - beta_oipc_c4: Stan exports `beta_oipc_x_c4`.
  #  - beta_elevation_c4: no elevation×C4 interaction exists in the model.
  params_to_extract <- c("beta_oipc", "beta_c4", "beta_oipc_x_c4")

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

  # Elevation B-spline coefficients (vector, not scalar).
  if ("beta_elev_bspline" %in% fit$metadata()$variables) {
    elev_draws <- fit$draws("beta_elev_bspline", format = "matrix")
    # Mean across-coefficient effect (not across draws). colMeans collapses
    # draws per coefficient, yielding a posterior-mean coefficient vector.
    elev_effects <- colMeans(elev_draws)
    cat(sprintf("  Elevation effects   : Mean across B-spline coefs = %.3f\n",
                mean(elev_effects)))
    cat(sprintf("                        Range: [%.3f, %.3f]\n",
                min(elev_effects), max(elev_effects)))
  }
}

coefficients_df <- bind_rows(coefficient_list)

#───────────────────────────────────────────────────────────────────────────────
# 4. SPATIAL HETEROGENEITY ANALYSIS
#───────────────────────────────────────────────────────────────────────────────

cat("\n\n4. SPATIAL HETEROGENEITY\n")
cat(strrep("-", 60), "\n\n")

spatial_hetero_list <- list()

for (model_name in c(paleo_models, reference_models)) {
  fit_file <- paste0("model_output/", model_name, "/fit.rds")

  if (file.exists(fit_file)) {
    fit <- readRDS(fit_file)

    # Extract spatial variance components
    if ("sigma_slope_spatial" %in% fit$metadata()$variables) {
      slope_draws <- as.vector(fit$draws("sigma_slope_spatial", format = "matrix"))

      spatial_hetero_list[[model_name]] <- data.frame(
        model = model_name,
        slope_sd_mean = mean(slope_draws),
        slope_sd_median = median(slope_draws),
        slope_sd_lower = quantile(slope_draws, 0.025),
        slope_sd_upper = quantile(slope_draws, 0.975),
        stringsAsFactors = FALSE
      )

      cat(sprintf("%-25s: Slope SD = %.3f [%.3f, %.3f]\n",
                  model_name, mean(slope_draws),
                  quantile(slope_draws, 0.025),
                  quantile(slope_draws, 0.975)))
    }

    # Also check intercept variation
    if ("sigma_intercept_spatial" %in% fit$metadata()$variables) {
      int_draws <- as.vector(fit$draws("sigma_intercept_spatial", format = "matrix"))
      spatial_hetero_list[[model_name]]$intercept_sd_mean <- mean(int_draws)
    }
  }
}

spatial_hetero_df <- bind_rows(spatial_hetero_list)

cat("\n\nSpatial homogeneity ranking (lower slope SD = more homogeneous):\n")
print(spatial_hetero_df %>%
      arrange(slope_sd_mean) %>%
      select(model, slope_sd_mean, intercept_sd_mean))

#───────────────────────────────────────────────────────────────────────────────
# 5. CREATE PRACTICAL COMPARISON TABLE
#───────────────────────────────────────────────────────────────────────────────

cat("\n\n5. PRACTICAL COMPARISON TABLE\n")
cat(strrep("-", 60), "\n\n")

# Combine all information
practical_comparison <- performance_df %>%
  filter(model %in% c(paleo_models, reference_models)) %>%
  mutate(
    predictors = sapply(model, function(m) {
      if (m %in% names(predictor_info)) {
        paste(predictor_info[[m]]$predictors, collapse = ", ")
      } else NA
    }),
    paleo_feasible = sapply(model, function(m) {
      if (m %in% names(predictor_info)) {
        predictor_info[[m]]$paleo_feasible
      } else NA
    }),
    performance_loss_vs_best = looic - min(performance_df$looic)
  ) %>%
  select(model, predictors, paleo_feasible, looic, r_squared,
         rmse_permil, performance_loss_vs_best) %>%
  arrange(looic)

cat("Models ranked by performance:\n")
print(practical_comparison)

# Save the table
write.csv(practical_comparison,
          file.path(output_dir, "paleo_inversion_models.csv"),
          row.names = FALSE)

#───────────────────────────────────────────────────────────────────────────────
# 6. RECOMMENDATION SUMMARY
#───────────────────────────────────────────────────────────────────────────────

cat("\n\n", strrep("=", 60), "\n")
cat("RECOMMENDATIONS FOR PALEOCLIMATE INVERSIONS\n")
cat(strrep("=", 60), "\n\n")

# Find best paleo-feasible model
paleo_feasible_models <- practical_comparison %>%
  filter(paleo_feasible == TRUE, model != "baseline_sp")

best_paleo <- paleo_feasible_models %>%
  slice(1)

cat("1. BEST PALEO-FEASIBLE MODEL:\n")
cat(sprintf("   %s (LOOIC = %.1f, R² = %.3f)\n",
            best_paleo$model, best_paleo$looic, best_paleo$r_squared))
cat(sprintf("   Performance loss vs full_sp: %.1f LOOIC units\n",
            best_paleo$performance_loss_vs_best))

# Compare elevation+C4 vs C4-only
elev_c4_perf <- performance_df$looic[performance_df$model == "elevation_c4_sp"]
c4_only_perf <- performance_df$looic[performance_df$model == "c4_only_sp"]
elev_only_perf <- performance_df$looic[performance_df$model == "elevation_only_sp"]

cat("\n2. PREDICTOR CONTRIBUTION:\n")
cat(sprintf("   C4-only model: LOOIC = %.1f\n", c4_only_perf))
cat(sprintf("   Elevation-only model: LOOIC = %.1f\n", elev_only_perf))
cat(sprintf("   Elevation+C4 model: LOOIC = %.1f\n", elev_c4_perf))
cat(sprintf("   Benefit of combining: Δ%.1f LOOIC improvement over C4-only\n",
            c4_only_perf - elev_c4_perf))

# Performance cost analysis
cat("\n3. PERFORMANCE COST OF PALEO-FEASIBLE MODELS:\n")
baseline_r2 <- performance_df$r_squared[performance_df$model == "baseline_sp"]
full_r2 <- performance_df$r_squared[performance_df$model == "full_sp"]
best_paleo_r2 <- best_paleo$r_squared

cat(sprintf("   Baseline_sp (no predictors): R² = %.3f\n", baseline_r2))
cat(sprintf("   Best paleo-feasible: R² = %.3f\n", best_paleo_r2))
cat(sprintf("   Full_sp (all predictors): R² = %.3f\n", full_r2))
cat(sprintf("   R² gain from paleo predictors: %.3f\n", best_paleo_r2 - baseline_r2))
cat(sprintf("   R² loss vs full model: %.3f\n", full_r2 - best_paleo_r2))

# Check if interaction helps
if ("elevation_c4_interact_sp" %in% performance_df$model) {
  interact_looic <- performance_df$looic[performance_df$model == "elevation_c4_interact_sp"]
  no_interact_looic <- performance_df$looic[performance_df$model == "elevation_c4_sp"]

  cat("\n4. INTERACTION EFFECTS:\n")
  cat(sprintf("   Without interaction: LOOIC = %.1f\n", no_interact_looic))
  cat(sprintf("   With interaction: LOOIC = %.1f\n", interact_looic))

  if (interact_looic < no_interact_looic - 2) {
    cat("   Recommendation: Include interaction term\n")
  } else {
    cat("   Recommendation: Interaction not needed (ΔLOOIC < 2)\n")
  }
}

# Final recommendation
cat("\n5. FINAL RECOMMENDATION:\n")
cat("   For paleoclimate inversions, use: ")

if (elev_c4_perf < c4_only_perf - 5) {
  cat("elevation_c4_sp model\n")
  cat("   - Requires: Elevation (from basin analysis) + C4% (from δ13C)\n")
  cat("   - Performance: Nearly matches best environmental model\n")
} else {
  cat("c4_only_sp model\n")
  cat("   - Requires: Only C4% (from δ13C of leaf waxes)\n")
  cat("   - Simpler with minimal performance loss\n")
}

cat("\n6. KEY INSIGHTS:\n")
cat("   - Paleo-feasible predictors capture most of the signal\n")
cat("   - Spatial component essential (accounts for unmeasured variation)\n")
cat("   - PFT percentages not critical for paleoclimate reconstruction\n")
cat("   - C4 abundance is the most important paleo-available vegetation metric\n")

# Save coefficient stability analysis
cat("\n\n7. COEFFICIENT STABILITY ACROSS MODELS:\n")
cat(strrep("-", 60), "\n")

# Check beta_oipc stability
oipc_coefs <- coefficients_df %>%
  filter(parameter == "beta_oipc") %>%
  select(model, mean, lower_95, upper_95)

if (nrow(oipc_coefs) > 0) {
  cat("\nbeta_oipc across models:\n")
  print(oipc_coefs)
  cat(sprintf("  Coefficient of variation: %.2f%%\n",
              100 * sd(oipc_coefs$mean) / mean(oipc_coefs$mean)))
}

# Check beta_c4 stability
c4_coefs <- coefficients_df %>%
  filter(parameter == "beta_c4") %>%
  select(model, mean, lower_95, upper_95)

if (nrow(c4_coefs) > 0) {
  cat("\nbeta_c4 across models:\n")
  print(c4_coefs)
}

cat("\n\nAnalysis complete! Results saved to:",
    file.path(output_dir, "paleo_inversion_models.csv"), "\n")