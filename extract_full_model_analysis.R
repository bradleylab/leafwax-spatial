# Comprehensive analysis of full models for Section 3.4
library(tidyverse)
library(posterior)

source("scripts/posterior_helpers.R")

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

# Models to compare — full 14-model roster per manuscript/TABLES.md (Table 1).
# Order matches config.yaml for stable table row ordering.
all_models <- c(
  "baseline", "baseline_sp",
  "baseline_env", "baseline_env_sp",
  "baseline_veg", "baseline_veg_sp",
  "full", "full_sp",
  "full_interact", "full_interact_sp",
  "elevation_only_sp", "elevation_c4_sp",
  "c4_only_sp", "elevation_c4_interact_sp"
)

# Full-model subset used for per-coefficient extraction (section 2).
# Section 1 runs over all 14 models; section 2 only over the four "full" ones.
full_models <- c("full", "full_sp", "full_interact", "full_interact_sp")

performance_metrics <- list()

for (model_name in all_models) {
  # Helpers raise a clear error if any piece is missing; wrap in tryCatch
  # so a missing model doesn't abort the whole loop.
  bundle <- tryCatch({
    list(draws     = load_draws(model_name),
         loo       = load_loo(model_name),
         stan_data = load_stan_data(model_name))
  }, error = function(e) {
    cat(sprintf("%-20s: %s\n", model_name, conditionMessage(e)))
    NULL
  })
  if (is.null(bundle)) next

  loo_result <- bundle$loo
  stan_data  <- bundle$stan_data

  # Posterior predictive means: y_rep is a draws matrix with one column per
  # observation; colMeans gives the per-obs posterior mean.
  y_rep       <- as_draws_matrix(subset_draws(bundle$draws, variable = "d2H_rep"))
  y_pred_mean <- colMeans(y_rep)
  y_obs       <- stan_data$d2H_wax

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

  draws <- tryCatch(load_draws(model_name),
                    error = function(e) { cat("  skip:", conditionMessage(e), "\n"); NULL })
  if (is.null(draws)) next
  config <- load_config(model_name)
  vars_present <- variables(draws)

  # Parameters to extract. Names must match 4d_leaf_wax_spatial_model.stan exactly.
  # Stan exports `beta_precip` (not `beta_precip_amount`) and `beta_oipc_x_*`
  # (not `beta_oipc_*`). There are no `beta_precip_*` interaction terms in the
  # model; those extraction requests were stale and have been removed.
  params_to_extract <- c("beta_oipc", "beta_c4", "beta_tree", "beta_shrub",
                        "beta_grass", "beta_precip")

  # Add OIPC interaction terms if applicable (no precip interactions in Stan)
  if (grepl("interact", model_name)) {
    params_to_extract <- c(params_to_extract,
                          "beta_oipc_x_c4",
                          "beta_oipc_x_tree", "beta_oipc_x_shrub", "beta_oipc_x_grass")
  }

  for (param in params_to_extract) {
    if (!(param %in% vars_present)) next
    d <- as.numeric(as_draws_matrix(subset_draws(draws, variable = param)))

    coef_stats <- data.frame(
      model = model_name,
      parameter = param,
      mean = mean(d),
      median = median(d),
      sd = sd(d),
      lower_95 = quantile(d, 0.025),
      upper_95 = quantile(d, 0.975),
      lower_50 = quantile(d, 0.25),
      upper_50 = quantile(d, 0.75),
      significant = (quantile(d, 0.025) > 0) | (quantile(d, 0.975) < 0),
      stringsAsFactors = FALSE
    )

    coefficient_list[[paste(model_name, param, sep = "_")]] <- coef_stats

    sig_marker <- if(coef_stats$significant) "*" else " "
    cat(sprintf("  %-20s: %7.3f [%7.3f, %7.3f]%s\n",
                param, coef_stats$mean, coef_stats$lower_95,
                coef_stats$upper_95, sig_marker))
  }

  # Elevation enters via B-spline coefficients `beta_elev_bspline` (a vector),
  # not a scalar `beta_elevation`. Summarize mean magnitude across coefficients.
  # beta_elev_bspline is NOT in draws_to_save (Phase 5 W1 widened list); pull
  # its summary from diagnostics.rds where all parameters are summarized.
  summ <- load_summaries(model_name)
  elev_rows <- summ[grepl("^beta_elev_bspline", summ$variable), ]
  if (nrow(elev_rows) > 0) {
    # Use posterior-mean coefficient values; across-coefficient stats only.
    elev_mean_abs <- mean(abs(elev_rows$mean))
    elev_sd       <- sd(elev_rows$mean)
    cat(sprintf("  Elevation (|β| avg) : %7.3f (SD = %.3f)\n", elev_mean_abs, elev_sd))
  }
}

coefficients_df <- bind_rows(coefficient_list)

#───────────────────────────────────────────────────────────────────────────────
# 3. VARIANCE DECOMPOSITION FOR SPATIAL MODELS
#───────────────────────────────────────────────────────────────────────────────

cat("\n\n3. VARIANCE DECOMPOSITION ANALYSIS\n")
cat(strrep("-", 60), "\n")

# All 9 spatial models per manuscript/TABLES.md (Table 3 variance decomposition).
spatial_models <- c(
  "baseline_sp",
  "baseline_env_sp",
  "baseline_veg_sp",
  "full_sp",
  "full_interact_sp",
  "elevation_only_sp",
  "elevation_c4_sp",
  "c4_only_sp",
  "elevation_c4_interact_sp"
)

variance_components <- list()

for (model_name in spatial_models) {
  summ <- tryCatch(load_summaries(model_name),
                   error = function(e) { cat(sprintf("%-20s skip: %s\n", model_name, conditionMessage(e))); NULL })
  if (is.null(summ)) next

  # Use Stan's own variance decomposition from generated quantities.
  # 4d_leaf_wax_spatial_model.stan:504-511 defines var_spatial_intercept,
  # var_spatial_slope, var_residual on a consistent scale, and
  # prop_variance_spatial / prop_variance_residual as the final proportions.
  # The prior ad-hoc calculation mixed sigma_intercept_spatial (original ‰)
  # with sigma (standardized) and inflated the spatial share to 100%.
  get_mean <- function(var) {
    row <- summ[summ$variable == var, , drop = FALSE]
    if (nrow(row) != 1) return(NA_real_)
    row$mean
  }

  var_int      <- get_mean("var_spatial_intercept")
  var_slope    <- get_mean("var_spatial_slope")
  var_spatial  <- get_mean("var_total_spatial")
  var_resid    <- get_mean("var_residual")
  var_total    <- get_mean("var_total")
  prop_spatial <- get_mean("prop_variance_spatial")
  prop_resid   <- get_mean("prop_variance_residual")

  # Additional reporting: posterior-mean residual SD in original ‰.
  sigma_resid_orig <- get_mean("sigma_residual_original")

  variance_components[[model_name]] <- data.frame(
    model = model_name,
    var_spatial_intercept = var_int,
    var_spatial_slope     = var_slope,
    var_total_spatial     = var_spatial,
    var_residual          = var_resid,
    var_total             = var_total,
    var_intercept_pct     = 100 * var_int   / var_total,
    var_slope_pct         = 100 * var_slope / var_total,
    var_spatial_total_pct = 100 * prop_spatial,
    var_obs_pct           = 100 * prop_resid,
    sigma_residual_permil = sigma_resid_orig,
    stringsAsFactors = FALSE
  )

  cat(sprintf("\n%s:\n", model_name))
  cat(sprintf("  Spatial variance (total): %.1f%%\n",
              variance_components[[model_name]]$var_spatial_total_pct))
  cat(sprintf("    - Intercept: %.1f%%\n",
              variance_components[[model_name]]$var_intercept_pct))
  cat(sprintf("    - Slope:     %.1f%%\n",
              variance_components[[model_name]]$var_slope_pct))
  cat(sprintf("  Residual variance:        %.1f%% (σ = %.2f ‰)\n",
              variance_components[[model_name]]$var_obs_pct, sigma_resid_orig))
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

# Build a per-obs predictor matrix for multicollinearity analysis.
# Environmental covariates in stan_data are 1124×9 matrices (obs × spatial
# scale); use the narrowest-scale column (index 1) as the representative
# point value. Scalar per-obs variables come from the sediment frame.
sd_full <- load_stan_data("full")
sed     <- load_sediment()

.col1 <- function(m) if (is.matrix(m)) m[, 1] else as.numeric(m)

predictors <- data.frame(
  oipc      = .col1(sd_full$oipc_values),
  elevation = .col1(sd_full$elevation_values)
)
# c4 in stan_data: try weighted first, fall back to sediment scalar
if (!is.null(sd_full$c4_values) && is.matrix(sd_full$c4_values)) {
  predictors$c4 <- .col1(sd_full$c4_values)
} else {
  predictors$c4 <- sed$c4_mean_filled
}
if (!is.null(sd_full$precip_values) && is.matrix(sd_full$precip_values)) {
  predictors$precip <- .col1(sd_full$precip_values)
} else if (!is.null(sed$annual_precip)) {
  predictors$precip <- sed$annual_precip
}
if (!is.null(sd_full$pft_tree) && is.matrix(sd_full$pft_tree)) {
  predictors$tree  <- .col1(sd_full$pft_tree)
  predictors$shrub <- .col1(sd_full$pft_shrub)
  predictors$grass <- .col1(sd_full$pft_grass)
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