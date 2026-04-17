# Extract vegetation coefficients from all models
library(tidyverse)
library(posterior)

source("scripts/posterior_helpers.R")

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

  bundle <- tryCatch(
    list(draws = load_draws(model_name),
         summ  = load_summaries(model_name),
         config = load_config(model_name)),
    error = function(e) { cat("  skip:", conditionMessage(e), "\n"); NULL })
  if (is.null(bundle)) next
  draws <- bundle$draws; summ <- bundle$summ; config <- bundle$config
  vars_present <- variables(draws)

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
    if (!(param %in% vars_present)) {
      cat("\n  ", param, ": Not found in model\n")
      next
    }
    d <- as.numeric(as_draws_matrix(subset_draws(draws, variable = param)))

    # Calculate statistics
    param_stats <- data.frame(
      model = model_name,
      parameter = param,
      mean = mean(d),
      median = median(d),
      sd = sd(d),
      lower_95 = quantile(d, 0.025),
      upper_95 = quantile(d, 0.975),
      lower_50 = quantile(d, 0.25),
      upper_50 = quantile(d, 0.75),
      n_eff = NA_real_,
      rhat  = NA_real_,
      stringsAsFactors = FALSE
    )
    param_stats$significant <- (param_stats$lower_95 > 0) | (param_stats$upper_95 < 0)

    # Convergence diagnostics from diagnostics.rds summary
    s_row <- summ[summ$variable == param, , drop = FALSE]
    if (nrow(s_row) == 1) {
      param_stats$n_eff <- s_row$ess_bulk
      param_stats$rhat  <- s_row$rhat
    }

    all_coefficients[[paste(model_name, param, sep = "_")]] <- param_stats

    cat("\n  ", param, ":\n")
    cat("    Mean (SD):", sprintf("%.3f (%.3f)", param_stats$mean, param_stats$sd), "\n")
    cat("    Median:", sprintf("%.3f", param_stats$median), "\n")
    cat("    95% CI: [", sprintf("%.3f, %.3f", param_stats$lower_95, param_stats$upper_95), "]\n")
    cat("    Significant:", param_stats$significant, "\n")
    if (!is.na(param_stats$rhat)) cat("    Rhat:", sprintf("%.3f", param_stats$rhat), "\n")
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

# Variance decomposition uses Stan's own generated quantities
# (prop_variance_spatial, prop_variance_residual; see extract_full_model_analysis.R
# for the rationale — the old sigma_int^2 + sigma^2 calculation mixed units
# and inflated the spatial share to 100%%).
.variance_bundle <- function(model_name) {
  summ <- tryCatch(load_summaries(model_name),
                   error = function(e) { cat("  skip:", conditionMessage(e), "\n"); NULL })
  if (is.null(summ)) return(NULL)
  g <- function(v) {
    row <- summ[summ$variable == v, , drop = FALSE]
    if (nrow(row) != 1) NA_real_ else row$mean
  }
  list(
    sigma_resid_permil = g("sigma_residual_original"),
    prop_spatial       = g("prop_variance_spatial"),
    prop_residual      = g("prop_variance_residual"),
    var_intercept      = g("var_spatial_intercept"),
    var_slope          = g("var_spatial_slope"),
    var_total          = g("var_total")
  )
}

for (model_name in spatial_veg_models) {
  vc <- .variance_bundle(model_name)
  if (is.null(vc) || is.na(vc$prop_spatial)) next
  cat("\n", model_name, ":\n")
  cat(sprintf("  Spatial (GP) variance: %.1f%%\n", 100 * vc$prop_spatial))
  cat(sprintf("    Intercept component: %.1f%%\n",
              100 * vc$var_intercept / vc$var_total))
  cat(sprintf("    Slope component:     %.1f%%\n",
              100 * vc$var_slope / vc$var_total))
  cat(sprintf("  Residual variance:     %.1f%% (σ = %.2f ‰)\n",
              100 * vc$prop_residual, vc$sigma_resid_permil))
  variance_results[[model_name]] <- list(
    sigma_intercept = NA_real_,           # kept for downstream field names
    sigma_slope     = NA_real_,
    sigma_obs       = vc$sigma_resid_permil,
    prop_gp         = vc$prop_spatial
  )
}

# Compare to baseline_sp (no vegetation) using the same Stan-based decomposition.
baseline_vc <- .variance_bundle("baseline_sp")
if (!is.null(baseline_vc) && !is.na(baseline_vc$prop_spatial)) {
  gp_prop_base <- baseline_vc$prop_spatial
  cat("\n\nComparison to baseline_sp (no vegetation):\n")
  cat(sprintf("  Baseline GP variance: %.1f%%\n", 100 * gp_prop_base))
  for (model_name in names(variance_results)) {
    reduction <- 100 * (gp_prop_base - variance_results[[model_name]]$prop_gp) / gp_prop_base
    cat(sprintf("  %s reduces GP variance by %.1f%%\n", model_name, reduction))
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