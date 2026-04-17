# Extract spatial variance components for vegetation models
library(tidyverse)
library(posterior)

source("scripts/posterior_helpers.R")

cat("SPATIAL VARIANCE COMPONENTS FOR VEGETATION MODELS\n")
cat("==================================================\n\n")

# Models to check
spatial_veg_models <- c("baseline_veg_sp", "full_sp", "full_interact_sp")

# Stan-correct variance decomposition. Pulls posterior means of:
#   prop_variance_spatial  = var_total_spatial / var_total
#   prop_variance_residual = var_residual / var_total
#   var_spatial_intercept, var_spatial_slope, var_total
#   sigma_residual_original (‰)
# from diagnostics.rds$all_params_summary.
# See 4d_leaf_wax_spatial_model.stan:504-511 for the generated-quantity
# definitions. The prior sigma_int^2 + sigma_slope^2 + sigma^2 scheme
# mixed original-scale GP SDs with a standardized residual σ.
.variance_bundle <- function(model_name) {
  summ <- tryCatch(load_summaries(model_name),
                   error = function(e) { cat("  skip:", conditionMessage(e), "\n"); NULL })
  if (is.null(summ)) return(NULL)
  g <- function(v) {
    row <- summ[summ$variable == v, , drop = FALSE]
    if (nrow(row) != 1) NA_real_ else row$mean
  }
  list(
    model             = model_name,
    sigma_resid_permil = g("sigma_residual_original"),
    prop_spatial      = g("prop_variance_spatial"),
    prop_residual     = g("prop_variance_residual"),
    var_intercept     = g("var_spatial_intercept"),
    var_slope         = g("var_spatial_slope"),
    var_total_spatial = g("var_total_spatial"),
    var_residual      = g("var_residual"),
    var_total         = g("var_total")
  )
}

variance_results <- list()

for (model_name in spatial_veg_models) {
  cat("\nProcessing:", model_name, "\n")
  cat(strrep("-", 40), "\n")

  vc <- .variance_bundle(model_name)
  if (is.null(vc) || is.na(vc$prop_spatial)) {
    cat("  Spatial variance parameters not found\n")
    next
  }

  cat(sprintf("  Spatial (GP) variance: %.1f%%\n", 100 * vc$prop_spatial))
  cat(sprintf("    Intercept component: %.1f%%\n", 100 * vc$var_intercept / vc$var_total))
  cat(sprintf("    Slope component:     %.1f%%\n", 100 * vc$var_slope     / vc$var_total))
  cat(sprintf("  Residual variance:     %.1f%% (σ = %.2f ‰)\n",
              100 * vc$prop_residual, vc$sigma_resid_permil))

  variance_results[[model_name]] <- data.frame(
    model              = model_name,
    sigma_resid_permil = vc$sigma_resid_permil,
    var_intercept      = vc$var_intercept,
    var_slope          = vc$var_slope,
    var_total_spatial  = vc$var_total_spatial,
    var_residual       = vc$var_residual,
    var_total          = vc$var_total,
    prop_intercept     = vc$var_intercept / vc$var_total,
    prop_slope         = vc$var_slope     / vc$var_total,
    prop_residual      = vc$prop_residual,
    prop_spatial_total = vc$prop_spatial
  )
}

# Compare to baseline_sp if available
cat("\n\nComparison to baseline_sp (no vegetation):\n")
cat(strrep("-", 40), "\n")

vc_base <- .variance_bundle("baseline_sp")
if (!is.null(vc_base) && !is.na(vc_base$prop_spatial)) {
  cat(sprintf("  Baseline_sp GP variance: %.1f%% (σ_resid=%.2f ‰)\n",
              100 * vc_base$prop_spatial, vc_base$sigma_resid_permil))
  for (model_name in names(variance_results)) {
    reduction <- 100 * (vc_base$prop_spatial - variance_results[[model_name]]$prop_spatial_total) /
                       vc_base$prop_spatial
    cat(sprintf("  %s reduces GP variance by %.1f%%\n", model_name, reduction))
  }
}

# Save results
if (length(variance_results) > 0) {
  variance_df <- bind_rows(variance_results)
  out_dir <- "model_analysis/tables"
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(variance_df, file.path(out_dir, "vegetation_variance_components.csv"),
            row.names = FALSE)
  cat("\n\nVariance components saved to:",
      file.path(out_dir, "vegetation_variance_components.csv"), "\n")

  # Summary table (posterior-mean proportions and residual σ in ‰)
  cat("\nSummary Table:\n")
  cat(strrep("-", 72), "\n")
  cat(sprintf("%-20s %10s %10s %8s %8s\n",
              "Model", "σ_resid (‰)", "prop_int", "prop_sl", "prop_GP"))
  cat(strrep("-", 72), "\n")
  for (i in seq_len(nrow(variance_df))) {
    cat(sprintf("%-20s %10.2f %10.3f %8.3f %8.3f\n",
                variance_df$model[i],
                variance_df$sigma_resid_permil[i],
                variance_df$prop_intercept[i],
                variance_df$prop_slope[i],
                variance_df$prop_spatial_total[i]))
  }
}

cat("\nAnalysis complete!\n")
