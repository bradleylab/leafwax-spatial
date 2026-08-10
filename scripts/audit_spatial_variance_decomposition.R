#!/usr/bin/env Rscript

# Audit the spatial intercept/slope variance decomposition used by Figure 2
# and the manuscript text.
#
# This script compares three calculations from the saved posterior draws:
#   1. The legacy reporting convention, which evaluates the spatial slope
#      contribution against the standardized `oipc_d2h20` field. Despite its
#      historical name, this field is the unweighted mean of the OIPC pixels
#      extracted within the preprocessing radius, not the model's fitted
#      scale-weighted predictor.
#   2. The same component-sum convention evaluated against each posterior
#      draw's fitted, scale-weighted OIPC predictor.
#   3. The realized fitted spatial contribution, including covariance between
#      the spatial intercept and spatial slope contributions.
#
# Run from the repository root against the authoritative chordal mirror:
#   LEAFWAX_RUN_DIR=results/c2_run_20260728_chordal/model_output \
#     Rscript scripts/audit_spatial_variance_decomposition.R
#
# The calculation uses every retained posterior draw and is deterministic.
# The output CSV is tracked with the script so the audit numbers do not exist
# only in an interactive session.
#
# No Stan change or refit is required. The fitted likelihood already uses
# `oipc_weighted` and the local spatial slope. The Stan generated-quantities
# diagnostic named `var_total_spatial` simply adds the two marginal component
# variances and omits their covariance; it is not used for the manuscript
# comparison or for the realized-spatial diagnostic reported by this audit.

suppressPackageStartupMessages(library(posterior))

source("scripts/posterior_helpers.R")
source("scripts/spatial_variance_helpers.R")

SPATIAL_MODELS <- c(
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

OUTPUT_CSV <- Sys.getenv(
  "LEAFWAX_VARIANCE_AUDIT_OUT",
  unset = "scripts/reference_outputs/spatial_variance_decomposition_audit.csv"
)
SEDIMENT <- load_sediment()

audit_model <- function(model) {
  message("Auditing ", model, "...")

  draws <- load_draws(model)
  stan_data <- load_stan_data(model)
  config <- load_config(model)
  sediment <- SEDIMENT

  n_draws <- ndraws(draws)
  n_obs <- stan_data$N
  n_scales <- stan_data$n_scales
  if (!identical(dim(stan_data$oipc_values), c(n_obs, n_scales))) {
    stop("OIPC predictor matrix does not align with posterior arrays for ", model)
  }
  if (nrow(sediment) != n_obs) {
    stop("Sediment data do not align with posterior arrays for ", model)
  }
  fitted_oipc <- fitted_oipc_draw_matrix(draws, stan_data, model)
  scale_weight_max_sum_error <- attr(
    fitted_oipc, "scale_weight_max_sum_error"
  )
  legacy_oipc <-
    (sediment$oipc_d2h20 - config$scaling_params$oipc_mean) /
    config$scaling_params$oipc_sd
  if (any(!is.finite(legacy_oipc))) {
    stop("Legacy OIPC field contains non-finite values for ", model)
  }

  legacy <- summarize_spatial_components(draws, legacy_oipc, model)
  fitted <- summarize_spatial_components(draws, fitted_oipc, model)

  legacy_component_sum <-
    legacy$marginal_intercept_variance_std2 +
    legacy$marginal_slope_variance_std2
  fitted_component_sum <-
    fitted$marginal_intercept_variance_std2 +
    fitted$marginal_slope_variance_std2

  data.frame(
    model = model,
    n_draws = n_draws,
    n_obs = n_obs,
    n_scales = n_scales,
    scale_weight_max_sum_error_before_normalization =
      scale_weight_max_sum_error,
    legacy_intercept_share_pct =
      legacy$marginal_intercept_share_pct,
    legacy_spatial_share_component_sum_pct =
      100 * legacy_component_sum /
      (legacy_component_sum + legacy$residual_variance_std2),
    fitted_intercept_share_pct =
      fitted$marginal_intercept_share_pct,
    fitted_slope_share_pct =
      fitted$marginal_slope_share_pct,
    fitted_spatial_share_component_sum_pct =
      100 * fitted_component_sum /
      (fitted_component_sum + fitted$residual_variance_std2),
    fitted_spatial_share_realized_pct =
      fitted$realized_spatial_share_of_spatial_plus_residual_pct,
    fitted_covariance_twice_mean = fitted$twice_covariance_std2,
    fitted_covariance_share_of_realized_spatial_pct =
      100 * fitted$twice_covariance_std2 /
      fitted$realized_spatial_variance_std2,
    fitted_var_intercept_mean = fitted$marginal_intercept_variance_std2,
    fitted_var_slope_mean = fitted$marginal_slope_variance_std2,
    fitted_var_spatial_realized_mean = fitted$realized_spatial_variance_std2,
    fitted_var_residual_mean = fitted$residual_variance_std2,
    fitted_identity_max_abs_error = fitted$identity_max_abs_error,
    legacy_identity_max_abs_error = legacy$identity_max_abs_error,
    check.names = FALSE
  )
}

audit <- do.call(rbind, lapply(SPATIAL_MODELS, audit_model))

dir.create(dirname(OUTPUT_CSV), recursive = TRUE, showWarnings = FALSE)
write.csv(audit, OUTPUT_CSV, row.names = FALSE)

cat("\nWrote", OUTPUT_CSV, "\n\n")
print(audit[, c(
  "model",
  "legacy_intercept_share_pct",
  "fitted_intercept_share_pct",
  "fitted_spatial_share_component_sum_pct",
  "fitted_spatial_share_realized_pct",
  "fitted_covariance_share_of_realized_spatial_pct"
)], row.names = FALSE)

cat(sprintf(
  paste0(
    "\nFitted-predictor intercept share of component variance: %.1f--%.1f%%\n",
    "Fitted-predictor spatial share (component-sum convention): %.1f--%.1f%%\n",
    "Fitted-predictor spatial share (including covariance): %.1f--%.1f%%\n"
  ),
  min(audit$fitted_intercept_share_pct),
  max(audit$fitted_intercept_share_pct),
  min(audit$fitted_spatial_share_component_sum_pct),
  max(audit$fitted_spatial_share_component_sum_pct),
  min(audit$fitted_spatial_share_realized_pct),
  max(audit$fitted_spatial_share_realized_pct)
))
