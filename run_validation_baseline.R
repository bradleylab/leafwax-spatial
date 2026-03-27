# run_validation_baseline.R
# Phase 2: Re-fit baseline_veg_sp with fixed OIPC range thresholds
# This script reuses existing prepared Stan data but recompiles the updated model
#
# Usage: Rscript run_validation_baseline.R
# Expected runtime: ~2-3 hrs on r6i.8xlarge (32 vCPUs, 8 chains)

library(cmdstanr)
library(posterior)
library(loo)

cat("=== VALIDATION RUN: baseline_veg_sp with fixed range thresholds ===\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

setwd("/home/shared/leafwax_spatial/spatial_leafwax_model")

# Load the SAME prepared data that was used for the Sep 2025 run
stan_data <- readRDS("prepared_data_consolidated/stan_data_baseline_veg_sp.rds")
config <- readRDS("prepared_data_consolidated/config_baseline_veg_sp.rds")

cat("N observations:", stan_data$N, "\n")
cat("N scales:", stan_data$n_scales, "\n")
cat("N knots:", stan_data$n_pp_knots, "\n")
cat("Config: chains =", config$chains, ", iter =", config$iter, "\n")
cat("Config: adapt_delta =", config$adapt_delta, ", max_treedepth =", config$max_treedepth, "\n\n")

# Compile the UPDATED Stan model (with fixed range thresholds)
cat("Compiling updated Stan model...\n")
compile_start <- Sys.time()
model <- cmdstan_model("4d_leaf_wax_spatial_model.stan")
compile_time <- as.numeric(difftime(Sys.time(), compile_start, units = "mins"))
cat("Model compiled in", round(compile_time, 1), "minutes\n\n")

# Fit the model
cat("Starting MCMC sampling...\n")
fit_start <- Sys.time()

fit <- model$sample(
  data = stan_data,
  seed = ifelse(is.null(config$scaling_params$stan_seed), 314, config$scaling_params$stan_seed),
  chains = config$chains,
  parallel_chains = config$chains,  # Use all 8 chains in parallel
  iter_warmup = as.integer(config$iter * 0.5),
  iter_sampling = as.integer(config$iter * 0.5),
  adapt_delta = config$adapt_delta,
  max_treedepth = config$max_treedepth,
  refresh = 100
)

fit_time <- as.numeric(difftime(Sys.time(), fit_start, units = "mins"))
cat("\nModel fitting completed in", round(fit_time, 1), "minutes\n\n")

# Create output directory
outdir <- "model_output/baseline_veg_sp_fixed_range"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Save fit object
fit$save_object(file.path(outdir, "fit.rds"))

# Extract and save posterior draws (key parameters only, like the original pipeline)
draws <- fit$draws(format = "df")
key_params <- c("beta_0", "beta_oipc", "sigma", "lambda_decay", "effective_scale_km",
                 "beta_c4", "beta_oipc_x_c4", "beta_tree", "beta_shrub", "beta_grass",
                 "beta_oipc_x_tree", "beta_oipc_x_shrub", "beta_oipc_x_grass",
                 "ls_intercept_km", "ls_slope_km",
                 "sigma_intercept_spatial", "sigma_slope_spatial",
                 "min_oipc_slope", "max_oipc_slope", "sd_oipc_slope",
                 "var_spatial_intercept", "var_spatial_slope",
                 "prop_variance_spatial", "prop_variance_residual",
                 "rmse", "r_squared",
                 "range_threshold_low", "range_threshold_high")
available_params <- intersect(key_params, names(draws))
saveRDS(draws[, c(available_params, ".chain", ".iteration", ".draw")],
        file.path(outdir, "posterior_draws.rds"))

# Save diagnostics
diag_summary <- fit$diagnostic_summary()
saveRDS(diag_summary, file.path(outdir, "diagnostics.rds"))

# Compute and save LOO
loo_result <- fit$loo()
saveRDS(loo_result, file.path(outdir, "loo.rds"))

# Save tau diagnostics (these are in generated quantities now)
tau_slope_cols <- grep("^tau_final_slope", names(draws), value = TRUE)
tau_int_cols <- grep("^tau_final_intercept", names(draws), value = TRUE)
range_cols <- grep("^knot_oipc_ranges", names(draws), value = TRUE)
if (length(tau_slope_cols) > 0) {
  tau_draws <- draws[, c(tau_slope_cols, tau_int_cols, range_cols)]
  saveRDS(tau_draws, file.path(outdir, "tau_diagnostics.rds"))
}

# Save runtime info
runtime_info <- list(
  model = "baseline_veg_sp_fixed_range",
  start_time = fit_start,
  end_time = Sys.time(),
  elapsed_mins = fit_time,
  compile_mins = compile_time,
  config = config,
  git_commit = system("git rev-parse HEAD", intern = TRUE),
  description = "baseline_veg_sp with fixed OIPC range factor thresholds (data-relative)"
)
saveRDS(runtime_info, file.path(outdir, "runtime_info.rds"))

# Print summary
cat("\n=== RESULTS SUMMARY ===\n")
cat("beta_oipc: mean =", mean(draws$beta_oipc), 
    " 95% CI = [", quantile(draws$beta_oipc, 0.025), ",", quantile(draws$beta_oipc, 0.975), "]\n")
cat("sigma_slope_spatial: mean =", mean(draws$sigma_slope_spatial),
    " 95% CI = [", quantile(draws$sigma_slope_spatial, 0.025), ",", quantile(draws$sigma_slope_spatial, 0.975), "]\n")
cat("sigma_intercept_spatial: mean =", mean(draws$sigma_intercept_spatial), "\n")
cat("ls_km: mean =", mean(draws$ls_intercept_km), "\n")
cat("RMSE: mean =", mean(draws$rmse), "\n")
cat("R-squared: mean =", mean(draws$r_squared), "\n")

if ("range_threshold_low" %in% names(draws)) {
  cat("\nRange thresholds (standardized):\n")
  cat("  Low (25%):", mean(draws$range_threshold_low), "\n")
  cat("  High (60%):", mean(draws$range_threshold_high), "\n")
}

# Compare to original
cat("\n=== COMPARISON TO ORIGINAL baseline_veg_sp ===\n")
orig <- readRDS("model_output/baseline_veg_sp/posterior_draws.rds")
cat("ORIGINAL beta_oipc: mean =", mean(orig$beta_oipc),
    " 95% CI = [", quantile(orig$beta_oipc, 0.025), ",", quantile(orig$beta_oipc, 0.975), "]\n")
cat("FIXED    beta_oipc: mean =", mean(draws$beta_oipc),
    " 95% CI = [", quantile(draws$beta_oipc, 0.025), ",", quantile(draws$beta_oipc, 0.975), "]\n")

cat("\nDone! Results saved to:", outdir, "\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
