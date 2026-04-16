# run_sensitivity.R
# Phase 4: Prior/Hyperparameter Sensitivity Analysis
#
# Varies one prior at a time and re-fits baseline_veg_sp on real data.
# Uses the fixed range_factor thresholds.
#
# Usage: Rscript run_sensitivity.R [experiment_id]
#   experiment_id: "4a2", "4a3", "4a4", "4b2", "4b3", "4c2", "4c3", or "all"
#   Baselines (4a1, 4b1, 4c1) are just the fixed-range run from Phase 2.
#
# Expected runtime: ~2-3 hrs per experiment on r6i.8xlarge

library(cmdstanr)
library(posterior)
library(loo)

args <- commandArgs(trailingOnly = TRUE)
run_ids <- if (length(args) > 0 && args[1] != "all") args else "all"

cat("=== SENSITIVITY ANALYSIS ===\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

setwd("/home/shared/leafwax_spatial/spatial_leafwax_model")

# Define experiments (non-baseline only)
experiments <- list(
  "4a2" = list(
    name = "beta_oipc_prior_wider",
    description = "beta_oipc ~ Normal(0.8, 1.0)",
    type = "stan_model",
    find_line = "beta_oipc ~ normal(0.8, 0.3);",
    replace_line = "beta_oipc ~ normal(0.8, 1.0);"
  ),
  "4a3" = list(
    name = "beta_oipc_prior_shifted",
    description = "beta_oipc ~ Normal(0.5, 1.0)",
    type = "stan_model",
    find_line = "beta_oipc ~ normal(0.8, 0.3);",
    replace_line = "beta_oipc ~ normal(0.5, 1.0);"
  ),
  "4a4" = list(
    name = "beta_oipc_prior_uninformative",
    description = "beta_oipc ~ Normal(0, 2.0)",
    type = "stan_model",
    find_line = "beta_oipc ~ normal(0.8, 0.3);",
    replace_line = "beta_oipc ~ normal(0, 2.0);"
  ),
  "4b2" = list(
    name = "pc_slope_relaxed",
    description = "PC prior: P(sigma_slope > 0.5) = 0.05",
    type = "stan_data",
    modifications = list(pc_prior_slope_u = 0.5, pc_prior_slope_alpha = 0.05)
  ),
  "4b3" = list(
    name = "pc_slope_very_relaxed",
    description = "PC prior: P(sigma_slope > 1.0) = 0.10",
    type = "stan_data",
    modifications = list(pc_prior_slope_u = 1.0, pc_prior_slope_alpha = 0.10)
  ),
  "4c2" = list(
    name = "ls_longer",
    description = "log_ls_spatial ~ Normal(-0.5, 0.4)",
    type = "stan_model",
    find_line = "log_ls_spatial_raw[1] ~ normal(-1.0, 0.4);",
    replace_line = "log_ls_spatial_raw[1] ~ normal(-0.5, 0.4);"
  ),
  "4c3" = list(
    name = "ls_shorter",
    description = "log_ls_spatial ~ Normal(-1.5, 0.4)",
    type = "stan_model",
    find_line = "log_ls_spatial_raw[1] ~ normal(-1.0, 0.4);",
    replace_line = "log_ls_spatial_raw[1] ~ normal(-1.5, 0.4);"
  )
)

# Determine which experiments to run
if (run_ids[1] == "all") {
  to_run <- names(experiments)
} else {
  to_run <- run_ids
}

cat("Running experiments:", paste(to_run, collapse = ", "), "\n\n")

# Load base data
stan_data_base <- readRDS("prepared_data_consolidated/stan_data_baseline_veg_sp.rds")
config <- readRDS("prepared_data_consolidated/config_baseline_veg_sp.rds")

# Read base Stan model
stan_code_base <- readLines("4d_leaf_wax_spatial_model.stan")

# Run each experiment
results_summary <- list()

for (exp_id in to_run) {
  exp <- experiments[[exp_id]]
  if (is.null(exp)) {
    cat("Unknown experiment:", exp_id, "- skipping\n")
    next
  }
  cat("--- Experiment", exp_id, ":", exp$description, "---\n")

  outdir <- file.path("model_output", paste0("sensitivity_", exp$name))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  stan_data_exp <- stan_data_base

  if (exp$type == "stan_data") {
    for (field in names(exp$modifications)) {
      stan_data_exp[[field]] <- exp$modifications[[field]]
    }
    model <- cmdstan_model("4d_leaf_wax_spatial_model.stan")

  } else if (exp$type == "stan_model") {
    stan_code_mod <- stan_code_base
    line_idx <- grep(exp$find_line, stan_code_mod, fixed = TRUE)
    if (length(line_idx) == 0) {
      cat("  WARNING: Could not find prior line '", exp$find_line, "' - skipping\n\n")
      next
    }
    leading_ws <- sub("^(\\s*).*", "\\1", stan_code_mod[line_idx[1]])
    stan_code_mod[line_idx[1]] <- paste0(leading_ws, exp$replace_line,
                                          "  // SENSITIVITY ", exp_id)
    tmp_stan <- file.path(outdir, "model_modified.stan")
    writeLines(stan_code_mod, tmp_stan)
    model <- cmdstan_model(tmp_stan)
  }

  saveRDS(exp, file.path(outdir, "experiment_config.rds"))

  fit_start <- Sys.time()
  fit <- model$sample(
    data = stan_data_exp,
    seed = 314,
    chains = config$chains,
    parallel_chains = config$chains,
    iter_warmup = as.integer(config$iter * 0.5),
    iter_sampling = as.integer(config$iter * 0.5),
    adapt_delta = config$adapt_delta,
    max_treedepth = config$max_treedepth,
    refresh = 200
  )
  fit_time <- as.numeric(difftime(Sys.time(), fit_start, units = "mins"))
  cat("  Completed in", round(fit_time, 1), "minutes\n")

  draws <- fit$draws(format = "df")
  key_params <- c("beta_0", "beta_oipc", "sigma", "lambda_decay",
                   "sigma_intercept_spatial", "sigma_slope_spatial",
                   "ls_intercept_km", "rmse", "r_squared",
                   "var_spatial_intercept", "var_spatial_slope",
                   "min_oipc_slope", "max_oipc_slope", "sd_oipc_slope")
  available <- intersect(key_params, names(draws))
  saveRDS(draws[, c(available, ".chain", ".iteration", ".draw")],
          file.path(outdir, "posterior_draws.rds"))
  saveRDS(fit$diagnostic_summary(), file.path(outdir, "diagnostics.rds"))
  saveRDS(fit$loo(), file.path(outdir, "loo.rds"))
  saveRDS(list(experiment = exp_id, elapsed_mins = fit_time,
               git_commit = system("git rev-parse HEAD", intern = TRUE)),
          file.path(outdir, "runtime_info.rds"))

  results_summary[[exp_id]] <- data.frame(
    experiment = exp_id,
    name = exp$name,
    beta_oipc_mean = mean(draws$beta_oipc),
    beta_oipc_lo = quantile(draws$beta_oipc, 0.025),
    beta_oipc_hi = quantile(draws$beta_oipc, 0.975),
    sigma_slope = mean(draws$sigma_slope_spatial),
    sigma_intercept = mean(draws$sigma_intercept_spatial),
    looic = fit$loo()$estimates["looic", "Estimate"],
    runtime_mins = fit_time,
    stringsAsFactors = FALSE
  )

  cat("  beta_oipc:", round(mean(draws$beta_oipc), 3),
      "[", round(quantile(draws$beta_oipc, 0.025), 3), ",",
      round(quantile(draws$beta_oipc, 0.975), 3), "]\n\n")
}

# Print summary table
cat("\n=== SENSITIVITY ANALYSIS SUMMARY ===\n\n")
summary_df <- do.call(rbind, results_summary)
rownames(summary_df) <- NULL
print(summary_df, digits = 3)

saveRDS(summary_df, "model_output/sensitivity_summary.rds")
cat("\nEnd time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
