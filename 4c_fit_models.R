#!/usr/bin/env Rscript
#───────────────────────────────────────────────────────────────────────────────
# 4c_fit_models.R
#
# Fit Bayesian spatial models using CmdStanR
# Runs multiple model configurations with different covariate combinations
# Handles model compilation, sampling, and diagnostics
#
# Input: results/stan_data_*.rds (prepared Stan data)
# Output: results/model_fits/*.rds (fitted model objects)
#───────────────────────────────────────────────────────────────────────────────

library(cmdstanr)
library(posterior)
library(tidyverse)

# Load configuration
source("0_load_config.R")

# Check CmdStan installation
tryCatch({
  cmdstan_version()
}, error = function(e) {
  stop("CmdStan not found. Please install with: install_cmdstan()")
})

# Define output directories from config
for (dir_name in CONFIG$output_dirs) {
  dir.create(dir_name, showWarnings = FALSE)
}

# Check if Stan file exists
stan_file <- config$scripts$stan_model
if (!file.exists(stan_file)) {
  stop("Stan model file not found: ", stan_file)
}

cat("Using Stan model:", stan_file, "\n")

# Compile Stan model
cat("Compiling Stan model...\n")
compile_start <- Sys.time()
mod <- cmdstan_model(stan_file, compile = TRUE)
compile_time <- difftime(Sys.time(), compile_start, units = "mins")
cat("✓ Model compiled in", round(compile_time, 1), "minutes\n\n")

# Get list of prepared datasets
prepared_data_dir <- CONFIG$output_dirs$prepared_data
prepared_data_files <- list.files(prepared_data_dir, pattern = "^stan_data_.*\\.rds$", full.names = TRUE)
if (length(prepared_data_files) == 0) {
  stop("No prepared data files found. Please run 4b_stan_prep.R first.")
}

# Check if specific models were requested (for parallel execution)
if (!exists("model_names")) {
  # No specific models requested - run all
  # Use gsub to remove both prefix and suffix
  model_names <- gsub("^stan_data_", "", basename(prepared_data_files))
  model_names <- gsub("\\.rds$", "", model_names)
} else {
  # Verify requested models exist
  available_models <- gsub("^stan_data_", "", basename(prepared_data_files))
  available_models <- gsub("\\.rds$", "", available_models)
  invalid_models <- setdiff(model_names, available_models)
  if (length(invalid_models) > 0) {
    stop("Requested models not found: ", paste(invalid_models, collapse = ", "))
  }
}

cat("Will fit", length(model_names), "models:", paste(model_names, collapse = ", "), "\n\n")

# Create summary data frame for results
fit_summary <- data.frame(
  model = character(),
  status = character(),
  runtime_mins = numeric(),
  divergences = integer(),
  max_treedepth = integer(),
  low_ebfmi = integer(),
  max_rhat = numeric(),
  min_ess_bulk = numeric(),
  stringsAsFactors = FALSE
)

# Fit each model
for (model_name in model_names) {
  cat("\n", strrep("=", 65), "\n")
  cat("Fitting model:", model_name, "\n")
  cat(strrep("=", 65), "\n\n")

  elapsed_time <- NA  # Initialize in case we load existing model
  status <- "unknown"  # Must be set before tryCatch
  fit <- NULL

  stan_data_file <- file.path(CONFIG$output_dirs$prepared_data, paste0("stan_data_", model_name, ".rds"))
  config_file <- file.path(CONFIG$output_dirs$prepared_data, paste0("config_", model_name, ".rds"))
  output_dir <- file.path(CONFIG$output_dirs$model_output, model_name)
  
  # Check if model already fitted
  fit_file <- file.path(output_dir, "fit.rds")
  if (file.exists(fit_file)) {
    cat("✓ Model already fit - loading existing results\n")
    fit <- readRDS(fit_file)
    status <- "completed"
    # Load runtime info if available
    runtime_file <- file.path(output_dir, "runtime_info.rds")
    if (file.exists(runtime_file)) {
      runtime_info <- readRDS(runtime_file)
      elapsed_time <- runtime_info$elapsed_mins
    }
    # Still compute diagnostics for summary
  } else {
    # Create output directory
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Load data and config
    stan_data <- readRDS(stan_data_file)
    config <- readRDS(config_file)
    
    # Validate stan_data
    cat("Validating Stan data...\n")
    cat("  N =", stan_data$N, "\n")
    cat("  n_scales =", stan_data$n_scales, "\n")
    cat("  include_c4 =", stan_data$include_c4, "\n")
    cat("  include_pft =", stan_data$include_pft, "\n")
    cat("  include_gp =", stan_data$include_gp, "\n")
    cat("  include_elevation =", stan_data$include_elevation, "\n")
    if (stan_data$include_gp == 1) {
      cat("  n_pp_knots =", stan_data$n_pp_knots, "\n")
      cat("  Coordinates provided for internal kernel computation\n")
    }
    cat("\n")
    
    # Set sampling parameters
    cat("Sampling configuration:\n")
    cat("  Chains:", config$chains, "\n")
    cat("  Total iterations:", config$iter, "\n")
    cat("  Warmup:", config$iter / 2, "\n")
    cat("  Sampling:", config$iter / 2, "\n\n")
    
    # Run Stan
    cat("Running Stan...\n")
    start_time <- Sys.time()
    
    tryCatch({
      fit <- mod$sample(
        data = stan_data,
        chains = config$chains,
        iter_warmup = floor(config$iter * CONFIG$warmup_ratio),
        iter_sampling = floor(config$iter * (1 - CONFIG$warmup_ratio)),
        parallel_chains = config$chains,
        seed = CONFIG$stan_seed,
        refresh = CONFIG$refresh,
        output_dir = output_dir,
        save_warmup = FALSE,
        adapt_delta = config$adapt_delta,
        max_treedepth = config$max_treedepth
      )
      
      elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
      
      # Save fit object
      cat("\nSaving fitted model...\n")
      saveRDS(fit, fit_file)
      cat("✓ Fit saved to:", fit_file, "\n")
      
      # Save runtime info
      runtime_info <- list(
        model = model_name,
        start_time = start_time,
        end_time = Sys.time(),
        elapsed_mins = elapsed_time,
        config = config
      )
      saveRDS(runtime_info, file.path(output_dir, "runtime_info.rds"))
      
      status <- "completed"
      
    }, error = function(e) {
      cat("\n✗ ERROR during sampling:", e$message, "\n")
      elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
      
      # Save error info
      error_info <- list(
        model = model_name,
        error = e$message,
        traceback = traceback(),
        timestamp = Sys.time()
      )
      saveRDS(error_info, file.path(output_dir, "error_info.rds"))
      
      status <<- "failed"
      fit <<- NULL
    })
  }
  
  # Compute diagnostics if fit exists
  if (!is.null(fit) && exists("fit")) {
    cat("\nComputing diagnostics...\n")
    
    # Basic diagnostics
    diagnostics <- fit$diagnostic_summary()
    
    # Parameter summary for key parameters
    params_to_check <- c("beta_0", "beta_oipc", "sigma", "lambda_decay", "effective_scale_km")
    if (stan_data$include_c4 == 1) {
      params_to_check <- c(params_to_check, "beta_c4")
      if (stan_data$include_pft == 1) {
        params_to_check <- c(params_to_check, "beta_oipc_x_c4")
      }
    }
    if (stan_data$include_pft == 1) {
      params_to_check <- c(params_to_check, "beta_tree", "beta_shrub", "beta_grass")
      # Add interaction parameters
      params_to_check <- c(params_to_check, "beta_oipc_x_tree", "beta_oipc_x_shrub", "beta_oipc_x_grass")
    }
    if (stan_data$include_gp == 1) {
  	   params_to_check <- c(params_to_check,
                       "ls_intercept_km", "ls_slope_km", 
                       "sigma_intercept_spatial", "sigma_slope_spatial")
	}
    
    param_summary <- fit$summary(variables = params_to_check)
    
    # Get overall convergence metrics
    all_params_summary <- fit$summary()
    max_rhat <- max(all_params_summary$rhat, na.rm = TRUE)
    min_ess_bulk <- min(all_params_summary$ess_bulk, na.rm = TRUE)
    
    # Print diagnostics
    cat("\nDiagnostic summary:\n")
    cat("  Divergent transitions:", sum(diagnostics$num_divergent), "\n")
    cat("  Max treedepth reached:", sum(diagnostics$num_max_treedepth), "\n")
    cat("  E-BFMI warnings:", sum(diagnostics$ebfmi < 0.2), "\n")
    cat("  Max R-hat:", round(max_rhat, 3), "\n")
    cat("  Min ESS bulk:", round(min_ess_bulk, 0), "\n")
    
    cat("\nKey parameter summary:\n")
    print(param_summary, n = 25)
    
    # Save detailed diagnostics
    diagnostics_full <- list(
      summary = diagnostics,
      param_summary = param_summary,
      all_params_summary = all_params_summary,
      max_rhat = max_rhat,
      min_ess_bulk = min_ess_bulk
    )
    saveRDS(diagnostics_full, file.path(output_dir, "diagnostics.rds"))
    
    # Add to summary data frame
    fit_summary <- rbind(fit_summary, data.frame(
      model = model_name,
      # Override status if fit is NULL despite status not being "failed"
      status = ifelse(is.null(fit) && status != "failed", "failed", status),
      runtime_mins = ifelse(!is.na(elapsed_time), round(elapsed_time, 1), NA),
      divergences = sum(diagnostics$num_divergent),
      max_treedepth = sum(diagnostics$num_max_treedepth),
      low_ebfmi = sum(diagnostics$ebfmi < 0.2),
      max_rhat = round(max_rhat, 3),
      min_ess_bulk = round(min_ess_bulk, 0),
      stringsAsFactors = FALSE
    ))
    
    # Extract and save posterior draws for key parameters
    cat("\nExtracting posterior draws...\n")
    draws <- fit$draws(variables = params_to_check, format = "draws_df")
    saveRDS(draws, file.path(output_dir, "posterior_draws.rds"))
    
    if (!is.na(elapsed_time)) {
      cat("\n✓ Model fitting completed in", round(elapsed_time, 1), "minutes\n")
    }
  }
}

# Save overall summary
cat("\n", strrep("=", 65), "\n")
cat("ALL MODELS COMPLETED\n")
cat(strrep("=", 65), "\n\n")

print(fit_summary)

results_dir <- CONFIG$output_dirs$results
saveRDS(fit_summary, file.path(results_dir, "model_fit_summary.rds"))
write.csv(fit_summary, file.path(results_dir, "model_fit_summary.csv"), row.names = FALSE)

cat("\nSummary saved to:\n")
cat("  -", file.path(results_dir, "model_fit_summary.rds"), "\n")
cat("  -", file.path(results_dir, "model_fit_summary.csv"), "\n")

# Completion guard: only write success marker if ALL models completed.
# The launcher script uses pipeline_4c_complete.rds as the auto-shutdown gate.
n_completed <- sum(fit_summary$status == "completed")
n_failed <- sum(fit_summary$status == "failed")
n_expected <- length(model_names)

# Cross-check: count error_info.rds files independently of fit_summary
error_files <- list.files(CONFIG$output_dirs$model_output,
                          pattern = "error_info\\.rds$",
                          recursive = TRUE, full.names = TRUE)
n_error_files <- length(error_files)

completion_info <- list(
  completed_at = Sys.time(),
  models_fitted = fit_summary$model[fit_summary$status == "completed"],
  models_failed = fit_summary$model[fit_summary$status == "failed"],
  n_expected = n_expected,
  n_completed = n_completed,
  n_failed = n_failed,
  n_error_files = n_error_files,
  total_runtime_mins = sum(fit_summary$runtime_mins, na.rm = TRUE)
)

if (n_completed == n_expected && n_error_files == 0) {
  saveRDS(completion_info, file.path(results_dir, "pipeline_4c_complete.rds"))
  cat("\n✓ All", n_expected, "models completed successfully.\n")
  cat("  Completion marker written: pipeline_4c_complete.rds\n")
} else {
  saveRDS(completion_info, file.path(results_dir, "pipeline_4c_INCOMPLETE.rds"))
  cat("\n✗ INCOMPLETE:", n_completed, "of", n_expected, "models completed,",
      n_failed, "failed,", n_error_files, "error files found.\n")
  if (n_failed > 0) {
    cat("  Failed models:", paste(fit_summary$model[fit_summary$status == "failed"],
                                  collapse = ", "), "\n")
  }
  if (n_error_files > 0) {
    cat("  Error files:", paste(error_files, collapse = "\n              "), "\n")
  }
  cat("  INCOMPLETE marker written: pipeline_4c_INCOMPLETE.rds\n")
  cat("  The launcher will NOT auto-shutdown.\n")
}