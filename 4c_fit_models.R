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
library(loo)

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

  # Load stan_data and config eagerly — needed by both the resume-fit path
  # (post-sampling diagnostics reference stan_data$include_gp etc.) and the
  # fresh-fit path. Both files are small.
  stan_data <- readRDS(stan_data_file)
  config <- readRDS(config_file)

  # Skip-load only if fit.rds is newer than its stan_data input. A bare
  # file.exists() check would silently reuse a stale fit after stan_data
  # was rebuilt against a new compilation.
  fit_file <- file.path(output_dir, "fit.rds")
  fit_is_current <- file.exists(fit_file) &&
    file.mtime(fit_file) > file.mtime(stan_data_file)
  if (fit_is_current) {
    cat("✓ Model already fit (newer than stan_data) - loading existing results\n")
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
    
    # Scalar parameters for the diagnostic summary (printed to stdout).
    # `beta_precip` is gated by include_precip in the Stan model
    # (4d_leaf_wax_spatial_model.stan:297).
    params_to_check <- c("beta_0", "beta_oipc", "sigma",
                         "lambda_decay", "effective_scale_km")
    if (isTRUE(stan_data$include_precip == 1)) {
      params_to_check <- c(params_to_check, "beta_precip")
    }
    if (isTRUE(stan_data$include_elevation == 1)) {
      # Retain the elevation B-spline coefficients (and their smoothing SD) in the
      # saved draws so the released posteriors carry the fitted elevation response.
      # These were previously dropped here at the save step (not at export), which
      # is why has_elevation came back FALSE downstream. The package derives
      # has_elevation from column presence (load_posteriors.R), so retaining them
      # flips it automatically and makes Fig S6 reproducible.
      params_to_check <- c(params_to_check, "beta_elev_bspline", "tau_elev_bspline")
    }
    if (stan_data$include_c4 == 1) {
      params_to_check <- c(params_to_check, "beta_c4")
      # The C4 interaction is active whenever include_veg_interactions && include_c4
      # (see 4d_leaf_wax_spatial_model.stan:306), independent of PFT. Gating this
      # on include_pft skipped saving beta_oipc_x_c4 for C4-interaction models
      # without PFT (e.g. elevation_c4_interact_sp). Gate on the interaction flag.
      if (stan_data$include_veg_interactions == 1) {
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
    
    # Posterior draws to save. The widened set includes every variable that
    # downstream analysis reads per-draw. Indexed arrays are requested by base
    # name; cmdstanr expands them internally. Stored as draws_array for lower
    # RAM / smaller file size vs draws_df; readers convert locally if needed.
    #
    # Always-present (response-level):
    draws_to_save <- c(params_to_check,
                       "mu", "d2H_rep", "log_lik", "scale_weights")
    # GP-only (per-observation + per-knot spatial effects):
    if (stan_data$include_gp == 1) {
      draws_to_save <- c(draws_to_save,
                         "alpha_spatial",
                         "beta_oipc_spatial",
                         "z_intercept_spatial",
                         "z_slope_spatial")
    }

    cat("\nExtracting posterior draws (", length(draws_to_save),
        "variables requested)...\n", sep = "")
    draws <- fit$draws(variables = draws_to_save, format = "draws_array")
    saveRDS(draws, file.path(output_dir, "posterior_draws.rds"))

    # Emit loo.rds so analysis never needs fit.rds or chain CSVs for model
    # comparison (Table 1 LOOIC, SE, p_eff, Pareto-k counts). Computing from
    # the draws we just saved keeps the data source consistent with the
    # widened posterior_draws.rds.
    cat("\nComputing loo...\n")
    log_lik_array <- posterior::subset_draws(draws, variable = "log_lik")
    loo_cores <- min(4L, max(1L, parallel::detectCores() - 1L))
    loo_result <- tryCatch(
      loo::loo(log_lik_array, cores = loo_cores),
      error = function(e) {
        message("loo::loo failed for ", model_name, ": ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(loo_result)) {
      saveRDS(loo_result, file.path(output_dir, "loo.rds"))
      cat("  elpd_loo =", round(loo_result$estimates["elpd_loo", "Estimate"], 1),
          "(SE =", round(loo_result$estimates["elpd_loo", "SE"], 1), ")\n")
    }

    if (!is.na(elapsed_time)) {
      cat("\n✓ Model fitting completed in", round(elapsed_time, 1), "minutes\n")
    }
  }
}

# Save overall summary
cat("\n", strrep("=", 65), "\n")
cat("MODEL FITTING FINISHED — verifying per-model completeness\n")
cat(strrep("=", 65), "\n\n")

print(fit_summary)

results_dir <- CONFIG$output_dirs$results

# ── Per-model completeness + status artifacts (concurrency-safe) ────────────────
# When this process fits a single model (the SLURM array case — one model per
# task), many array tasks run 4c concurrently against a SHARED results/ dir.
# Writing model_fit_summary.* / pipeline_4c_complete.rds there would have every
# task clobber the others (last-writer-wins), and each single-model task would
# write a "pipeline complete" marker meaning only "my 1 model finished". Instead
# every task writes per-model artifacts under model_output/<model>/ (no shared
# path, no race): fit_status.rds always, and a DONE marker ONLY when the full
# core artifact set is on disk. A separate aggregator
# (scripts/aggregate_chordal_status.R) reads these to judge the whole batch.
core_artifacts_complete <- function(model) {
  mdir <- file.path(CONFIG$output_dirs$model_output, model)
  all(file.exists(file.path(mdir, c("fit.rds", "diagnostics.rds",
                                    "posterior_draws.rds"))))
}
loo_present <- function(model) {
  file.exists(file.path(CONFIG$output_dirs$model_output, model, "loo.rds"))
}

# Iterate the REQUESTED models, not fit_summary rows: a fully-failed fit adds no
# fit_summary row (its diagnostics block is skipped when fit is NULL), yet it
# still needs a fit_status.rds and, critically, removal of any stale DONE.
sget <- function(col, i, default = NA) if (!is.na(i)) fit_summary[[col]][i] else default
for (m in model_names) {
  i       <- match(m, fit_summary$model)   # NA if this model produced no row (failed)
  mdir    <- file.path(CONFIG$output_dirs$model_output, m)
  dir.create(mdir, showWarnings = FALSE, recursive = TRUE)
  core_ok <- core_artifacts_complete(m)
  loo_ok  <- loo_present(m)
  m_status <- if (!is.na(i)) fit_summary$status[i] else "failed"
  is_ok   <- identical(m_status, "completed") && core_ok
  status_row <- list(
    model = m,
    status = m_status,
    core_complete = core_ok,      # fit.rds + diagnostics.rds + posterior_draws.rds
    loo_present = loo_ok,         # best-effort; loo::loo can legitimately fail
    divergences = sget("divergences", i),
    max_treedepth = sget("max_treedepth", i),
    low_ebfmi = sget("low_ebfmi", i),
    max_rhat = sget("max_rhat", i),
    min_ess_bulk = sget("min_ess_bulk", i),
    runtime_mins = sget("runtime_mins", i),
    completed_at = Sys.time()
  )
  saveRDS(status_row, file.path(mdir, "fit_status.rds"))
  done_mark <- file.path(mdir, "DONE")
  running_mark <- file.path(mdir, "RUNNING")
  if (is_ok) {
    writeLines(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), done_mark)
    unlink(running_mark)                       # completed: clear the in-progress flag
    # A successful (re)fit supersedes any prior failure: drop this model's stale
    # error_info.rds so a recovered model is not held failed by an old record.
    unlink(file.path(mdir, "error_info.rds"))
  } else if (file.exists(done_mark)) {
    # A prior run's stale DONE must not survive a re-run that came up incomplete.
    file.remove(done_mark)
  }
  if (!loo_ok) {
    cat("  WARNING:", m, "has no loo.rds (loo::loo failed?) — LOOIC/Table-1",
        "will be unavailable for this model.\n")
  }
}

# Aggregate outcome for THIS process's requested models.
n_expected  <- length(model_names)
n_completed <- sum(vapply(model_names, function(m)
  identical(fit_summary$status[match(m, fit_summary$model)], "completed") &&
    core_artifacts_complete(m), logical(1)))
n_failed    <- n_expected - n_completed
# Scope the error scan to the REQUESTED models only. A global recursive scan makes
# one failed array task fail unrelated successful tasks, and a stale error file
# from another model's prior failure would wrongly hold this task failed. (This
# model's own stale error_info.rds was already removed above on success.)
error_files <- file.path(CONFIG$output_dirs$model_output, model_names, "error_info.rds")
error_files <- error_files[file.exists(error_files)]

completion_info <- list(
  completed_at = Sys.time(),
  requested = model_names,
  models_complete = model_names[vapply(model_names, function(m)
    identical(fit_summary$status[match(m, fit_summary$model)], "completed") &&
      core_artifacts_complete(m), logical(1))],
  n_expected = n_expected,
  n_completed = n_completed,
  n_failed = n_failed,
  n_error_files = length(error_files),
  total_runtime_mins = sum(fit_summary$runtime_mins, na.rm = TRUE)
)

if (n_expected == 1) {
  # Array / single-model mode: model-scoped summary only (no shared clobber).
  m <- model_names[[1]]
  write.csv(fit_summary, file.path(results_dir, paste0("model_fit_summary_", m, ".csv")),
            row.names = FALSE)
  saveRDS(completion_info, file.path(results_dir, paste0("fit_completion_", m, ".rds")))
} else {
  # Full-set / legacy launcher mode: shared summary + auto-shutdown gate marker.
  saveRDS(fit_summary, file.path(results_dir, "model_fit_summary.rds"))
  write.csv(fit_summary, file.path(results_dir, "model_fit_summary.csv"), row.names = FALSE)
  marker <- if (n_completed == n_expected && length(error_files) == 0)
    "pipeline_4c_complete.rds" else "pipeline_4c_INCOMPLETE.rds"
  saveRDS(completion_info, file.path(results_dir, marker))
  cat("\nCompletion marker written:", marker, "\n")
}

if (n_completed == n_expected && length(error_files) == 0) {
  cat("\n✓", n_completed, "of", n_expected, "requested model(s) completed with a",
      "full artifact set.\n")
} else {
  cat("\n✗ INCOMPLETE:", n_completed, "of", n_expected, "requested model(s) complete;",
      n_failed, "failed/incomplete,", length(error_files), "error file(s) found.\n")
  if (length(error_files) > 0) {
    cat("  Error files:", paste(error_files, collapse = "\n              "), "\n")
  }
  # Nonzero exit so the SLURM task (and the fit wrapper's `set -e`) register the
  # failure instead of false-succeeding. This is the backstop to the shell-level
  # complete-artifact-set resume guard in job_fit_chordal.sh.
  quit(status = 1)
}