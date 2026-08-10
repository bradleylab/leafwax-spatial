# posterior_helpers.R — helpers for reading a selected fitted-model run.
#
# Model-run selection (in priority order):
#   1. Environment variable LEAFWAX_RUN_DIR (absolute or relative path)
#   2. results/c2_run_20260728_chordal/model_output (reported analysis)
#
# Override at the shell:
#   LEAFWAX_RUN_DIR=<run>/model_output Rscript 5a_model_validation.R

library(posterior)

DEFAULT_RUN_NAME <- file.path("c2_run_20260728_chordal", "model_output")

.resolve_run_dir <- function() {
  env_override <- Sys.getenv("LEAFWAX_RUN_DIR", unset = "")
  if (nzchar(env_override)) {
    return(normalizePath(env_override, mustWork = FALSE))
  }
  # Try ../results/<run> first (sourced from scripts/), then results/<run>
  candidate <- normalizePath(
    file.path("..", "results", DEFAULT_RUN_NAME),
    mustWork = FALSE
  )
  if (!dir.exists(candidate)) {
    candidate <- normalizePath(
      file.path("results", DEFAULT_RUN_NAME),
      mustWork = FALSE
    )
  }
  candidate
}

MODEL_RUN_DIR <- .resolve_run_dir()
RUN_ID <- basename(dirname(MODEL_RUN_DIR))
PREPARED_DATA <- file.path(MODEL_RUN_DIR, "_prepared_data")

load_draws <- function(model) {
  path <- file.path(MODEL_RUN_DIR, model, "posterior_draws.rds")
  if (!file.exists(path)) {
    stop("posterior_draws.rds missing for '", model, "' at ", path,
         "\nRun 4c_fit_models.R on C2 against this model's output dir, ",
         "or set LEAFWAX_RUN_DIR to the correct mirror.")
  }
  readRDS(path)
}

load_summaries <- function(model) {
  path <- file.path(MODEL_RUN_DIR, model, "diagnostics.rds")
  if (!file.exists(path)) stop("diagnostics.rds missing at ", path)
  readRDS(path)$all_params_summary
}

load_sampler_diag <- function(model) {
  path <- file.path(MODEL_RUN_DIR, model, "diagnostics.rds")
  if (!file.exists(path)) stop("diagnostics.rds missing at ", path)
  readRDS(path)$summary
}

load_loo <- function(model) {
  path <- file.path(MODEL_RUN_DIR, model, "loo.rds")
  if (!file.exists(path)) {
    stop("loo.rds missing for '", model, "' at ", path,
         "\nRun 4c_fit_models.R on C2 against this model's output dir, ",
         "or set LEAFWAX_RUN_DIR to the correct mirror.")
  }
  readRDS(path)
}

load_stan_data <- function(model) {
  path <- file.path(PREPARED_DATA, paste0("stan_data_", model, ".rds"))
  if (!file.exists(path)) stop("stan_data missing at ", path)
  readRDS(path)
}

load_config <- function(model) {
  path <- file.path(PREPARED_DATA, paste0("config_", model, ".rds"))
  if (!file.exists(path)) stop("model config missing at ", path)
  readRDS(path)
}

load_sediment <- function() {
  path <- file.path(MODEL_RUN_DIR, "3_sediment_ready_for_modeling.rds")
  if (!file.exists(path)) stop("sediment data missing at ", path)
  readRDS(path)
}

# Convert the Stan coefficient (SD wax per SD precipitation) to the physical
# calibration slope (per mil wax per per mil precipitation), and back again.
# Spatial slope perturbations use the same coefficient scale as beta_oipc.
slope_model_to_physical <- function(slope, scaling_params) {
  required <- c("d2H_sd", "oipc_sd")
  if (!all(required %in% names(scaling_params))) {
    stop("scaling_params must contain d2H_sd and oipc_sd")
  }
  slope * scaling_params$d2H_sd / scaling_params$oipc_sd
}

slope_physical_to_model <- function(slope, scaling_params) {
  required <- c("d2H_sd", "oipc_sd")
  if (!all(required %in% names(scaling_params))) {
    stop("scaling_params must contain d2H_sd and oipc_sd")
  }
  slope * scaling_params$oipc_sd / scaling_params$d2H_sd
}
