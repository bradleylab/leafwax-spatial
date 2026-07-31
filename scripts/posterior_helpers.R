# posterior_helpers.R — root-level twin of manuscript/figure_code/posterior_helpers.R
# with an APRIL_RUN path resolved from the repo root (one level up from scripts/).
#
# Sourced by 5*_*.R and 7_paleo_inversion.R. Keep this file and the manuscript
# twin in sync.
#
# APRIL_RUN selection (in priority order):
#   1. Environment variable LEAFWAX_RUN_DIR (absolute or relative path)
#   2. results/c2_run_20260626/model_output   (default frozen refit)
#
# The variable is still named APRIL_RUN for backward compatibility with 5a-5e
# callers; semantically it is "the C2 run mirror to read posteriors from".
# Override at the shell:
#   LEAFWAX_RUN_DIR=results/c2_run_20260626/model_output Rscript 5a_model_validation.R

library(posterior)

DEFAULT_RUN_NAME <- file.path("c2_run_20260626", "model_output")

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

APRIL_RUN <- .resolve_run_dir()
PREPARED_DATA <- file.path(APRIL_RUN, "_prepared_data")

load_draws <- function(model) {
  path <- file.path(APRIL_RUN, model, "posterior_draws.rds")
  if (!file.exists(path)) {
    stop("posterior_draws.rds missing for '", model, "' at ", path,
         "\nRun 4c_fit_models.R on C2 against this model's output dir, ",
         "or set LEAFWAX_RUN_DIR to the correct mirror.")
  }
  readRDS(path)
}

load_summaries <- function(model) {
  path <- file.path(APRIL_RUN, model, "diagnostics.rds")
  if (!file.exists(path)) stop("diagnostics.rds missing at ", path)
  readRDS(path)$all_params_summary
}

load_sampler_diag <- function(model) {
  path <- file.path(APRIL_RUN, model, "diagnostics.rds")
  if (!file.exists(path)) stop("diagnostics.rds missing at ", path)
  readRDS(path)$summary
}

load_loo <- function(model) {
  path <- file.path(APRIL_RUN, model, "loo.rds")
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
  path <- file.path(APRIL_RUN, "3_sediment_ready_for_modeling.rds")
  if (!file.exists(path)) stop("sediment data missing at ", path)
  readRDS(path)
}
