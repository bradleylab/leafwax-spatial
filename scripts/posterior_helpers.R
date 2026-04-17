# posterior_helpers.R — root-level twin of manuscript/figure_code/posterior_helpers.R
# with an APRIL_RUN path resolved from the repo root (one level up from scripts/).
#
# Used by extract_*.R, analyze_*.R, and 5*_*.R when they live at the repo root.
# Keep the two helper files in sync.

library(posterior)

APRIL_RUN <- normalizePath(
  file.path("..", "results", "c2_run_20260414"),
  mustWork = FALSE
)
# When sourced from the repo root (not from scripts/), fall back.
if (!dir.exists(APRIL_RUN)) {
  APRIL_RUN <- normalizePath(
    file.path("results", "c2_run_20260414"),
    mustWork = FALSE
  )
}
PREPARED_DATA <- file.path(APRIL_RUN, "_prepared_data")

load_draws <- function(model) {
  path <- file.path(APRIL_RUN, model, "posterior_draws.rds")
  if (!file.exists(path)) {
    stop("posterior_draws.rds missing for '", model, "' at ", path,
         "\nRun 4c_fit_models.R on C2 against this model's output dir ",
         "(Phase 5 W2).")
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
         "\nRun 4c_fit_models.R on C2 against this model's output dir ",
         "(Phase 5 W2).")
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
