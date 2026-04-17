# posterior_helpers.R — thin wrappers over the April 2026 run's rds outputs.
#
# After Phase 5 W1 and W2, every <model>/ directory under
# results/c2_run_20260414/ contains a self-contained rds bundle:
#   posterior_draws.rds   — draws_array with all variables downstream reads
#   diagnostics.rds       — summaries (rhat/ess/quantiles) for every parameter
#   loo.rds               — psis_loo object for model comparison
#
# Analysis code must read ONLY those files (and _prepared_data/*.rds for
# stan_data / config). No fit.rds, no chain CSVs. See CLAUDE.md for the
# reproducibility contract.

library(posterior)

# Repo root is two levels up from manuscript/figure_code/.
APRIL_RUN <- normalizePath(
  file.path("..", "..", "results", "c2_run_20260414"),
  mustWork = FALSE
)
PREPARED_DATA <- file.path(APRIL_RUN, "_prepared_data")

load_draws <- function(model) {
  path <- file.path(APRIL_RUN, model, "posterior_draws.rds")
  if (!file.exists(path)) {
    stop("posterior_draws.rds missing for model '", model, "' at ", path,
         "\nRun 4c_fit_models.R on C2 against this model's output dir ",
         "(Phase 5 W2) to regenerate the widened draws.")
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
