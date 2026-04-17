#!/bin/bash
#SBATCH --job-name=leafwax_pp
#SBATCH -A compute2-alexander.s.bradley
#SBATCH --partition=general-cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --array=0-13
#SBATCH --output=logs/pp_%A_%a.out
#SBATCH --error=logs/pp_%A_%a.err

set -eo pipefail
module load apptainer

WORKDIR=/scratch2/fs1/alexander.s.bradley/leafwax_run
SIF=${WORKDIR}/leafwax-spatial.sif
cd "${WORKDIR}"

# 14 models — same order/names as slurm/job_fit.sh so array index maps cleanly.
MODELS=(
    baseline baseline_veg baseline_env full full_interact
    baseline_sp baseline_veg_sp baseline_env_sp c4_only_sp elevation_only_sp
    elevation_c4_sp elevation_c4_interact_sp full_sp full_interact_sp
)
MODEL="${MODELS[${SLURM_ARRAY_TASK_ID}]}"

echo "=== Postprocess model: ${MODEL} ==="
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Start: $(date)"

# Per-model R wrapper. Rereads existing fit.rds, extracts widened draws and
# emits posterior_draws.rds (overwrite) + new loo.rds. Matches Phase 5 W1
# logic in 4c_fit_models.R on master.
RUNNER=/tmp/pp_${MODEL}_${SLURM_JOB_ID}.R
cat > "${RUNNER}" <<'REOF'
suppressPackageStartupMessages({
  library(cmdstanr); library(posterior); library(loo)
})
cmdstanr::set_cmdstan_path("/root/.cmdstan/cmdstan-2.36.0")

args <- commandArgs(trailingOnly = TRUE)
model_name <- args[1]

model_dir <- file.path("model_output", model_name)
prep_dir  <- "prepared_data"
stan_data <- readRDS(file.path(prep_dir, paste0("stan_data_", model_name, ".rds")))

fit_file <- file.path(model_dir, "fit.rds")
if (!file.exists(fit_file)) stop("fit.rds missing: ", fit_file)
cat("Loading fit.rds ...\n"); fit <- readRDS(fit_file)

params_to_check <- c("beta_0", "beta_oipc", "sigma",
                     "lambda_decay", "effective_scale_km")
if (isTRUE(stan_data$include_precip == 1))
  params_to_check <- c(params_to_check, "beta_precip")
if (isTRUE(stan_data$include_c4 == 1)) {
  params_to_check <- c(params_to_check, "beta_c4")
  if (isTRUE(stan_data$include_pft == 1))
    params_to_check <- c(params_to_check, "beta_oipc_x_c4")
}
if (isTRUE(stan_data$include_pft == 1)) {
  params_to_check <- c(params_to_check,
                       "beta_tree", "beta_shrub", "beta_grass",
                       "beta_oipc_x_tree", "beta_oipc_x_shrub", "beta_oipc_x_grass")
}
if (isTRUE(stan_data$include_gp == 1)) {
  params_to_check <- c(params_to_check,
                       "ls_intercept_km", "ls_slope_km",
                       "sigma_intercept_spatial", "sigma_slope_spatial")
}

draws_to_save <- c(params_to_check,
                   "mu", "d2H_rep", "log_lik", "scale_weights")
if (isTRUE(stan_data$include_gp == 1)) {
  draws_to_save <- c(draws_to_save,
                     "alpha_spatial", "beta_oipc_spatial",
                     "z_intercept_spatial", "z_slope_spatial")
}

cat("Extracting ", length(draws_to_save), " variable groups as draws_array...\n", sep = "")
draws <- fit$draws(variables = draws_to_save, format = "draws_array")
cat("  dim(draws) = ", paste(dim(draws), collapse = "x"), "\n", sep = "")
saveRDS(draws, file.path(model_dir, "posterior_draws.rds"))
cat("Wrote posterior_draws.rds\n")

cat("Computing loo...\n")
log_lik_array <- posterior::subset_draws(draws, variable = "log_lik")
loo_result <- tryCatch(
  loo::loo(log_lik_array, cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 4))),
  error = function(e) { message("loo failed: ", conditionMessage(e)); NULL }
)
if (!is.null(loo_result)) {
  saveRDS(loo_result, file.path(model_dir, "loo.rds"))
  cat("Wrote loo.rds: elpd_loo = ",
      round(loo_result$estimates["elpd_loo", "Estimate"], 1),
      " (SE = ", round(loo_result$estimates["elpd_loo", "SE"], 1), ")\n", sep = "")
}
cat("Done: ", model_name, "\n")
REOF

apptainer exec --no-home --containall \
    --bind "${WORKDIR}:${WORKDIR}" \
    --bind /tmp:/tmp \
    --pwd "${WORKDIR}" \
    --env CMDSTAN=/root/.cmdstan/cmdstan-2.36.0 \
    "${SIF}" \
    Rscript "${RUNNER}" "${MODEL}"

rm -f "${RUNNER}"
echo "End: $(date)"
