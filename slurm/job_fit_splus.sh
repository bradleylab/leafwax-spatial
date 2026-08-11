#!/bin/bash
#SBATCH --job-name=leafwax_splus
#SBATCH -A compute2-alexander.s.bradley
#SBATCH --partition=general-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=96:00:00
#SBATCH --array=0-3
#SBATCH --output=logs/splus_%A_%a.out
#SBATCH --error=logs/splus_%A_%a.err

# Spatial+ refits (Dupont, Wood & Augustin 2022) for the four spatial
# models that ground the slope claims in §3 / §4 of the manuscript.
# Each task uses the same compiled Stan model as the main fits; only
# the predictor matrices in stan_data are residualized against
# coordinates (per 3i_spatial_plus_residuals.R). Walltime mirrors the main
# chordal fitting array; full_interact_sp_splus is the slowest expected element.

set -eo pipefail

WORKDIR=/scratch2/fs1/alexander.s.bradley/leafwax_run
SIF=${WORKDIR}/leafwax-spatial.sif

# Order matches the four base spatial models that ground the headline
# slope claim. Subset, not all 14 — see manuscript discussion of
# Spatial+ as a confounding-bias bound, not a wholesale model swap.
MODELS=(
    baseline_sp_splus
    baseline_env_sp_splus
    full_sp_splus
    full_interact_sp_splus
)

MODEL=${MODELS[$SLURM_ARRAY_TASK_ID]}

module load ris
module load apptainer/1.3.6   # bumped from 1.3.4 after the 2026-07-15 C2 upgrade

cd "${WORKDIR}"
mkdir -p logs model_output/${MODEL}

# Resume guard mirrors the main fitting array: only skip if fit.rds is newer than
# its stan_data, otherwise a freshly residualized stan_data would be
# silently shadowed by a stale fit.
FIT_RDS="model_output/${MODEL}/fit.rds"
STAN_DATA="prepared_data/stan_data_${MODEL}.rds"
if [ -f "$FIT_RDS" ] && [ -f "$STAN_DATA" ] && [ "$FIT_RDS" -nt "$STAN_DATA" ]; then
    echo "${MODEL}: fit.rds newer than stan_data, skipping"
    exit 0
fi
if [ -f "$FIT_RDS" ]; then
    echo "${MODEL}: fit.rds exists but older than (or has no) stan_data — refitting"
fi
if [ ! -f "$STAN_DATA" ]; then
    echo "${MODEL}: ERROR — stan_data missing at $STAN_DATA."
    echo "  Run 3i_spatial_plus_residuals.R first to generate splus stan_data."
    exit 1
fi

echo "=== Fitting Spatial+ model: ${MODEL} ==="
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Start: $(date)"

# Per-task R wrapper. Set model_names before sourcing 4c so the auto-discovery
# before sourcing 4c so the auto-discovery short-circuits to just this
# one model.
RUNNER=/tmp/fit_${MODEL}_${SLURM_JOB_ID}.R
cat > "${RUNNER}" << REOF
cmdstanr::set_cmdstan_path("/root/.cmdstan/cmdstan-2.36.0")
source("0_load_config.R")
model_names <- c("${MODEL}")
source(CONFIG\$scripts\$fit_script)
REOF

apptainer exec --no-home --containall \
    --bind "${WORKDIR}:${WORKDIR}" \
    --bind /tmp:/tmp \
    --pwd "${WORKDIR}" \
    --env CMDSTAN=/root/.cmdstan/cmdstan-2.36.0 \
    "${SIF}" \
    Rscript "${RUNNER}"

rm -f "${RUNNER}"
echo "Done: $(date)"
