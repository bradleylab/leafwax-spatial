#!/bin/bash
#SBATCH --job-name=leafwax_fit
#SBATCH -A compute2-alexander.s.bradley
#SBATCH --partition=general-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=96:00:00
#SBATCH --array=0-13
#SBATCH --output=logs/fit_%A_%a.out
#SBATCH --error=logs/fit_%A_%a.err

set -eo pipefail

WORKDIR=/scratch2/fs1/alexander.s.bradley/leafwax_run
SIF=${WORKDIR}/leafwax-spatial.sif

# 14 models — order matches config.yaml
MODELS=(
    baseline
    baseline_veg
    baseline_env
    full
    full_interact
    baseline_sp
    baseline_veg_sp
    baseline_env_sp
    c4_only_sp
    elevation_only_sp
    elevation_c4_sp
    elevation_c4_interact_sp
    full_sp
    full_interact_sp
)

MODEL=${MODELS[$SLURM_ARRAY_TASK_ID]}

module load ris
module load apptainer/1.3.4

cd "${WORKDIR}"
mkdir -p model_output/${MODEL}

# Resume guard: skip only if fit.rds is newer than its stan_data input,
# so a stale fit can't silently shadow a freshly rebuilt stan_data.
FIT_RDS="model_output/${MODEL}/fit.rds"
STAN_DATA="prepared_data/stan_data_${MODEL}.rds"
if [ -f "$FIT_RDS" ] && [ -f "$STAN_DATA" ] && [ "$FIT_RDS" -nt "$STAN_DATA" ]; then
    echo "${MODEL}: fit.rds newer than stan_data, skipping"
    exit 0
fi
if [ -f "$FIT_RDS" ]; then
    echo "${MODEL}: fit.rds exists but is older than (or has no) stan_data — refitting"
fi
if [ ! -f "$STAN_DATA" ]; then
    echo "${MODEL}: ERROR — stan_data missing at $STAN_DATA. Run prep first."
    exit 1
fi

echo "=== Fitting model: ${MODEL} ==="
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Start: $(date)"

# Write per-model R wrapper
RUNNER=/tmp/fit_${MODEL}_${SLURM_JOB_ID}.R
cat > "${RUNNER}" << REOF
cmdstanr::set_cmdstan_path("/root/.cmdstan/cmdstan-2.36.0")
source("0_load_config.R")
model_names <- c("${MODEL}")
source(CONFIG\$scripts\$fit_script)
REOF

# --no-home prevents Apptainer from bind-mounting the user's home over /root,
# which would hide /root/.cmdstan/ (where CmdStan is installed in the image).
# --containall isolates the container filesystem fully; we bind-mount only
# the working directory and /tmp (for the runner script).
apptainer exec --no-home --containall \
    --bind "${WORKDIR}:${WORKDIR}" \
    --bind /tmp:/tmp \
    --pwd "${WORKDIR}" \
    --env CMDSTAN=/root/.cmdstan/cmdstan-2.36.0 \
    "${SIF}" \
    Rscript "${RUNNER}"

rm -f "${RUNNER}"
echo "Done: $(date)"
