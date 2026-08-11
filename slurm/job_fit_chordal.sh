#!/bin/bash
#SBATCH --job-name=leafwax_fit_chordal
#SBATCH -A compute2-alexander.s.bradley
#SBATCH --partition=general-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=96:00:00
# Node isolation: each model owns its node. Co-scheduled array tasks each spawn
# eight CmdStan chains plus compilation processes and can exceed the per-user
# process limit even when memory use is modest.
#SBATCH --exclusive
#SBATCH --array=0-16
#SBATCH --output=logs/fit_chordal_%A_%a.out
#SBATCH --error=logs/fit_chordal_%A_%a.err

# Complete fitting batch: 9 spatial models using the chordal metric, 5
# non-spatial models, and 3 range-factor sensitivity variants. Every task writes
# a complete posterior artifact set from the prepared data for that model.
#
# Compute2 requirements: `-A compute2-alexander.s.bradley` and apptainer/1.3.6.
# Run prep first (job_prep.sh) into a fresh prepared_data directory so the
# chordal stan_data is not shadowed by a stale cache. After all fits complete,
# snapshot model_output -> results/c2_run_<date>_chordal/ and build the run
# manifest (scripts/build_run_manifest.R).

set -eo pipefail

WORKDIR=/scratch2/fs1/alexander.s.bradley/leafwax_run
SIF=${WORKDIR}/leafwax-spatial.sif

# 17 models. Index 0 is baseline_sp so the pilot (`sbatch --array=0`) is unchanged.
# Order: 9 chordal spatial, then 3 range_factor-off sensitivities, then 5 non-spatial.
MODELS=(
    baseline_sp
    baseline_veg_sp
    baseline_env_sp
    c4_only_sp
    elevation_only_sp
    elevation_c4_sp
    elevation_c4_interact_sp
    full_sp
    full_interact_sp
    baseline_sp_rfoff
    baseline_env_sp_rfoff
    full_interact_sp_rfoff
    baseline
    baseline_veg
    baseline_env
    full
    full_interact
)

MODEL=${MODELS[$SLURM_ARRAY_TASK_ID]}

module load ris
module load apptainer/1.3.6   # 2026-07-15 upgrade: 1.3.4 retired, 1.3.6 is current

cd "${WORKDIR}"
mkdir -p model_output/${MODEL}

# Per-model artifact paths (used by the pre-fit provenance gate and the resume
# guard below).
STAN_DATA="prepared_data/stan_data_${MODEL}.rds"
FIT_RDS="model_output/${MODEL}/fit.rds"
DIAG_RDS="model_output/${MODEL}/diagnostics.rds"
DRAWS_RDS="model_output/${MODEL}/posterior_draws.rds"
DONE_MARK="model_output/${MODEL}/DONE"

if [ ! -f "$STAN_DATA" ]; then
    echo "${MODEL}: ERROR — stan_data missing at $STAN_DATA. Run prep first."
    exit 1
fi

# Pre-fit gate — runs BEFORE the resume/skip guard, so even a task with apparently
# complete artifacts cannot exit success without a completed prep, a valid+current
# fit-time provenance record, and no prep→fit drift. The in-SIF check verifies the
# shared code/config/inputs AND this model's own stan_data/config hashes AND the
# current SIF md5 (computed host-side) against what prep captured.
if [ ! -f results/prep_complete.flag ]; then
    echo "${MODEL}: ERROR — results/prep_complete.flag missing. Run job_prep.sh first."
    exit 1
fi
if [ ! -f results/fit_provenance.rds ]; then
    echo "${MODEL}: ERROR — results/fit_provenance.rds missing. Prep must capture provenance first."
    exit 1
fi
CUR_SIF_MD5=$(md5sum "${SIF}" | awk '{print $1}')
apptainer exec --no-home --containall \
    --bind "${WORKDIR}:${WORKDIR}" --bind /tmp:/tmp --pwd "${WORKDIR}" "${SIF}" \
    Rscript scripts/check_provenance_current.R \
      --prov results/fit_provenance.rds --model "${MODEL}" --current-sif-md5 "${CUR_SIF_MD5}"

# Resume guard: skip ONLY if the COMPLETE per-model artifact set is present and
# newer than its stan_data input. fit.rds alone is not enough — 4c writes fit.rds
# before diagnostics/draws, so a crash between them would otherwise let a re-run
# false-succeed. Runs AFTER the provenance gate above.
STATUS_RDS="model_output/${MODEL}/fit_status.rds"
set_is_complete_and_current() {
    [ -f "$DONE_MARK" ] || return 1
    [ -f "$STATUS_RDS" ] || return 1
    for f in "$FIT_RDS" "$DIAG_RDS" "$DRAWS_RDS"; do
        [ -f "$f" ] || return 1
        [ "$f" -nt "$STAN_DATA" ] || return 1
    done
    return 0
}
if set_is_complete_and_current; then
    echo "${MODEL}: complete artifact set newer than stan_data (provenance verified), skipping"
    exit 0
fi
if [ -f "$FIT_RDS" ]; then
    echo "${MODEL}: fit.rds present but artifact set incomplete/stale — 4c will resume from fit.rds and regenerate downstream artifacts"
fi

# Entering the (re)fit path: invalidate any prior completion state IMMEDIATELY, so
# a crash before 4c's end block cannot leave a stale DONE/fit_status that the
# aggregator (or a later skip) would read as complete. 4c writes a fresh
# fit_status.rds + DONE only after the full artifact set is on disk. The RUNNING
# marker records that a (re)fit was started but not yet completed.
rm -f "$DONE_MARK" "model_output/${MODEL}/fit_status.rds"
date +"%Y-%m-%dT%H:%M:%S%z" > "model_output/${MODEL}/RUNNING"

echo "=== Fitting model (chordal): ${MODEL} ==="
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
