#!/bin/bash
#SBATCH --job-name=leafwax_6b
#SBATCH -A compute2-alexander.s.bradley
#SBATCH --partition=general-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --array=0-3
#SBATCH --output=logs/6b_%A_%a.out
#SBATCH --error=logs/6b_%A_%a.err

set -eo pipefail

WORKDIR=/scratch2/fs1/alexander.s.bradley/leafwax_run
SIF=${WORKDIR}/leafwax-spatial.sif

module load ris
module load apptainer/1.3.4

cd "${WORKDIR}"
mkdir -p logs rtmp

SCENARIOS=(rho00 rho03 rho05 empirical)
SCEN=${SCENARIOS[$SLURM_ARRAY_TASK_ID]}

echo "=== 6b spatial confounding sim: scenario ${SCEN} ==="
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo "Start: $(date)"

apptainer exec --no-home --containall \
    --bind "${WORKDIR}:${WORKDIR}" \
    --bind "${WORKDIR}/rtmp:/tmp" \
    --pwd "${WORKDIR}" \
    --env TMPDIR=/tmp \
    --env CMDSTAN=/root/.cmdstan/cmdstan-2.36.0 \
    "${SIF}" \
    Rscript 6b_spatial_confounding_simulation.R "${SCEN}"

echo "Done: $(date)"
