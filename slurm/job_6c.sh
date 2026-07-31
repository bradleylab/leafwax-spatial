#!/bin/bash
#SBATCH --job-name=leafwax_6c
#SBATCH -A compute2-alexander.s.bradley
#SBATCH --partition=general-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --array=0-6
#SBATCH --output=logs/6c_%A_%a.out
#SBATCH --error=logs/6c_%A_%a.err

set -eo pipefail

WORKDIR=/scratch2/fs1/alexander.s.bradley/leafwax_run
SIF=${WORKDIR}/leafwax-spatial.sif

module load ris
module load apptainer/1.3.6

cd "${WORKDIR}"
mkdir -p logs rtmp

EXPERIMENTS=(4a2 4a3 4a4 4b2 4b3 4c2 4c3)
EXP=${EXPERIMENTS[$SLURM_ARRAY_TASK_ID]}

echo "=== 6c prior sensitivity: experiment ${EXP} ==="
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo "Start: $(date)"

apptainer exec --no-home --containall \
    --bind "${WORKDIR}:${WORKDIR}" \
    --bind "${WORKDIR}/rtmp:/tmp" \
    --pwd "${WORKDIR}" \
    --env TMPDIR=/tmp \
    --env CMDSTAN=/root/.cmdstan/cmdstan-2.36.0 \
    --env OPENBLAS_CORETYPE=Haswell \
    "${SIF}" \
    Rscript 6c_prior_sensitivity.R "${EXP}"

echo "Done: $(date)"
