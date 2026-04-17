#!/bin/bash
# Submit the full pipeline: prep job → 14 fitting array jobs
# Run from: /scratch2/fs1/alexander.s.bradley/leafwax_run/
#
# Usage:
#   bash slurm/submit.sh           # full pipeline (prep + 14 models)
#   bash slurm/submit.sh --fit-only  # skip prep, just submit fitting array
#   bash slurm/submit.sh --test 6    # single model test (array index 6 = baseline_veg_sp)

set -eo pipefail
cd /scratch2/fs1/alexander.s.bradley/leafwax_run

if [ "$1" = "--fit-only" ]; then
    echo "Submitting fitting array only (no prep)..."
    JOB_B=$(sbatch --parsable slurm/job_fit.sh)
    echo "Fitting array: ${JOB_B}"
elif [ "$1" = "--test" ]; then
    IDX=${2:-6}
    echo "Submitting single test model (array index ${IDX})..."
    JOB_T=$(sbatch --parsable --array=${IDX} slurm/job_fit.sh)
    echo "Test job: ${JOB_T}"
else
    echo "Submitting prep job..."
    JOB_A=$(sbatch --parsable slurm/job_prep.sh)
    echo "Prep job: ${JOB_A}"

    echo "Submitting fitting array (depends on prep)..."
    JOB_B=$(sbatch --parsable --dependency=afterok:${JOB_A} slurm/job_fit.sh)
    echo "Fitting array: ${JOB_B}"
fi

echo ""
echo "Monitor:"
echo "  squeue -u alexander.s.bradley"
echo "  tail -f logs/fit_*_*.out"
