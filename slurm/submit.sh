#!/bin/bash
# DEPRECATED for the chordal run. This launches the OLD standardized-metric
# 14-model workflow (job_fit.sh). The chordal refit (spec v2 §3) uses a DIFFERENT
# path — see slurm/README_chordal_run.md:
#     sbatch slurm/job_prep.sh
#     sbatch --array=0 slurm/job_fit_chordal.sh    # pilot gate first
#     sbatch slurm/job_fit_chordal.sh              # full 12-fit array
#     Rscript scripts/aggregate_chordal_status.R   # batch completion gate
# Running this script for the chordal run would fit the wrong (standardized)
# models. It is retained only to reproduce the frozen standardized run; it
# refuses to run without an explicit acknowledgment.
#
# Original usage (legacy standardized workflow):
#   bash slurm/submit.sh --legacy             # full pipeline (prep + 14 models)
#   bash slurm/submit.sh --legacy --fit-only  # skip prep, just submit fitting array
#   bash slurm/submit.sh --legacy --test 6    # single model test (array index 6)

set -eo pipefail

if [ "$1" != "--legacy" ]; then
    echo "ERROR: submit.sh is the DEPRECATED standardized-metric workflow." >&2
    echo "For the chordal run use job_prep.sh + job_fit_chordal.sh (see" >&2
    echo "slurm/README_chordal_run.md). To intentionally run the legacy" >&2
    echo "standardized 14-model workflow, pass --legacy as the first argument." >&2
    exit 1
fi
shift   # consume --legacy; remaining args ($1/$2) keep their original meaning

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
