#!/bin/bash
#SBATCH --job-name=leafwax_prep
#SBATCH -A compute2-alexander.s.bradley
#SBATCH --partition=general-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=logs/prep_%j.out
#SBATCH --error=logs/prep_%j.err

set -eo pipefail

WORKDIR=/scratch2/fs1/alexander.s.bradley/leafwax_run
SIF=${WORKDIR}/leafwax-spatial.sif

module load apptainer

cd "${WORKDIR}"
mkdir -p logs prepared_data model_output results

# If step-3 RDS already exists (transferred from EC2), skip steps 1+3
if [ -f "results/3_sediment_ready_for_modeling.rds" ]; then
    echo "Step-3 RDS found — skipping steps 1+3"
else
    echo "Running step 1: C4 raster extraction"
    apptainer exec --no-home --containall --bind "${WORKDIR}:${WORKDIR}" --bind /tmp:/tmp --pwd "${WORKDIR}" "${SIF}" \
        Rscript 1_extract_c4_raster.R

    echo "Running step 3: data preparation (includes resample fix)"
    apptainer exec --no-home --containall --bind "${WORKDIR}:${WORKDIR}" --bind /tmp:/tmp --pwd "${WORKDIR}" "${SIF}" \
        Rscript 3_prep_data.R
fi

echo "Running step 4b: Stan data preparation for all models"
apptainer exec --no-home --containall --bind "${WORKDIR}:${WORKDIR}" --bind /tmp:/tmp --pwd "${WORKDIR}" "${SIF}" \
    Rscript 4b_stan_prep.R

echo "Verifying prepared data:"
ls -lh prepared_data/stan_data_*.rds | wc -l
ls -lh prepared_data/stan_data_*.rds

touch results/prep_complete.flag
echo "Prep complete at: $(date)"
