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

INPUT_CSV="input_data/global_data_c29.csv"
PREP_RDS="results/3_sediment_ready_for_modeling.rds"
PREP_SIDECAR="results/3_sediment_ready_for_modeling.input_md5"

# ── Prerequisite raster checks ─────────────────────────────────────────────
# 3_prep_data.R reads C4, MODIS-PFT, and TerraClimate rasters. C4 is built by
# step 1 below. MODIS-PFT (2d_) and TerraClimate (2f_) are produced from large
# downloads on EC2/local and must be pre-staged into ${WORKDIR}/results/.
required_rasters=(
    "results/2d_MODIS_PFT_3classes_Downsampled.tif"
    "results/2f_TerraClimate_ppt_mean_2001_2019.tif"
    "results/2f_TerraClimate_soil_mean_2001_2019.tif"
    "results/2f_TerraClimate_tmax_mean_2001_2019.tif"
    "results/2f_TerraClimate_vpd_mean_2001_2019.tif"
)
missing_rasters=()
for r in "${required_rasters[@]}"; do
    [ -f "$r" ] || missing_rasters+=("$r")
done
if [ ${#missing_rasters[@]} -gt 0 ]; then
    echo "ERROR: Required env rasters missing in ${WORKDIR}/results/:"
    printf '  - %s\n' "${missing_rasters[@]}"
    echo ""
    echo "Pre-stage them from local input_data/ via scp before running prep."
    echo "See README and ec2_rasters_location memory for raster sources."
    exit 1
fi
echo "✓ All required env rasters present"

# ── Input CSV md5-aware skip guard ─────────────────────────────────────────
# Skip prep only if the cached prep RDS was produced from the *same* input
# CSV. Previously this only checked file existence, which silently reused v8
# prep when v10 input was loaded into a stale WORKDIR.
if [ ! -f "$INPUT_CSV" ]; then
    echo "ERROR: Input CSV missing: ${WORKDIR}/${INPUT_CSV}"
    exit 1
fi
current_md5=$(md5sum "$INPUT_CSV" | awk '{print $1}')

prep_is_current=false
if [ -f "$PREP_RDS" ] && [ -f "$PREP_SIDECAR" ]; then
    cached_md5=$(awk -F': ' '/^input_md5:/ {print $2}' "$PREP_SIDECAR" | tr -d '[:space:]')
    if [ "$cached_md5" = "$current_md5" ]; then
        prep_is_current=true
    fi
fi

if $prep_is_current; then
    echo "Step-3 RDS up to date (input md5 matches) — skipping steps 1+3"
    echo "  cached md5:  $cached_md5"
else
    if [ -f "$PREP_RDS" ] && [ ! -f "$PREP_SIDECAR" ]; then
        echo "WARNING: Step-3 RDS exists but no sidecar — assuming stale. Re-running prep."
    elif [ -f "$PREP_RDS" ]; then
        echo "Step-3 RDS sidecar md5 mismatch — re-running prep."
        echo "  cached:  ${cached_md5:-<none>}"
        echo "  current: $current_md5"
    fi

    echo "Running step 1: C4 raster extraction"
    apptainer exec --no-home --containall --bind "${WORKDIR}:${WORKDIR}" --bind /tmp:/tmp --pwd "${WORKDIR}" "${SIF}" \
        Rscript 1_extract_c4_raster.R

    echo "Running step 3: data preparation (includes resample fix)"
    apptainer exec --no-home --containall --bind "${WORKDIR}:${WORKDIR}" --bind /tmp:/tmp --pwd "${WORKDIR}" "${SIF}" \
        Rscript 3_prep_data.R

    # Verify 3_prep_data.R wrote the sidecar
    if [ ! -f "$PREP_SIDECAR" ]; then
        echo "ERROR: 3_prep_data.R did not write input md5 sidecar at $PREP_SIDECAR"
        echo "       Downstream cache guards will be unsafe — failing fast."
        exit 1
    fi
fi

echo "Running step 4b: Stan data preparation for all models"
apptainer exec --no-home --containall --bind "${WORKDIR}:${WORKDIR}" --bind /tmp:/tmp --pwd "${WORKDIR}" "${SIF}" \
    Rscript 4b_stan_prep.R

echo "Verifying prepared data:"
ls -lh prepared_data/stan_data_*.rds | wc -l
ls -lh prepared_data/stan_data_*.rds

touch results/prep_complete.flag
echo "Prep complete at: $(date)"
