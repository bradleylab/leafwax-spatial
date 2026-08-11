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

module load ris
module load apptainer/1.3.6

cd "${WORKDIR}"
mkdir -p logs prepared_data model_output results

# Clear any stale success marker at the START of a prep attempt, so a failed
# rerun can never leave a previous run's prep_complete.flag standing. It is
# re-written only after the completeness assertions pass at the end.
rm -f results/prep_complete.flag

# 3_prep_data.R reads the audited calibration (not the source compilation);
# the skip guard must watch that same file so a stale prep can't be reused.
INPUT_CSV="input_data/leafwax_d2h_c29_calibration_v1.csv"
PREP_RDS="results/3_sediment_ready_for_modeling.rds"
PREP_SIDECAR="results/3_sediment_ready_for_modeling.input_md5"

# ── Input CSV md5-aware skip guard ─────────────────────────────────────────
# Skip prep only if the cached prep RDS was produced from the *same* input
# CSV rather than relying on file existence alone.
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

    # ── Prerequisite raster checks (only when steps 1+3 will actually run) ──────
    # 3_prep_data.R reads C4, MODIS-PFT, and TerraClimate rasters. This guard sits
    # inside the rebuild branch: when the cached step-3 RDS is current (input md5
    # matches), steps 1+3 are skipped and these rasters are never read, so requiring
    # them then would abort a valid reuse of the cached model-ready data for inputs it
    # does not touch. C4 is built by step 1 below; MODIS-PFT (2d_) and TerraClimate
    # (2f_) are produced from large downloads on EC2/local and must be pre-staged
    # into ${WORKDIR}/results/ before a rebuild.
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
        echo "See README.md and slurm/README_chordal_run.md for raster sources and staging."
        exit 1
    fi
    echo "✓ All required env rasters present"

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
# 4b_stan_prep.R exits nonzero if any model fails to prepare or a stan_data file
# is missing; `set -e` above aborts here before the completeness assertion below.
apptainer exec --no-home --containall --bind "${WORKDIR}:${WORKDIR}" --bind /tmp:/tmp --pwd "${WORKDIR}" "${SIF}" \
    Rscript 4b_stan_prep.R

echo "Verifying prepared data:"
ls -lh prepared_data/stan_data_*.rds

# ── Completeness assertion (fail closed) ────────────────────────────────────
# Do not trust a bare count. Assert a stan_data_<model>.rds AND a config_<model>.rds
# exists for EVERY model in the expected set 4b wrote, and that no per-model prep
# error files were left behind. Only then is the flag safe to write.
#
# NB: counts use `find`, not `ls glob | wc -l`. Under `set -eo pipefail` the
# no-match case makes `ls error_*.rds` exit nonzero, and pipefail propagates that
# through the command substitution, aborting the (successful) zero-error path
# before the flag is written. `find` returns 0 rows and exit 0 when nothing matches.
EXPECTED_LIST="prepared_data/expected_models.txt"
if [ ! -f "$EXPECTED_LIST" ]; then
    echo "ERROR: $EXPECTED_LIST missing — 4b_stan_prep.R did not record the expected model set."
    exit 1
fi
missing=()
while IFS= read -r m; do
    [ -z "$m" ] && continue
    if [ ! -f "prepared_data/stan_data_${m}.rds" ] || [ ! -f "prepared_data/config_${m}.rds" ]; then
        missing+=("$m")
    fi
done < "$EXPECTED_LIST"
n_expected=$(grep -cve '^[[:space:]]*$' "$EXPECTED_LIST")
n_present=$(find prepared_data -maxdepth 1 -name 'stan_data_*.rds' | wc -l | tr -d ' ')
err_files=$(find prepared_data -maxdepth 1 -name 'error_*.rds' | wc -l | tr -d ' ')
if [ ${#missing[@]} -gt 0 ] || [ "$err_files" -gt 0 ]; then
    echo "ERROR: prep incomplete — NOT writing prep_complete.flag."
    [ ${#missing[@]} -gt 0 ] && printf '  missing stan_data/config: %s\n' "${missing[@]}"
    [ "$err_files" -gt 0 ] && echo "  $err_files per-model error file(s) in prepared_data/"
    exit 1
fi
echo "✓ All $n_expected expected datasets present ($n_present stan_data files, 0 errors)"

# ── Capture fit-time provenance from INSIDE the SIF (before the fit array) ───────
# The authoritative manifest requires this record; it captures the SIF's actual
# software versions + code/input hashes, which cannot be reconstructed confidently
# after the workspace changes. SIF checksum + git state are host-side facts.
echo "Capturing fit-time provenance from inside the SIF..."
SIF_MD5=$(md5sum "${SIF}" | awk '{print $1}')
GIT_HEAD=$(git -C "${WORKDIR}" rev-parse HEAD 2>/dev/null || echo NA)
if git -C "${WORKDIR}" rev-parse --git-dir >/dev/null 2>&1; then
    [ -n "$(git -C "${WORKDIR}" status --porcelain 2>/dev/null)" ] && GIT_DIRTY=1 || GIT_DIRTY=0
else
    GIT_DIRTY=NA
fi
# --env CMDSTAN so cmdstan_version() resolves inside --containall (matches the fit
# runner). capture_fit_provenance.R fails closed on any missing field, so `set -e`
# aborts here (before the flag) if the record would be incomplete.
apptainer exec --no-home --containall \
    --bind "${WORKDIR}:${WORKDIR}" --bind /tmp:/tmp --pwd "${WORKDIR}" \
    --env CMDSTAN=/root/.cmdstan/cmdstan-2.36.0 "${SIF}" \
    Rscript scripts/capture_fit_provenance.R \
      --out results/fit_provenance.rds \
      --sif-md5 "$SIF_MD5" --sif-path "${SIF}" \
      --git-head "$GIT_HEAD" --git-dirty "$GIT_DIRTY"
if [ ! -f results/fit_provenance.rds ]; then
    echo "ERROR: fit provenance capture failed — NOT writing prep_complete.flag."
    exit 1
fi
echo "✓ Fit-time provenance captured (results/fit_provenance.rds)"

touch results/prep_complete.flag
echo "Prep complete at: $(date)"
