#!/bin/bash
#───────────────────────────────────────────────────────────────────────────────
# launch_pipeline.sh
#
# Launches the leafwax spatial modeling pipeline with parallel model fitting
# and safe auto-shutdown on EC2. Auto-shutdown requires ALL models to
# complete with no errors; partial failure keeps the instance running.
#
# Usage:
#   ./launch_pipeline.sh               # full pipeline (all models)
#   ./launch_pipeline.sh baseline_veg_sp  # single model (validation)
#───────────────────────────────────────────────────────────────────────────────

set -eo pipefail

# Ensure yq is available in non-login shells (EC2 installs to ~/bin)
export PATH="/home/ubuntu/bin:$PATH"

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_DIR="logs/run_${TIMESTAMP}"
mkdir -p "$LOG_DIR" prepared_data model_output results

echo "======================================================================"
echo "LEAF WAX d2H SPATIAL MODELING PIPELINE"
echo "======================================================================"
echo "Started at: $(date)"
echo "Log directory: $LOG_DIR"
echo

#───────────────────────────────────────────────────────────────────────────────
# Config parsing
#───────────────────────────────────────────────────────────────────────────────
if [ ! -f "config.yaml" ]; then
    echo "ERROR: config.yaml not found." >&2
    exit 1
fi

# Requires yq for YAML parsing (nested keys like scripts.prep_script)
if ! command -v yq &>/dev/null; then
    echo "ERROR: yq not found. Install yq or add to PATH." >&2
    exit 1
fi

get_config() {
    yq ".${1}" config.yaml 2>/dev/null | tr -d '"'
}

PREP_SCRIPT=$(get_config 'scripts.prep_script')
FIT_SCRIPT=$(get_config 'scripts.fit_script')
INPUT_FILE=$(get_config 'input_data')
MAX_PARALLEL=$(get_config 'max_parallel_models')
MAX_PARALLEL=${MAX_PARALLEL:-4}

if [ -z "$PREP_SCRIPT" ] || [ -z "$FIT_SCRIPT" ]; then
    echo "ERROR: Could not read scripts from config.yaml" >&2
    exit 1
fi

#───────────────────────────────────────────────────────────────────────────────
# Determine which models to run
#───────────────────────────────────────────────────────────────────────────────
if [ $# -gt 0 ]; then
    # Single model mode (for validation runs)
    MODELS=("$@")
    echo "Single-model mode: ${MODELS[*]}"
else
    # All models from config (yq required, checked above)
    MODELS=($(yq '.model_configs | keys | .[]' config.yaml 2>/dev/null))
    echo "Full pipeline: ${#MODELS[@]} models"
fi

EXPECTED_MODELS=${#MODELS[@]}
echo "Expected models: $EXPECTED_MODELS"
echo "Models: ${MODELS[*]}"
echo

#───────────────────────────────────────────────────────────────────────────────
# Step 1: Data preparation (sequential, sources 4a_spatial_functions.R)
#───────────────────────────────────────────────────────────────────────────────
echo "======================================================================"
echo "STEP 1: DATA PREPARATION"
echo "======================================================================"

PREP_LOG="$LOG_DIR/data_prep.log"

# Process substitution preserves Rscript exit code through tee
Rscript "$PREP_SCRIPT" > >(tee "$PREP_LOG") 2>&1
echo "✓ Data preparation completed (log: $PREP_LOG)"
echo

#───────────────────────────────────────────────────────────────────────────────
# Step 2: Model fitting (parallel, one Rscript per model)
#───────────────────────────────────────────────────────────────────────────────
echo "======================================================================"
echo "STEP 2: MODEL FITTING ($EXPECTED_MODELS models, max $MAX_PARALLEL parallel)"
echo "======================================================================"

# Generate per-model runner scripts
for model in "${MODELS[@]}"; do
    cat > "run_${model}.R" << REOF
source("0_load_config.R")
model_names <- c("$model")
source(CONFIG\$scripts\$fit_script)
REOF
done

# Launch models with controlled parallelism
PIDS=()
for model in "${MODELS[@]}"; do
    # Skip if already fitted
    if [ -f "model_output/${model}/fit.rds" ]; then
        echo "  ✓ $model already completed — skipping"
        continue
    fi

    # Wait for a slot
    while [ "$(jobs -rp | wc -l)" -ge "$MAX_PARALLEL" ]; do
        sleep 30
    done

    MODEL_LOG="$LOG_DIR/${model}.log"
    echo "  Starting $model (log: $MODEL_LOG)"
    Rscript "run_${model}.R" > >(tee "$MODEL_LOG") 2>&1 &
    PIDS+=($!)
done

# Wait for all background jobs
echo
echo "Waiting for all models to finish..."
FAILED=0
for pid in "${PIDS[@]}"; do
    if ! wait "$pid"; then
        FAILED=$((FAILED + 1))
    fi
done

echo
if [ "$FAILED" -gt 0 ]; then
    echo "✗ $FAILED model job(s) returned non-zero exit code"
else
    echo "✓ All model jobs returned successfully"
fi

# Clean up runner scripts
rm -f run_*.R fit_single_model.R

#───────────────────────────────────────────────────────────────────────────────
# Step 3: Check completion and decide on shutdown
#───────────────────────────────────────────────────────────────────────────────
echo
echo "======================================================================"
echo "STEP 3: COMPLETION CHECK"
echo "======================================================================"

COMPLETE_MARKER="results/pipeline_4c_complete.rds"
INCOMPLETE_MARKER="results/pipeline_4c_INCOMPLETE.rds"
ERROR_FILES=$(find model_output -name "error_info.rds" 2>/dev/null | wc -l)
FIT_FILES=$(find model_output -name "fit.rds" 2>/dev/null | wc -l)

echo "  Expected models: $EXPECTED_MODELS"
echo "  Fit files found: $FIT_FILES"
echo "  Error files found: $ERROR_FILES"

if [ -f "$COMPLETE_MARKER" ] && [ "$ERROR_FILES" -eq 0 ]; then
    echo
    echo "✓ ALL MODELS COMPLETED SUCCESSFULLY"
    echo "  Pipeline finished at: $(date)"
    echo "  Total fit files: $FIT_FILES"
    echo
    echo "  Auto-shutdown in 5 minutes..."
    echo "  (Cancel with: sudo shutdown -c)"
    sudo shutdown -h +5 "Pipeline complete — auto-shutdown"
elif [ -f "$INCOMPLETE_MARKER" ]; then
    echo
    echo "✗ INCOMPLETE: Some models failed."
    echo "  Check: $INCOMPLETE_MARKER"
    echo "  Error files:"
    find model_output -name "error_info.rds" -exec echo "    {}" \;
    echo
    echo "  Instance will NOT auto-shutdown. Investigate and re-run."
else
    echo
    echo "✗ No completion marker found."
    echo "  This may mean the pipeline was interrupted before the fitting"
    echo "  loop finished writing the summary."
    echo
    echo "  Instance will NOT auto-shutdown. Investigate and re-run."
fi

echo
echo "Pipeline script finished at: $(date)"
