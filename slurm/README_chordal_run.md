# Spatial-model fitting workflow

This runbook describes the cluster workflow for the chordal-distance models.
The batch contains the 14 manuscript comparison models and three
`range_factor` sensitivity fits. Every model is fit from prepared calibration
and environmental data.

Cluster account, partition, filesystem, and container paths in the SLURM scripts
must be adapted to the target system. Model fitting must run on compute nodes,
not on a login node.

## 1. Static checks

Run these checks from the repository root before preparing the batch:

```bash
Rscript scripts/test_chordal_acceptance.R
Rscript scripts/validate_chordal_prep.R
```

The first script checks chordal geometry, range-factor logic, and metric
dispatch. The second exercises the complete Stan-data preparation path against
the audited calibration input.

## 2. Prepare model inputs

Prepare all configured models in a fresh `prepared_data/` directory:

```bash
sbatch slurm/job_prep.sh
```

`4b_stan_prep.R` records a configuration and source-code signature so that a
changed metric or model specification invalidates an older prepared-data cache.

## 3. Pilot fit and validation

Fit `baseline_sp` before launching the full array:

```bash
sbatch --array=0 slurm/job_fit_chordal.sh
```

After the pilot completes, check convergence and the fitted length-scale bounds:

```bash
Rscript scripts/check_pilot.R \
  model_output/baseline_sp \
  prepared_data/stan_data_baseline_sp.rds
```

The pilot check requires finite diagnostics, R-hat below 1.01, effective sample
sizes above 400, no divergent transitions or maximum-treedepth events, adequate
E-BFMI, and no posterior pile-up at the configured length-scale bounds.

Verify that the R implementation of the spatial prediction reproduces the Stan
generated quantities:

```bash
apptainer exec --no-home --containall \
  --bind <project-parent>:/project \
  --pwd /project/leafwax-spatial \
  --env OPENBLAS_CORETYPE=Haswell \
  <container.sif> \
  Rscript scripts/check_stan_r_parity.R \
    model_output/baseline_sp \
    prepared_data/stan_data_baseline_sp.rds
```

The parity check sources the spatial interpolation implementation from the
companion `leafwax` package, so the bind mount must include both repositories.

### OpenBLAS compatibility

The production container uses an OpenBLAS version that predates Intel Sapphire
Rapids. On those processors, the automatically selected AVX-512 kernel can
produce incorrect Cholesky and triangular-solve results for the knot covariance
matrix. Set `OPENBLAS_CORETYPE=Haswell` before R starts for every containerized R
step that performs Gaussian-process prediction. Stan fitting is unaffected
because Stan uses Eigen rather than the R-linked OpenBLAS library.

## 4. Fit the complete batch

After the pilot gates pass, launch all configured models:

```bash
sbatch slurm/job_fit_chordal.sh
```

Each task writes `fit.rds`, `diagnostics.rds`, `posterior_draws.rds`, and a
`fit_status.rds` record under `model_output/<model>/`. A `DONE` marker is written
only when the complete artifact set exists.

Check batch completeness with:

```bash
Rscript scripts/aggregate_chordal_status.R model_output prepared_data
```

The aggregator exits nonzero unless every configured model has a complete,
fresh artifact set.

## 5. Preserve provenance

The preparation job runs `scripts/capture_fit_provenance.R` inside the container.
Preserve its output with the fitted models and prepared data. Build the run
manifest after copying the completed run to its archival directory:

```bash
Rscript scripts/build_run_manifest.R \
  --chordal-output results/<chordal-run>/model_output \
  --chordal-prepared results/<chordal-run>/_prepared_data \
  --fit-provenance results/<chordal-run>/fit_provenance.rds \
  --out results/<chordal-run>/RUN_MANIFEST_chordal
```

## 6. Regenerate reported quantities

Select the completed run with `LEAFWAX_RUN_DIR`, then run the deterministic
post-fit scripts listed in `NUMBERS.md`. Posterior files distributed through the
companion data archive must retain
`attr(posterior, "spatial_metric") = "chordal"`.
