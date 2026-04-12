# Leaf wax hydrogen isotope spatial calibration

Bayesian spatial model for calibrating n-alkane (n-C29) hydrogen isotope ratios
(d2H) in lake surface sediments against precipitation isotope composition and
environmental covariates. Companion code to Bradley (2026), submitted to
Geochimica et Cosmochimica Acta.

## What this does

The model relates leaf wax d2H in surface sediments to precipitation d2H (from
the OIPC), with spatially varying intercept and slope, estimated via a Gaussian
process with a Matern 3/2 kernel. Environmental covariates (C4 vegetation
fraction, plant functional type composition, elevation, precipitation amount)
are averaged across 9 spatial scales using distance-weighted kernels, and the
optimal scale is estimated jointly with all other parameters.

Fourteen model variants are compared via leave-one-out cross-validation to
evaluate the contribution of each covariate group and the spatial random
effects.

## Requirements

Everything runs inside a container. Pull it from GHCR:

```bash
apptainer pull docker://ghcr.io/bradleylab/leafwax-spatial:latest
```

Or build locally:

```bash
docker build -t leafwax-spatial .
```

The container includes R 4.4.1, CmdStan 2.36.0, terra, sf, tidyverse, and all
other dependencies. Nothing else needs to be installed.

## Input data

The pipeline expects the following in the working directory:

- `input_data/global_data_c29.csv` -- sediment d2H compilation (1131 sites)
- `input_data/GlobalPrecip/d2h_MA.tif` -- OIPC mean annual d2H of precipitation
- `input_data/GlobalPrecip/d2h_se_MA.tif` -- OIPC standard error
- `input_data/elevation_5KMmn_GMTEDmn.tif` -- GMTED2010 elevation
- `input_data/C4_distribution_NUS_v2.2.nc` -- C4 grass fraction (NUS v2.2)
- `results/2d_MODIS_PFT_3classes_Downsampled.tif` -- MODIS PFT (tree/shrub/grass)
- `results/2f_TerraClimate_*.tif` -- TerraClimate annual means (ppt, soil, tmax, vpd)

Raster preprocessing scripts (`1_extract_c4_raster.R`, `2c_reproject_modis.R`,
`2d_downsample_modis.R`, `2f_process_terraclimate.R`) produce the derived
rasters from the raw source data.

## Running the pipeline

### On a single machine (EC2 or local)

```bash
# Step 1: Extract C4 raster from NetCDF
Rscript 1_extract_c4_raster.R

# Step 3: Data preparation (raster extraction, spatial averaging)
Rscript 3_prep_data.R

# Steps 4a-c: Stan data assembly and model fitting
./launch_pipeline.sh                     # all 14 models
./launch_pipeline.sh baseline_veg_sp     # single model
```

### On SLURM (e.g., WashU Compute2)

```bash
# Submit prep job, then 14 parallel fitting jobs
bash slurm/submit.sh

# Or step by step:
sbatch slurm/job_prep.sh
sbatch --dependency=afterok:<prep_job_id> --array=0-13 slurm/job_fit.sh
```

Each fitting job requests 8 CPUs (one per MCMC chain) and 120 GB RAM. All 14
models run simultaneously. Wall time is determined by the slowest spatial model,
typically 3-6 hours.

## Pipeline structure

| Step | Script | What it does |
|------|--------|-------------|
| 0 | `0_load_config.R` | Load `config.yaml` |
| 1 | `1_extract_c4_raster.R` | C4 fraction from NUS NetCDF |
| 2c | `2c_reproject_modis.R` | Reproject MODIS land cover to WGS84 |
| 2d | `2d_downsample_modis.R` | Aggregate MODIS to PFT percentages |
| 2f | `2f_process_terraclimate.R` | Compute TerraClimate annual means |
| 3 | `3_prep_data.R` | Raster extraction, C4 imputation, grid alignment |
| 4a | `4a_spatial_functions.R` | Spatial weighting, interactions, B-spline basis |
| 4b | `4b_stan_prep.R` | Standardization, Stan data assembly |
| 4c | `4c_fit_models.R` | Model fitting loop (CmdStanR) |
| 4d | `4d_leaf_wax_spatial_model.stan` | Bayesian model specification |
| 5a-e | `5a_*.R` through `5e_*.R` | Validation, diagnostics, spatial analysis |

## Post-hoc analysis

Scripts for sensitivity analysis, simulated data recovery, and spatial
confounding checks:

- `run_simulated_recovery.R` -- parameter recovery from synthetic data
- `run_sensitivity.R` -- prior sensitivity analysis
- `run_confounding_test_v2.R` -- spatial confounding simulation (Paciorek 2010)
- `analyze_paleo_inversion_models.R` -- paleoclimate inversion framework

These scripts require completed model fits in `model_output/`.

## Configuration

All hyperparameters, model definitions, and MCMC settings are in `config.yaml`.

## License

MIT. See [LICENSE](LICENSE).

## Citation

> Bradley, A.S. (2026). [Title]. Geochimica et Cosmochimica Acta. [DOI pending]
