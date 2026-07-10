# Leaf wax hydrogen isotope spatial calibration

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20172576.svg)](https://doi.org/10.5281/zenodo.20172575)

Companion code to Bradley (2026), submitted to *Geochimica et Cosmochimica Acta*.

Hierarchical Bayesian spatial model that calibrates sedimentary leaf-wax
n-C₂₉ hydrogen isotope ratios (δ²H<sub>wax</sub>) against precipitation
isotope composition and environmental covariates. Spatial autocorrelation
is captured by a Gaussian process with a Matérn 3/2 kernel; intercept and
slope vary across Earth's surface. Fourteen model variants are compared by
leave-one-out cross-validation.

## What is in this repo

| Path | Contents |
|---|---|
| `0_*.R` | Project configuration loader |
| `1_*.R`, `2a_*` … `2f_*` | Environmental raster extraction (C₄ NUS, MODIS PFT, TerraClimate) |
| `2g_*.R` … `2i_*.R` | Calibration data audit, sample-level archive classification, frozen-dataset build |
| `3_prep_data.R`, `3b_*` … `3h_*` | Data preparation + diagnostics |
| `4a_*.R` … `4c_*.R`, `4d_*.stan` | Stan data assembly + model fitting |
| `5a_*.R` … `5e_*_weighted.R` | Post-fit validation + diagnostics |
| `6a_simulated_recovery.R` | Parameter recovery from synthetic data |
| `6b_spatial_confounding_simulation.R` | Spatial-confounding simulation (Paciorek 2010) |
| `6c_prior_sensitivity.R` | Prior / hyperparameter sensitivity |
| `7_paleo_inversion.R` | Paleoclimate inversion application |
| `scripts/posterior_helpers.R`, `scripts/table_helpers.R` | Helpers sourced by the numbered scripts |
| `slurm/` | SLURM submission scripts (one production example) |
| `Dockerfile`, `.github/workflows/build-container.yml` | Containerized environment + CI |
| `config.yaml` | MCMC settings + 14 model specifications |
| `input_data/global_data_c29.csv` | Compiled δ²H<sub>wax</sub> dataset (1,136 rows from 73 publications; the modeling pipeline retains 1,128 after filtering, the n reported in Bradley 2026) |

The repo does **not** include manuscript drafts, figure/table generation
code, or the upstream pipeline that builds the input compilation.
Everything needed to reproduce the model fits and the diagnostics from the
input CSV is in the numbered scripts.

## Quick start (containerized)

```bash
git clone https://github.com/bradleylab/leafwax-spatial.git
cd leafwax-spatial

# (1) Pull the container
apptainer pull docker://ghcr.io/bradleylab/leafwax-calibration:latest
#   or:  docker pull ghcr.io/bradleylab/leafwax-calibration:latest

# (2) Place the four environmental rasters in input_data/  (see below)

# (3) Run the pipeline (single machine; SLURM example in slurm/)
apptainer exec leafwax-calibration_latest.sif Rscript 1_extract_c4_raster.R
apptainer exec leafwax-calibration_latest.sif Rscript 2c_reproject_modis.R
# … 2d, 2f, 3, then 4b/4c/4d via slurm or directly
```

The container ships R 4.4.1, CmdStan 2.36.0, terra, sf, tidyverse, posterior,
loo, bayesplot, and cmdstanr. No additional installs are needed. Re-build
locally with `docker build -t leafwax-spatial .` if you prefer.

## Obtaining input rasters

`input_data/global_data_c29.csv` is tracked in this repo. The four
environmental rasters that the pipeline reads at runtime (~308 MB total)
are public datasets and must be downloaded separately:

| File | Source | Where to obtain |
|---|---|---|
| `input_data/C4_distribution_NUS_v2.2.nc` (~241 MB) | Luo et al. 2024 — global C₄ vegetation distribution | National University of Singapore data repository (or contact the authors) |
| `input_data/GlobalPrecip/d2h_MA.tif` (~13 MB) | Bowen 2018 — Online Isotopes in Precipitation Calculator (OIPC), mean-annual δ²H | https://wateriso.utah.edu/waterisotopes/pages/data_access/oipc.html |
| `input_data/GlobalPrecip/d2h_se_MA.tif` (~14 MB) | Bowen 2018 — OIPC mean-annual δ²H standard error | (same as above) |
| `input_data/elevation_5KMmn_GMTEDmn.tif` (~30 MB) | Amatulli et al. 2018 — GMTED2010 elevation, 5 km mean | https://www.earthenv.org/topography (variable: `elevation_5KMmn_GMTEDmn.tif`) |

Re-projection and downsampling of MODIS land cover and TerraClimate are
handled by `2c_*`, `2d_*`, `2e_*`, and `2f_*` directly from the public
source APIs — no manual download required for those.

## Running on a single machine

After the rasters are in place:

```bash
Rscript 1_extract_c4_raster.R          # C4 fraction from the NUS NetCDF
python3 2a_download_modis.py           # MODIS PFT (downloads from MODIS server)
Rscript 2c_reproject_modis.R
Rscript 2d_downsample_modis.R
python3 2e_download_terraclimate.py    # TerraClimate annual means
Rscript 2f_process_terraclimate.R

# Calibration data audit + freeze (produces data/frozen/leafwax_d2h_c29_calibration_v1)
Rscript 2g_data_audit.R                # duplicate detection + coordinate archive split
Rscript 2h_archive_overrides.R         # sample-level archive classes (fetches Gensel PANGAEA)
Rscript 2i_freeze_calibration.R        # apply decisions -> frozen calibration dataset

Rscript 3_prep_data.R                  # joins all covariates + spatial averaging

# Model fitting (loops 14 model variants)
Rscript 4b_stan_prep.R                 # Stan data assembly per model
Rscript 4c_fit_models.R                # MCMC sampling via cmdstanr
```

On a single machine expect ~3–6 h per spatial model with 8 chains.
Validation / diagnostic scripts (`5a_*` through `5e_*`) and simulation
studies (`6a–6c`) consume the resulting `fit.rds` / `posterior_draws.rds`
artifacts and can be run independently.

## Running on SLURM

Production fits are most efficient on a cluster. `slurm/submit.sh` chains
the prep job and 14 parallel fitting jobs (one per model variant); each
fitting job requests 8 CPUs and ≈120 GB RAM. Wall time is set by the
slowest spatial model (~3–6 h).

```bash
sbatch slurm/job_prep.sh
sbatch --dependency=afterok:<prep_job_id> --array=0-13 slurm/job_fit.sh
sbatch --dependency=afterok:<fit_job_id> slurm/postprocess_fits.sh
# or:
bash slurm/submit.sh
```

Scripts target WashU Compute2 conventions (account flag, partition); adapt
to your cluster as needed.

## Configuration

`config.yaml` defines MCMC settings, the 14 model variants, raster paths,
and output directories. Every numbered script begins by sourcing
`0_load_config.R`, which reads `config.yaml`.

## Repository scope (what is intentionally not here)

- **Compilation building.** The upstream Python / R pipeline that produces
  `input_data/global_data_c29.csv` from primary sources (PANGAEA datasets,
  Ladd 2021 compilation, individual paper supplements, etc.) is not part of
  this repo. The shipped artifact is the CSV; the build scripts and
  intermediate spreadsheets stay with the author.
- **Manuscript.** Figure code, table code, narrative drafts, and rendered
  artifacts are not part of this repo.
- **Model fits.** `results/` and `model_output/` are git-ignored. Producing
  fresh fits from a clean clone takes the wall-clock time noted above. The
  published-version posterior draws are archived at Zenodo:
  [10.5281/zenodo.20085465](https://doi.org/10.5281/zenodo.20085465)
  (`leafwax v10 model posteriors`, CC-BY-4.0). The leafwax R package
  downloads these posteriors on first use; manual download is not required.

## License

MIT — see [LICENSE](LICENSE).

## Citation

Cite both the software archive and the related manuscript:

```bibtex
@software{bradley_leafwax_spatial_2026,
  author  = {Bradley, Alexander S.},
  title   = {leafwax-spatial: hierarchical Bayesian spatial calibration
             of leaf-wax hydrogen isotopes},
  year    = {2026},
  doi     = {10.5281/zenodo.20172575},
  url     = {https://doi.org/10.5281/zenodo.20172575}
}

@unpublished{bradley_leafwax_paper_2026,
  author = {Bradley, Alexander S.},
  title  = {Spatial modeling improves the calibration of leaf wax
            hydrogen isotopes to precipitation},
  year   = {2026},
  note   = {Manuscript in preparation}
}
```

The `@software` DOI is the concept DOI — it always resolves to the
latest version. To cite a specific release, replace it with that
release's version DOI from the Zenodo deposit page.
