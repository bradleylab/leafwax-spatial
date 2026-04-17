# METHODS.md — Leaf Wax Spatial Calibration

Cross-reference of manuscript methods against the code that produced the
April 2026 model run. Line numbers refer to the canonical pipeline in
`leafwax_gca_working/` as of the Apr 2026 fit (repo: `bradleylab/leafwax-spatial`,
commit current at this file's Git blame).

Every quantitative claim below cites `script.R:LINE`. Where code and prior
manuscript text disagree, the discrepancy is flagged at the end of each
section and collected in §9.

---

## 1. Data compilation & preprocessing

### 1.1 Sediment data

- 1131 n-C29 alkane δ²H measurements from lake and river surface sediments.
- Filter: `chain == 29`, non-NA for `d2H_wax`, `latitude`, `longitude`
  (`3_prep_data.R:55-60`). Input file is `input_data/global_data_c29.csv`
  (`3_prep_data.R:40`), built by the Python pipeline in
  `data/compilation/` (driver: `build_final_dataset.py`).

### 1.2 Measurement uncertainty

| Reported `d2H_wax_err` | Assigned value (‰) |
|---|---|
| NA | 3 |
| 0 | 1 |
| otherwise | as reported |

Code: `3_prep_data.R:49-53`.

### 1.3 Coordinate handling

- CRS: WGS84 (EPSG:4326).
- Gaussian jitter on duplicate coordinates: SD = 0.0001°, seed = 12345.
  Code: `4b_stan_prep.R:70-79`; jitter magnitude from `config.yaml:72`
  (`coordinate_jitter: 0.0001`).

### 1.4 Environmental covariates

| Variable | Source raster | Native grid | Loaded at |
|---|---|---|---|
| δ²H precip (mean annual) | OIPC `d2h_MA.tif` | 2083 × 4320, 0.0833° | `3_prep_data.R:73` |
| δ²H precip SE | OIPC `d2h_se_MA.tif` | 2083 × 4320, 0.0833° | `3_prep_data.R:74` |
| Elevation (GMTED) | `elevation_5KMmn_GMTEDmn.tif` | ~5 km | `3_prep_data.R:75` |
| PFT (tree/shrub/grass) | `2d_MODIS_PFT_3classes_Downsampled.tif` | 2186 × 4371, ~0.0824° | `3_prep_data.R:76` |
| C4 fraction | NUS v2.2, 2001–2019 mean | 360 × 720, 0.5° | `1_extract_c4_raster.R` → `results/1_C4_total_mean.tif` (`3_prep_data.R:72`) |
| Annual precipitation | TerraClimate ppt, 2001–2019 mean | ~4 km | `3_prep_data.R:79` |

Soil moisture, max temperature, and VPD are extracted
(`3_prep_data.R:80-82`) but excluded from every model (see §2.3). Sites are
**not** filtered for missing values in these unused climate variables
(`3_prep_data.R:591-600`).

### 1.4b Raster grid harmonization

Interactions (OIPC × C4, OIPC × PFT) are computed as pixel-level products
inside each site's extraction radius. This requires C4 and PFT rasters on
the OIPC grid. C4 and PFT are therefore resampled to the OIPC reference
grid before extraction:

| Raster | Native | Resampled to | Method | Code |
|---|---|---|---|---|
| C4 | 360 × 720, 0.5° | 2083 × 4320, 0.0833° | bilinear | `3_prep_data.R:99` |
| PFT | 2186 × 4371 | 2083 × 4320, 0.0833° | bilinear → clamp[0,1] → per-pixel renormalize so Tree+Shrub+Grass = 1 | `3_prep_data.R:102-111` |

Ocean / no-data PFT pixels (row-sum = 0) are left as 0/0/0
(`3_prep_data.R:107-109`). Grid equality is asserted
(`3_prep_data.R:113-118`). PFT renormalization is verified against a
tolerance of 1e-6 (`3_prep_data.R:120-130`).

OIPC, OIPC SE, elevation, and TerraClimate rasters are used at native
resolution.

### 1.5 C4 imputation

Applied to NA C4 pixels after resampling:

- `|latitude| >= 50°` → C4 = 0
- mean GMTED elevation ≥ 1500 m → C4 = 0
- otherwise: left as NA

Thresholds: `lat_threshold = 50` (`3_prep_data.R:25`),
`elev_threshold = 1500` (`3_prep_data.R:26`). Imputation function:
`apply_c4_ecological_rules_gmted` at `3_prep_data.R:509-536`; invoked at
`3_prep_data.R:539-549`.

### 1.6 Multi-scale spatial averaging

Covariates are averaged at 9 distance scales using a Gaussian distance
kernel: `w = exp(-d² / (2 r²))` where `d` is in km and `r` is the scale.

- Scales (km): `{1, 3, 5, 10, 20, 40, 70, 100, 150}` — `config.yaml:29`.
- Extraction radius: 5° (`3_prep_data.R:24`, `max_radius_deg`).
- Kernel: `gaussian_weights` at `4a_spatial_functions.R:18-20`;
  distance-weighted aggregation in `compute_weighted_mean` at
  `4a_spatial_functions.R:33-…` (degree → km via `deg_to_km`,
  `4a_spatial_functions.R:23-30`).

---

## 2. Exploratory data analysis

### 2.1 Spatial autocorrelation

Moran's I is computed on OLS residuals (`3e_spatial_clustering_analysis.R`).

### 2.2 Collinearity screening

- VIF threshold = 5 (`3d_analyze_collinearity.R`).
- Elastic-net regularization: α = 0.5, 10-fold CV
  (`3d_analyze_collinearity.R`).

### 2.3 Variable selection outcome

- **Retained:** δ²H precip (OIPC), elevation, C4, PFT (tree/shrub/grass),
  precipitation.
- **Excluded (collinearity):** VPD, max temperature, soil moisture.

Exclusion is hard-coded at `4b_stan_prep.R:200-202`:

```r
include_temp = FALSE,   # Always FALSE due to collinearity
include_vpd  = FALSE,   # Always FALSE due to collinearity
include_soil = FALSE,   # Always FALSE due to collinearity
```

---

## 3. Statistical model

### 3.1 Full model

```
d2H_wax,i ~ N(mu_i, sigma_total^2)

mu_i = beta_0(s_i)
     + beta_OIPC(s_i) * OIPC_i
     + f_elev(elev_i)
     + beta_C4 * C4_i
     + beta_precip * precip_i
     + beta_tree * tree_i + beta_shrub * shrub_i + beta_grass * grass_i
     + beta_OIPC_x_C4 * (OIPC_i * C4_i)
     + beta_OIPC_x_tree * (OIPC_i * tree_i)
     + beta_OIPC_x_shrub * (OIPC_i * shrub_i)
     + beta_OIPC_x_grass * (OIPC_i * grass_i)
```

Stan linear-predictor assembly: `4d_leaf_wax_spatial_model.stan:340-386`.
Interactions are pixel-level products within each site's extraction radius,
distance-weighted, then aggregated across 9 scales (§1.6, §3.2).

### 3.2 Multi-scale weighting (Stan side)

Predictors enter the likelihood as `X_i = Σ_s w_s · X_{i,s}` with weights

```
w_s = exp(-r_s / lambda) / Σ exp(-r_s / lambda)
```

Code: `4d_leaf_wax_spatial_model.stan:267-271`. Lambda is estimated when
`estimate_lambda == 1` (`config.yaml:32-36`):

- Prior: `lambda_decay ~ Lognormal(2.5, 0.5)` — median ≈ 12 km,
  approximate 95% CI 4–28 km (`4d_leaf_wax_spatial_model.stan:434`).
- Bounds: `<lower=1, upper=400>` km (`4d_leaf_wax_spatial_model.stan:247`).

### 3.3 Elevation B-spline

- 9 interior knots, degree 3 (cubic) → 13 total basis functions
  (`config.yaml:59-60`: `n_elevation_knots: 9`,
  `elevation_spline_degree: 3`).
- Interior knots at equally spaced quantiles within the data range.
- Random-walk smoothing prior on coefficients (`k = 2, …, 13`):

  ```
  beta_elev_bspline[k] ~ Normal(beta_elev_bspline[k-1], tau_elev_bspline)
  ```

  — `4d_leaf_wax_spatial_model.stan:398-400`.
- `tau_elev_bspline_raw ~ Normal⁺(0, 1)`
  (`4d_leaf_wax_spatial_model.stan:396`).
- First coefficient: `beta_elev_bspline[1] ~ Normal(0, 2)`
  (`4d_leaf_wax_spatial_model.stan:397`). See §9 Discrepancy 2.

### 3.4 Three-component observation variance

Observation-level variance combines measurement error, propagated OIPC
error, and residual variance:

```
sigma_total^2 = d2H_wax_err^2
              + (beta_OIPC(s_i) * oipc_se_weighted)^2
              + sigma^2
```

The OIPC SE term is scaled by the local OIPC slope so predictor-unit
uncertainty is converted to response-unit uncertainty (standard Gaussian
error-in-variables propagation).

Code (likelihood): `4d_leaf_wax_spatial_model.stan:457-461`:

```stan
real total_var = square(d2H_wax_err[n])
               + square(beta_oipc_spatial[n] * oipc_se_weighted[n])
               + square(sigma);
real total_sd = sqrt(total_var);
d2H_wax[n] ~ normal(mu[n], total_sd);
```

The same variance formula is used in `generated quantities` for `log_lik`
(LOO-CV) and posterior predictive checks.

---

## 4. Spatial framework

### 4.1 Gaussian processes

Two independent GPs for spatially-varying coefficients:

- Intercept: `beta_0(s) = beta_0_global + GP_intercept(s)`
- OIPC slope: `beta_OIPC(s) = beta_OIPC_global + GP_slope(s)`

Construction: `4d_leaf_wax_spatial_model.stan:315-337`. Both GPs share a
single length scale (§4.7).

### 4.2 Matérn 3/2 covariance

```
k(s_i, s_j) = eta^2 * (1 + sqrt(3) * d / rho) * exp(-sqrt(3) * d / rho)
```

- `eta^2 = 1.0` (fixed amplitude).
- `rho = exp(log_ls_spatial)` (estimated, in standardized coordinate units).

Kernel functions: `4d_leaf_wax_spatial_model.stan:4-41`
(`cov_matern32`, `cov_matern32_cross`).

### 4.3 Predictive-process approximation

- `m = 125` knots for all spatial models (`config.yaml:186, 200, 214, …`,
  `n_pp_knots: 125` on every `_sp` model).
- Projection: `K_cross * K_knots^{-1} * w_knots`
  (`4d_leaf_wax_spatial_model.stan:325-336`).
- Diagonal jitter: `1e-4` (`config.yaml:71`, applied at
  `4d_leaf_wax_spatial_model.stan:327`).

### 4.4 Knot placement: spherical Fibonacci lattice

Deterministic equal-area global grid via golden-angle algorithm
(`4a_spatial_functions.R:258-…`, invoked when
`CONFIG$knot_selection_method == "regular_global"`,
`config.yaml:55`).

### 4.5 Spatially adaptive regularization

Knot-level density is counted within a radius of 0.2 standardized units
(`config.yaml:40`, `density_radius_std: 0.2`) with a density-scaling
constant of 10 (`config.yaml:41`).

Density-based τ tiers (`4d_leaf_wax_spatial_model.stan:185-194`):

- `d == 0`: τ = 0.50 (both intercept and slope)
- `0 < d < 10`: τ = 0.50 + 0.30 · (d/10)
- `d ≥ 10`: τ = 0.80

Slope τ is then multiplied by a data-relative OIPC range factor
(`4d_leaf_wax_spatial_model.stan:200-209`):

- knot OIPC range < 0.25 × global max range → factor = 0.2
- 0.25 × ≤ range < 0.60 × → factor = 0.5
- ≥ 0.60 × → factor = 1.0

Thresholds are relative to the maximum OIPC range across all knots
(`4d_leaf_wax_spatial_model.stan:179`, `max_oipc_range = max(oipc_range_at_knots)`).

### 4.6 PC priors for spatial SDs

- Intercept GP: `P(sigma_intercept > 20 ‰) = 0.05`
  (`config.yaml:44-47`).
- Slope GP: `P(sigma_slope > 0.3) = 0.05`
  (`config.yaml:49-52`).

Implementation via exponential rate parameter
`lambda = -log(alpha) / u` (`4d_leaf_wax_spatial_model.stan:218-219`);
applied at `4d_leaf_wax_spatial_model.stan:443-444`.

### 4.7 GP length-scale prior

```
log_ls_spatial ~ Normal(-1.0, 0.4),   bounds <lower=-2, upper=0>
```

Code: `4d_leaf_wax_spatial_model.stan:250` (bounds), `:440` (prior).

---

## 5. Complete prior table

| Parameter | Prior | Stan line |
|---|---|---|
| `beta_0` | Normal(0, 5) | 391 |
| `beta_oipc` | Normal(0.8, 0.3) | 392 |
| `tau_elev_bspline_raw` | Normal⁺(0, 1) | 396 |
| `beta_elev_bspline[1]` | Normal(0, 2) | 397 |
| `beta_elev_bspline[k>1]` | Normal(`beta_elev[k-1]`, `tau_elev_bspline`) | 399 |
| `beta_c4` | Normal(0, 2) | 405 |
| `beta_precip` | Normal(0, 0.5) | 410 |
| `beta_tree`, `beta_shrub`, `beta_grass` | Normal(0, 2) | 415–417 |
| `beta_oipc_x_c4` | Normal(0, 0.5) | 423 |
| `beta_oipc_x_tree`, `_shrub`, `_grass` | Normal(0, 0.5) | 426–428 |
| `lambda_decay` | Lognormal(2.5, 0.5), [1, 400] | 434, bounds at 247 |
| `log_ls_spatial` | Normal(−1.0, 0.4), [−2, 0] | 440, bounds at 250 |
| `sigma_intercept` | Exponential(PC prior) | 443 |
| `sigma_slope` | Exponential(PC prior) | 444 |
| `z_intercept_spatial[k]` | Normal(0, `tau_spatial_intercept[k]`) | 448 |
| `z_slope_spatial[k]` | Normal(0, `tau_spatial_slope[k]`) | 449 |
| `sigma` (residual) | Normal⁺(0, 2) | 454 |

All line numbers refer to `4d_leaf_wax_spatial_model.stan`.

---

## 6. Model configurations

### 6.1 Manuscript models (14)

Defined in `config.yaml:105-302`. Five non-spatial and nine spatial
variants compare combinations of covariates (OIPC, C4, PFT, elevation,
precipitation), vegetation interactions, and spatial varying coefficients
(GP):

| Model | C4 | PFT | elev | precip | veg int. | GP |
|---|---|---|---|---|---|---|
| `baseline` (`:107`) | – | – | – | – | – | – |
| `baseline_veg` (`:121`) | ✓ | ✓ | – | – | ✓ | – |
| `baseline_env` (`:135`) | – | – | ✓ | ✓ | – | – |
| `full` (`:149`) | ✓ | ✓ | ✓ | ✓ | – | – |
| `full_interact` (`:163`) | ✓ | ✓ | ✓ | ✓ | ✓ | – |
| `baseline_sp` (`:178`) | – | – | – | – | – | ✓ |
| `baseline_veg_sp` (`:192`) | ✓ | ✓ | – | – | ✓ | ✓ |
| `baseline_env_sp` (`:206`) | – | – | ✓ | ✓ | – | ✓ |
| `c4_only_sp` (`:220`) | ✓ | – | – | – | – | ✓ |
| `elevation_only_sp` (`:234`) | – | – | ✓ | – | – | ✓ |
| `elevation_c4_sp` (`:248`) | ✓ | – | ✓ | – | – | ✓ |
| `elevation_c4_interact_sp` (`:262`) | ✓ | – | ✓ | – | ✓ | ✓ |
| `full_sp` (`:276`) | ✓ | ✓ | ✓ | ✓ | – | ✓ |
| `full_interact_sp` (`:290`) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |

### 6.2 MCMC settings

Per-model settings in `config.yaml:105-302`; sampler call at
`4c_fit_models.R:146-158`:

| Model family | Chains | Iter (total) | Warmup | `adapt_delta` | `max_treedepth` |
|---|---|---|---|---|---|
| Non-spatial (`baseline`, `baseline_veg`, `baseline_env`, `full`, `full_interact`) | 8 | 2000 | 50% | 0.95 | 15 |
| `baseline_veg_sp` | 8 | 3000 | 50% | 0.95 | 14 |
| `baseline_sp`, `baseline_env_sp`, `c4_only_sp`, `elevation_only_sp`, `elevation_c4_sp` | 8 | 3000 | 50% | 0.99 | 12–15 |
| `elevation_c4_interact_sp`, `full_sp`, `full_interact_sp` | 8 | 4000 | 50% | 0.99 | 12 |

Global defaults: `stan_seed = 314`, `warmup_ratio = 0.5`,
`refresh = 100`, `max_parallel_models = 4`
(`config.yaml:77-85`).

### 6.3 Convergence criteria

- R-hat < 1.01
- ESS (bulk) > 900
- Zero divergent transitions
- E-BFMI > 0.3

Diagnostics pulled via `fit$diagnostic_summary()` at
`4c_fit_models.R:202`.

### 6.1 Fit-step outputs

Each call to `4c_fit_models.R` writes the following into
`model_output/<model_name>/`:

| File | Contents | Purpose |
|---|---|---|
| `fit.rds` | `CmdStanFit` stub with absolute paths to the per-chain CSVs written during sampling. | Kept for HPC-local reanalysis. Not read by downstream analysis. |
| `posterior_draws.rds` | `draws_array` with the widened variable set: all scalars in `params_to_check`, plus `mu`, `d2H_rep`, `log_lik`, `scale_weights`, and (for GP models) `alpha_spatial`, `beta_oipc_spatial`, `z_intercept_spatial`, `z_slope_spatial`. | Per-draw posterior used by every figure and table producer. |
| `diagnostics.rds` | `list` with `summary` (per-chain divergences / E-BFMI), `param_summary` (key scalars), `all_params_summary` (mean/sd/q5/q95/rhat/ess for every parameter). | Summary-level reads when per-draw is overkill. |
| `loo.rds` | `psis_loo` from `loo::loo(log_lik_array)`. | Table 1 LOOIC / SE / p_eff / n_hi_k. |
| `runtime_info.rds` | Elapsed sampling time. | Runtime reporting. |

Downstream analysis reads only `posterior_draws.rds`, `diagnostics.rds`,
and `loo.rds`. It never reads `fit.rds` or chain CSVs, which keeps every
figure and table reproducible on any machine with R and the rds files.

---

## 7. Standardization

| Variable | Transformation | Code |
|---|---|---|
| `d2H_wax` | (x − mean) / sd | `4a_spatial_functions.R:614` |
| `d2H_wax_err` | x / sd(d2H_wax) | `4a_spatial_functions.R:615` |
| OIPC δ²H | (x − mean) / sd | `4a_spatial_functions.R:861` |
| OIPC SE | x / sd(OIPC) | `4a_spatial_functions.R:862` |
| Elevation | → km, then (x_km − mean_km) / sd_km | `4a_spatial_functions.R:877-878` |
| C4 | (x − 20) / 25 (fixed constants) | `config.yaml:63-65`, applied at `4a_spatial_functions.R:868` |
| Precipitation | (x − mean) / sd | `4a_spatial_functions.R:911` |
| PFT | proportions 0–1, not standardized | — |
| OIPC × C4 | x / (sd_OIPC · sd_C4) | `4a_spatial_functions.R:970` |
| OIPC × PFT | x / sd_OIPC | `4a_spatial_functions.R:977-979` |
| Coordinates (GP) | per-dimension (x − mean) / sd | `4b_stan_prep.R` / downstream |

Scaling parameters (d2H mean/sd, OIPC mean/sd, elevation mean/sd, C4
from config) are assembled in `SCALING_PARAMS` at `4b_stan_prep.R:90-99`.

### 7.1 NA handling in aggregated matrices

After distance-weighted averaging at 9 scales, some site-scale cells can
be NA (coastal sites, missing covariates). Each NA is filled with the
site's mean across other scales, or 0 on the standardized scale
(= population mean) if all scales are NA. PFT fractions default to 0.33
(equal split). OIPC SE defaults to 0.1.

Filling helper: `fill_na_row` at `4a_spatial_functions.R:989`; applied
per-matrix at `4a_spatial_functions.R:1021-1042+`.

---

## 8. Post-fit validation

- **LOO-CV and WAIC**: `5a_model_validation.R:46-91`. Comparison via
  `loo_compare` at `5a_model_validation.R:85-88`; outputs saved as
  `results/loo_results.rds`, `results/loo_comparison.rds`.
- **Spatial cross-validation** (held-out regions + knot-density bins):
  `spatial_cv_analysis` at `5b_overfitting_diagnostics.R:30-129`.
  Regions are defined geographically
  (`5b_overfitting_diagnostics.R:59-66`): Americas, Europe, Africa,
  Asia, Oceania.
- **Spatial-pattern diagnostics** (maps, variograms, Moran's I on
  posterior residuals): `5c_spatial_patterns.R`.
- **Spatial-confounding check** (collinearity between spatial effects and
  δ²H precip, across all `_sp` models): `5d_spatial_confounding_check.R`.
- **Simulated recovery, sensitivity, paleo-inversion**:
  `run_simulated_recovery.R`, `run_sensitivity.R`,
  `run_confounding_test_v2.R`, `analyze_paleo_inversion_models.R`.

---

## 9. Discrepancies between code and previously-drafted manuscript text

These flags are carried forward from the December 2025 drafts so the
revision can resolve them.

### Discrepancy 1 — multi-scale weighting kernel

- **Manuscript (Supp Eq 8):** exponential, `w = exp(-d / r)`.
- **Code (`4a_spatial_functions.R:19`):** Gaussian,
  `w = exp(-d² / (2 r²))`.

The Stan-side scale weighting (`exp(-r_s / lambda)`,
`4d_leaf_wax_spatial_model.stan:269`) is exponential and consistent
between code and manuscript.

### Discrepancy 2 — B-spline first coefficient prior

- **Table S1 (as drafted):** `beta_elev[1] ~ Normal(0, tau_elev)`.
- **Code (`4d_leaf_wax_spatial_model.stan:397`):**
  `beta_elev_bspline[1] ~ Normal(0, 2)` — fixed SD of 2.

### Discrepancy 3 — adaptive regularization formulas

- **Table S1 and code agree.**
- **Supplement Eqs 19–21 (as drafted):** stale (different starting
  values, different ramp, exponential vs. linear).

---

## 10. Changes since the December 2025 draft

| # | Change | Dec 2025 | Apr 2026 |
|---|---|---|---|
| A | Raster grid harmonization | not documented | C4 and PFT resampled to OIPC grid (§1.4b) |
| B | OIPC SE error propagation | `square(oipc_se_weighted[n])` | `square(beta_oipc_spatial[n] * oipc_se_weighted[n])` (§3.4) |
| C | `beta_precip` prior | Normal(0, 0.02) | Normal(0, 0.5) (§5) |
| D | Sample size / filtering | 818 sites, filtered on unused soil/temp/VPD | 1131 sites, no filter on unused climate (§1.1, §2.3) |
| E | Model count | 6 | 14 (§6.1) |
| F | Adaptive regularization range factor | fixed thresholds (5, 15 std units) | 0.25 / 0.60 × max OIPC range at knots (§4.5) |

Item E differs from an earlier note claiming "14 manuscript + 1 variant
(15 total)"; the April `config.yaml:105-302` contains exactly 14 model
entries.
