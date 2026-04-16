# METHODS.md — Leaf Wax Spatial Calibration Model (Updated 2026-04-10)

> Cross-reference of manuscript methods and code implementation.
> Updated against code on branch `fix/pre-run-patches` (commit `6bb6fdd`), repo `bradleylab/leafwax-spatial`.
> Original METHODS.md written Dec 2025 against `leafwax-spatial-master/` snapshot.
> This version incorporates all pre-flight code review fixes (P0-1 through P1-6).

---

## Changes from December 2025 METHODS.md

| ID | Change | Old | New | Section |
|----|--------|-----|-----|---------|
| A | Raster grid harmonization added | Not documented | C4 and PFT resampled to OIPC grid | 1.4b (new) |
| B | OIPC SE error propagation corrected | `square(oipc_se_weighted[n])` | `square(beta_oipc_spatial[n] * oipc_se_weighted[n])` | 3.4 |
| C | beta_precip prior widened | `Normal(0, 0.02)` | `Normal(0, 0.5)` | 5 (prior table) |
| D | Site count updated, unused filter removed | 818 sites; filtered on soil/temp/VPD NAs | 1131 sites; no filter on unused climate variables | 1.1, 2.3 |
| E | Model count updated | 6 models | 15 models (14 manuscript + 1 variant) | 6.2 |
| F | Range factor thresholds now data-relative | Fixed thresholds (5, 15 std units) | 0.25 and 0.60 of max OIPC range | 4.5 |

---

## 1. Data Compilation & Preprocessing

### 1.1 Sediment Data

- 1131 n-C29 alkane d2H measurements from lake and river surface sediments
- Compiled from published datasets; 129 sites added from Perplexity-compiled sources (cross-checked against original supplements)
- Filtered: `chain == 29`, non-NA for d2H_wax, latitude, longitude
- **Code**: `3_prep_data.R:40-63`

### 1.2 Measurement Uncertainty

| Condition | d2H_wax_err (permil) |
|---|---|
| Unreported (NA) | 3 |
| Reported as 0 | 1 |
| Otherwise | as reported |

**Code**: `3_prep_data.R:49-53`

### 1.3 Coordinate Handling

- Coordinate system: WGS84 (EPSG:4326)
- Gaussian jitter for duplicate coordinates: SD = 0.0001 degrees, seed = 12345
- **Code**: `4b_stan_prep.R:88-97`, `config.yaml:71`

### 1.4 Environmental Covariates

| Variable | Source | Native Resolution | Code |
|---|---|---|---|
| d2H_precip | OIPC (d2h_MA.tif) | 2083 x 4320, 0.0833 deg | `3_prep_data.R:73` |
| d2H_precip SE | OIPC (d2h_se_MA.tif) | 2083 x 4320, 0.0833 deg | `3_prep_data.R:74` |
| Elevation | GMTED (elevation_5KMmn_GMTEDmn.tif) | ~5 km | `3_prep_data.R:75` |
| PFT (tree/shrub/grass) | MODIS 3-class (2d_MODIS_PFT_3classes_Downsampled.tif) | 2186 x 4371, ~0.0824 deg | `3_prep_data.R:76` |
| C4 fraction | NUS v2.2 (C4_distribution_NUS_v2.2.nc), 2001-2019 mean | 360 x 720, 0.5 deg | `1_extract_c4_raster.R` |
| Annual precipitation | TerraClimate (ppt, 2001-2019 mean) | ~4 km | `3_prep_data.R:79` |

Note: Soil moisture, max temperature, and VPD are extracted from TerraClimate but excluded from all model configurations due to collinearity (see Section 2.3). Sites are NOT filtered for missing values in these unused variables.

### 1.4b Raster Grid Harmonization (P0-1 fix)

Interaction computations require that OIPC and covariate rasters produce equal-length pixel vectors within each site's extraction radius. Because the source rasters have different native grids, C4 and PFT rasters are resampled to the OIPC reference grid before extraction:

| Raster | Native Grid | Resampled To | Method |
|---|---|---|---|
| C4 | 360 x 720, 0.5 deg | 2083 x 4320, 0.0833 deg | bilinear |
| PFT | 2186 x 4371, ~0.0824 deg | 2083 x 4320, 0.0833 deg | bilinear + clamp + renormalize |

For PFT: bilinear interpolation of fractional layers (Tree, Shrub, Grass) can produce values slightly outside [0, 1] near class boundaries. Each layer is clamped to [0, 1], then per-pixel values are renormalized so Tree + Shrub + Grass = 1.0. Pixels with rowsum = 0 (ocean/no-data) are left as 0/0/0.

OIPC, OIPC SE, elevation, and TerraClimate rasters are used at native resolution (no resampling needed for main effects; interactions only involve OIPC with C4 and PFT).

**Code**: `3_prep_data.R:86-163`

### 1.5 C4 Imputation

Rules applied to NA C4 pixels (after resampling):
- |latitude| >= 50 degrees -> C4 = 0
- Mean GMTED elevation >= 1500 m -> C4 = 0
- Otherwise: left as NA

**Code**: `3_prep_data.R:24-26` (thresholds), `3_prep_data.R:548-575` (imputation function)

### 1.6 Multi-Scale Spatial Averaging

Covariates averaged at 9 distance scales using a Gaussian weighting kernel.

Scales: {1, 3, 5, 10, 20, 40, 70, 100, 150} km
**Code**: `config.yaml:27`

Extraction radius: 5 degrees (`3_prep_data.R:24`)

#### [DISCREPANCY 1] Weighting Kernel

- **Manuscript** (Supp Eq 8): Exponential decay: `w = exp(-d / r)`
- **Code** (`4a_spatial_functions.R:18-19`): Gaussian kernel: `w = exp(-d^2 / (2 * r^2))`

**Impact**: Gaussian drops off faster at short distances but slower at long distances than exponential. The Stan-side scale weighting (Supp Eq 10) uses `exp(-r / lambda)` and IS consistent between manuscript and code.

---

## 2. Exploratory Data Analysis

### 2.1 Spatial Autocorrelation

- Moran's I = 0.614 on OLS residuals
- **Code**: `3e_spatial_clustering_analysis.R`

### 2.2 Collinearity Screening

- VIF threshold = 5
- Elastic net regularization: alpha = 0.5, 10-fold CV
- **Code**: `3d_analyze_collinearity.R:113-125` (VIF), `3d_analyze_collinearity.R:189-215` (elastic net)

### 2.3 Variable Selection Outcome

**Retained**: d2H_precip (OIPC), elevation, C4, PFT (tree/shrub/grass), precipitation
**Excluded**: VPD, max temperature, soil moisture (due to collinearity)

**Code**: `4b_stan_prep.R:219-221` — always set to FALSE:
```r
include_temp = FALSE,   # Always FALSE due to collinearity
include_vpd = FALSE,    # Always FALSE due to collinearity
include_soil = FALSE,   # Always FALSE due to collinearity
```

Note: The data prep filter on these three variables has been removed. Sites missing unused climate variables are no longer excluded from modeling.

---

## 3. Statistical Model Specification

### 3.1 Full Model (Main Eq 1-2)

```
d2H_wax,i ~ N(mu_i, sigma_total^2)

mu_i = beta_0(s_i) + beta_OIPC(s_i) * OIPC_i
       + f_elev(elev_i)
       + beta_C4 * C4_i
       + beta_precip * precip_i
       + beta_tree * tree_i + beta_shrub * shrub_i + beta_grass * grass_i
       + beta_OIPC_x_C4 * (OIPC_i * C4_i)
       + beta_OIPC_x_tree * (OIPC_i * tree_i)
       + beta_OIPC_x_shrub * (OIPC_i * shrub_i)
       + beta_OIPC_x_grass * (OIPC_i * grass_i)
```

Interactions (OIPC x C4, OIPC x PFT) are computed as pixel-level products within the extraction radius, then distance-weighted and aggregated. This requires aligned grids (see Section 1.4b).

**Code**: `4d_leaf_wax_spatial_model.stan:332-380`

### 3.2 Multi-Scale Weighting (Supp Eq 10)

Each predictor is a weighted average across 9 spatial scales:
```
X_weighted = sum_s( w_s * X_s )
w_s = exp(-r_s / lambda) / sum( exp(-r_s / lambda) )
```

Lambda prior: `lambda ~ lognormal(2.5, 0.5)` — median ~12 km, 95% CI ~4-28 km
Lambda bounds: `<lower=1, upper=400>` km

**Code**: `4d_leaf_wax_spatial_model.stan:259-264`, `config.yaml:30-34`

### 3.3 B-Spline for Elevation

- 9 interior knots, degree 3 (cubic), 13 total basis functions
- Interior knots placed at equally spaced quantiles within data range
- Random walk prior on coefficients:
  ```
  beta_elev[k] ~ Normal(beta_elev[k-1], tau_elev)  for k = 2, ..., 13
  ```

#### [DISCREPANCY 2] First B-Spline Coefficient Prior

- **Table S1**: `beta_elev[1] ~ Normal(0, tau_elev)`
- **Code** (`4d_leaf_wax_spatial_model.stan:396`): `beta_elev_bspline[1] ~ normal(0, 2)` — fixed SD = 2

### 3.4 Three Uncertainty Components (error-in-variables model)

The observation-level variance combines three sources:

```
sigma_total^2 = sigma_measurement^2 + (beta_OIPC(s_i) * sigma_OIPC)^2 + sigma_residual^2
```

- `sigma_measurement` = d2H_wax_err[n] (per-observation, standardized by d2H_wax SD)
- `sigma_OIPC` = oipc_se_weighted[n] (per-observation, standardized by OIPC SD)
- `sigma_residual` = sigma (estimated parameter)

The OIPC SE term is scaled by the local OIPC slope `beta_oipc_spatial[n]` to convert from predictor-unit uncertainty to response-unit uncertainty. This is standard Gaussian error propagation for error-in-variables models: for y = beta * x + ..., Var(y from x error) = beta^2 * Var(x).

**Code**: `4d_leaf_wax_spatial_model.stan:462-466`
```stan
// Likelihood with OIPC measurement error (error-in-variables propagation)
real total_var = square(d2H_wax_err[n]) + square(beta_oipc_spatial[n] * oipc_se_weighted[n]) + square(sigma);
real total_sd = sqrt(total_var);
d2H_wax[n] ~ normal(mu[n], total_sd);
```

The same variance formula is used in the generated quantities block for log_lik (LOO-CV) and posterior predictive checks.

---

## 4. Spatial Framework

### 4.1 Gaussian Processes

Two independent GPs for spatially-varying coefficients:
- Intercept: `beta_0(s) = beta_0_global + GP_intercept(s)`
- OIPC slope: `beta_OIPC(s) = beta_OIPC_global + GP_slope(s)`

**Code**: `4d_leaf_wax_spatial_model.stan:308-330`

### 4.2 Matern 3/2 Covariance Kernel

```
k(s_i, s_j) = eta^2 * (1 + sqrt(3) * d / rho) * exp(-sqrt(3) * d / rho)
```

- eta^2 = 1.0 (fixed amplitude)
- rho = exp(log_ls_spatial) (estimated length scale in standardized coordinate units)
- Both GPs share a single length scale

**Code**: `4d_leaf_wax_spatial_model.stan:6-23`

### 4.3 Predictive Process Approximation

- m = 125 knots
- Projection: `K_cross * K_knots^{-1} * w_knots`
- Jitter for numerical stability: 1e-4

**Code**: `4d_leaf_wax_spatial_model.stan:319-329`

### 4.4 Knot Placement: Spherical Fibonacci Lattice

Deterministic equal-area global grid using the golden angle algorithm.

**Code**: `4a_spatial_functions.R:270-289`, `config.yaml:54`

### 4.5 Spatially Adaptive Regularization

Knot-level density computed within radius = 0.2 standardized units.
Density scaling factor = 10.

Regularization tiers:
- density = 0: tau = 0.50 (both intercept and slope)
- density < 10: tau = 0.50 + 0.30 * (density/10) (linear ramp)
- density >= 10: tau = 0.80 (cap)

Slope tau is further multiplied by a data-relative OIPC range factor:
- range < 0.25 * max_global_range: factor = 0.2
- range 0.25-0.60 * max_global_range: factor = 0.5
- range >= 0.60 * max_global_range: factor = 1.0

(Thresholds are relative to the maximum OIPC range across all knots, not fixed values in standardized units.)

**Code**: `4d_leaf_wax_spatial_model.stan:181-208`

#### [DISCREPANCY 3] Adaptive Regularization Formulas

Table S1 and code AGREE. Supplement Eqs 19-21 have stale formulas (different starting values, different ramp, exponential vs linear).

### 4.6 PC Priors for Spatial SDs

- Intercept GP: P(sigma_intercept > 20 permil) = 0.05
- Slope GP: P(sigma_slope > 0.3) = 0.05

**Code**: `config.yaml:42-51`, `4d_leaf_wax_spatial_model.stan:441-442`

### 4.7 GP Length Scale Prior

```
log_ls_spatial ~ Normal(-1.0, 0.4)
```
Bounds: `<lower=-2, upper=0>`

**Code**: `4d_leaf_wax_spatial_model.stan:243,438`

---

## 5. Complete Prior Table

| Parameter | Prior (Code) | Stan Line | Table S1 | Match? | Note |
|---|---|---|---|---|---|
| beta_0 (intercept) | Normal(0, 5) | 390 | Normal(0, 5) | YES | |
| beta_oipc (OIPC slope) | Normal(0.8, 0.3) | 391 | Normal(0.8, 0.3) | YES | |
| beta_c4 | Normal(0, 2) | 404 | Normal(0, 2) | YES | |
| **beta_precip** | **Normal(0, 0.5)** | **412** | **Normal(0, 0.5)** | **YES** | **Changed from 0.02** |
| beta_tree | Normal(0, 2) | 413 | Normal(0, 2) | YES | |
| beta_shrub | Normal(0, 2) | 414 | Normal(0, 2) | YES | |
| beta_grass | Normal(0, 2) | 415 | Normal(0, 2) | YES | |
| beta_oipc_x_c4 | Normal(0, 0.5) | 421 | Normal(0, 0.5) | YES | |
| beta_oipc_x_tree | Normal(0, 0.5) | 424 | Normal(0, 0.5) | YES | |
| beta_oipc_x_shrub | Normal(0, 0.5) | 425 | Normal(0, 0.5) | YES | |
| beta_oipc_x_grass | Normal(0, 0.5) | 426 | Normal(0, 0.5) | YES | |
| tau_elev_bspline | Normal+(0, 1) | 395 | Normal+(0, 1) | YES | |
| beta_elev[1] | Normal(0, 2) | 396 | Normal(0, tau_elev) | **NO** | [DISCREPANCY 2] |
| beta_elev[k>1] | Normal(beta_elev[k-1], tau_elev) | 398 | same | YES | |
| sigma (residual) | Normal+(0, 2) | 452 | Normal+(0, 2) | YES | |
| lambda_decay | Lognormal(2.5, 0.5) | 432 | Lognormal(2.5, 0.5) | YES | |
| log_ls_spatial | Normal(-1.0, 0.4) | 438 | Normal(-1.0, 0.4) | YES | |
| sigma_intercept | Exponential(PC prior) | 441 | Exponential(PC prior) | YES | |
| sigma_slope | Exponential(PC prior) | 442 | Exponential(PC prior) | YES | |
| z_intercept[k] | Normal(0, tau_int[k]) | 446 | same | YES | |
| z_slope[k] | Normal(0, tau_slope[k]) | 447 | same | YES | |

---

## 6. Model Configurations

### 6.1 Manuscript Models

15 model variants (14 manuscript comparison + 1 fixed-range variant) defined in `config.yaml:108-306`. Models evaluate combinations of covariates (OIPC, C4, PFT, elevation, precipitation), vegetation interactions (OIPC x C4, OIPC x PFT), and spatial varying coefficients (GP).

### 6.2 MCMC Settings

MCMC settings verified against actual chain CSV headers from the most recent successful run (Aug-Sep 2025):

| Setting | Non-spatial models | Spatial _veg models | Spatial _env/full models |
|---|---|---|---|
| Chains | 8 | 8 | 8 |
| Iterations (total) | 2000 | 3000 | 3000-4000 |
| Warmup ratio | 50% | 50% | 50% |
| adapt_delta | 0.95 | 0.95 | 0.99 |
| max_treedepth | 15 | 14 | 12 |

Seed: 314. Max parallel models: 4.

### 6.3 Convergence Criteria

- R-hat < 1.01
- ESS > 900
- No divergent transitions
- E-BFMI > 0.3

---

## 7. Data Standardization

| Variable | Standardization | Code |
|---|---|---|
| d2H_wax | (value - mean) / sd | `4b_stan_prep.R:628` |
| d2H_wax_err | value / sd(d2H_wax) | `4b_stan_prep.R:629` |
| OIPC d2H | (value - mean) / sd | `4a_spatial_functions.R:874` |
| OIPC SE | value / sd(OIPC) | `4a_spatial_functions.R:876` |
| Elevation | Convert to km, then (value - mean_km) / sd_km | `4a_spatial_functions.R:885-892` |
| C4 | (value - 20) / 25 (fixed constants from config) | `config.yaml:62-64` |
| Precipitation | (value - mean) / sd | `4a_spatial_functions.R:924-925` |
| PFT | Proportions (0-1), NOT standardized | |
| OIPC x C4 interaction | value / (sd_OIPC * sd_C4) | `4a_spatial_functions.R:984` |
| OIPC x PFT interactions | value / sd_OIPC | `4a_spatial_functions.R:991-993` |
| Coordinates (for GP) | (value - mean) / sd per dimension | |

### 7.1 NA Handling in Aggregated Matrices

After computing distance-weighted averages at 9 scales, some site-scale cells may have NA values (e.g., coastal sites where the extraction radius extends over ocean). These are filled with the site's mean across other scales, or 0 on the standardized scale (= population mean) if all scales are NA. PFT fractions default to 0.33 (equal split). OIPC SE defaults to 0.1 (small positive value for numerical stability).

The pipeline logs how many site-scale cells per matrix are NA-filled before proceeding to Stan.

**Code**: `4a_spatial_functions.R:1000-1075` (fill_na_row helper and loop)

---

## Discrepancy Summary

### DISCREPANCY 1: Multi-Scale Spatial Averaging Kernel

- **Manuscript** (Supp Eq 8): `w = exp(-d / r)` (exponential)
- **Code** (`4a_spatial_functions.R:19`): `w = exp(-d^2 / (2*r^2))` (Gaussian)
- Stan-side scale weighting uses exponential and IS consistent.

### DISCREPANCY 2: B-Spline First Coefficient Prior

- **Table S1**: `beta_elev[1] ~ Normal(0, tau_elev)`
- **Code** (`4d_leaf_wax_spatial_model.stan:396`): `beta_elev[1] ~ normal(0, 2)` — fixed SD

### DISCREPANCY 3: Adaptive Regularization Formulas

- **Code and Table S1**: AGREE
- **Supplement Eqs 19-21**: Stale formulas (different starting values, exponential vs linear ramp)
