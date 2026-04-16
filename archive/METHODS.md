# METHODS.md — Leaf Wax Spatial Calibration Model

> Cross-reference of manuscript methods and code implementation.
> Sources: `leafwax_gca_text_2025-12-04.docx`, `leafwax_gca_supplement_2025-12-02.docx`, `Table_S1_priors.pdf`, and GitHub repo `leafwax-spatial-master/` (Dec 2, 2025).
> All code paths are relative to `leafwax-spatial-master/`.

---

## 1. Data Compilation & Preprocessing

**Manuscript**: Section 2.1, Supplement S2.1

### 1.1 Sediment Data

- 818 n-C29 alkane d2H measurements from lake surface sediments
- Filtered: `chain == 29`, non-NA for d2H_wax, latitude, longitude
- **Code**: `3_prep_data.R:40-63`

### 1.2 Measurement Uncertainty (Supp S2.1.2)

| Condition | d2H_wax_err (permil) |
|---|---|
| Unreported (NA) | 3 |
| Reported as 0 | 1 |
| Otherwise | as reported |

**Code**: `3_prep_data.R:49-53`
```r
d2H_wax_err = case_when(
  is.na(d2H_wax_err) ~ 3,
  d2H_wax_err == 0 ~ 1,
  TRUE ~ d2H_wax_err
)
```

### 1.3 Coordinate Handling (Supp S2.1.3)

- Coordinate system: WGS84 (EPSG:4326)
- Gaussian jitter for duplicate coordinates: SD = 0.0001 degrees
- Jitter seed: 12345
- **Code**: `4b_stan_prep.R:88-97`, `config.yaml:71` (`coordinate_jitter: 0.0001`)

### 1.4 Environmental Covariates

| Variable | Source | Code |
|---|---|---|
| d2H_precip | OIPC (d2h_MA.tif) | `3_prep_data.R:73-74` |
| d2H_precip SE | OIPC (d2h_se_MA.tif) | `3_prep_data.R:74` |
| Elevation | GMTED (elevation_5KMmn_GMTEDmn.tif) | `3_prep_data.R:75` |
| PFT (tree/shrub/grass) | MODIS 3-class (2d_MODIS_PFT_3classes_Downsampled.tif) | `3_prep_data.R:76` |
| C4 fraction | NUS v2.2 (C4_distribution_NUS_v2.2.nc), 2001-2019 mean | `1_extract_c4_raster.R` |
| Annual precipitation | TerraClimate (ppt, 2001-2019 mean) | `3_prep_data.R:79` |
| Soil moisture | TerraClimate (soil, 2001-2019 mean) | `3_prep_data.R:80` |
| Max temperature | TerraClimate (tmax, 2001-2019 mean) | `3_prep_data.R:81` |
| VPD | TerraClimate (vpd, 2001-2019 mean), converted hPa to kPa | `3_prep_data.R:82,574` |

### 1.5 C4 Imputation (Supp S2.1.1)

Rules applied to NA C4 pixels:
- |latitude| >= 50 degrees -> C4 = 0
- Mean GMTED elevation >= 1500 m -> C4 = 0
- Otherwise: left as NA

**Code**: `3_prep_data.R:24-26` (thresholds), `3_prep_data.R:470-497` (imputation function)

### 1.6 Multi-Scale Spatial Averaging

**Manuscript** (Supp Eq 7-9): Covariates averaged at 9 distance scales using a weighting kernel.

Scales: {1, 3, 5, 10, 20, 40, 70, 100, 150} km
**Code**: `config.yaml:27`

Extraction radius: 5 degrees (`3_prep_data.R:24`)

#### **[DISCREPANCY 1]** Weighting Kernel

- **Manuscript** (Supp Eq 8): Exponential decay: `w = exp(-d / r)`
- **Code** (`4a_spatial_functions.R:18-19`): Gaussian kernel: `w = exp(-d^2 / (2 * r^2))`

Used throughout `compute_weighted_mean()` (`4a_spatial_functions.R:33-56`) and `compute_weighted_interaction()` (`4a_spatial_functions.R:59-101`):
```r
weights <- exp(-(dists_km^2) / (2 * scale_km^2))
```

**Impact**: Gaussian drops off faster at short distances but slower at long distances than exponential.

**Note**: The Stan-side scale weighting (Supp Eq 10) uses `exp(-r / lambda)` and IS consistent between manuscript and code (`4d_leaf_wax_spatial_model.stan:261-263`):
```stan
scale_weights[i] = exp(-distance_scales[i] / lambda_decay);
```

---

## 2. Exploratory Data Analysis

**Manuscript**: Section 2.2, Supplement S2.2

### 2.1 Spatial Autocorrelation

- Moran's I = 0.614 on OLS residuals (Section 2.2)
- **Code**: `3e_spatial_clustering_analysis.R`

### 2.2 Collinearity Screening

- VIF threshold = 5 (Section 2.2)
- Elastic net regularization: alpha = 0.5, 10-fold CV (Supp S2.2.2)
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

---

## 3. Statistical Model Specification

**Manuscript**: Section 2.3, Main Eq 1-2, Supplement Eq 12-13

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

**Code**: `4d_leaf_wax_spatial_model.stan:332-380` (linear predictor construction)

### 3.2 Multi-Scale Weighting (Supp Eq 10)

Each predictor is a weighted average across 9 spatial scales:
```
X_weighted = sum_s( w_s * X_s )
w_s = exp(-r_s / lambda) / sum( exp(-r_s / lambda) )
```

where lambda is the estimated or fixed decay parameter (in km).

**Code**: `4d_leaf_wax_spatial_model.stan:259-264`
```stan
scale_weights[i] = exp(-distance_scales[i] / lambda_decay);
scale_weights = scale_weights / sum(scale_weights);
```

Lambda prior: `lambda ~ lognormal(2.5, 0.5)` — i.e., median ~12 km, 95% CI ~4-28 km
Lambda bounds: `<lower=1, upper=400>` km

**Code**: `config.yaml:30-34`, `4d_leaf_wax_spatial_model.stan:240,426`

### 3.3 B-Spline for Elevation (Supp S2.3.2, Eq 14)

- 9 interior knots
- Degree 3 (cubic)
- Total basis functions: 9 + 3 + 1 = 13
- Interior knots placed at equally spaced quantiles within data range
- Random walk prior on coefficients (Supp Eq 14):
  ```
  beta_elev[k] ~ Normal(beta_elev[k-1], tau_elev)  for k = 2, ..., 13
  ```

#### **[DISCREPANCY 2]** First B-Spline Coefficient Prior

- **Table S1**: `beta_elev[1] ~ Normal(0, tau_elev)` — tied to estimated smoothness parameter
- **Code** (`4d_leaf_wax_spatial_model.stan:390`): `beta_elev_bspline[1] ~ normal(0, 2)` — fixed SD = 2

**Impact**: Code uses a fixed weakly informative prior for the first coefficient; Table S1 implies it's tied to the same smoothness parameter that governs the random walk.

**Code**: `config.yaml:58-59`, `4a_spatial_functions.R:134-153` (basis generation), `4d_leaf_wax_spatial_model.stan:388-393` (priors)

### 3.4 Three Uncertainty Components (Supp Eq 15)

```
sigma_total^2 = sigma_measurement^2 + sigma_OIPC^2 + sigma_residual^2
```

- `sigma_measurement` = d2H_wax_err (per-observation, standardized)
- `sigma_OIPC` = OIPC SE weighted across scales (per-observation, standardized)
- `sigma_residual` = sigma (estimated parameter)

**Code**: `4d_leaf_wax_spatial_model.stan:449-453`
```stan
real total_var = square(d2H_wax_err[n]) + square(oipc_se_weighted[n]) + square(sigma);
real total_sd = sqrt(total_var);
d2H_wax[n] ~ normal(mu[n], total_sd);
```

---

## 4. Spatial Framework

**Manuscript**: Section 2.4, Supplement S2.4

### 4.1 Gaussian Processes (Supp Eq 16)

Two independent GPs for spatially-varying coefficients:
- Intercept: `beta_0(s) = beta_0_global + GP_intercept(s)`
- OIPC slope: `beta_OIPC(s) = beta_OIPC_global + GP_slope(s)`

**Code**: `4d_leaf_wax_spatial_model.stan:308-330`

### 4.2 Matern 3/2 Covariance Kernel (Supp Eq 18)

```
k(s_i, s_j) = eta^2 * (1 + sqrt(3) * d / rho) * exp(-sqrt(3) * d / rho)
```

where:
- eta^2 = 1.0 (fixed amplitude) **[UNSPECIFIED IN MANUSCRIPT]**
- rho = exp(log_ls_spatial) (estimated length scale in standardized coordinate units)
- d = Euclidean distance in standardized coordinate space

**Code**: `4d_leaf_wax_spatial_model.stan:6-23` (kernel function), line 318:
```stan
matrix[n_pp_knots, n_pp_knots] K_knots = cov_matern32(knot_coords, 1.0, ls_spatial);
```

Both GPs share a single length scale.

### 4.3 Predictive Process Approximation (Supp S2.4.2, Eq 17)

- m = 125 knots
- Projection: `K_cross * K_knots^{-1} * w_knots`
- Jitter for numerical stability: 1e-4 added to diagonal

**Code**: `4d_leaf_wax_spatial_model.stan:319-329`
```stan
matrix[n_pp_knots, n_pp_knots] K_knots_jitter = K_knots +
    diag_matrix(rep_vector(1e-4, n_pp_knots));
alpha_spatial += K_cross * mdivide_left_spd(K_knots_jitter, knot_intercepts);
beta_oipc_spatial += K_cross * mdivide_left_spd(K_knots_jitter, knot_slopes);
```

### 4.4 Knot Placement: Spherical Fibonacci Lattice (Supp S2.4.4)

Deterministic equal-area global grid using the golden angle algorithm:
```
theta = golden_angle * (i - 1)
z = 1 - 2 * (i - 0.5) / n_knots
lat = asin(z) * 180 / pi
lon = (theta mod 2*pi) * 180 / pi - 180
```

**Code**: `4a_spatial_functions.R:270-289`, `config.yaml:54` (`knot_selection_method: "regular_global"`)

### 4.5 Spatially Adaptive Regularization (Supp S2.4.5)

Knot-level density computed within radius = 0.2 standardized units.
Density scaling factor = 10.

**Code**: `4a_spatial_functions.R:188-223` (density calculation), `config.yaml:37-39`

#### **[DISCREPANCY 3]** Regularization Formulas

**Supplement Eqs 19-21** (3-tier with exponential ramp):
- density = 0: tau_slope = 0.40, tau_intercept = 0.50
- density 1-10: tau_slope = 0.40 + 0.20*(n/10), tau_intercept = 0.50 + 0.20*(n/10)
- density > 10: tau_slope = 0.60 + 0.40*(1 - exp(-(n-10)/30)), tau_intercept = 0.70 + 0.30*(1 - exp(-(n-10)/30))

**Table S1**: `tau in [0.5, 0.8] * range_factor (slopes only)`

**Code** (`4d_leaf_wax_spatial_model.stan:178-202`):
- density = 0: tau = 0.50 (both slope and intercept)
- density < 10: tau = 0.50 + 0.30 * (d/10) (both)
- density >= 10: tau = 0.80 (both)
- Then slope tau multiplied by OIPC range factor:
  - range < 5: factor = 0.2
  - range 5-15: factor = 0.5
  - range >= 15: factor = 1.0

**Summary**: Table S1 and code agree. Supplement Eqs 19-21 have stale/different formulas (different starting values, different ramp, exponential vs linear).

OIPC range computed within 1000 km search radius (`4d_leaf_wax_spatial_model.stan:154-175`), minimum 5 data points required.

### 4.6 PC Priors for Spatial SDs (Supp S2.4.6)

Penalized complexity priors:
- Intercept GP: P(sigma_intercept > 20 permil) = 0.05 -> Exponential(lambda) where lambda = -log(0.05) / (20/d2H_sd)
- Slope GP: P(sigma_slope > 0.3) = 0.05 -> Exponential(lambda) where lambda = -log(0.05) / 0.3

Note: The intercept threshold is scaled by d2H_wax_sd for standardized space.

**Code**: `config.yaml:42-51`, `4b_stan_prep.R:1398-1404`, `4d_leaf_wax_spatial_model.stan:211-212,435-436`

### 4.7 GP Length Scale Prior

```
log_ls_spatial ~ Normal(-1.0, 0.4)
```
Bounds: `<lower=-2, upper=0>`

Conversion to km: `ls_km = exp(log_ls_spatial) * coord_scale_km`
where `coord_scale_km = mean(coord_scaling) * 111.0`

**Code**: `4d_leaf_wax_spatial_model.stan:243,432`, `4d_leaf_wax_spatial_model.stan:142`

---

## 5. Complete Prior Table

All priors cross-referenced between Table S1 and Stan code (`4d_leaf_wax_spatial_model.stan`).

| Parameter | Prior (Code) | Code Line | Table S1 | Match? | Justification |
|---|---|---|---|---|---|
| beta_0 (intercept) | Normal(0, 5) | 384 | Normal(0, 5) | YES | Weakly informative on standardized scale |
| beta_oipc (OIPC slope) | Normal(0.8, 0.3) | 385 | Normal(0.8, 0.3) | YES | Prior centered on expected apparent fractionation |
| beta_c4 | Normal(0, 2) | 398 | Normal(0, 2) | YES | Weakly informative |
| beta_precip | Normal(0, 0.02) | 402 | Normal(0, 0.02) | YES | Small effect expected relative to OIPC |
| beta_tree | Normal(0, 2) | 407 | Normal(0, 2) | YES | Weakly informative |
| beta_shrub | Normal(0, 2) | 408 | Normal(0, 2) | YES | Weakly informative |
| beta_grass | Normal(0, 2) | 409 | Normal(0, 2) | YES | Weakly informative |
| beta_oipc_x_c4 | Normal(0, 0.5) | 415 | Normal(0, 0.5) | YES | Moderate interaction expected |
| beta_oipc_x_tree | Normal(0, 0.5) | 418 | Normal(0, 0.5) | YES | Moderate interaction expected |
| beta_oipc_x_shrub | Normal(0, 0.5) | 419 | Normal(0, 0.5) | YES | Moderate interaction expected |
| beta_oipc_x_grass | Normal(0, 0.5) | 420 | Normal(0, 0.5) | YES | Moderate interaction expected |
| tau_elev_bspline | |Normal|(0, 1) (half-normal via <lower=0>) | 389 | Normal+(0, 1) | YES | Smoothness of elevation spline |
| beta_elev[1] | Normal(0, 2) | 390 | Normal(0, tau_elev) | **NO** | **[DISCREPANCY 2]** |
| beta_elev[k>1] | Normal(beta_elev[k-1], tau_elev) | 391-392 | Normal(beta_elev[k-1], tau_elev) | YES | Random walk |
| sigma (residual) | Normal+(0, 2) (half-normal via <lower=0>) | 446 | Normal+(0, 2) | YES | Weakly informative |
| lambda_decay | Lognormal(2.5, 0.5) | 426 | Lognormal(2.5, 0.5) | YES | Median ~12 km |
| log_ls_spatial | Normal(-1.0, 0.4) | 432 | Normal(-1.0, 0.4) | YES | Length scale in log-standardized units |
| sigma_intercept | Exponential(PC_lambda_int) | 435 | Exponential(PC prior) | YES | PC prior: P(sigma > 20 permil) = 0.05 |
| sigma_slope | Exponential(PC_lambda_slope) | 436 | Exponential(PC prior) | YES | PC prior: P(sigma > 0.3) = 0.05 |
| z_intercept_spatial[k] | Normal(0, tau_intercept[k]) | 440 | Normal(0, tau[k]) | YES | Density-based regularization |
| z_slope_spatial[k] | Normal(0, tau_slope[k]) | 441 | Normal(0, tau[k] * range_factor) | YES | Density + OIPC range factor |

---

## 6. Model Configurations

### 6.1 Manuscript Models

14 model variants described in manuscript (Supp S2.3.7), evaluating combinations of covariates and spatial effects.

### 6.2 Code Models

`config.yaml` defines 6 models (a subset rerun):

| Model | C4 | PFT | Elevation | Precip | Interactions | GP | Knots | Chains | Iter | adapt_delta | max_treedepth |
|---|---|---|---|---|---|---|---|---|---|---|---|
| baseline_veg | Y | Y | N | N | Y | N | 0 | 8 | 2000 | 0.95 | 15 |
| baseline_veg_sp | Y | Y | N | N | Y | Y | 125 | 8 | 3000 | 0.95 | 14 |
| full | Y | Y | Y | Y | N | N | 0 | 8 | 2000 | 0.95 | 15 |
| full_sp | Y | Y | Y | Y | N | Y | 125 | 8 | 4000 | 0.99 | 12 |
| full_interact | Y | Y | Y | Y | Y | N | 0 | 8 | 2000 | 0.95 | 15 |
| full_interact_sp | Y | Y | Y | Y | Y | Y | 125 | 8 | 4000 | 0.99 | 12 |

**Code**: `config.yaml:97-180`

### 6.3 MCMC Settings

| Setting | Value | Source |
|---|---|---|
| Chains | 8 | `config.yaml` per model |
| Iterations | 2000-4000 | `config.yaml` per model |
| Warmup ratio | 50% | `config.yaml:77` |
| adapt_delta | 0.95 (non-spatial) / 0.99 (spatial) | `config.yaml` per model **[UNSPECIFIED IN MANUSCRIPT]** |
| max_treedepth | 12-15 | `config.yaml` per model **[UNSPECIFIED IN MANUSCRIPT]** |
| Seed | 314 | `config.yaml:76` **[UNSPECIFIED IN MANUSCRIPT]** |
| Refresh | 100 | `config.yaml:78` |
| Max parallel models | 4 | `config.yaml:84` |

### 6.4 Convergence Criteria

- R-hat < 1.01
- ESS > 900
- No divergent transitions
- E-BFMI > 0.3

---

## 7. Data Standardization

**Code**: `4b_stan_prep.R:601-631` (initial), `4b_stan_prep.R:871-998` (aggregated values)

### 7.1 Standardization Rules

| Variable | Standardization | Code |
|---|---|---|
| d2H_wax | (value - mean) / sd | `4b_stan_prep.R:628` |
| d2H_wax_err | value / sd(d2H_wax) | `4b_stan_prep.R:629` |
| OIPC d2H | (value - mean) / sd | `4b_stan_prep.R:875` |
| OIPC SE | value / sd(OIPC) | `4b_stan_prep.R:876` |
| Elevation | Convert to km, then (value - mean_km) / sd_km | `4b_stan_prep.R:889-892` |
| C4 | (value - 20) / 25 (fixed constants) | `config.yaml:62-64`, `4b_stan_prep.R:880-882` |
| Precipitation | (value - mean) / sd (computed from data) | `4b_stan_prep.R:924-925` |
| PFT (tree/shrub/grass) | Proportions (0-1), NOT standardized | `4b_stan_prep.R:1346-1348` |
| OIPC x C4 interaction | value / (sd_OIPC * sd_C4) | `4b_stan_prep.R:984` |
| OIPC x PFT interactions | value / sd_OIPC | `4b_stan_prep.R:991-993` |
| Coordinates (for GP) | (value - mean) / sd per dimension | `4b_stan_prep.R:1279-1282` |

### 7.2 Scaling Parameters

Computed from the full dataset before model fitting (`4b_stan_prep.R:109-118`):
```r
SCALING_PARAMS <- list(
  d2H_mean = mean(sediment$d2H_wax),
  d2H_sd = sd(sediment$d2H_wax),
  oipc_mean = mean(sediment$oipc_d2h20),
  oipc_sd = sd(sediment$oipc_d2h20),
  elev_mean = mean(sediment$elevation_gmted),
  elev_sd = sd(sediment$elevation_gmted),
  c4_mean = 20,  # Fixed from config
  c4_sd = 25     # Fixed from config
)
```

### 7.3 Back-Transformation

**Code**: `4d_leaf_wax_spatial_model.stan:468-470`
```stan
real intercept_original = beta_0 * d2H_wax_sd_original + d2H_wax_mean_original;
real sigma_residual_original = sigma * d2H_wax_sd_original;
```

---

## 8. Model Validation

**Manuscript**: Section 2.5, Supplement S2.5

### 8.1 Leave-One-Region-Out Cross-Validation

6 geographic regions defined by coordinate boundaries:

| Region | Definition |
|---|---|
| Americas | lon < -30 |
| Europe | -30 <= lon < 60 AND lat > 35 |
| Africa | -30 <= lon < 60 AND lat <= 35 |
| Asia | 60 <= lon < 140 AND lat > -10 |
| Oceania | lon >= 140 OR lat < -10 |

**Code**: `5b_overfitting_diagnostics.R:59-65`

Note: The continent definitions in `4a_spatial_functions.R:316-337` use slightly different boundaries (e.g., separating North/South America, including Australia separately). The validation code uses the 5-region scheme above.

### 8.2 Model Comparison Metrics

- LOO-IC (leave-one-out information criterion via PSIS-LOO)
- WAIC
- Bayesian R-squared
- RMSE (in standardized and original scale)
- **Code**: `5a_model_validation.R`, `config.yaml:87-89`

### 8.3 Overfitting Diagnostics

- p_eff/n ratio (effective parameters / observations)
- Length scale vs mean nearest-neighbor distance
- Pareto k diagnostics (k > 0.7 flagged)
- Prediction interval width analysis
- **Code**: `5b_overfitting_diagnostics.R`

### 8.4 Convergence Criteria (applied in code)

- R-hat < 1.01
- ESS > 900
- No divergent transitions
- E-BFMI > 0.3

---

## 9. Bayesian Inversion

**Manuscript**: Supplement S4, Eq 22-24

### 9.1 Inversion Framework

Given observed d2H_wax at a downcore location, invert the forward model to estimate d2H_precip:

```
d2H_precip_inv = (d2H_wax - beta_0(s) - other_terms) / beta_OIPC(s)
```

### 9.2 Four Uncertainty Sources Propagated

1. **Analytical uncertainty**: measurement error on d2H_wax
2. **Model residual uncertainty**: sigma (residual SD)
3. **Parameter uncertainty**: posterior distributions of all coefficients
4. **GP prediction uncertainty**: spatial interpolation uncertainty at new locations

### 9.3 Key Assumptions

- **Stationarity**: Modern calibration relationships assumed to hold in the past
- **Ergodicity**: Spatial variability used as proxy for temporal variability

---

## Appendix A: Unspecified Items

Values that appear only in code, not in manuscript or Table S1:

| Item | Value | Code Location |
|---|---|---|
| adapt_delta (non-spatial) | 0.95 | `config.yaml` per model |
| adapt_delta (spatial) | 0.99 | `config.yaml` per model |
| max_treedepth | 12-15 depending on model | `config.yaml` per model |
| seed | 314 | `config.yaml:76` |
| kernel_jitter | 1e-4 | `config.yaml:70`, `4d_leaf_wax_spatial_model.stan:320` |
| coordinate_jitter SD | 0.0001 degrees | `config.yaml:71` |
| coordinate_jitter seed | 12345 | `4b_stan_prep.R:89` |
| density_radius | 0.2 standardized units | `config.yaml:38` |
| density_scaling | 10 | `config.yaml:39` |
| OIPC range factor breakpoints | <5, 5-15, >=15 | `4d_leaf_wax_spatial_model.stan:194-199` |
| OIPC range search radius | 1000 km | `4d_leaf_wax_spatial_model.stan:154` |
| OIPC range min count | >5 points required | `4d_leaf_wax_spatial_model.stan:173` |
| GP amplitude (eta^2) | 1.0 (fixed) | `4d_leaf_wax_spatial_model.stan:318` |
| lambda_decay bounds | 1-400 km | `4d_leaf_wax_spatial_model.stan:240` |
| log_ls_spatial bounds | -2 to 0 | `4d_leaf_wax_spatial_model.stan:243` |
| max_condition_number | 1e6 | `config.yaml:72` |
| min_eigenvalue threshold | -1e-6 | `config.yaml:73` |
| B-spline basis centering | center=TRUE, scale=FALSE | `4b_stan_prep.R:1240` |
| Raster extraction radius | 5 degrees | `3_prep_data.R:24` |

---

## Appendix B: Discrepancy Summary

### DISCREPANCY 1: Multi-Scale Spatial Averaging Kernel

- **Manuscript** (Supp Eq 8): `w = exp(-d / r)` — exponential decay
- **Code** (`4a_spatial_functions.R:19`): `w = exp(-d^2 / (2*r^2))` — Gaussian kernel
- **Impact**: Gaussian drops off faster at short distances, slower at long distances
- **Scope**: Affects all pre-computed weighted averages of covariates

### DISCREPANCY 2: B-Spline First Coefficient Prior

- **Table S1**: `beta_elev[1] ~ Normal(0, tau_elev)` — tied to estimated smoothness
- **Code** (`4d_leaf_wax_spatial_model.stan:390`): `beta_elev[1] ~ normal(0, 2)` — fixed SD
- **Impact**: Code uses a fixed weakly informative prior; Table S1 implies adaptivity

### DISCREPANCY 3: Adaptive Regularization Formulas

- **Supplement Eqs 19-21**: 3-tier with different starting values and exponential ramp
  - density=0: tau_slope=0.40, tau_intercept=0.50
  - density 1-10: tau_slope=0.40+0.20(n/10), tau_intercept=0.50+0.20(n/10)
  - density>10: tau_slope=0.60+0.40(1-exp(-(n-10)/30)), tau_intercept=0.70+0.30(1-exp(-(n-10)/30))
- **Table S1**: `tau in [0.5, 0.8] * range_factor (slopes only)`
- **Code** (`4d_leaf_wax_spatial_model.stan:178-202`):
  - density=0: tau=0.50 (both)
  - density<10: tau=0.50+0.30*(d/10) (both)
  - density>=10: tau=0.80 (both)
  - slope tau *= range_factor: <5->0.2, 5-15->0.5, >=15->1.0
- **Summary**: Table S1 and code AGREE. Supplement Eqs 19-21 have stale formulas.

---

## Appendix C: Code File Index

| File | Purpose |
|---|---|
| `0_load_config.R` | Load config.yaml |
| `1_extract_c4_raster.R` | C4 vegetation fraction from NetCDF |
| `3_prep_data.R` | Data compilation, raster extraction, C4 imputation |
| `3d_analyze_collinearity.R` | VIF, PCA, elastic net for variable selection |
| `3e_spatial_clustering_analysis.R` | Moran's I, spatial autocorrelation |
| `4a_spatial_functions.R` | Spatial weighting, B-spline basis, knot selection, validation |
| `4b_stan_prep.R` | Standardization, Stan data assembly |
| `4c_fit_models.R` | Model fitting orchestration |
| `4d_leaf_wax_spatial_model.stan` | Full Bayesian model (Stan) |
| `5a_model_validation.R` | LOO-CV, WAIC, posterior predictive checks |
| `5b_overfitting_diagnostics.R` | Overfitting diagnostics, spatial CV |
| `config.yaml` | All hyperparameters, model configurations |
