# Diagnostic 01: OIPC-Elevation Correlation & Slope GP Regularization

**Date**: 2026-03-27
**Purpose**: Understand the two-stage slope attenuation pattern and check slope GP regularization behavior.

## OIPC-Elevation Correlation

Data from `stan_data_full_sp.rds` (818 observations, standardized).

| Scale (km) | r(OIPC, Elevation) |
|------------|-------------------|
| 1 | -0.379 |
| 3 | -0.383 |
| 5 | -0.383 |
| 10 | -0.384 |
| 20 | -0.382 |
| 40 | -0.374 |
| 70 | -0.350 |
| 100 | -0.322 |
| 150 | -0.285 |

OIPC and elevation are moderately negatively correlated (higher elevation → more depleted precipitation δ²H). OLS R² = 0.144 (elevation explains ~14% of OIPC variance).

## Partial Correlations

| Relationship | r |
|---|---|
| d2Hwax ~ OIPC | 0.839 |
| d2Hwax ~ Elevation | -0.242 |
| d2Hwax ~ OIPC \| Elevation (partial) | 0.832 |

The partial correlation barely changes (0.839 → 0.832), meaning elevation does NOT substantially mediate the linear d2Hwax-OIPC relationship in OLS. The 0.11 slope drop in the Bayesian spatial models (from 0.61 with GP alone to 0.50 with GP+elevation) is therefore NOT explained by simple linear collinearity.

## OLS Slopes

| Model | OIPC slope | R² |
|---|---|---|
| d2Hwax ~ OIPC | 0.824 | 0.703 |
| d2Hwax ~ OIPC + Elevation | 0.858 | 0.710 |
| d2Hwax ~ OIPC + Elevation + Elevation² | 0.860 | 0.711 |

**In OLS, adding elevation slightly INCREASES the OIPC slope** (0.824 → 0.858). This is the opposite of what happens in the Bayesian spatial models. The elevation-induced slope drop from 0.61 to 0.50 is therefore a spatial model phenomenon, not a simple confounding effect.

**Interpretation**: The elevation B-spline in the spatial model can capture non-linear, spatially-structured relationships between elevation and d2Hwax that interact with the GP intercept. When both the GP and elevation spline are present, they jointly explain more of the variance that was previously attributed to the OIPC slope.

## CRITICAL FINDING: Slope GP is Effectively Disabled

### OIPC Range at Knots

The adaptive regularization multiplies `tau_slope` by a range_factor based on OIPC range within 1000 km of each knot:
- OIPC range < 5 (standardized): factor = 0.2
- OIPC range 5-15: factor = 0.5
- OIPC range >= 15: factor = 1.0

**Results**:
- Global OIPC range in data: 4.97 standardized units (194 permil)
- Threshold for factor=0.2: 5 standardized units = 195 permil
- **The threshold exceeds the entire global range of the data**
- **125/125 knots have range_factor = 0.2**
- **125/125 knots have final tau_slope < 0.20**

### Tau Slope Values
```
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
 0.1000  0.1000  0.1000  0.1112  0.1060  0.1600
```

All knots have `tau_slope` between 0.10 and 0.16. This means `z_slope_spatial[k] ~ Normal(0, 0.10-0.16)` — the slope random effects are constrained to be essentially zero.

### Data Density at Knots
```
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
   0.00    0.00    0.00    4.92    1.00  100.00
```
- 91/125 knots have zero nearby data points
- Most knots are in ocean or data-sparse regions (Fibonacci sphere is global)

### Consequence

The slope GP (`sigma_slope_spatial`) has a posterior mean of 0.69 but 95% CI [0.007, 2.05] — it is unidentified because the regularization prevents it from doing anything. The observed slope attenuation (0.83 → 0.45) is **entirely driven by the intercept GP and elevation spline**, not by the slope GP.

This is consistent with the manuscript's interpretation that the intercept GP captures spatially-structured variance. The slope GP was intended to allow spatially-varying OIPC relationships, but the regularization effectively prevents this.

## Implication for Validation Plan

The OIPC range thresholds (5/15 standardized units) were set in a regime where they are unreachable given the actual data range. This is likely an oversight rather than an intentional design choice. Two options:

1. **Accept**: The slope GP is effectively off. The model is really: global OIPC slope + spatially-varying intercept + covariates. This is a simpler and more interpretable model. Report it honestly.

2. **Investigate**: Relax the range thresholds (e.g., use percentiles of the actual data range rather than fixed values) and see if allowing spatial slope variation changes the results or model fit.

Option 2 should be part of the sensitivity analysis regardless, since reviewers may ask about it.
