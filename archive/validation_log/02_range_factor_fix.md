# Validation 02: Range Factor Threshold Fix

**Date**: 2026-03-28
**Git commit**: `fe6b9ba` (fix), `22d5350` (scripts)
**Model**: baseline_veg_sp with data-relative OIPC range thresholds
**Runtime**: ~6 hrs on r6i.8xlarge (some chains slower than others)

## Change Made

Replaced fixed OIPC range thresholds (5/15 standardized units — unreachable) with
data-relative thresholds (25%/60% of max observed OIPC range at knots).

Stan model `4d_leaf_wax_spatial_model.stan` lines 193-199:
- Old: `if (oipc_range < 5)` → factor=0.2 for ALL 125 knots
- New: `if (oipc_range < 0.25 * max_range)` → factor=0.2 for 49 knots, 0.5 for 32, 1.0 for 44

New thresholds (standardized): low = 0.60, high = 1.44

## Results Comparison

| Parameter | Original (frozen GP) | Fixed thresholds | Change |
|---|---|---|---|
| beta_oipc mean | 0.569 | 0.555 | -0.014 (negligible) |
| beta_oipc 95% CI | [0.461, 0.682] | [0.415, 0.686] | Slightly wider |
| sigma_slope_spatial mean | 0.688 | 0.547 | Better identified |
| sigma_slope_spatial 95% CI | [0.007, 2.05] | [0.216, 1.08] | Much tighter |
| sigma_intercept_spatial | 45.08 | 48.33 | +3.25 permil |
| Length scale (km) | 4,247 | 4,309 | +62 km (negligible) |
| RMSE | — | 0.398 | — |
| R² | — | 0.842 | — |

## Interpretation

**The slope attenuation from 0.83 (non-spatial) to ~0.56 (spatial) is robust to the range factor fix.**

Unfreezing the slope GP had minimal effect on beta_oipc (-0.014, well within posterior uncertainty). The slope GP is now better identified (sigma_slope CI narrowed from [0.007, 2.05] to [0.22, 1.08]), confirming that it has *some* spatial variation to express, but this variation does not substantially change the global slope estimate.

The intercept GP remains the primary driver of slope attenuation. This is consistent with the spatial confounding interpretation: OIPC has strong spatial structure, and the intercept GP absorbs spatially-correlated variance that OLS attributes to the OIPC slope.

## Implication

The fixed range thresholds should be kept (they are more principled than the old unreachable thresholds), but the manuscript's core finding — that spatial models reduce the apparent OIPC slope by ~30-45% — is not an artifact of the slope GP regularization.
