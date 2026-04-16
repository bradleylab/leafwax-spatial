# Baseline Snapshot — Pre-Validation State

**Date**: 2026-03-27
**Source**: EC2 instance i-053cbbe7833ef0b97 (r6i.large)
**Working directory**: `/home/shared/leafwax_spatial/spatial_leafwax_model`
**Git commit**: `08a4def` (branch: main, repo: bradleylab/leafwax-spatial)
**Stan model version**: `updated_bspline_matern` (Matern 3/2 kernel)
**N observations**: 818

## EC2 Code vs GCA Local Snapshot

The EC2 `config.yaml` currently defines **16 models** (systematic factorial design), while the GCA local snapshot has 6 models. The Sep 24-25, 2025 model runs used a 6-model config matching the GCA submission (saved in `prepared_data_consolidated/config_*.rds`). The `config.yaml` was later updated to the 16-model format.

Key differences between EC2 config.yaml (current) and GCA config (used for Sep 2025 runs):
- EC2 current: 12 spatial scales [1,3,5,10,20,40,70,100,150,220,300,400] vs GCA: 9 scales [1,3,5,10,20,40,70,100,150]
- EC2 current: 120 knots vs GCA run: 125 knots
- EC2 current: elevation_spline_degree=2 vs GCA run: degree=3
- EC2 current: density_scaling=50 vs GCA run: 10
- EC2 current: seed=123 vs GCA run: 314
- EC2 current: no PC prior config section vs GCA: explicit PC priors

**The Sep 2025 model runs used the GCA-style config** (125 knots, 9 scales, PC priors, etc.), which is preserved in the `config_*.rds` files.

## Model Results — beta_oipc (OIPC Slope) Across All Models

### Non-spatial models (no GP)
| Model | beta_oipc mean | 95% CI | Runtime |
|-------|---------------|--------|---------|
| baseline | 0.824 | [0.789, 0.861] | <1 min |
| baseline_veg | 0.814 | [0.772, 0.859] | 5 min |
| baseline_env | 0.839 | [0.798, 0.881] | 59 min |
| full | 0.834 | [0.785, 0.883] | 109 min |
| full_interact | 0.833 | [0.785, 0.882] | 139 min |

### Spatial models (with GP, no elevation)
| Model | beta_oipc mean | 95% CI | Runtime |
|-------|---------------|--------|---------|
| baseline_sp | 0.612 | [0.515, 0.717] | 311 min |
| c4_only_sp | 0.605 | [0.512, 0.704] | 302 min |
| baseline_veg_sp | 0.569 | [0.461, 0.682] | 945 min |

### Spatial models (with GP + elevation)
| Model | beta_oipc mean | 95% CI | Runtime |
|-------|---------------|--------|---------|
| elevation_only_sp | 0.502 | [0.368, 0.630] | 2000 min |
| elevation_c4_sp | 0.504 | [0.369, 0.630] | 1854 min |
| elevation_c4_interact_sp | 0.505 | [0.374, 0.630] | 2207 min |
| baseline_env_sp | 0.498 | [0.364, 0.621] | 1842 min |

### Spatial models (with GP + elevation + all covariates)
| Model | beta_oipc mean | 95% CI | Runtime |
|-------|---------------|--------|---------|
| full_sp | 0.448 | [0.307, 0.580] | 1608 min |
| full_interact_sp | 0.446 | [0.306, 0.580] | 1798 min |

## Key Observations

1. **Non-spatial models**: beta_oipc consistently ~0.82-0.84 regardless of covariates
2. **GP alone**: drops slope to ~0.61 (25% reduction)
3. **GP + elevation**: drops to ~0.50 (39% reduction)
4. **GP + elevation + veg + precip**: drops to ~0.45 (45% reduction)
5. **Pattern**: The slope reduction occurs in two stages — GP accounts for ~0.22 drop, elevation for another ~0.11, other covariates for another ~0.05

## GP Parameters (baseline_veg_sp)

- **Length scale**: mean = 4,247 km (shared between intercept and slope GPs)
- **sigma_intercept_spatial**: mean = 45.08 permil
- **sigma_slope_spatial**: mean = 0.69 (very wide CI: [0.007, 2.05])

## Priors Used (from Stan model on EC2)

| Parameter | Prior | Notes |
|-----------|-------|-------|
| beta_0 | Normal(0, 5) | |
| beta_oipc | Normal(0.8, 0.3) | Centered on expected fractionation |
| beta_c4 | Normal(0, 2) | |
| beta_tree/shrub/grass | Normal(0, 2) | |
| beta_oipc_x_* | Normal(0, 0.5) | Interaction terms |
| tau_elev_bspline | Normal+(0, 1) | Half-normal |
| beta_elev[1] | Normal(0, 2) | Fixed SD (DISCREPANCY 2 from METHODS.md) |
| sigma | Normal+(0, 2) | Half-normal |
| lambda_decay | Lognormal(2.5, 0.5) | Bounds [1, 400] km |
| log_ls_spatial | Normal(-1.0, 0.4) | Bounds [-2, 0] |
| sigma_intercept | Exponential(PC_lambda) | PC: P(sigma > 20‰) = 0.05 |
| sigma_slope | Exponential(PC_lambda) | PC: P(sigma > 0.3) = 0.05 |
| z_intercept_spatial[k] | Normal(0, tau_intercept[k]) | Density-based |
| z_slope_spatial[k] | Normal(0, tau_slope[k] * range_factor) | Density + OIPC range |

## Adaptive Regularization (from Stan model)

tau values: density-dependent in [0.5, 0.8], then slope tau multiplied by OIPC range factor:
- OIPC range < 5 (standardized): factor = 0.2
- OIPC range 5-15: factor = 0.5
- OIPC range >= 15: factor = 1.0

## Scaling Parameters (for back-transformation)

- d2H_wax: mean = -179.25, sd = 38.74
- OIPC d2H: mean = -65.28, sd = 39.04
- Elevation: mean = 1265.86 m, sd = 1073.55 m
- C4: mean = 20, sd = 25 (fixed)
- Precipitation: mean = 790.87 mm, sd = 552.91 mm
