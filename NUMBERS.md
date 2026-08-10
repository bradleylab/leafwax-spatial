# Reproducibility reference for the leaf-wax δ²H analysis

This file records the quantities regenerated from the fitted models used in the
associated manuscript. The spatial Gaussian process uses chordal distance
(3-D Euclidean distance on a sphere with radius 6,371 km), for which the
Matérn 3/2 covariance is positive definite by construction. Non-spatial
quantities do not depend on this distance metric.

## Provenance

- **Run:** `results/c2_run_20260728_chordal/model_output/` (14 models plus
  simulation and prior-sensitivity outputs).
  Authoritative manifest: `results/c2_run_20260728_chordal/RUN_MANIFEST_chordal.rds`.
- **Fit-time calibration snapshot:** `data/frozen/leafwax_d2h_c29_calibration_v1.csv`,
  md5 `bb52649130a02b9b7d897325899b4f5a`.
  The fit-time file labeled 27 Río Bermejo records with the Repasch et al. (2021)
  article DOI. The public compilation corrects those provenance fields to the 2020
  PANGAEA deposit (10.1594/PANGAEA.925616) and associates the records with Dosch
  et al. (2024; 10.5194/esurf-12-907-2024). Coordinates, isotope measurements,
  archive classes, and modeled values are identical. The
  corrected public calibration CSV has md5 `81d8276e28e4c05dc7784604e1120baf`
  and sha256 `34d90e3207f81056e486501e1da50498eaee6f58a4ddeffcd26a194272171ad5`.
- **Regenerated:** via `scripts/regen_tables.R` + `scripts/regen_manuscript_numbers.R`
  + `scripts/regen_simulation_numbers.R` with `LEAFWAX_RUN_DIR` pointed at the chordal run.
- Sediment for classical stats: `results/3_sediment_ready_for_modeling.rds` (1,128 rows).

## How to regenerate (replication)

All commands run from the analysis-repo root with the chordal run selected via
`LEAFWAX_RUN_DIR`. `MO` below denotes
`results/c2_run_20260728_chordal/model_output`. Generated files are written to
`model_analysis/reported_outputs/` by default. Set `LEAFWAX_OUTPUT_DIR` to use a
different directory; only the output path changes, never a value.

| Artifact | Command | Output |
|---|---|---|
| Core tables + reported scalars (slopes, variance, RMSE, CIs, C₃/C₄, vegetation coefficients) | `LEAFWAX_RUN_DIR=MO Rscript scripts/regen_tables.R` + `scripts/regen_manuscript_numbers.R` | `model_analysis/reported_outputs/` |
| Spatial intercept/slope component audit | `LEAFWAX_RUN_DIR=MO Rscript scripts/audit_spatial_variance_decomposition.R` | `scripts/reference_outputs/spatial_variance_decomposition_audit.csv` |
| Simulation and threshold quantities not in the table scripts (recovery + confounding intervals, intercept–OIPC correlations, vegetation effects, detection thresholds) | `LEAFWAX_RUN_DIR=MO Rscript scripts/regen_simulation_numbers.R` | `model_analysis/reported_outputs/SIMULATION_NUMBERS.md` plus CSV sources for Tables S4 and S5 |
| Precip-scale reconstruction SD + detection thresholds | `LEAFWAX_RUN_DIR=MO Rscript scripts/regen_precip_space.R` | `REGENERATED_NUMBERS_precip_space.md` |
| Table S3 (regional in-sample RMSE, 14×6) | `LEAFWAX_RUN_DIR=MO Rscript 5b_overfitting_diagnostics.R` | `model_analysis/reported_outputs/Table_S2_regional_performance.csv` |
| Table S5 (prior sensitivity, 7 variants) | fits of `baseline_veg_sp` via `6c_prior_sensitivity.R` | physical means 0.717–0.752 |
| Positive-definiteness verification | `LEAFWAX_RUN_DIR=MO Rscript scripts/verify_pd_knots.R` | `model_analysis/reported_outputs/PD_VERIFICATION.md` (minimum eigenvalue 2.8e-7) |

## Non-spatial results

- OLS Fig 1: β₀ **−121.90 ± 1.43**, β_OIPC **0.833 ± 0.019**, R² **0.643**, RMSE **23.0‰**,
  σ̂ 23.07‰, n 1,128, 95% PI half-width at x̄ 45.3‰.
- Moran's I (OLS resid, k=8) **0.583**; raw δ²H_wax **0.845**. Banded: 100 km 0.564 → 5,000 km 0.060.
- Mean nearest-neighbour **33.8 km** (median 3.9; 62.5% < 10 km; **92.2%** < 100 km).
- Dataset: n **1,128**; 49 unique DOIs; 904/1,128 (80.1%) measured (not raster) elevation.
  Continents: Asia 551, Americas 373, Africa 142, Oceania 32, Europe 30 = 1,128.

## Spatial-model results

- Non-spatial physical slope **0.81** (baseline 0.809, 0.778–0.841).
- **Spatial-adjusted physical-slope range 0.63–0.73** (full_sp 0.628 → baseline_veg_sp 0.726).
- Marginal spatial-component variance: the intercept represents **86–91%**
  (86.1–90.7) and the slope contribution **9–14%** (9.3–13.9) of their
  summed marginal variances. The slope contribution is evaluated with each
  posterior draw's fitted, scale-weighted OIPC predictor. These percentages
  exclude intercept–slope covariance and are not a partition of total
  δ²H_wax variance.
- RMSE **21 → 16** (baseline 21.2 → spatial 15.7–16.3).
- GP length scales: **3,457–3,699 km** (rounded to 3,460–3,700 km).
- Intercept–OIPC correlation **r = 0.28–0.41** (r² 0.08–0.17); baseline_sp 0.414 anchors the confounding test.
- Spatial intercept SD **12.2–17.2‰**; obs-level range −58 to +36‰ (full_sp widest).
- Spatial slope field in physical units: mean **0.61–0.71**, SD **≈0.13**,
  IQR **≈0.16**, range 0.18–1.02.
- Residual σ (spatial) **≈16.0‰**; 95% PI width **66.0‰** (spatial) vs **80.2‰** (non-spatial); reduction 17.7%.

## Per-model slopes (table2_global_params.csv), 95% CI

| model | β_OIPC | 95% CI | GP scale (km) |
|---|---:|---|---:|
| baseline (non-sp) | 0.809 | 0.778–0.841 | — |
| baseline_sp | 0.674 | 0.554–0.797 | 3,589 |
| baseline_env_sp | 0.641 | 0.483–0.795 | 3,699 |
| baseline_veg_sp | 0.726 | 0.570–0.887 | 3,457 |
| full_sp | 0.628 | 0.459–0.794 | 3,652 |
| full_interact_sp | 0.683 | 0.502–0.867 | 3,520 |
| elevation_only_sp | 0.692 | 0.545–0.834 | 3,558 |
| elevation_c4_sp | 0.690 | 0.537–0.839 | 3,577 |
| c4_only_sp | 0.671 | 0.546–0.794 | 3,596 |
| elevation_c4_interact_sp | 0.679 | 0.527–0.831 | 3,581 |

## Model comparison (table1_model_performance.csv, LOO)

- **Best model: full_sp** (LOOIC 1,305.6). Spatial variants span only ~54 LOOIC (1,305.6–1,359.8).
- baseline_sp ΔLOOIC **54.3 [SE 11.6]**; the five non-spatial models trail by ΔLOOIC **384–557 [SE 41–47]**.
- Effective parameters p_eff: spatial 44–58, non-spatial 3–18.

## C₃/C₄ and vegetation

**Main effects, full_sp (95% CI, table4_environmental.csv):**
- shrub **+0.264 [+0.019, +0.509]** resolved positive
- tree **−0.263 [−0.477, −0.056]** resolved negative
- grass **−0.193 [−0.362, −0.022]** resolved negative
- C4 +0.043 [−0.052, +0.136] not resolved

**δ²H_precip × PFT interactions, full_interact_sp (table5):**
- ×shrub **+0.265 [+0.078, +0.444]** resolved (steeper slope in shrubs)
- ×tree **−0.161 [−0.346, +0.021]** not resolved
- ×C4 +0.035 [−0.069, +0.141] not resolved
- ×grass −0.048 [−0.177, +0.079] not resolved
- Only the shrub interaction is resolved in this model.

## Convergence and divergences

- All models: max R̂ ≤ **1.008** (baseline_sp / baseline_env_sp); min bulk ESS **1,156** (full_interact),
  min spatial ESS **1,161** (elevation_only_sp); max-treedepth hits 0; min E-BFMI 0.711 (elevation_c4_sp).
- Divergences: **full_interact_sp 4, elevation_c4_sp 4, full_sp 1, elevation_c4_interact_sp 1**,
  all other 10 models **0**. Excluding divergent draws changes no reported quantity by more than 0.0009 SD.
- PD: chordal Matérn PD by construction; Cholesky no-jitter succeeded on all 1,800 draws/matrix in all 9
  spatial models; smallest eigenvalue 2.8e-7; Stan uses 1e-4 inference jitter.

## Uncertainty propagation / detection thresholds

- σ_residual (wax) median **16.0‰**; combined σ_total = √(15.98² + 3²) = **16.26‰** (analytical SD 3‰).
- Detection threshold (**δ²H_wax scale**), two-sample: ρ=0 **45.1**, ρ=0.5 **32.4**, ρ=0.8 **21.5**, ρ=0.9 **16.3‰**. The autocorrelation factor applies only to residual variance; independent analytical error is not reduced by ρ.
- **δ²H_precip scale** (per-draw propagation, `regen_precip_space.R`, reconstruction model **baseline_env_sp**):
  - single-sample reconstruction SD **25.5‰ [20.4, 33.7]**.
  - two-sample threshold: ρ=0 **70.6‰**, ρ=0.5 **50.8**, ρ=0.8 **33.7**, ρ=0.9 **25.5‰**.
  - Table S9 grid (fixed residual variance, wax-scale constant 45‰): slope 0.63→72,
    0.64→70 (fitted baseline_env_sp), 0.73→62, 0.81→56‰.
- Fig 6 endpoints (**δ²H_precip scale**, 95% confidence): spatial ρ=0.9 → **25.5‰**; non-spatial ρ=0 → **72.6‰** (each uses its model's physical global slope; `regen_precip_space.R`).

## Confounding simulation (Table S4 and Figure 3)

The knob ρ_c and the achieved correlation differ because two continental-scale smooth fields are
incidentally correlated over a finite site set (achieved r ≈ 0.40 already at ρ_c = 0). The empirical
scenario's knob (ρ_c = 0.011) is calibrated so its achieved r (0.41) matches the real-data
intercept–δ²H_precip correlation. Ordered by achieved correlation (GP posterior median, N = 1,128):

| achieved r | knob ρ_c | physical true slope | GP mean | 95% CI | bias | covers input |
|---:|---:|---:|---:|---|---:|:--:|
| 0.404 | 0.0 | 0.727 | 0.836 | [0.731, 0.949] | +0.109 | no |
| 0.414 (empirical) | 0.011 | 0.727 | 0.843 | [0.738, 0.956] | +0.115 | no |
| 0.622 | 0.3 | 0.727 | 1.081 | [0.979, 1.193] | +0.353 | no |
| 0.736 | 0.5 | 0.727 | 1.241 | [1.143, 1.347] | +0.514 | no |

- Empirical scenario: OLS **1.075**, GP **0.843**, true **0.727**; the spatial GP closes
  (1.075−0.843)/(1.075−0.727) ≈ **67%** of the OLS-to-truth gap, leaving **+0.115 (~33%)** residual.

## Simulated recovery — 95% CI (quantiles 0.025 and 0.975)

- 3a uniform 0.727 → **0.679 [0.572, 0.784]** (CI covers truth).
- 3b spatially varying mean 0.675 → **0.663 [0.556, 0.769]** (CI covers truth).
- 3c uniform 0.727, injected intercept–predictor r = 0.96 → **1.002 [0.915, 1.088]**
  (over-recovers and excludes the simulated value).

## Prior and hyperparameter sensitivity (`6c_prior_sensitivity.R`)

Reference plus seven alternative-prior variants of `baseline_veg_sp`;
physical β_OIPC posterior **mean** [95% CI]:

| variant | mean | 95% CI |
|---|---:|---|
| Reference | 0.726 | 0.570–0.887 |
| β_OIPC wider N(0.8, 1.0) | 0.721 | 0.556–0.890 |
| β_OIPC shifted N(0.5, 1.0) | 0.717 | 0.553–0.888 |
| β_OIPC uninformative N(0, 2.0) | 0.719 | 0.558–0.890 |
| σ_slope relaxed | 0.748 | 0.576–0.922 |
| σ_slope very relaxed | 0.752 | 0.560–0.939 |
| GP length scale longer N(8.27, 0.4) | 0.724 | 0.566–0.887 |
| GP length scale shorter N(7.27, 0.4) | 0.729 | 0.569–0.893 |

**Range of means 0.717–0.752**, all CIs overlapping — slope is prior-robust.

## Supplementary-table outputs

- **Table S3 (regional in-sample RMSE, 14 models × 5 regions + Overall)** — `5b_overfitting_diagnostics.R`
  → `model_analysis/reported_outputs/Table_S2_regional_performance.csv` (Overall n 1,128; Africa 142, Americas 373,
  Asia 551, Europe 30, Oceania 32).
- **Table S8 (global parameters)** — `regen_tables.R` → `table2_global_params_body.tex`
  (baseline_sp physical slope 0.674, GP scale 3,589 km).
- **Tables S2 / S9** — prior descriptor physical-slope range 0.63–0.73; detection-threshold grid above.
