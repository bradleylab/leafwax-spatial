# NUMBERS.md — authoritative chordal-run quantities for the leaf-wax δ²H manuscript

**Single source of truth for every quotable number.** All values below are
regenerated from the chordal-metric refit; the manuscript (main text + supplement)
must draw from this file, not from older drafts, the frozen great-circle run, or
the GCA numbers.

The model's spatial GP metric changed from great-circle (standardized coordinate
units) to **chordal** (3-D Euclidean distance on the sphere, R = 6,371 km), so the
Matérn 3/2 kernel is positive definite by construction. Every spatial quantity below
is the chordal value; non-spatial quantities (OLS, Moran, nearest-neighbour) are
metric-independent and unchanged.

## Provenance

- **Run:** `results/c2_run_20260728_chordal/model_output/` (14 models, chordal metric;
  Tier-C simulations confounding_v2 + simrecovery + prior sensitivity refit chordal).
  Authoritative manifest: `results/c2_run_20260728_chordal/RUN_MANIFEST_chordal.rds`.
- **Calibration input:** `data/frozen/leafwax_d2h_c29_calibration_v1.csv`, md5
  `bb52649130a02b9b7d897325899b4f5a` (unchanged; the metric change does not touch the data).
- **Regenerated:** via `scripts/regen_tables_v10.R` + `scripts/regen_manuscript_numbers.R`
  + `scripts/phase_c_numbers.R` with `LEAFWAX_RUN_DIR` pointed at the chordal run.
- Sediment for classical stats: `results/3_sediment_ready_for_modeling.rds` (1,128 rows).

## How to regenerate (replication)

All commands run from the analysis-repo root with the chordal run selected via
`LEAFWAX_RUN_DIR`. `MO` below = `results/c2_run_20260728_chordal/model_output`,
`RUN` = `results/c2_run_20260728_chordal`. Set `LEAFWAX_RETRACE_OUT_DIR` to redirect
outputs into a sandbox instead of overwriting the manuscript; only the output path
changes, never a value.

| Artifact | Command | Output |
|---|---|---|
| Core tables + main-text scalars (slopes, variance, RMSE, CIs, C₃/C₄, veg coeffs) | `LEAFWAX_RUN_DIR=MO Rscript scripts/regen_tables_v10.R` + `scripts/regen_manuscript_numbers.R` | `*/tables/*.csv`, `REGENERATED_NUMBERS_v10.md` |
| Chordal-era numbers not in the table scripts (recovery + confounding intervals, intercept–OIPC correlations, veg effects, detection thresholds) | `LEAFWAX_RUN_DIR=MO Rscript scripts/phase_c_numbers.R` | `retrace/chordal/PHASE_C_NUMBERS.md` |
| Precip-scale reconstruction SD + detection thresholds | `LEAFWAX_RUN_DIR=MO Rscript scripts/regen_precip_space.R` | `REGENERATED_NUMBERS_precip_space_v10.md` |
| Table S3 (regional in-sample RMSE, 14×6) | `LEAFWAX_RUN_DIR=MO Rscript 5b_overfitting_diagnostics.R` | `manuscript/tables/Table_S2_regional_performance.csv` |
| Table S4 + Fig 3 (confounding stress test) | `LEAFWAX_RUN_DIR=MO Rscript manuscript/drafts/comms_ee/analysis/make_figures.R NG3` | prints "Figure NG3 source values"; `figures/Figure_3.pdf` |
| Table S5 (prior sensitivity, 7 variants) | chordal refits of `baseline_veg_sp` via `6c_prior_sensitivity.R` | means 0.690–0.724 |
| PD verification (chordal Matérn PD by construction) | `LEAFWAX_RUN_DIR=MO Rscript scripts/verify_pd_knots.R` | `retrace/chordal/PD_VERIFICATION.md` (min eig 2.8e-7) |

## Key chordal changes vs the frozen great-circle run

The claim-affecting shifts (see `retrace/map_main.md` / `retrace/map_supp.md` for the full ~195-number map):

- **Spatial slope range 0.53–0.62 → 0.60–0.70** (still well below the non-spatial 0.78).
- **Best-LOO model flips** full_interact_sp → **full_sp**.
- **Two credible-interval resolution flips (C4/veg):** shrub main effect resolves *positive*
  (was unresolved); the δ²H_p×tree interaction becomes *unresolved* (was resolved).
- **GP length scale 3,610–3,930 → 3,460–3,700 km.**
- **PD argument reframed:** chordal Matérn is PD by construction (Banerjee 2005; Gneiting 2013);
  the empirical eigenvalue check is now confirmation, not the primary defence.
- Detection thresholds ~10% lower; confounding recovery re-anchored to the achieved correlation.

## Non-spatial / metric-independent quantities (unchanged from frozen)

- OLS Fig 1: β₀ **−121.90 ± 1.43**, β_OIPC **0.833 ± 0.019**, R² **0.643**, RMSE **23.0‰**,
  σ̂ 23.07‰, n 1,128, 95% PI half-width at x̄ 45.3‰.
- Moran's I (OLS resid, k=8) **0.583**; raw δ²H_wax **0.845**. Banded: 100 km 0.564 → 5,000 km 0.060.
- Mean nearest-neighbour **33.8 km** (median 3.9; 62.5% < 10 km; **92.2%** < 100 km).
- Dataset: n **1,128**; 49 unique DOIs; 904/1,128 (80.1%) measured (not raster) elevation.
  Continents: Asia 551, Americas 373, Africa 142, Oceania 32, Europe 30 = 1,128.

## Chordal spatial results (headline)

- Non-spatial slope **0.78** (baseline 0.779, 0.749–0.809).
- **Spatial-adjusted slope range 0.60–0.70** (full_sp 0.604 → baseline_veg_sp 0.699).
- Spatial variance **~42–60%** (41.7–59.6); intercept share of spatial variance **~89–94%** (88.7–93.9).
- RMSE **21 → 16** (baseline 21.2 → spatial 15.7–16.3).
- GP length scales: chordal range **3,457–3,699 km** (quote "3,460–3,700 km").
- Intercept–OIPC correlation **r = 0.28–0.41** (r² 0.08–0.17); baseline_sp 0.414 anchors the confounding test.
- Spatial intercept SD **12.2–17.2‰**; obs-level range −58 to +36‰ (full_sp widest).
- Spatial slope field: mean **0.58–0.68**, SD **≈0.13**, IQR **≈0.16**, range 0.17–0.98.
- Residual σ (spatial) **≈16.0‰**; 95% PI width **66.0‰** (spatial) vs **80.2‰** (non-spatial); reduction 17.7%.

## Per-model slopes (table2_global_params.csv), 95% CI

| model | β_OIPC | 95% CI | GP scale (km) |
|---|---:|---|---:|
| baseline (non-sp) | 0.779 | 0.749–0.809 | — |
| baseline_sp | 0.649 | 0.533–0.767 | 3,589 |
| baseline_env_sp | 0.617 | 0.465–0.765 | 3,699 |
| baseline_veg_sp | 0.699 | 0.549–0.854 | 3,457 |
| full_sp | 0.604 | 0.442–0.764 | 3,652 |
| full_interact_sp | 0.657 | 0.483–0.835 | 3,520 |
| elevation_only_sp | 0.666 | 0.524–0.803 | 3,558 |
| elevation_c4_sp | 0.664 | 0.517–0.808 | 3,577 |
| c4_only_sp | 0.646 | 0.526–0.765 | 3,596 |
| elevation_c4_interact_sp | 0.654 | 0.508–0.800 | 3,581 |

## Model comparison (table1_model_performance.csv, LOO)

- **Best model: full_sp** (LOOIC 1,305.6). Spatial variants span only ~54 LOOIC (1,305.6–1,359.8).
- baseline_sp ΔLOOIC **54.3 [SE 11.6]**; the five non-spatial models trail by ΔLOOIC **384–557 [SE 41–47]**.
- Effective parameters p_eff: spatial 44–58, non-spatial 3–18.

## C₃/C₄ + vegetation (reframe-critical)

**Main effects, full_sp (95% CI, table4_environmental.csv):**
- shrub **+0.264 [+0.019, +0.509]** resolved **positive** ← flip vs frozen (was unresolved)
- tree **−0.263 [−0.477, −0.056]** resolved negative
- grass **−0.193 [−0.362, −0.022]** resolved negative
- C4 +0.043 [−0.052, +0.136] not resolved

**δ²H_precip × PFT interactions, full_interact_sp (table5):**
- ×shrub **+0.265 [+0.078, +0.444]** resolved (steeper slope in shrubs)
- ×tree **−0.161 [−0.346, +0.021]** NOT resolved ← flip vs frozen (was resolved −0.241)
- ×C4 +0.035 [−0.069, +0.141] not resolved
- ×grass −0.048 [−0.177, +0.079] not resolved
- ⇒ PFT modulates the slope across woody/herbaceous lines; the tree-interaction claim no longer resolves.

## Convergence / divergences (chordal, adapt_delta 0.999 on the four models that needed it)

- All models: max R̂ ≤ **1.008** (baseline_sp / baseline_env_sp); min bulk ESS **1,156** (full_interact),
  min spatial ESS **1,161** (elevation_only_sp); max-treedepth hits 0; min E-BFMI 0.711 (elevation_c4_sp).
- Divergences: **full_interact_sp 4, elevation_c4_sp 4, full_sp 1, elevation_c4_interact_sp 1**,
  all other 10 models **0**. Divergence-sensitivity (Tier B): max effect on any reported quantity 0.0009 SD.
- PD: chordal Matérn PD by construction; Cholesky no-jitter succeeded on all 1,800 draws/matrix in all 9
  spatial models; smallest eigenvalue 2.8e-7; Stan uses 1e-4 inference jitter.

## Uncertainty propagation / detection thresholds

- σ_residual (wax) median **16.0‰**; combined σ_total = √(15.98² + 3²) = **16.26‰** (analytical SD 3‰).
- Detection threshold (**δ²H_wax scale**), two-sample: ρ=0 **45.1**, ρ=0.5 **31.9**, ρ=0.8 **20.2**, ρ=0.9 **14.3‰**.
- **δ²H_precip scale** (per-draw propagation, `regen_precip_space.R`, reconstruction model **baseline_env_sp**):
  - single-sample reconstruction SD **26.5‰ [21.2, 35.0]** → quote "≈26‰ (21–35‰)".
  - two-sample threshold: ρ=0 **73.4‰**, ρ=0.5 **52.8**, ρ=0.8 **35.0**, ρ=0.9 **26.5‰** → headline "≈73‰".
  - Table S9 grid (fixed residual variance, wax-scale constant 45‰): slope 0.50→90, 0.60→75,
    0.62→73 (fitted baseline_env_sp), 0.70→64, 0.78→58‰.
- Fig 6 endpoints (wax scale): spatial ρ=0.9 → **14.3‰**; non-spatial ρ=0 → **58.7‰** (σ_baseline 20.96‰).

## Confounding simulation (Table S4 / Fig 3, chordal confounding_v2) — verified

The knob ρ_c and the achieved correlation differ because two continental-scale smooth fields are
incidentally correlated over a finite site set (achieved r ≈ 0.40 already at ρ_c = 0). The empirical
scenario's knob (ρ_c = 0.011) is calibrated so its achieved r (0.41) matches the real-data
intercept–δ²H_precip correlation. Ordered by achieved correlation (GP posterior median, N = 1,128):

| achieved r | knob ρ_c | input β (true, re-std) | GP median | 95% CI | bias | covers input |
|---:|---:|---:|---:|---|---:|:--:|
| 0.404 | 0.0 | 0.506 | 0.582 | [0.509, 0.661] | +0.076 | no |
| 0.414 (empirical) | 0.011 | 0.504 | 0.584 | [0.511, 0.662] | +0.080 | no |
| 0.622 | 0.3 | 0.447 | 0.665 | [0.602, 0.734] | +0.217 | no |
| 0.736 | 0.5 | 0.422 | 0.720 | [0.663, 0.782] | +0.298 | no |

- Empirical headline (main text): OLS **0.748**, GP **0.584**, true **0.504**; the spatial GP closes
  (0.748−0.584)/(0.748−0.504) ≈ **67%** of the OLS-to-truth gap, leaving **+0.08 (~33%)** residual.
- The originally reported anchor ρ_c = 0.45 corresponds under this calibration to an achieved
  correlation of **0.71** (1.7× the observed 0.41) and is superseded.

## Simulated recovery (simrecovery, chordal) — 95% CI (quantile 0.025/0.975)

- 3a uniform 0.70 → **0.654 [0.551, 0.755]** (CI covers truth).
- 3b spatially varying mean 0.65 → **0.638 [0.535, 0.740]** (CI covers truth).
- 3c uniform 0.70, injected intercept–predictor r = 0.96 → **0.964 [0.880, 1.047]** (over-recovers,
  excludes truth) — the near-worst-case pilot motivating the graded confounding test.

## Prior / hyperparameter sensitivity (chordal, `6c_prior_sensitivity.R`)

7 variants of `baseline_veg_sp` refit chordal; β_OIPC posterior **mean** [95% CI]:

| variant | mean | 95% CI |
|---|---:|---|
| Reference | 0.698 | 0.549–0.854 |
| β_OIPC wider N(0.8, 1.0) | 0.694 | 0.535–0.857 |
| β_OIPC shifted N(0.5, 1.0) | 0.690 | 0.532–0.854 |
| β_OIPC uninformative N(0, 2.0) | 0.692 | 0.537–0.856 |
| σ_slope relaxed | 0.720 | 0.555–0.887 |
| σ_slope very relaxed | 0.724 | 0.539–0.903 |
| GP length scale longer N(8.27, 0.4) | 0.697 | 0.545–0.854 |
| GP length scale shorter N(7.27, 0.4) | 0.702 | 0.548–0.860 |

**Range of means 0.690–0.724**, all CIs overlapping — slope is prior-robust.

## Supplement tables regenerated on the chordal run

- **Table S3 (regional in-sample RMSE, 14 models × 5 regions + Overall)** — `5b_overfitting_diagnostics.R`
  → `manuscript/tables/Table_S2_regional_performance.csv` (Overall n 1,128; Africa 142, Americas 373,
  Asia 551, Europe 30, Oceania 32).
- **Table S8 (global parameters)** — `regen_tables_v10.R` → `table2_global_params_body.tex`
  (baseline_sp slope 0.649, GP scale 3,589 km).
- **Tables S2 / S9** — prior descriptor slope range 0.60–0.70; detection-threshold grid above.

## OPEN

- `manuscript/METHODS.md` is git-untracked; its chordal §4.2 edit (d = chordal km; PD by
  construction) is present and self-verified but not part of any committed diff.

## Provenance notes

- Deposit: the 14 posteriors in `bradleylab/leafwax-data` are rebuilt chordal on branch
  `chordal-deposit` (v3.0.0, spatial_metric stamped, elevation coefficients retained); the frozen
  great-circle deposit (v2.0.0 = zenodo.21286445) is preserved at leafwax-data commit 5a477e3.
- The frozen great-circle run is preserved unchanged at `results/c2_run_20260626/` for comparison.
