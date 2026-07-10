# NUMBERS.md — authoritative frozen-run quantities for the leaf-wax δ²H manuscript

**Single source of truth for every quotable number.** All values below are
regenerated from the frozen-data refit; the rewrite (main text + supplement)
must draw from this file, not from older drafts or the GCA numbers.

## Provenance

- **Run:** `results/c2_run_20260626/model_output/` (14 models fit 2026-06-23 – 06-29;
  new analyses confounding_v2 + simrecovery 2026-07-02).
- **Calibration input:** `data/frozen/leafwax_d2h_c29_calibration_v1.csv`, md5
  `bb52649130a02b9b7d897325899b4f5a` (verified against the C2 authoritative copy).
- **Regenerated:** 2026-07-03 via `scripts/regen_tables_v10.R` +
  `scripts/regen_manuscript_numbers.R` with `LEAFWAX_RUN_DIR` pointed at the frozen run
  (→ `model_analysis/tables/*.csv`, `manuscript/drafts/REGENERATED_NUMBERS_v10.md`).
- Sediment for classical stats: `results/3_sediment_ready_for_modeling.rds` (frozen, 1128 rows).

## How to regenerate (replication)

All commands run from the analysis-repo root with the frozen run selected via
`LEAFWAX_RUN_DIR` (models) / `LEAFWAX_CONFOUNDING_DIR` (confounding sim). `MO` below
= `results/c2_run_20260626/model_output`, `RUN` = `results/c2_run_20260626`.

| Artifact | Command | Output |
|---|---|---|
| Core tables + main-text scalars (slopes, variance, RMSE, CIs, C₃/C₄, veg coeffs) | `LEAFWAX_RUN_DIR=MO Rscript scripts/regen_tables_v10.R` + `scripts/regen_manuscript_numbers.R` | `model_analysis/tables/*.csv` |
| Precip-scale reconstruction SD + detection thresholds | `LEAFWAX_RUN_DIR=MO Rscript scripts/regen_precip_space.R` | `manuscript/drafts/REGENERATED_NUMBERS_precip_space_v10.md` |
| Table S3 (regional in-sample RMSE, 14×6) | `LEAFWAX_RUN_DIR=MO Rscript 5b_overfitting_diagnostics.R` | `manuscript/tables/Table_S2_regional_performance.csv` |
| Table S4 + Fig 3 (confounding stress test) | `LEAFWAX_CONFOUNDING_DIR=RUN Rscript manuscript/drafts/comms_ee/analysis/make_figures.R NG3` | prints "Figure NG3 source values"; `figures/Figure_3.pdf` |
| Table S4 caption correlation (real-data ρ) | `5d_spatial_confounding_check.R` logic: `cor(colMeans(alpha_spatial), oipc_values[,1])` on `baseline_env_sp` | 0.4555 |
| Table S5 (prior sensitivity, 7 variants) | C2 SLURM 1930423 refits of `baseline_veg_sp` | medians 0.614–0.639 |
| Table S7 + within/between decomposition | `Rscript manuscript/drafts/comms_ee/analysis/within_between_decomposition.R` | `manuscript/drafts/comms_ee/tables/{regional_slopes,within_between_decomposition}.csv` |
| Overfitting PI widths by density stratum (S2.5.3) | `LEAFWAX_RUN_DIR=MO Rscript manuscript/drafts/comms_ee/analysis/regen_predictive_interval_widths.R` | `manuscript/drafts/comms_ee/tables/predictive_interval_widths.csv`; baseline_env_sp 64.3 dense→65.0 sparse, all 9 spatial models keep ordering (NN-distance terciles). Replaces earlier ad-hoc 63.8/66.1. |
| Fig S4 panels (integration OLS) | baseline `stan_data$oipc_values` (scales 1,3,5,10,20,40,70,100,150 km): A=`oipc_mean`, B=col@10 km, C=equal-weight, D=Bayesian effective scale | slopes 0.833 / 0.797 / 0.815 / 0.779 |
| Fig 5 (vegetation main effects + interactions) | `LEAFWAX_RUN_DIR=MO Rscript manuscript/drafts/comms_ee/analysis/make_figures.R VEG` | `figures/Figure_5.pdf` |

## ⚠️ STALE numbers to fix in both manuscripts

| quantity | current draft | frozen truth | source |
|---|---|---|---|
| **sample size** | 1,129 | **1,128** | frozen sediment nrow; sed_min.csv |
| **Africa count** | 143 | **142** | continent re-sum (dropped row was African) |
| `baseline_sp` slope CI upper | 0.671 | **0.675** | table2 |
| min bulk ESS | 1,098 | **1,097** (elevation_c4_sp) | REGEN §per-model diag |
| `full_sp` interpolated range | −47/+52 → already fixed to −48/+53 | −48.4/+52.8 | interpolated_field_ranges.csv |
| **`full_sp` tree main effect** | "grass and tree negative" | **tree now −0.176 [−0.386,+0.036], crosses 0 (unresolved); grass stays −0.252 resolved** | vegetation_coefficients.csv |

Continent breakdown (frozen): Asia 551, Americas 373, Africa **142**, Oceania 32, Europe 30 = **1128**.

## ✅ HOLDS (verified unchanged on frozen data)

- Non-spatial slope **0.78** (baseline 0.779, 0.749–0.809).
- Spatial-adjusted slope range **0.53–0.62** (full_sp 0.529 → baseline_veg_sp 0.623).
- Spatial variance **48–57%** (47.8–56.8); intercept share **96–98%** (96.5–97.7).
- RMSE **21→16** (baseline 21.2 → spatial 15.5–16.0); ~25% reduction.
- GP length scales: frozen range **3,611–3,931 km** (quote as "~3,610–3,930 km"; the draft's "3,600–3,950" pads the upper bound). *(→ supplement one-liner per reframe.)*
- Moran's I (OLS resid, k=8) **0.583**; raw δ²H_wax 0.845.
- Mean nearest-neighbor **33.8 km** (median 3.9; 62.5% <10 km; 92.2% <100 km).
- OLS Fig 1: β₀ **−121.9 ± 1.4**, β_OIPC **0.833 ± 0.019**, R² **0.643**, RMSE **23.0‰**.
- Intercept–OIPC correlation **r = 0.39–0.51** (r² 0.15–0.26).
- Spatial intercept SD **14.2–17.2‰**; obs-level range −41 to +48‰.
- Spatial slope field: mean **0.52–0.62**, SD **~0.12–0.15**, range 0.23–0.97.
- Residual σ (spatial) **≈15.6‰**; 95% PI width 64.6 (spatial) vs 80.2‰ (non-spatial).

## Per-model slopes (table2_global_params.csv), 95% CI

| model | β_OIPC | 95% CI |
|---|---:|---|
| baseline (non-sp) | 0.779 | 0.749–0.809 |
| baseline_sp | 0.575 | 0.473–0.675 |
| baseline_env_sp | 0.544 | 0.423–0.662 |
| baseline_veg_sp | 0.623 | 0.491–0.750 |
| full_sp | 0.529 | 0.409–0.644 |
| full_interact_sp | 0.589 | 0.442–0.733 |
| elevation_only_sp | 0.585 | 0.465–0.694 |
| elevation_c4_sp | 0.584 | 0.465–0.697 |
| c4_only_sp | 0.573 | 0.472–0.670 |
| elevation_c4_interact_sp | 0.561 | 0.434–0.681 |

## C₃/C₄ + vegetation (reframe-critical)

**Main effects, full_sp (95% CI, vegetation_coefficients.csv):**
- grass **−0.252 [−0.419, −0.085]** resolved negative
- tree **−0.176 [−0.386, +0.036]** NOT resolved (weakened vs GCA)
- shrub 0.060 [−0.183, 0.300] not resolved
- C4 0.044 [−0.045, 0.133] not resolved

**δ²H_precip × PFT interactions, full_interact_sp (table5):** — HEADLINE, all hold
- ×tree **−0.241 [−0.414, −0.063]** resolved (shallower slope in trees)
- ×shrub **+0.268 [+0.110, +0.431]** resolved (steeper in shrubs)
- ×C4 +0.021 [−0.075, 0.120] not resolved
- ×grass −0.027 [−0.140, 0.088] not resolved
- ⇒ PFT modulates the slope, cutting across woody/herbaceous, **not along C₃/C₄ lines**.

## Convergence / divergences (frozen, adapt_delta 0.999 for full_sp/full_interact_sp)

- All models: max R̂ ≤ **1.0082** (elevation_c4_sp); min bulk ESS **1,097** (elevation_c4_sp); max-treedepth hits 0; min E-BFMI 0.743.
- Divergences: **full_sp 2**, **full_interact_sp 4**, all other 12 models **0**.
- (Divergence immateriality already established: `results/divergence_sensitivity.csv`.)

## Uncertainty propagation / detection thresholds

- σ_residual (wax) median **15.6‰**; combined σ_total = √(15.64²+3²) = **15.93‰**.
- Detection threshold (**δ²H_wax scale**), two-sample: ρ=0 **44.1**, ρ=0.5 **31.2**, ρ=0.8 **19.7**, ρ=0.9 **14.0‰**.
- **δ²H_precip scale** — CONFIRMED by per-draw propagation (`scripts/regen_precip_space.R`,
  `REGENERATED_NUMBERS_precip_space_v10.md`), reconstruction model **baseline_env_sp**:
  - single-sample reconstruction SD **29.3‰ [24.1, 37.6]** → quote "≈29‰ (24–38‰)" ✓ matches draft.
  - two-sample threshold: ρ=0 **81.1‰**, ρ=0.5 **58.4**, ρ=0.8 **38.8**, ρ=0.9 **29.5‰** → "≈81‰" ✓.
  - (full_sp gives SD 30.0, ρ=0 threshold 83.2; baseline_sp 28.2 / 78.0 — model-dependent; the draft's 29/81 = baseline_env_sp.)
- Fig 5 endpoints (wax scale): spatial ρ=0.9 → **14.0‰**; non-spatial ρ=0 → **58.7‰** (σ_baseline 20.96‰).

## Confounding simulation (Fig 3, confounding_v2, 2026-07-02) — verified

| scenario | achieved ρ | true.std | GP post.mean | OLS (oipc[,1]) | % bias absorbed |
|---|---:|---:|---:|---:|---:|
| rho00 | 0.000 | 0.474 | 0.511 | 0.734 | 85.8 |
| rho03 | 0.639 | 0.414 | 0.613 | 0.798 | 48.1 |
| empirical | 0.726 | 0.394 | 0.664 | 0.826 | 37.5 |
| rho05 | 0.752 | 0.389 | 0.682 | 0.835 | 34.2 |

**Note on OLS definition:** the table's OLS column uses `oipc_values[,1]` (1 km point).
The manuscript (Fig 3 / `make_figures.R` NG3 and supplement S2.6.2–3) uses the multi-scale
**rowMeans(oipc)** predictor, giving OLS = **0.735 / 0.800 / 0.828 / 0.837** (rho00/rho03/
empirical/rho05) and GP posterior **median** 0.510 / 0.612 / 0.662 / 0.681 (vs the means above).
Table S4 and the S2.6.2 prose use these rowMeans/median values (empirical OLS 0.828, GP 0.662).

Empirical headline (main text): true **0.39**, OLS **0.83**, GP **0.66**, ~40% removed / ~60% remains.

**Supplement Table S4 form** (nominal ρ_c on the x-axis of Fig 3; GP posterior median + 95% CI;
extracted from `results/c2_run_20260626/confounding_v2_3c_*` via `make_figures.R` NG3, N=1128):

| ρ_c (nominal) | input β (true.std) | GP median | 95% CI | bias | covers input |
|---:|---:|---:|---|---:|:--:|
| 0.0 | 0.474 | 0.510 | [0.424, 0.607] | +0.036 | yes |
| 0.3 | 0.414 | 0.612 | [0.538, 0.697] | +0.198 | no |
| 0.45 (empirical) | 0.394 | 0.662 | [0.593, 0.743] | +0.268 | no |
| 0.5 | 0.389 | 0.681 | [0.613, 0.759] | +0.292 | no |

## Simulated recovery (simrecovery, 2026-07-02) — verified

- 3a uniform 0.70 → **0.68 [0.61, 0.75]** (CI covers truth).
- 3b spatially varying mean 0.65 → **0.67 [0.60, 0.75]** (CI covers truth).

## Prior / hyperparameter sensitivity (frozen re-run, SLURM 1930423, 2026-07-03)

7 variants of `baseline_veg_sp` refit on frozen data; β_OIPC posterior median [95% CI]:

| variant | median | 95% CI |
|---|---:|---|
| beta_oipc_prior_wider | 0.615 | 0.476–0.750 |
| beta_oipc_prior_shifted | 0.614 | 0.476–0.748 |
| beta_oipc_prior_uninformative | 0.615 | 0.482–0.744 |
| pc_slope_relaxed | 0.629 | 0.494–0.762 |
| pc_slope_very_relaxed | 0.639 | 0.500–0.779 |
| ls_longer | 0.627 | 0.497–0.755 |
| ls_shorter | 0.619 | 0.491–0.746 |

**Range of medians 0.614–0.639**, all CIs overlapping — slope is prior-robust
(GCA quoted 0.614–0.638; holds on frozen data).

## Supplement tables regenerated on frozen run (2026-07-03)

- **Table S3 (regional in-sample RMSE, 14 models × 5 regions + Overall)** — regenerated by
  `5b_overfitting_diagnostics.R` with `LEAFWAX_RUN_DIR=results/c2_run_20260626/model_output`
  → `manuscript/tables/Table_S2_regional_performance.csv` (Overall n=1128, Africa n=142,
  Americas 373, Asia 551, Europe 30, Oceania 32). This CSV is the authority for the Table S3 body.
- **Table S7 (region-specific OLS slopes)** — regenerated by `within_between_decomposition.R`
  on the frozen 1128-row sediment → `manuscript/drafts/comms_ee/tables/regional_slopes.csv`.
  Africa n 142, slope −0.058 (quote −0.06); other regions unchanged at stated precision.
  Within/between slopes: between (n-weighted) 0.810, within 0.842, global pooled 0.833.
- **Fig S4 (OLS under 4 spatial-integration choices)** — panels recomputed on frozen data
  (predictors from baseline `stan_data$oipc_values`; scales 1,3,5,10,20,40,70,100,150 km):
  A (oipc_mean, =Fig 1) **0.833** / R² 0.643 / RMSE 23.0; B (10 km column) 0.797 / 0.686 / 21.6
  (unchanged); C (equal-weight) 0.815 / 0.693 / 21.4 (unchanged); D (Bayesian effective scale)
  slope **0.779** (baseline posterior median) / R² 0.697 / RMSE 21.2. Only A and D slopes shifted
  (0.832→0.833, 0.778→0.779) vs the pre-freeze figure; the scatter is visually unchanged.

## OPEN — none. Phase 1 + supplement-table numbers are locked.

## RESOLVED

- ✅ Full manuscript number audit (main + supplement vs frozen run, codex + independent, 2026-07-03).
  Fixed: main R2 slopes (0.779/0.575/0.53–0.62), Fig 1c CI span (−0.38 to 1.84), Fig 2 span (0.53–0.84),
  supplement simrecovery 3a/3b (0.68/0.67), confounding biases (+0.20/+0.27/+0.29), 92.2%, R² 0.697.
- ✅ Supplement diagnostics verified on frozen run (subagent, 2026-07-03): OLS inverse example, banded
  Moran (100 km band 0.57→0.56), elastic-net OIPC coef 0.66→0.65, VIF max-temp 19.8→19.4 / soil 5.6→5.5,
  density Medium count 431→430 (225+430+473=1128), PI widths 63.9/66.2, PD min eigenvalue 8.5e-5→2.0e-4,
  divergence-sensitivity ranges (exact match). Elastic-net ordering reworded (max temp is 2nd, excluded for VIF).
- ✅ Supplement Tables S3, S5, S7 + Fig S4 regenerated on frozen run and synced to supplement.tex (2026-07-03).
- ✅ Supplement Table S4 (confounding stress test) synced to frozen confounding_v2 (2026-07-03).
- ✅ Table S4 caption note verified: real-data `baseline_env_sp` spatial-intercept ↔ δ²H_precip
  correlation = **0.4555** on frozen 1128 data (alpha_spatial posterior means vs `oipc_values[,1]`,
  the `5d_spatial_confounding_check.R` definition) — unchanged from the 1129 value, so "≈0.45" holds.
- ✅ Prior-sensitivity re-run on frozen data (medians 0.614–0.639; table above).
- ✅ Precip-scale reconstruction SD (29‰) + two-sample thresholds (81‰) confirmed via per-draw propagation on frozen `baseline_env_sp` (see detection-thresholds section).
- Prior tables backed up at `results/c2_run_20260626/tables_prev_20260703/` (pre-frozen, for reference).
