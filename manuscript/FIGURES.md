# Figure provenance — leafwax GCA manuscript

Phase 4 audit (2026-04-16). Maps each figure script in
`manuscript/figure_code/` and `manuscript/figures/supplement_figs/` to its
input `fit.rds` files and rendered PDF(s). All hardcoded paths updated to
point at the April 2026 run (`results/c2_run_20260414/`).

## Main text figures

| Figure (.R)                                   | Input `fit.rds` (April run)                       | Other inputs                                                                                  | Output PDF (manuscript/figures/main_figs/)       | Regen? | Notes                                                                                           |
| --------------------------------------------- | ------------------------------------------------- | --------------------------------------------------------------------------------------------- | ------------------------------------------------ | :----: | ----------------------------------------------------------------------------------------------- |
| `Figure_01_ols_regression.R`                  | `baseline/posterior_draws.rds`                     | `3_sediment_ready_for_modeling.rds`, `_prepared_data/stan_data_baseline.rds`, `common_functions.R` | Figure 1 (OLS regression panels)                 |   Y    | Sources `common_functions.R` (in `figure_code/`). Emitted under `main_figs/` via `make figures`. |
| `Figure_02_all_environmental_variables_noborders.R` | none (analytical)                          | `3_sediment_ready_for_modeling.rds`, `input_data/GlobalPrecip/d2h_MA.tif`, other rasters        | 8-panel environmental Figure 2                    |   Y    | CSV-based inputs replaced with rds read.                                                        |
| `figure_03_spatial_confounding_revised.R`     | `baseline/posterior_draws.rds`, `baseline_sp/posterior_draws.rds` | `3_sediment_ready_for_modeling.rds`                                        | Figure 3 (spatial confounding)                   |   Y    | Fixed units bug (Phase 5 W3 corrections log below).                                              |
| `Figure_04_spatial_maps.R`                    | `full_sp/posterior_draws.rds`                      | `_prepared_data/stan_data_full_sp.rds`                                                         | Figure 4 (spatial intercepts/slopes)              |   Y    | Script now emits `Figure_04.{pdf,png}` (W7 cleanup fixed the legacy `Figure_03` name).          |
| `Figure_05_detection_thresholds.R`            | `full_sp/posterior_draws.rds`, `baseline/posterior_draws.rds` | `_prepared_data/stan_data_full_sp.rds`, `_prepared_data/stan_data_baseline.rds` | Figure 5 (detection thresholds) |   Y    | σ no longer hardcoded — reads Stan `sigma_residual_original` (fallback to `sigma × d2H_sd`).    |

## Supplement figures

| Figure (.pdf)                           | Producing script (if known)                         | Regen?  | Notes                                                                                |
| --------------------------------------- | --------------------------------------------------- | :-----: | ------------------------------------------------------------------------------------ |
| `Figure_S1_03_spatial_weighting.pdf`    | unknown — likely `3_prep_data.R` weighting diagnostics | ?     | No `.R` in `supplement_figs/` other than `merge_figs.R`; trace to archive if needed. |
| `Figure_S2_pairwise_correlations.pdf`   | likely `3d_analyze_collinearity.R`                   | Y      | Pairwise VIF/PCA panel.                                                              |
| `Figure_S3_residuals.pdf`               | likely `5a_model_validation.R` or `5b_overfitting_diagnostics.R` | Y | Residual diagnostics.                                                      |
| `Figure_S4_ols_spatial_scale.pdf`       | `archive/Ω_old_drafts/.../Figure_S4_ols_spatial_scale.R` | Y  | Producer lives only in archive; copy into `figure_code/` before regen.               |
| `Figure_S5_elevation.pdf`               | unknown — trace                                     | ?       |                                                                                     |
| `supplementary_figures.pdf` (merged)    | `supplement_figs/merge_figs.R`                      | Y       | Wraps the individual S1–S5 PDFs into one file.                                        |

## Per-artifact regeneration checklist

All main figures read from `results/c2_run_20260414/` via the
`posterior_helpers.R` API — no script reads `fit.rds` or chain CSVs.
Regeneration is a single command: `make figures` from the repo root
(W7 Makefile). Outputs land in `manuscript/figures/main_figs/`; scratch
copies under `manuscript/figure_code/` are gitignored.

Supplement figure producers still require tracing through the archive.
After tracing, lift supplement `.R` files into `figure_code/` (or a
dedicated subdir) and add Makefile rules.

## Corrections log

**2026-04-16 — Figure 3 spatial intercept back-transform corrected (Phase 5 W3).**
Prior versions of `figure_03_spatial_confounding_revised.R` multiplied
`alpha_spatial[i]` by `sigma_intercept_spatial` to obtain per-mille spatial
effects. This is a units error: `alpha_spatial` is defined on the
*standardized* response scale (4d_leaf_wax_spatial_model.stan:316,
331-335) and already includes the global intercept `beta_0`, while
`sigma_intercept_spatial` is the GP intercept SD in original per-mille units
(Stan:488). Multiplying the two mixes scales.

The corrected transform isolates the GP residual contribution in ‰:
```
alpha_spatial_contribution_‰ = (alpha_spatial − beta_0) × d2H_sd
```
where `d2H_sd = stan_data$scaling_params$d2H_sd` is the standardization factor
used in prep (4a_spatial_functions.R:614). This matches Stan's own variance
accounting at `4d_leaf_wax_spatial_model.stan:504` (`variance(alpha_spatial − beta_0)`)
and Stan's original-scale derivation at line 488.

The same corrupted transform existed in `5e_spatial_intercept_correlations.R`
and `5e_spatial_intercept_correlations_weighted.R`; both fixed in the same
commit. Downstream regenerated figures should be compared against the prior
"±60 ‰" narrative — the magnitude range may shift.

## Known blockers

1. ~~**cmdstanr chain-CSV paths.**~~ **RESOLVED (Phase 5 W1–W5, 2026-04-16).**
   `posterior_draws.rds` widened at fit time to include `mu`, `d2H_rep`,
   `log_lik`, `alpha_spatial`, `beta_oipc_spatial`, `z_intercept_spatial`,
   `z_slope_spatial`, `scale_weights` (plus existing scalars); `loo.rds`
   emitted in the same step. Every figure and analysis script migrated to
   the `posterior_helpers.R` API. No script reads `fit.rds` anymore.
2. ~~`Figure_05_detection_thresholds.R` hardcodes σ from the Dec-2025 run.~~
   **RESOLVED (Phase 5 W5).** Reads `sigma_residual_original` per draw
   from `full_sp/posterior_draws.rds` and `baseline/posterior_draws.rds`.
3. Supplement S1, S5 producers unidentified — need to trace or mark MANUAL.
4. `common_functions.R` did not exist in `figure_code/`; copied + re-rooted
   for the manuscript/ layout this phase.

## Smoke test (Figure_01, 2026-04-16)

`head -n 25 Figure_01_ols_regression.R | Rscript -` from `manuscript/figure_code/`
completes without error. Confirms:
- `common_functions.R` sources cleanly (with lazy `gt`/`fields` load).
- `load_project_config()` finds `0_load_config.R` via `REPO_ROOT`.
- All three `readRDS()` calls resolve (sediment, stan_data_baseline, baseline/fit.rds).

First error appears at line 26 (`fit_baseline$summary(...)`) due to blocker (1).

## Acceptance

- All main figure scripts have their input paths updated to
  `results/c2_run_20260414/` (or repo-root-relative paths that resolve there).
- Supplement producers partially traced; open items flagged.
- Figure 5 σ-refresh + S1/S5 producer tracing carried forward as new tasks.
