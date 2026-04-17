# Table provenance — leafwax GCA manuscript

Phase 4 audit (2026-04-16). Maps each `manuscript/tables/*.tex` to the
pipeline script that produces its numeric content.

**Important finding:** No R script in the pipeline writes LaTeX directly.
Every `.tex` file is typeset by hand from CSV or console output. Regeneration
against the April 2026 run therefore requires (a) re-running the producer
script to refresh the CSV/console numbers, then (b) manually updating the
`.tex` source. Mark `regenerate = Y` on tables whose numeric values depend on
fits; those values are currently stale (reflect the Dec 2025 run).

## Main text tables

| Table (.tex)                       | Producer (CSV / console)                                              | Input fits (April run)                    | Regen? | Notes                                                                                         |
| ---------------------------------- | --------------------------------------------------------------------- | ----------------------------------------- | :----: | --------------------------------------------------------------------------------------------- |
| `table1_model_performance.tex`     | `extract_full_model_analysis.R:363` → `full_model_performance.csv`    | all 14 models                             |   Y    | Columns: Rhat, ESS, LOOIC, SE, ΔLOOIC, p_eff, n_hi_k. Hand-assembled from CSV.                |
| `table2_global_params.tex`         | `analyze_paleo_inversion_models.R:293` → `practical_comparison.csv`   | all sp models (RMSE, R², β, λ, GP scale)  |   Y    | Also pulls `beta_oipc`, `effective_scale_km`, knot SDs from `extract_full_model_analysis.R`.  |
| `table3_variance_decomp.tex`       | `extract_full_model_analysis.R:358` → `variance_df.csv` (spatial/residual/intercept/slope) | all `*_sp` models |   Y    | 5c_spatial_patterns.R also reports the same decomposition.                                    |
| `table4_environmental.tex`         | `extract_vegetation_coefficients.R:111` → `vegetation_coefficients.csv` + `extract_full_model_analysis.R:353` → `coefficients_df.csv` | models with env/veg predictors |   Y    | Combines β_C4, β_tree, β_shrub, β_grass, β_precip with 90% CIs.                                |

## Supplement tables

| Table (.tex)                                    | Producer                                                | Input fits (April run)          | Regen? | Notes                                                                                       |
| ----------------------------------------------- | ------------------------------------------------------- | ------------------------------- | :----: | ------------------------------------------------------------------------------------------- |
| `Table_S1_priors.tex`                           | **MANUAL** — documents prior specifications             | none (derived from `config.yaml` + `4d_leaf_wax_spatial_model.stan`) |   N*   | *Update only if prior specs changed. Dec→Apr: `beta_precip` prior widened (see config.yaml:105-302). |
| `Table_S2_regional_performance.tex`             | `5b_overfitting_diagnostics.R:69-81` → `regional_performance` (console `print`) | `full_sp` + comparators |   Y    | No CSV output; transcribed from console. Consider adding `write.csv` in 5b.                 |
| `table_S3_compilations-UNRELIABLE.tex`          | **MANUAL** (marked UNRELIABLE)                          | n/a                             |   —    | Name flags for deletion or rebuild against `data/HrenBrandon/` cross-check.                 |

## Markdown table drafts

| File                                     | Status                                                          |
| ---------------------------------------- | --------------------------------------------------------------- |
| `d2H_wax_regression_tables.md`           | Draft prose + tables predating April run. Superseded by main .tex. |
| `d2H_wax_regressions_UPDATED.md`         | Later draft still pre-April. Do not carry into submission.      |

## Producing-script `fit.rds` input convention (April run)

All extract/analyze scripts expect fits at `model_output/<model>/fit.rds`
relative to the repo root. April run lives at
`results/c2_run_20260414/<model>/fit.rds`. Before regeneration, either:

1. Symlink `model_output -> results/c2_run_20260414` at the repo root, or
2. Patch the extract scripts to read `results/c2_run_20260414/<model>/fit.rds`.

Option (1) is minimal-change and leaves the scripts untouched.

## Acceptance

- Every `.tex` is mapped to a producer or flagged MANUAL.
- `table_S3_compilations-UNRELIABLE.tex` flagged as candidate for deletion.
- Producer scripts identified — no `.tex`-writing code exists in the repo;
  `.tex` files will need hand edits after CSV regeneration.
