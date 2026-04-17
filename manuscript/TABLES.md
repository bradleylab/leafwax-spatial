# Table provenance — leafwax GCA manuscript

Phase 5 W6 (2026-04-17). Every main-text and Table S2 cell is produced
by an R script that reads only `results/c2_run_20260414/<model>/*.rds`
and writes both a `.csv` source-of-truth and a `.tex` (or `.tex` body
fragment for tables with landscape longtable wrappers). Regeneration
is one command: `make tables` from the repo root (W7 Makefile).

## Main text tables

| Table (.tex)                       | Producer → CSV                                                        | Input rds (April run)                     | Regen? | Notes                                                                                          |
| ---------------------------------- | --------------------------------------------------------------------- | ----------------------------------------- | :----: | ---------------------------------------------------------------------------------------------- |
| `table1_model_performance.tex`     | `extract_full_model_analysis.R` §8 → `table1_model_performance.csv`   | 14 × `loo.rds` + `diagnostics.rds`        |   Y    | ΔLOOIC via `loo::loo_compare`; Rhat / ESS from diagnostics; n_hi_k = count of Pareto-k > 0.7. |
| `table2_global_params.tex`         | wrapper `\input{table2_global_params_body.tex}`                       | 14 × `posterior_draws.rds`                |   Y    | Body produced by `extract_full_model_analysis.R` §11. Per-draw stats from `mu`; β₀ on ‰ scale.|
| `table3_variance_decomp.tex`       | `extract_full_model_analysis.R` §9 → `table3_variance_decomp.csv`     | 9 spatial `diagnostics.rds` (GQs)         |   Y    | Uses Stan `prop_variance_*` + `var_spatial_*` generated quantities.                            |
| `table4_environmental.tex`         | `extract_full_model_analysis.R` §10 → `table4_environmental.csv`      | 14 × `posterior_draws.rds` (β_* scalars)  |   Y    | 90% CIs from per-draw quantiles; missing params render as `-`.                                 |

## Supplement tables

| Table (.tex)                                    | Producer → CSV                                                         | Input rds (April run)           | Regen? | Notes                                                                                       |
| ----------------------------------------------- | ---------------------------------------------------------------------- | ------------------------------- | :----: | ------------------------------------------------------------------------------------------- |
| `Table_S1_priors.tex`                           | **MANUAL** — documents prior specifications                            | none (from `config.yaml` + Stan) |   N*   | *Update only if prior specs change. Apr: `beta_precip` prior widened vs Dec config.        |
| `Table_S2_regional_performance.tex`             | wrapper `\input{Table_S2_regional_performance_body.tex}`               | 14 × `posterior_draws.rds`       |   Y    | Body produced by `5b_overfitting_diagnostics.R` tail. Companion CSV lives in `manuscript/tables/`. |
| `table_S3_compilations-UNRELIABLE.tex`          | **MANUAL** (marked UNRELIABLE)                                         | n/a                             |   —    | Name flags for deletion or rebuild against `data/HrenBrandon/`.                             |

## Markdown table drafts

| File                                     | Status                                                          |
| ---------------------------------------- | --------------------------------------------------------------- |
| `d2H_wax_regression_tables.md`           | Pre-April draft. Superseded by main `.tex`.                     |
| `d2H_wax_regressions_UPDATED.md`         | Later pre-April draft. Do not carry into submission.            |

## Numeric audit

`manuscript/table_code/check_tables_numeric.R` re-opens each `.tex` and
verifies that every numeric token aligns (within tolerance) with its
companion CSV. Runs at the tail of `make pipeline` (target `audit`).
Fails non-zero on any drift — catches hand-edits and producer/CSV
desynchronization.

## RMSE / R² estimator choice for Table 2 vs narrative numbers

The manuscript narrative (Section 3.x) may cite RMSE / R² values from
the first console section of `extract_full_model_analysis.R`:

- `full_sp` RMSE = **15.7 ‰**, R² = **0.834**
- `baseline_sp` RMSE = 16.3 ‰, R² = 0.817

Table 2 reports the same models' per-draw posterior summaries:

- `full_sp` RMSE = **16.1 ‰** [15.9, 16.4], R² = **0.824** [0.818, 0.829]
- `baseline_sp` RMSE = 16.5 ‰ [16.2, 16.7], R² = 0.817 [0.811, 0.822]

These are two valid but different estimators:

- **Point estimate** (Section 1 / narrative): compute the single
  posterior-mean prediction `colMeans(mu)` once, then evaluate RMSE
  or R² against the observations. Produces one scalar per model.
- **Per-draw summary** (Table 2 / its CI): evaluate RMSE and R²
  separately for each draw's `mu[s, :]`, then summarize with mean and
  95% quantiles. Produces a posterior over the metric itself.

The per-draw *mean* is always greater than (RMSE) or less than (R²)
the point estimate by Jensen's inequality — averaging predictions
before the non-linear metric reduces apparent error. The 0.4 ‰ / 0.01
gap in `full_sp` is typical, not a bug.

**Consistency rule for the manuscript**: if narrative prose cites the
point-estimate number (15.7 / 0.834), that citation must use the same
estimator every time it appears. Table 2's cells stay on the per-draw
summary because a CI requires a posterior over the metric. Caption
now states the convention explicitly. Do not replace Table 2's mean
column with the point estimate to "match narrative" — that would
disconnect the mean from the CI.

If a reviewer flags the discrepancy, switch the narrative to cite
Table 2's values (or provide both with the estimator named).

## Acceptance (W6 final)

- Every `.tex` emits deterministically from rds + code; `make tables`
  followed by `check_tables_numeric.R` must pass on a clean checkout
  plus populated `results/c2_run_20260414/`.
- CSV siblings live alongside each producer's `.tex`; never edit the
  `.tex` body by hand.
- `Table_S1_priors.tex` remains hand-maintained; flag any prior edit
  against `config.yaml`.
- `table_S3_compilations-UNRELIABLE.tex` still pending delete / rebuild
  decision (cleanup pass).
