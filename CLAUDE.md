# CLAUDE.md — repo-specific standards for leafwax-spatial

This file encodes project-specific rules for any agent working on this repository. It applies on top of, not in place of, the operator's global CLAUDE.md. When a conflict arises, the stricter rule wins.

This is a Bayesian spatial calibration of n-C29 alkane δ²H against precipitation δ²H (OIPC) across 1131 sites, fit in Stan via `cmdstanr` on WashU Compute2. The output is a peer-reviewed GCA manuscript. Standards below are what make the manuscript defensible.

## Scientific integrity — non-negotiable

1. **No imputed, guessed, interpolated, or fabricated numbers.** Every numeric value that appears in a figure, table, or prose claim must trace to (a) a draw from the chain CSVs of the April 2026 Stan run (under `results/c2_run_20260414/`), (b) a deterministic transform of those draws, or (c) a documented input in `input_data/` or `data/`. If a number cannot be sourced this way, it cannot appear in the manuscript. State explicitly when a value is unavailable.
2. **Same-draws invariant.** The April 2026 posterior is the ground truth. Analysis code may read more or fewer columns from the chain CSVs, summarize differently, or re-run `loo()` — but it may not resample, bootstrap, perturb, or otherwise replace the existing draws with something else. Rerunning HMC would produce different draws (seed-dependent) and is not equivalent to re-analyzing the same posterior.
3. **Stan parameterization is authoritative.** When analysis code reads a parameter name (`alpha_spatial`, `sigma_intercept_spatial`, `beta_oipc_x_c4`, etc.), that name must appear in `4d_leaf_wax_spatial_model.stan`'s parameters or generated quantities block. If a name in analysis code does not match Stan, it is a bug, not a convention choice. Grep-verify before committing:
   ```bash
   grep -nE "^\s*(real|vector|matrix)\s+\w" 4d_leaf_wax_spatial_model.stan  # Stan exports
   grep -rn "<name>" *.R manuscript/                                        # usage
   ```
4. **Unit / scale must be reasoned about, not assumed.** Parameters in the Stan model live on a standardized scale. Figures and tables presenting original-scale values (δ²H in ‰) require an explicit back-transform using the scaling parameters saved with `stan_data`. When code multiplies, subtracts, or scales, the comment above it must state the intended units on both sides.
5. **All numeric outputs verifiable.** Every emitted `.tex` table cell and every summary number in prose must be reproducible by running the pipeline. Hand-typed numbers are prohibited except where explicitly flagged (e.g., hand-maintained priors table wrappers).

## Reproducibility contract

1. **One command rebuilds.** `make pipeline` from a clean clone, after `results/c2_run_20260414/` is populated, regenerates every main figure, every main table, and Table_S2 without manual intervention.
2. **Fit-step self-contained.** Running `4c_fit_models.R` on a fresh config + input set produces every rds file downstream analysis needs. There is no followup script, no separate `loo` step, no sync step, no hidden state. If someone reruns fitting from scratch, they get everything the first time.
3. **Analysis never reads `fit.rds` or chain CSVs.** Every figure, table, and `5*_*.R` analysis script reads only `posterior_draws.rds`, `diagnostics.rds`, `loo.rds`, `stan_data_<model>.rds`, `config_<model>.rds`, and `3_sediment_ready_for_modeling.rds`. This is what makes the analysis layer portable to any machine with R.
4. **No networked resources in user-facing flow.** No `aws s3 sync`, no `scp`, no HTTP fetches in figure/table scripts or their prereqs. Operators may move rds files onto and off of HPC as they see fit; that is not part of the pipeline spec.
5. **Git-tracked state is the manuscript.** Results (`results/`) are git-ignored. Code, methods, and per-artifact provenance (`manuscript/METHODS.md`, `manuscript/FIGURES.md`, `manuscript/TABLES.md`, `manuscript/PHASE5_PLAN.md`) are tracked and authoritative.

## Code-change standards for analysis scripts

1. **Verify against Stan before editing.** Before adding or renaming a parameter reference in an analysis script, grep the Stan model and confirm the name exists and is on the expected scale.
2. **No ad-hoc scale fixes.** If a number looks wrong, trace it back through the Stan generated quantities or the standardization in `4a_spatial_functions.R` / `stan_data$scaling_params`. Do not divide by `sigma` or multiply by a random SD to "make it match the figure."
3. **Flag, don't silently fix.** If a figure or table script has a pre-existing numerical error, document it in the commit message and in `manuscript/FIGURES.md` / `TABLES.md`. Do not bury a fix in an unrelated commit.
4. **Single source of truth for variable subsets.** The list of draws extracted into `posterior_draws.rds` lives in `4c_fit_models.R`. Analysis scripts should assume any variable in that list is available per-draw; indexed arrays outside that list are only available as summaries from `diagnostics.rds`.
5. **No new `fit$*()` calls.** Any new script touching posteriors must use the `posterior_helpers.R` API (see `manuscript/figure_code/posterior_helpers.R`). Calls to `fit$draws()` or `fit$summary()` are prohibited in analysis code.

## Peer-review defensibility checks

Before committing a figure/table change, the author (or agent) must be able to answer:
- Which parameter(s) are plotted / tabulated?
- Which `.stan` line defines each parameter?
- On what scale (standardized / original / log / probability)?
- Which transform (if any) is applied between the draws and the final number?
- Where does the number appear in the CSV/diagnostics lineage?

If any answer is "I'm not sure," the change is not ready.

## Canonical documents

| Document | Role |
|---|---|
| `manuscript/METHODS.md` | Authoritative methods doc, with line-level citations into the pipeline scripts. |
| `manuscript/FIGURES.md` | Per-figure provenance + regeneration checklist. |
| `manuscript/TABLES.md` | Per-table provenance + regeneration checklist. |
| `manuscript/PHASE5_PLAN.md` | Current plan for the figure/table regeneration refactor. |
| `4d_leaf_wax_spatial_model.stan` | Source of truth for all parameter names and scales. |

## Scope for agents

- Run `ruff`, type-checks, or project lints where they apply. For R: run `Rscript -e 'devtools::test()'` if tests exist; otherwise smoke-test by sourcing the affected script through the point of change.
- Follow the operator's global edit-safety and refactoring-discipline rules.
- When in doubt about a numerical operation, stop and ask. Fabrication risk is higher here than in most projects because the analysis chain is long and the pipeline already contains pre-existing bugs (see `manuscript/PHASE5_PLAN.md` open items).
