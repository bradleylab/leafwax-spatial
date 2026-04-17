# Phase 5 plan v4 — widen fit outputs; portable, peer-review-defensible analysis

**Date drafted:** 2026-04-16. v1, v2, v3 in git history.

## Principle

The fit step emits self-contained rds files that travel anywhere. Analysis reads only those rds files plus prepared inputs. Running `4c_fit_models.R` from scratch produces everything needed on the first pass. No separate loo step, no sync choreography, no hidden state.

The plan also fixes **pre-existing numerical bugs** uncovered by the Phase 4 audit and codex reviews:
- A double-scaling / wrong-units subtraction in `figure_03_spatial_confounding_revised.R` (and twin bugs in `5e_*_correlations.R` variants).
- Parameter-name drift across extract/analyze scripts that makes them read nonexistent Stan variables.
- Incomplete model coverage in extract scripts that feeds Tables 1 and 3.

These bugs are independent of the posterior (**no Stan refit required**) and must be fixed in Phase 5 — not deferred — because any regenerated figure or table that inherits them is not peer-review defensible.

## Data flow after Phase 5

```
fit step (4c_fit_models.R on HPC, reads CSV chains it just wrote)
  └─> <model>/posterior_draws.rds  (WIDE, ~500–700 MB per model)
  └─> <model>/diagnostics.rds       (summaries for all ~19,800 params)
  └─> <model>/loo.rds               (psis_loo object)
  └─> <model>/fit.rds               (kept for HPC-local reanalysis only)

analysis step (any machine with R, reads rds only)
  └─> reads posterior_draws.rds, diagnostics.rds, loo.rds,
           stan_data_<model>.rds, config_<model>.rds, 3_sediment_ready_for_modeling.rds
  └─> emits figures + .tex tables under manuscript/
```

## Key design decisions (forced by codex review of v3)

1. **Save as `draws_array`, not `draws_df`.** Same content, materially lower RAM, ~30% smaller on disk, no tibble conversion cost. Analysis code already uses `posterior::subset_draws()`, which handles arrays directly.
2. **No separate `4e` backfill script.** `4c_fit_models.R:101-110` already falls through to diagnostics + save when `fit.rds` exists. After W1, rerunning `4c_fit_models.R` on C2 against `model_output_20260414/` regenerates widened `posterior_draws.rds` + `loo.rds` for every existing fit automatically. One less script to maintain.
3. **Acknowledge stan_data / config prereqs.** The portable analysis layer needs `stan_data_<model>.rds` and `config_<model>.rds` from `_prepared_data/` — not just the model outputs. Makefile wires them explicitly. `3_prep_data.R` is the producer; it already runs pre-fit.
4. **All `fit.rds` readers must migrate.** Codex surfaced eight more readers I missed:
    - `extract_vegetation_coefficients.R`, `extract_vegetation_variance.R`
    - `5a_model_validation.R`, `5b_overfitting_diagnostics.R`, `5c_spatial_patterns.R`, `5d_spatial_confounding_check.R`
    - `5e_spatial_intercept_correlations.R`, `5e_spatial_intercept_correlations_weighted.R`
   Plus `launch_pipeline.sh:111,158` which keys off `fit.rds` at the orchestration layer (operator concern, noted but out of repo scope).
5. **Figure 3 fix is a scientific correction, not a cosmetic one.** Treat as its own workstream with Stan-verified unit tracking.
6. **`beta_precip_*` interactions are not open questions.** Grep of `4d_leaf_wax_spatial_model.stan` confirms they do not exist. Delete those extraction requests in W0. Also delete `beta_elevation` / `beta_elevation_c4` if same check confirms absence.

## Workstreams

### W0 — Parameter name drift + model coverage (1 commit)

**0a. Grep-verify what actually exists in Stan** (pre-work, not a commit):
```bash
grep -nE "^\s*(real|vector|matrix|array)" 4d_leaf_wax_spatial_model.stan
```
Confirm presence / absence of `beta_precip`, `beta_precip_*` interactions, `beta_elevation`, `beta_elevation_c4`, `sigma_intercept_spatial`, `sigma_slope_spatial`.

**0b. Fix the drift.** Required renames and deletions, all verified against Stan:

| Script | Action |
|---|---|
| `extract_full_model_analysis.R:103` | `beta_precip_amount` → `beta_precip` |
| `extract_full_model_analysis.R:107-109` | `beta_oipc_c4/tree/shrub/grass` → `beta_oipc_x_c4/tree/shrub/grass` |
| `extract_full_model_analysis.R:107-109` | `beta_precip_c4/tree/shrub/grass` → **delete** (not in Stan) |
| `extract_full_model_analysis.R:140` | `beta_elevation`, `beta_elevation_c4` → verify; delete if absent |
| `extract_vegetation_coefficients.R:204` | `sigma_gp_intercept`, `sigma_gp_slope` → `sigma_intercept_spatial`, `sigma_slope_spatial` |
| `analyze_paleo_inversion_models.R:171` | same family — cross-check and rename / delete |
| `5c_spatial_patterns.R:334` | `gp_intercept`, `gp_slope` → verify and rename |

**0c. Expand model coverage.** `extract_full_model_analysis.R:20-22` currently hits 8 of 14 models; `:159` hits 5 of 9 `_sp` models. Extend both loops to cover all applicable models so Tables 1 and 3 are complete.

**Acceptance:**
```bash
grep -rn 'beta_precip_amount\|beta_oipc_c4[^_x]\|beta_oipc_tree\|beta_oipc_shrub\|beta_oipc_grass\|beta_precip_c4\|beta_precip_tree\|beta_precip_shrub\|beta_precip_grass\|sigma_gp_intercept\|sigma_gp_slope' *.R manuscript/
# returns nothing
```
Plus a manual verification that every parameter reference in edited scripts has a matching declaration in `4d_leaf_wax_spatial_model.stan`.

### W1 — Widen `4c_fit_models.R` outputs + emit `loo.rds` (1 commit)

Edit `4c_fit_models.R:205-268`:

```r
# Scalars (same logic as current params_to_check)
draws_to_save <- c("beta_0", "beta_oipc", "beta_precip",  # beta_precip was missing
                   "sigma", "lambda_decay", "effective_scale_km")

if (stan_data$include_c4 == 1) {
  draws_to_save <- c(draws_to_save, "beta_c4")
  if (stan_data$include_pft == 1) {
    draws_to_save <- c(draws_to_save, "beta_oipc_x_c4")
  }
}
if (stan_data$include_pft == 1) {
  draws_to_save <- c(draws_to_save,
                     "beta_tree", "beta_shrub", "beta_grass",
                     "beta_oipc_x_tree", "beta_oipc_x_shrub", "beta_oipc_x_grass")
}
if (stan_data$include_gp == 1) {
  draws_to_save <- c(draws_to_save,
                     "ls_intercept_km", "ls_slope_km",
                     "sigma_intercept_spatial", "sigma_slope_spatial",
                     "alpha_spatial", "z_intercept_spatial", "z_slope_spatial")
}

# Always
draws_to_save <- c(draws_to_save, "mu", "d2H_rep", "log_lik", "scale_weights")

# Replace draws_df with draws_array for RAM/IO efficiency
draws <- fit$draws(variables = draws_to_save, format = "draws_array")
saveRDS(draws, file.path(output_dir, "posterior_draws.rds"))

# Emit loo
log_lik_draws <- subset_draws(draws, variable = "log_lik")
loo_result <- loo::loo(log_lik_draws, cores = getOption("mc.cores", 4))
saveRDS(loo_result, file.path(output_dir, "loo.rds"))
```

Update `manuscript/METHODS.md` to document the widened rds contents and automatic `loo.rds` emission.

**Acceptance:** run on a small test fit (locally or on a C2 sandbox): `readRDS(.../posterior_draws.rds)` is class `draws_array` with ~4k variables for a spatial model, ~1200 for a non-spatial. `readRDS(.../loo.rds)` is class `psis_loo`.

### W2 — Regenerate rds outputs for the April run (operator action, 0 commits)

On C2 `general-interactive` (5 days, 1 MIG) or a short `general-cpu` job, run:
```bash
cd /scratch2/fs1/alexander.s.bradley/leafwax_run
Rscript 4c_fit_models.R
```
4c skips Stan sampling because each model's `fit.rds` exists, falls through to diagnostics + save, and emits widened `posterior_draws.rds` + new `loo.rds` for every model. ~3 hours total CPU (loo is the bulk).

Sync the widened rds files to the canonical April mirror location. This is an **operator action**, not a git commit. It runs once. After this, `results/c2_run_20260414/` contains the data state every post-W1 commit assumes.

**Acceptance:** local `results/c2_run_20260414/` has 14 `loo.rds` files, class `psis_loo`; `posterior::variables(readRDS(".../full_sp/posterior_draws.rds"))` reports > 3500 variables.

### W3 — Figure 3 scientific correction (1 commit)

Independent of the draws-source refactor. This is a bug in the existing code that would produce a wrong figure regardless of which rds file it reads from.

Current state:
- `figure_03_spatial_confounding_revised.R:87-91` multiplies `alpha_spatial` by `sigma_intercept_spatial`.
- `:140` subtracts `alpha_spatial` directly from original-scale `d2H_wax`.

Stan definition:
- `alpha_spatial[i]` is already on the standardized intercept scale (`4d_leaf_wax_spatial_model.stan:316`), constructed from `z_intercept_spatial` × `sigma_intercept_spatial` internally (lines 331-335).
- Response is standardized in prep (`4a_spatial_functions.R:1280`).

Correct transform to original-scale spatial contribution:
```r
# alpha_spatial is already on the standardized intercept scale.
# To express it in per-mille original units:
alpha_spatial_orig <- (alpha_spatial - beta_0) * d2H_wax_sd_original
# d2H_wax_sd_original is in stan_data$scaling_params
```

**Required:**
1. Replace `alpha_spatial * sigma_intercept_spatial` with the correct back-transform in Figure_03.
2. Apply the same fix in `5e_spatial_intercept_correlations.R:86` and `5e_spatial_intercept_correlations_weighted.R:143` (same bug, copy-pasted).
3. Commit with diagnostic output: before/after comparison of the maximum spatial-contribution magnitude to confirm the fix doesn't break the downstream narrative ("±60 ‰ regional baseline fractionation"). If the magnitude changes materially, flag for author review before committing.
4. Update `manuscript/FIGURES.md` with a note documenting the correction.

**Acceptance:** Figure 3 produced from corrected code; spatial-contribution range reported in commit message; matches peer-review expectations (nondimensional analysis confirms the new expression is in ‰).

### W4 — Analysis helper layer (1 commit)

New file `manuscript/figure_code/posterior_helpers.R`:

```r
APRIL_RUN <- normalizePath(file.path("..", "..", "results", "c2_run_20260414"),
                           mustWork = FALSE)

load_draws <- function(model) {
  readRDS(file.path(APRIL_RUN, model, "posterior_draws.rds"))
}

load_summaries <- function(model) {
  readRDS(file.path(APRIL_RUN, model, "diagnostics.rds"))$all_params_summary
}

load_loo <- function(model) {
  path <- file.path(APRIL_RUN, model, "loo.rds")
  if (!file.exists(path)) {
    stop("loo.rds missing for ", model,
         " — rerun 4c_fit_models.R on C2 against the model's output dir.")
  }
  readRDS(path)
}

load_stan_data <- function(model) {
  readRDS(file.path(APRIL_RUN, "_prepared_data",
                    paste0("stan_data_", model, ".rds")))
}

load_config <- function(model) {
  readRDS(file.path(APRIL_RUN, "_prepared_data",
                    paste0("config_", model, ".rds")))
}
```

A matching helper for scripts that currently live at repo root (extract/analyze/5*-series): `scripts/posterior_helpers.R` (identical API, different `APRIL_RUN` relative path). Factor into a common file later if needed; duplicate fine for now.

**Acceptance:** `Rscript -e 'source("manuscript/figure_code/posterior_helpers.R"); x <- load_draws("baseline"); cat(class(x), length(posterior::variables(x)))'` prints a `draws_array` and a variable count.

### W5 — Migrate every `fit.rds` reader (10 commits, grouped)

Every script currently calling `readRDS(.../fit.rds)` or `fit$summary()` / `fit$draws()` must switch to the helper API.

| Script | Commit group |
|---|---|
| `Figure_01_ols_regression.R` | W5.1 (one commit) |
| `figure_03_spatial_confounding_revised.R` | combined with W3 or separate commit after |
| `Figure_04_spatial_maps.R` | W5.3 |
| `Figure_05_detection_thresholds.R` | W5.5 — use original-scale σ (see below) |
| `extract_full_model_analysis.R`, `extract_vegetation_coefficients.R`, `extract_vegetation_variance.R`, `analyze_paleo_inversion_models.R` | W5.6 (one commit per script, 4 commits) — after W0 renames |
| `5a_model_validation.R`, `5b_overfitting_diagnostics.R`, `5c_spatial_patterns.R`, `5d_spatial_confounding_check.R`, `5e_spatial_intercept_correlations*.R` | W5.10 (one commit per script, 6 commits) |

**Figure_05 correction.** Current code hardcodes σ = 15.2 / 21 ‰ but `sigma` in the fit is on the standardized scale. Use either `sigma_residual_original` if exported by Stan (`4d_leaf_wax_spatial_model.stan:479`), or back-transform: `sigma_orig <- mean(load_draws(m) |> subset_draws("sigma") |> as.numeric()) * stan_data$scaling_params$d2H_wax_sd_original`. Commit message documents which quantity was used and why.

**Acceptance per commit:** script runs end-to-end against `results/c2_run_20260414/` and emits its artifacts (PDF, CSV, or whatever it writes). No call to `fit$*()` remains.

### W6 — Table producers as pipeline outputs (5 commits)

**6a. `manuscript/table_code/table_helpers.R`** with two emitters:
- `emit_standalone_tex(df, out_path, caption, label, ...)` — for tables that are `\begin{table}...\end{table}` snippets.
- `emit_tabular_fragment(df, out_path, ...)` — for tables whose canonical `.tex` is a hand-maintained wrapper (landscape longtable, row coloring, grouped rows).

**6b-d. Producer edits (3 commits).**

| Table | Producer | Emitter | Notes |
|---|---|---|---|
| `table1_model_performance.tex` | `extract_full_model_analysis.R` | `emit_standalone_tex` | After W0+W5, uses `load_loo(model)` across all 14 models. |
| `table2_global_params.tex` | `analyze_paleo_inversion_models.R` + `extract_full_model_analysis.R` | `emit_tabular_fragment` → `table2_global_params_body.tex` | **Keep canonical filename as wrapper**; wrapper contains preamble + `\input{table2_global_params_body.tex}`. |
| `table3_variance_decomp.tex` | `extract_full_model_analysis.R` | `emit_standalone_tex` | After W0, all 9 spatial models. |
| `table4_environmental.tex` | `extract_vegetation_coefficients.R` + `extract_full_model_analysis.R` | `emit_standalone_tex` | β_C4, β_tree, β_shrub, β_grass, β_precip with 90% CIs. |

**6e. Table_S2 + numeric audit (1 commit).**

- `5b_overfitting_diagnostics.R` adds `write.csv(...)` + `emit_tabular_fragment(...)` for `Table_S2_regional_performance_body.tex`.
- `Table_S2_regional_performance.tex` becomes a hand-maintained wrapper with longtable / row coloring / Overall grouping + `\input{table_S2_regional_performance_body.tex}`.

New `manuscript/table_code/check_tables_numeric.R`:
- Parses each `.tex` for numeric cells in tabular rows (row/column aware — uses known header row to index columns).
- Loads companion CSV.
- For each (row, column), checks the `.tex` value matches the CSV value within a small tolerance (default `1e-3` for numbers, exact for integers).
- Errors with a structured diff: "Row X, column Y: tex has A.BC, csv has D.EF."

### W7 — Makefile (1 commit)

```make
APRIL  := results/c2_run_20260414
PREP   := $(APRIL)/_prepared_data

MODELS := baseline baseline_sp baseline_env baseline_env_sp \
          baseline_veg baseline_veg_sp full full_sp \
          full_interact full_interact_sp c4_only_sp \
          elevation_only_sp elevation_c4_sp elevation_c4_interact_sp

SED := $(APRIL)/3_sediment_ready_for_modeling.rds

draws    = $(APRIL)/$(1)/posterior_draws.rds
diag     = $(APRIL)/$(1)/diagnostics.rds
loo      = $(APRIL)/$(1)/loo.rds
stand    = $(PREP)/stan_data_$(1).rds
conf     = $(PREP)/config_$(1).rds

manuscript/tables/table1_model_performance.tex: \
    $(foreach m,$(MODELS),$(call loo,$(m))) \
    $(foreach m,$(MODELS),$(call diag,$(m))) \
    extract_full_model_analysis.R
	Rscript extract_full_model_analysis.R

manuscript/figures/main_figs/Figure_01.pdf: \
    $(call draws,baseline) $(call diag,baseline) $(call stand,baseline) $(SED) \
    manuscript/figure_code/Figure_01_ols_regression.R \
    manuscript/figure_code/common_functions.R \
    manuscript/figure_code/posterior_helpers.R
	cd manuscript/figure_code && Rscript Figure_01_ols_regression.R

# ... one recipe per figure and per table
# Figure 4 recipe depends on stan_data_full_sp.rds explicitly.

.PHONY: figures tables audit pipeline
figures: Figure_01.pdf Figure_02.pdf figure_03.pdf Figure_04.pdf Figure_05.pdf
tables:  table1.tex table2.tex table3.tex table4.tex Table_S2.tex
audit:   tables
	Rscript manuscript/table_code/check_tables_numeric.R
pipeline: figures tables audit
```

**Acceptance:** `make -n pipeline` lists every step. `make pipeline` from a clean clone + populated `results/c2_run_20260414/` runs to completion. `git diff manuscript/` after rebuild shows numeric deltas (expected) or nothing. `make audit` errors on any row/column mismatch.

### Not in this phase

- Supplement producer tracing (S1, S3, S5).
- `table_S3_compilations-UNRELIABLE.tex` decision.
- `launch_pipeline.sh` / `slurm/job_fit.sh` orchestration hardening (operator layer).
- CI smoke test.

## Commit ordering and dependencies (honest)

| # | Commit | Depends on | State of master after |
|---|---|---|---|
| 1 | W0 — param drift + coverage | none | buildable; analysis still fit.rds-based; tables empty for some cells until rerun |
| 2 | W1 — widen 4c + emit loo | W0 (for consistent names) | master buildable on HPC; April rds files not yet widened |
| **—** | **W2 — C2 rerun (operator action, NOT a commit)** | W1 code present on C2 | April mirror now holds widened rds + loo |
| 3 | W3 — Figure 3 bug fix | independent | Figure 3 correct |
| 4 | W4 — helper layer | W1 (for widened schema) | helper present, nothing yet uses it |
| 5a | W5.1 — Figure_01 migrate | W2 data, W4 helper | Figure 1 portable |
| 5b | W5 — figures 03/04/05 | W2 data, W4 helper, W3 | each figure portable |
| 6a | W5.6 — extract/analyze migrate (4 commits) | W0 + W2 + W4 | extract scripts portable |
| 6b | W5.10 — 5*-series migrate (6 commits) | W0 + W2 + W4 | 5*-series portable |
| 7 | W6 — table producers + audit (5 commits) | W5.6 ready | tables regenerable |
| 8 | W7 — Makefile | everything above | `make pipeline` works |

**Honest caveat:** commits 5a through 8 require W2 (the one-time C2 rerun) to have happened, and W2 is not a git commit. Between W1 landing and W2 completing, those later commits cannot be merged. Timeline: W1 → W2 on C2 (same day ideally) → rest of the train. If W2 gets delayed, either stall the rest of the train or keep a fallback path in the helpers temporarily.

**Total commits:** ~18 (1 W0 + 1 W1 + 1 W3 + 1 W4 + 4 W5 figs + 4 W5 extract + 6 W5 fives + 5 W6 + 1 W7).

## Acceptance criteria for Phase 5

1. **Pipeline invariant:** running `4c_fit_models.R` from scratch on a fresh stan_data set produces every rds file downstream analysis requires. Verified by running against one model end-to-end on C2.
2. **Analysis portability:** every figure, table producer, and `5*` analysis script runs on any machine with R + the rds files. Zero reads of `fit.rds`. Zero chain-CSV dependency.
3. **Parameter-name drift eliminated:** `grep` in W0 acceptance returns nothing. Manual verification: every Stan reference in analysis code has a matching declaration in `4d_leaf_wax_spatial_model.stan`.
4. **Figure 3 correctly represents spatial intercept contribution in original units:** corrected transform documented in the commit message; nondimensional analysis checks out.
5. **`make pipeline` reproduces** every main figure, tables 1–4, Table_S2 from a clean clone + populated April mirror.
6. **`make audit` passes** — every numeric cell in every regenerated `.tex` matches the companion CSV at its (row, column) position.
7. **FIGURES.md + TABLES.md updated** — every main artifact tagged `regenerate = pipeline` with its producer path; open items (S1, S3, S5, Table_S3) explicitly listed.

## Risks

| Risk | Mitigation |
|---|---|
| W2 C2 rerun hits a memory ceiling (widened log_lik draws) | Run models serially on `general-interactive`; request 64–128 GB. |
| Widened `posterior_draws.rds` too big to rsync to laptop | Document minimal subset (`full_sp`, `baseline`, `baseline_sp`) needed for lightweight analysis; keep others on external disk. |
| Figure 3 fix changes the manuscript's "±60 ‰" narrative | Check before committing; if magnitude shifts, stop and flag for author review. |
| Parameter renames in W0 reveal Stan doesn't export a quantity we thought it did | Delete the extraction, note in TABLES.md, flag whether the table column can be computed from a different parameter. |
| W6 wrapper filename mismatches break LaTeX build | Test compile the manuscript PDF locally before committing wrapper restructure. |
