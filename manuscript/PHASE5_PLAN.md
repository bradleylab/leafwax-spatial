# Phase 5 plan v3 — self-contained rds outputs + portable analysis

**Date drafted:** 2026-04-16 (third rewrite; v1/v2 in git history)

## Principle

The fit step produces **self-contained rds files** that travel anywhere. Analysis scripts read only those rds files — never `fit.rds`, never chain CSVs, never network resources. If someone reruns the fit step from scratch they get everything they need the first time; no out-of-band scripts, no sync steps, no hidden state.

## Why v1 and v2 got complicated

- v1 assumed `posterior_draws.rds` was a full posterior. It isn't — `4c_fit_models.R:205-220` saves only 5–17 scalar global parameters.
- v2 tried to paper over the gap with `diagnostics.rds` summaries (which have means + CIs for all 19,802 params but no per-draw samples) plus a one-off `generate_loo.R`.
- Both treated the shortfall as an analysis-layer problem. It's a fit-layer problem: the fit step saves too little.

## Fix

Widen the fit step's output so downstream analysis has per-draw access to everything it plots. Emit `loo.rds` at the same step. After Phase 5, the fit step's product per model is:

```
<model>/
  fit.rds             # cmdstanr stub, non-portable; kept for archival/on-HPC use
  posterior_draws.rds # WIDENED: all variables downstream analysis reads
  diagnostics.rds     # rhat/ess/sampler diagnostics (unchanged)
  loo.rds             # NEW: psis_loo object for model comparison
  runtime_info.rds    # (unchanged)
```

The analysis layer reads `posterior_draws.rds`, `diagnostics.rds`, `loo.rds`. That's it.

## Data sizes

Back-of-envelope for `full_sp` (1131 sites, 8 chains × 2000 = 16,000 draws, 8 bytes/double):

| Variable | Size per model |
|---|---|
| Current scalar set (17 vars) | ~1 MB |
| `log_lik[1..1131]` | 145 MB |
| `mu[1..1131]` | 145 MB |
| `d2H_rep[1..1131]` | 145 MB |
| `alpha_spatial[1..1131]` | 145 MB |
| `z_intercept_spatial[1..250]` + `z_slope_spatial[1..250]` | 64 MB |
| `scale_weights[1..9]` | 1 MB |
| **Widened total per model** | **~650 MB** |

14 models × ~650 MB ≈ **~9 GB** across the run. Movable by rsync/scp/external drive. Not small, not unreasonable for a Bayesian spatial calibration at this scale.

Non-spatial models (no `z_*`, no `alpha_spatial`) are ~430 MB each.

## Workstreams

### W0 — Parameter name drift (1 commit)

Several extract scripts reference names that never existed in the Stan model. Fix before anything else, because widening the draws set doesn't help if we're asking for the wrong names.

| Script | Old | New |
|---|---|---|
| `extract_full_model_analysis.R:103` | `beta_precip_amount` | `beta_precip` |
| `extract_full_model_analysis.R:107-109` | `beta_oipc_c4`, `beta_oipc_tree`, `beta_oipc_shrub`, `beta_oipc_grass` | `beta_oipc_x_c4`, `beta_oipc_x_tree`, `beta_oipc_x_shrub`, `beta_oipc_x_grass` |
| `extract_full_model_analysis.R:107-109` | `beta_precip_c4`, `beta_precip_tree`, `beta_precip_shrub`, `beta_precip_grass` | Verify against `4d_leaf_wax_spatial_model.stan` — if not exported, remove the extraction and note in TABLES.md |
| `extract_vegetation_coefficients.R:204` | `sigma_gp_intercept`, `sigma_gp_slope` | `sigma_intercept_spatial`, `sigma_slope_spatial` |
| `analyze_paleo_inversion_models.R:171` | same family as above | cross-check |

Also expand model coverage where the extract scripts silently cover a subset:
- `extract_full_model_analysis.R:20-22` currently covers 8 of 14 models → extend to all 14 for Table 1.
- `extract_full_model_analysis.R:159` variance loop covers 5 spatial models → extend to all 9 `_sp` models for Table 3.

**Acceptance:**
```bash
grep -rn 'beta_precip_amount\|beta_oipc_c4\|beta_oipc_tree\|beta_oipc_shrub\|beta_oipc_grass\|beta_precip_c4\|beta_precip_tree\|beta_precip_shrub\|beta_precip_grass\|sigma_gp_intercept\|sigma_gp_slope' *.R manuscript/
# returns nothing
```

### W1 — Widen `4c_fit_models.R` output + emit `loo.rds` (1 commit)

Edit `4c_fit_models.R:205-220` to build a `draws_to_save` variable list that includes **everything downstream analysis plots or computes with**:

```r
# Scalars (retained from current behavior)
draws_to_save <- params_to_check

# Observation-level (1131 sites)
draws_to_save <- c(draws_to_save, "mu", "d2H_rep", "log_lik")

# Spatial effects (only if GP is in the model)
if (stan_data$include_gp == 1) {
  draws_to_save <- c(draws_to_save, "alpha_spatial",
                     "z_intercept_spatial", "z_slope_spatial")
}

# Multi-scale weighting (always present)
draws_to_save <- c(draws_to_save, "scale_weights")
```

Replace the current save call at `4c_fit_models.R:267-268` with a single expanded call:

```r
draws <- fit$draws(variables = draws_to_save, format = "draws_df")
saveRDS(draws, file.path(output_dir, "posterior_draws.rds"))
```

Immediately after the draws save, add:

```r
cat("Computing loo...\n")
log_lik_array <- fit$draws("log_lik", format = "draws_array")
loo_result <- loo::loo(log_lik_array, cores = getOption("mc.cores", 4))
saveRDS(loo_result, file.path(output_dir, "loo.rds"))
```

Also update METHODS.md to document the expanded rds contents and the automatic `loo.rds` emission.

**Reproducibility invariant:** running `4c_fit_models.R` from scratch produces every rds needed for downstream analysis. No followup script. This is the test — someone cloning the repo and rerunning should need exactly one command, not two.

**Acceptance:** a dry-run of `4c_fit_models.R` for one small model (baseline) locally or on a test machine emits all four rds files. Inspect `posterior_draws.rds` — variable count matches the expected set. Inspect `loo.rds` — class is `psis_loo`.

### W2 — Backfill the April 2026 run (1 commit for the script; 1 separate C2 job for the execution)

The April fits pre-date W1, so their `posterior_draws.rds` is still narrow and there is no `loo.rds`. New script `4e_postprocess_existing_fits.R`:

```r
# Re-read existing fit.rds (works on the machine where chain CSVs live)
# Re-save a widened posterior_draws.rds using the same draws_to_save list as 4c
# Emit loo.rds
```

Run once on Compute2 (`general-interactive` partition or a short CPU batch job, ~3 hours total for 14 models) against `/scratch2/fs1/alexander.s.bradley/leafwax_run/model_output_20260414/`. If the chain CSVs are no longer on scratch, restore them first from `s3://bradleylab-public/tmp/leafwax_run/model_output_20260414/` — that's a one-liner, and irrelevant to anyone reproducing in the future (their fits will emit everything on the first pass via W1).

Sync the 14 refreshed `posterior_draws.rds` + 14 new `loo.rds` back to the canonical April mirror location. Commit the script; the resulting rds files live with the other run outputs (already git-ignored under `results/`).

**This is the only step that touches S3 and it is one-time. New fits under the W1-patched pipeline bypass it entirely.**

**Acceptance:** `ls results/c2_run_20260414/*/loo.rds` shows 14 files. `Rscript -e 'posterior::variables(readRDS("results/c2_run_20260414/full_sp/posterior_draws.rds")) |> length()'` reports ≥1200 (scalars + indexed arrays).

### W3 — Analysis helper layer (1 commit)

New file `manuscript/figure_code/posterior_helpers.R`. Three entry points:

```r
APRIL_RUN <- normalizePath(file.path("..", "..", "results", "c2_run_20260414"), mustWork = FALSE)

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
         " — rerun 4c_fit_models.R (post-W1) or 4e_postprocess_existing_fits.R (April run).")
  }
  readRDS(path)
}
```

No custom mean / quantile wrappers — callers use `posterior::summarise_draws()`, `posterior::subset_draws()`, or direct indexing on the `draws_df`. The `diagnostics.rds` summary tibble is already well-shaped.

**Acceptance:** `Rscript -e 'source("manuscript/figure_code/posterior_helpers.R"); length(load_draws("baseline"))'` runs, prints a count.

### W4 — Figure migrations (5 commits, one per script)

Each figure reads only `posterior_draws.rds` (widened) or `diagnostics.rds`. No `fit$*()` calls anywhere. Flag and/or fix bugs discovered during migration rather than carrying them forward silently.

| Script | Replace | With |
|---|---|---|
| `Figure_01_ols_regression.R` | all `fit_baseline$summary()` / `$draws()` | `draws <- load_draws("baseline")`; posterior-predictive ribbons from `subset_draws(draws, c("mu","d2H_rep","scale_weights"))` |
| `figure_03_spatial_confounding_revised.R` | all `fit$*()` on `baseline`, `baseline_sp` | `load_draws("baseline")`, `load_draws("baseline_sp")` with `alpha_spatial` array. **Fix likely double-scaling bug at lines 87-91** (`alpha_spatial` multiplied by `sigma_intercept_spatial`, but Stan at `4d_leaf_wax_spatial_model.stan:316,331-335` already has `alpha_spatial` on the intercept scale). Verify against Stan generated quantities before keeping the multiplication. |
| `Figure_04_spatial_maps.R` | iterated `fit$draws("z_intercept_spatial[k]")` | `load_draws("full_sp")` → `subset_draws(draws, variable = c("z_intercept_spatial","z_slope_spatial"), regex = FALSE)` |
| `Figure_05_detection_thresholds.R` | hardcoded σ = 15.2, 21 ‰ | `sigma_sp <- mean(load_draws("full_sp")$sigma)`; same for baseline |
| `Figure_02_all_environmental_variables_noborders.R` | data-only; no changes needed |  |

**Acceptance per script:** running the script end-to-end emits its PDF + PNG without error. Run it on the target machine; capture `Rscript <script>` stderr to confirm.

### W5 — Table producers as pipeline outputs (5 commits)

**W5a. `manuscript/table_code/table_helpers.R` (1 commit).** Two emitters:

- `emit_standalone_tex(df, out_path, caption, label, ...)` — full `\begin{table}...\end{table}` document for simple tables.
- `emit_tabular_fragment(df, out_path, ...)` — tabular body only, meant to be `\input`'d by a hand-maintained wrapper.

Both prepend a `% AUTOGENERATED <timestamp> from <script> — do not edit` banner.

**W5b-d. Producer edits (3 commits).**

| Table | Producer | Emitter | Notes |
|---|---|---|---|
| `table1_model_performance.tex` | `extract_full_model_analysis.R` | `emit_standalone_tex` | Needs `load_loo(model)` for LOOIC/SE/p_eff/n_hi_k across all 14 models (after W0 coverage expansion). |
| `table2_global_params.tex` | `analyze_paleo_inversion_models.R` + `extract_full_model_analysis.R` (both contribute per `TABLES.md:18-19`) | `emit_tabular_fragment` → `table2_global_params_body.tex` | **Keep `table2_global_params.tex` as the canonical filename and make it the wrapper** (contains the `\documentclass`, `longtable` preamble, `\input{table2_global_params_body.tex}`). Do NOT rename. |
| `table3_variance_decomp.tex` | `extract_full_model_analysis.R` | `emit_standalone_tex` | After W0, covers all 9 spatial models. |
| `table4_environmental.tex` | `extract_vegetation_coefficients.R` + `extract_full_model_analysis.R` | `emit_standalone_tex` | Merge of β_C4, β_tree, β_shrub, β_grass, β_precip with 90% CIs. |

**W5e. Table_S2 + numeric audit (1 commit).**

- `5b_overfitting_diagnostics.R` adds `write.csv(regional_performance, "manuscript/tables/table_S2_regional_performance.csv")` and `emit_tabular_fragment(..., "table_S2_regional_performance_body.tex")`.
- Canonical `Table_S2_regional_performance.tex` becomes the wrapper (keeps `longtable`, row coloring, "Overall" grouping by hand) with an `\input{table_S2_regional_performance_body.tex}` where the tabular body is today.

New script `manuscript/table_code/check_tables_numeric.R`:
- Parses each `.tex` for numeric cells in its tabular rows.
- Loads the companion CSV.
- For every numeric cell in the `.tex`, checks that the same value appears in the matching CSV row/column, **within the expected tolerance and in the expected position** — not a token-presence grep.
- Errors with a row diff if a cell is out of place.

### W6 — Makefile (1 commit)

Top-level `Makefile` with real prereqs. Every figure/table target depends on the specific model's `posterior_draws.rds`, `diagnostics.rds`, `loo.rds` (where used), and the sediment rds. Recipes `cd` into the directory the script expects so cwd-coupling is explicit.

```make
APRIL := results/c2_run_20260414
MODELS := baseline baseline_sp baseline_env baseline_env_sp \
          baseline_veg baseline_veg_sp full full_sp \
          full_interact full_interact_sp c4_only_sp \
          elevation_only_sp elevation_c4_sp elevation_c4_interact_sp

SED := $(APRIL)/3_sediment_ready_for_modeling.rds
draws = $(APRIL)/$(1)/posterior_draws.rds
diag  = $(APRIL)/$(1)/diagnostics.rds
loo   = $(APRIL)/$(1)/loo.rds

manuscript/tables/table1_model_performance.tex: \
    $(foreach m,$(MODELS),$(call loo,$(m))) \
    $(foreach m,$(MODELS),$(call diag,$(m))) \
    extract_full_model_analysis.R
	Rscript extract_full_model_analysis.R

manuscript/figures/main_figs/Figure_01.pdf: \
    $(call draws,baseline) $(call diag,baseline) $(SED) \
    manuscript/figure_code/Figure_01_ols_regression.R \
    manuscript/figure_code/common_functions.R \
    manuscript/figure_code/posterior_helpers.R
	cd manuscript/figure_code && Rscript Figure_01_ols_regression.R

# ... one recipe per figure + table

.PHONY: figures tables pipeline
figures: Figure_01.pdf Figure_02.pdf figure_03.pdf Figure_04.pdf Figure_05.pdf
tables:  table1.tex table2.tex table3.tex table4.tex Table_S2.tex
pipeline: figures tables
```

**Acceptance:** `make -n pipeline` from a clean clone + populated `results/c2_run_20260414/` prints every step needed. `make pipeline` runs to completion. `git diff manuscript/` shows only numeric deltas (expected) or nothing.

### Not in Phase 5

- Supplement S1 (`Figure_S1_03_spatial_weighting`) and S5 (`Figure_S5_elevation`) producer tracing — separate task. `FIGURES.md` carries them as open items.
- S3 (`Figure_S3_residuals`) — per `FIGURES.md`, producer is likely `5a_model_validation.R` or `5b_overfitting_diagnostics.R`; move to Phase 5 if trivial during W5, otherwise defer.
- `table_S3_compilations-UNRELIABLE.tex` — delete or rebuild from `data/HrenBrandon/`; decision deferred.
- `slurm/job_fit.sh:44` skip-existing logic and `launch_pipeline.sh:155` shutdown timing: noted as future hardening for the HPC orchestration. Not in scope because v3's reproducibility contract is "running `4c_fit_models.R` produces everything" — how SLURM is wired around that is an operator concern, not a repo concern.
- CI smoke test — Phase 6.

## Ordering and commit count

1. **W0** param name drift (1 commit)
2. **W1** widen `4c_fit_models.R` + emit `loo.rds` (1 commit)
3. **W2** `4e_postprocess_existing_fits.R` + one-time backfill (1 commit for the script; the C2 run is an operator action, not a git commit)
4. **W3** helper layer (1 commit)
5. **W4** figures (5 commits)
6. **W5** tables + audit (5 commits)
7. **W6** Makefile (1 commit)

**Total:** 15 commits. Each leaves `master` in a regenerable state as long as W2's backfill has run, which it will before any post-W3 work that depends on the widened draws.

## Acceptance criteria

1. **Pipeline invariant.** Running `4c_fit_models.R` from scratch produces every rds needed for Phase-5-downstream analysis. No followup script required. Verified by inspecting a fresh baseline-model run.
2. **Analysis portability.** Every figure/table script runs on any machine with R + the rds files. Zero reads of `fit.rds`. Zero network access. Zero chain-CSV dependency.
3. **`make pipeline` reproduces** every main text figure, tables 1–4, and Table_S2 from a clean clone + populated `results/c2_run_20260414/`.
4. **`check_tables_numeric.R` passes** row/column-aware — not token presence.
5. **Parameter name drift eliminated** — the full grep in W0 returns nothing.
6. **`manuscript/FIGURES.md` and `manuscript/TABLES.md` updated** — every main figure and every main-text/Table_S2 row tagged `regenerate = pipeline` with its producer path.

## Open questions for the author

- `beta_precip_*` interaction variants listed in the old extract scripts — do they exist in the Stan generated quantities, or are they leftover from an earlier model? Need to check `4d_leaf_wax_spatial_model.stan` before deciding whether to rename or delete in W0.
- Figure_03 suspected double-scaling bug (W4 note): confirm against the Stan model before "fixing" — may be correct as written if `alpha_spatial` in the fit is pre-unscaled.
