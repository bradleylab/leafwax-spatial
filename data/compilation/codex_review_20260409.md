# Codex Challenge Review — Leafwax Spatial Pipeline

**Date:** 2026-04-09
**Reviewer:** OpenAI Codex (gpt-5, reasoning=high), via `/codex challenge`
**Mode:** Adversarial challenge — actively try to break the pipeline
**Purpose:** Pre-flight check before burning ~5 days of EC2 compute on the final
14-model manuscript run.
**Snapshot reviewed:** `/tmp/leafwax_pipeline_for_review/` (git SHA `004bf54`),
pulled from EC2 `i-053cbbe7833ef0b97:/home/shared/leafwax_spatial/spatial_leafwax_model/`
on 2026-04-09 ~02:04 UTC.
**Context passed to codex:** full pipeline (13 files, 5,634 lines: R + Stan + bash + config)
plus `METHODS.md` (578-line manuscript-vs-code cross-reference, ~2 weeks stale).
**Token usage:** 4,263,523

**Gate: FAIL** — 4 P0 issues found. Do not start the run as-is.

---

## Claude evaluation summary (2026-04-09, post-codex)

After codex completed, Claude verified each finding against the actual code in the same snapshot and against raster metadata from EC2. The findings hold up broadly but severities shift. **Revised blocker list below.** The per-finding evaluation is interleaved under each section as "**Claude evaluation:**" notes.

### Revised blocker list (what actually blocks the 5-day run)

**P0 — must fix before launch:**

1. **P0-1 grid mismatch in interaction predictors.** VERIFIED with raster dimensions: OIPC 2083 × 4320 vs C4 360 × 720 vs PFT 2186 × 4371 — ~140× pixel-count mismatch in a 5° radius between OIPC and C4. `compute_weighted_interaction()` length-checks and returns `NA_real_` on every call; subsequent row-mean NA fill (`4a:1071-1099`) returns `fill_value=0` when the whole row is NA. **Every interaction predictor is uniformly zero across all sites and scales.** All 5 interaction models (`baseline_veg`, `full_interact`, `baseline_veg_sp`, `elevation_c4_interact_sp`, `full_interact_sp`) are structurally identical to their non-interaction counterparts. **This is worse than codex stated.** It also means every previous 818-row and 1002-row run had the same bug, so historical interaction-model results are noise.

2. **P0-4 partial-success shutdown.** VERIFIED. `launch_pipeline_4_parallel_w_shutdown.sh:482` requires only `>0` `fit.rds` files before shutdown. 13 of 14 would trigger shutdown with no warning.

**Should fix before launch (real impact, not pure P0):**

3. **P1-1 unused climate variable row dropping.** VERIFIED. `3_prep_data.R:560-562` drops sites missing `tc_soil_mean`, `tc_tmax_mean`, `tc_vpd_mean`, but `4b_stan_prep.R:219-221` hardcodes `include_temp/vpd/soil = FALSE` with comment "Always FALSE due to collinearity". Silent data loss on variables that never enter any model.

4. **P1-6 precip prior.** `beta_precip_raw[1] ~ normal(0, 0.02)` at `4d:408` is 100× tighter than peers (all `normal(0, 2)`). On a standardized predictor, 95% CI is ±0.04 — de facto exclusion of precipitation as a predictor. **Needs user decision:** intentional "test via near-zero prior" or leftover from debugging? Either way, should be justified in METHODS.md.

5. **P0-2 OIPC SE units (downgraded to P1).** Codex's dimensional-consistency argument is correct: `oipc_se_weighted[n]` is in (std OIPC)² units but gets added to response variance in (std response)² units. The correct term is `(beta_oipc_effective * oipc_se_weighted[n])²`. **But:** direction is conservative (over-inflates uncertainty by ~2–4×) and `oipc_sd ≈ d2H_wax_sd` keeps magnitudes accidentally in the ballpark. Model still runs correctly. Should be fixed for publication defensibility, not urgency.

**Can fix while editing but don't block:**

6. P0-3 pipefail. Trivial 1-line fix (`set -eo pipefail`). Defense-in-depth at `launch:144` catches total prep failure.
7. P1-3 YAML parser brittleness. yq is installed at `/home/ubuntu/bin/yq`, reachable in login shells only. Launcher has worked from interactive SSH before — no crisis. Fix by adding `export PATH="/home/ubuntu/bin:$PATH"` at top of launcher.
8. P1-4 PFT 0.33 imputation. Decision: drop sites with no PFT coverage OR document the uniform prior in METHODS.md.
9. P1-5 resume scope bug in `4c_fit_models.R`. Doesn't fire via launcher path (per-model R sessions).

**Pure cosmetic (P2):**
- P2-1 `include_veg_interactions` subset semantics — document in METHODS.md.
- P2-2 approximate deg-to-km — minor accuracy near poles.
- P2-3 `typical_n_obs <- 800` hardcoded in `3g_optimize_knots.R:146` — diagnostic-only impact.

### Bugs codex missed

1. **`4c_fit_models.R:189-190` scope bug in error handler.** Inside `tryCatch(..., error = function(e) { ...; status <- "failed"; fit <- NULL })`, the assignments are local to the error function. They do NOT modify the outer-scope `status` and `fit`. After a failed model in multi-model invocation, the outer `fit` retains the *previous* successful model's object and the diagnostics block at line 195 contaminates the failing model's `fit_summary` entry with stale diagnostics. **Doesn't fire via the launcher path** (each model runs in its own R process with only one model name), but fires on any direct invocation. Fix: use `fit <<- NULL` and `status <<- "failed"` or restructure to set these in the caller's environment.

2. **`count_running_models()` regex in launcher line 191** is brittle (`ps aux | grep -E "fit_single_model.R|run_.*\.R"`) but does not self-match the launcher or `run_final_pipeline.sh` (neither name ends in `.R`). Worth flagging but not actionable.

3. **Historical implication of P0-1:** every past 818-row and 1002-row manuscript run had zero-valued interaction predictors. LOO scores, parameter estimates, and model comparison rankings from past runs involving any interaction model are not trustworthy.

### Environmental note

`yq` IS installed on EC2 at `/home/ubuntu/bin/yq`. Reachable in login shells (`bash -lc`) but not in non-interactive ssh commands. Old manuscript-run logs in `Ω_archive/old_logs/` (e.g. `baseline_20250807_111712.log`) confirm the launcher has been run successfully in the past. Codex's initial P1-3 "launcher cannot start without yq" is incorrect in the common case; it applies only to non-interactive invocations.

---

## Verdict

Multiple P0 issues. Codex would **not** start the EC2 run as-is.

---

## P0 — DO NOT RUN

### P0-1. Interaction predictors computed on mismatched raster grids, silently zero-imputed

**Files/lines:**
- `3_prep_data.R:230` — C4 extraction
- `3_prep_data.R:244` — OIPC extraction
- `3_prep_data.R:283` — PFT extraction
- `4a_spatial_functions.R:59` — `compute_weighted_interaction()`
- `4a_spatial_functions.R:742,754,765,776` — interaction calls passing OIPC distances with C4/PFT values
- `4a_spatial_functions.R:1071` — NA zero-fill

**What's wrong:** C4, OIPC, and PFT are extracted on different raster grids
(different resolutions/projections). `compute_weighted_interaction()` requires
equal-length, aligned vectors. The interaction calls pass OIPC distances
together with C4 or PFT values "assuming aligned grids". When the grids don't
match, the mismatched-length inputs return NA, which is then zero-filled at
line 1071.

**Why it matters:** Corrupts **every interaction model**:
- `baseline_veg`
- `full_interact`
- `baseline_veg_sp`
- `elevation_c4_interact_sp`
- `full_interact_sp`

That's 5 of the 14 manuscript models, including `full_interact_sp` which is a
core comparison.

**Fix direction:** Store aligned pixel coordinates, or compute interactions
only after aligned site-scale aggregation (product of already-aggregated values
at each site).

---

### P0-2. OIPC uncertainty added to likelihood in wrong units / without slope scaling

**Files/lines:**
- `4a_spatial_functions.R:876` — OIPC SE scaled by `oipc_sd` (predictor units)
- `4d_leaf_wax_spatial_model.stan:456,468` — SE added directly to response variance

**What's wrong:** The OIPC standard error is in predictor units (precipitation
dD, roughly tens of permil). The Stan model adds it to the response variance
alongside `d2H_wax_err` and `sigma`. Propagating input uncertainty into the
response likelihood requires multiplying by the local slope: the correct term
is `(beta_precip * oipc_se)^2`, not `oipc_se^2`.

**Why it matters:** Dimensional inconsistency. Uncertainty is miscalibrated for
**all models** that include OIPC uncertainty propagation. Likely inflates
response variance by a factor of `1/beta_precip^2` relative to intent.

**Fix direction:** Multiply propagated SE by the local OIPC slope on the
response scale before adding to response variance.

---

### P0-3. Launcher can miss prep failures (missing `set -o pipefail`)

**Files/lines:**
- `launch_pipeline_4_parallel_w_shutdown.sh:7` — `set -e`
- `launch_pipeline_4_parallel_w_shutdown.sh:123` — `Rscript "$PREP_SCRIPT" 2>&1 | tee "$log_file"`
- `launch_pipeline_4_parallel_w_shutdown.sh:125` — exit code check

**What's wrong:** The script uses `set -e` but not `set -o pipefail`. When
`Rscript ... | tee ...` runs, `$?` is `tee`'s exit code (almost always 0),
not `Rscript`'s. A failed `4b_stan_prep.R` looks successful and the launcher
proceeds to fit models on broken/missing stan_data.

**Fix direction:** `set -eo pipefail` at top of launcher. Also check every
piped command's exit status explicitly.

---

### P0-4. Partial failure can still trigger auto-shutdown

**Files/lines:**
- `4c_fit_models.R:176` — failed-fit path
- `4c_fit_models.R:194` — only appends summary rows when `fit` exists
- `4c_fit_models.R:287` — writes completion marker unconditionally
- `launch_pipeline_4_parallel_w_shutdown.sh:410` — touches `pipeline_complete_*.txt`
- `launch_pipeline_4_parallel_w_shutdown.sh:482` — only requires `>0` `fit.rds` files before shutdown

**What's wrong:** If some models fail mid-run:
1. Failed models disappear from `fit_summary` (not appended)
2. Completion marker still gets written
3. Launcher shuts down the instance as long as **any** `fit.rds` exists

**Why it matters:** You can lose the instance with 1 of 14 models completed and
no obvious way to resume. Critical for a 5-day run where any individual model
could fail.

**Fix direction:** Shutdown trigger must require `fit_count == EXPECTED_MODELS`
AND zero `error_info.rds` files. Also persist a failure manifest separately.

---

## P1 — LIKELY BAD

### P1-1. Dropping rows on climate variables that are never used

**Files/lines:**
- `3_prep_data.R:551,560,561,562` — require soil/temp/VPD completeness in `sediment_ready`
- `4b_stan_prep.R:219,220,221` — hard-disable those covariates for every model

**What's wrong:** Prep step drops sites missing soil/temp/VPD, but `4b_stan_prep.R`
hard-disables those same covariates in every single model. Silent data loss for
variables that don't even enter the model.

**Why it matters:** Reduces effective N below 1131 with no documentation of
which sites were dropped or why. Biases the spatial distribution of the retained
dataset toward climate-data-rich regions.

**Fix direction:** Drop the completeness filter for variables not used by any
model in the config.

---

### P1-2. `3h_prior_checks.R` currently crashes

**Files/lines:**
- `3h_prior_checks.R:221` — defines `75/100/125/150` knot configs
- `3h_prior_checks.R:339` — indexes `all_densities[["50_knots"]]`

**What's wrong:** Index mismatch. `50_knots` isn't in the configs list.

**Why it matters:** Blocks the pre-pipeline gate. Easy fix.

---

### P1-3. Launcher YAML fallback parser broken without `yq`

**Files/lines:**
- `launch_pipeline_4_parallel_w_shutdown.sh:35` — grep fallback
- `launch_pipeline_4_parallel_w_shutdown.sh:53` — abort on empty

**What's wrong:** The grep fallback looks for literal keys like
`scripts.prep_script`, but the YAML only contains nested `prep_script:`. On a
box without `yq`, `PREP_SCRIPT` comes back empty and the launcher aborts.

**Why it matters:** Environment-dependent crash. Verify `yq` is installed on the
EC2 box OR fix the fallback parser.

---

### P1-4. Missing PFT silently imputed to `0.33/0.33/0.33`

**Files/lines:**
- `3_prep_data.R:293` — permits zero-PFT-pixel rows
- `3_prep_data.R:551` — does not filter them
- `4a_spatial_functions.R:1040,1144` — replace missing PFT scales with `0.33`

**What's wrong:** Sites with no PFT pixels in range get an undocumented uniform
vegetation prior (1/3 tree, 1/3 shrub, 1/3 grass) in every PFT model.

**Why it matters:** Silently injects a synthetic vegetation signal for data-sparse
sites. Could bias vegetation coefficients toward zero and inflate uncertainty.

**Fix direction:** Drop sites with no PFT coverage, or add an explicit
missingness indicator column.

---

### P1-5. Re-running `4c_fit_models.R` over existing `fit.rds` uses undefined `stan_data`

**Files/lines:**
- `4c_fit_models.R:97` — loads existing fit without re-loading its data
- `4c_fit_models.R:203` — branches on `stan_data$include_c4`/`include_pft`

**What's wrong:** When a model's `fit.rds` already exists, the script loads the
fit object but never loads the corresponding `stan_data_<model>.rds`. Then the
diagnostics block references `stan_data$include_c4` etc., which is undefined
(or worse, stale from a previous iteration).

**Why it matters:** Silent diagnostics corruption on resume, potentially a
crash depending on R's `exists()` semantics.

---

### P1-6. Precipitation prior is extremely restrictive

**Files/lines:**
- `4d_leaf_wax_spatial_model.stan:407` — `beta_precip ~ normal(0, 0.02)`
- `4a_spatial_functions.R:925` — precipitation standardization

**What's wrong:** A `normal(0, 0.02)` prior on a standardized predictor
effectively hard-codes a near-zero slope. Precipitation is the headline
covariate for a calibration model — pinning its slope to zero is a strong
choice that needs justification (or is a leftover bug from development).

**Fix direction:** Verify with the user whether this is intentional
(e.g., precipitation effect absorbed into the GP) or a leftover.

---

## P2 — WORTH KNOWING

### P2-1. `include_veg_interactions` is a narrow subset, not all pairwise

**File/line:** `4d_leaf_wax_spatial_model.stan:285`

**Details:** Includes only `OIPC×C4`, `OIPC×Tree`, `OIPC×Shrub`, `OIPC×Grass`.
No C4×PFT, no PFT×PFT, no elevation interactions.

**Impact:** Not a bug — a modeling choice. But it should be documented in
METHODS.md explicitly since "interactions" is ambiguous.

---

### P2-2. Spatial averaging is approximate, not geodesic

**Files/lines:**
- `3_prep_data.R:120` — Euclidean distance in degrees for extraction radius
- `4a_spatial_functions.R:23` — latitude-averaged deg-to-km conversion

**Impact:** Least trustworthy near the dateline and at high latitudes. For
global calibration this could matter at Arctic/Antarctic sites. Low priority
unless the manuscript cites this as a limitation.

---

### P2-3. `3g_optimize_knots.R` hardcodes `typical_n_obs <- 800`

**File/line:** `3g_optimize_knots.R:146`

**Impact:** Knot-scoring heuristic is stale for the 1131-site dataset. The
optimal knot count from this script may not match what's right for the new N.
Not a bug, just wrong input for the decision.

---

## Notes for downstream fixes

- Codex's verdict is gate FAIL based on 4 P0s. Recommend NOT launching until P0s are verified
  and fixed, and P1-2/P1-3 are resolved (blocking crashes), and P1-6 (precip prior)
  is confirmed intentional.
- P0-1 is the single biggest concern — corrupts 5 of 14 models silently. Must be
  verified by line-by-line reading of `4a_spatial_functions.R:59-100`,
  `4a_spatial_functions.R:742-780`, and `3_prep_data.R:200-290` before deciding
  on a fix.
- P0-2 also needs code verification — codex may have misread the Stan likelihood.
  Check `4d_leaf_wax_spatial_model.stan` data block, `transformed data`, and `model` block
  for how OIPC SE enters the response variance.
- METHODS.md is ~2 weeks stale and some discrepancies it documents may already
  be fixed in code. After P0/P1 triage, rewrite METHODS.md against the current
  snapshot at git SHA `004bf54`.

## Reproducing this review

```bash
# Pull snapshot
mkdir -p /tmp/leafwax_pipeline_for_review
cd /tmp/leafwax_pipeline_for_review
scp -i ~/Desktop/ai_help/aws/my-ec2-keypair.pem \
  ubuntu@<IP>:/home/shared/leafwax_spatial/spatial_leafwax_model/{0_load_config.R,1_extract_c4_raster.R,3_prep_data.R,3b_pre-flight_check.R,3d_analyze_collinearity.R,3g_optimize_knots.R,3h_prior_checks.R,4a_spatial_functions.R,4b_stan_prep.R,4c_fit_models.R,4d_leaf_wax_spatial_model.stan,config.yaml,launch_pipeline_4_parallel_w_shutdown.sh} \
  .
cp /Users/abradley/Desktop/proxy_uncertainty/leafwax_gca/METHODS.md .
git init -q && git add -A && git commit -q -m "snapshot"

# Invoke codex via /codex challenge with a custom prompt file, then:
codex exec "$(cat prompt.txt)" -C /tmp/leafwax_pipeline_for_review \
  -s read-only -c 'model_reasoning_effort="high"' \
  --enable web_search_cached --json 2>err.txt | json_parser.py
```
