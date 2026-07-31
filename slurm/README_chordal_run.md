# Chordal refit run — orchestration (spec v2 §3)

End-to-end procedure for the chordal-metric refit batch on Compute2. All fits run
on C2 (no local fits). The frozen run `results/c2_run_20260626` is preserved; the
chordal run becomes the new authoritative run only after re-trace.

## Batch composition (17 fits — ALL models refit fresh)

Revised 2026-07-18: refit **all 17 model_configs**; nothing is reused from the
frozen run. Why reuse was abandoned: **baseline_env's** frozen chain CSVs are
header-only (0 draw rows), so its `beta_elev` is unrecoverable and it MUST be
refit. (full/full_interact frozen CSVs *do* retain draws + `beta_elev`, but the
frozen `posterior_draws.rds` dropped `beta_elev` at save for all three, and a
uniform refit of all 17 is cleaner to audit than a hybrid.) Non-spatial models
run fast, so refitting all is simplest and self-consistent.

- **9 chordal spatial refits:** baseline_sp, baseline_veg_sp, baseline_env_sp,
  c4_only_sp, elevation_only_sp, elevation_c4_sp, elevation_c4_interact_sp,
  full_sp, full_interact_sp.
- **3 range_factor-off sensitivities:** baseline_sp_rfoff, baseline_env_sp_rfoff,
  full_interact_sp_rfoff (config `apply_range_factor: false`).
- **5 non-spatial refits:** baseline, baseline_veg, baseline_env, full,
  full_interact. No GP, so the chordal geometry is inert for them; they *should*
  target the same posterior as the frozen run, but that is verified post-fit by
  `scripts/verify_nonspatial_vs_frozen.R` (shared-parameter equivalence), not
  assumed. The 4c save-fix retains `beta_elev`, so baseline_env/full/full_interact
  carry elevation.

The frozen run `results/c2_run_20260626` is preserved and recorded for
COMPARISON only. Whole-file md5 cannot match (the refit elevation deposits carry
extra `beta_elev` columns); use `verify_nonspatial_vs_frozen.R` for the real
shared-parameter comparison.

## 0. C2 preflight (2026-07-15 upgrade — done 2026-07-18)

- Host keys refreshed on **pliny** — `ssh-keygen -R` takes ONE host per command
  (`ssh-keygen -R c2-login-001.ris.wustl.edu`, then `-002`, then `-003`; a brace
  expansion errors "Too many arguments"). ControlMaster re-established (password + 2FA).
- Account `-A compute2-alexander.s.bradley` (unchanged, verified).
- `apptainer/1.3.6` (SLURM scripts bumped from 1.3.4).
- See `~/.claude/rules/research-infrastructure.md` + vault
  `Dev_tools/Compute2 (WashU RIS) — Usage and July 2026 Upgrade.md`.

## 1. Pre-launch gate (before the batch — must pass)

Local static checks (already green):
- `Rscript scripts/test_chordal_acceptance.R` — chordal geometry, τ logic, metric
  dispatch (24/24).
- `Rscript scripts/validate_chordal_prep.R` — full `prepare_stan_data` path on the
  REAL frozen sediment (chordal coords, km reg, ls bounds, range_factor on/off).
- `stanc` syntax + R parse on the edited files.

On C2 (real data, converged — NOT a few-iter smoke):
1. **Prep** into a fresh dir (§2), then **pilot-fit baseline_sp only**:
   ```bash
   sbatch --array=0 slurm/job_fit_chordal.sh
   ```
2. **Convergence + length-scale bound pile-up gate** (pull the output back or run
   on C2):
   ```bash
   Rscript scripts/check_pilot.R model_output/baseline_sp \
     prepared_data/stan_data_baseline_sp.rds
   ```
   Gates: R-hat < 1.01, ESS > 400, 0 divergences, 0 max-treedepth hits, no low
   E-BFMI chains; the length-scale column present/finite, its two reported columns
   (ls_intercept_km == ls_slope_km) in agreement, and the posterior not piling at
   the length-scale bounds. **Bounds are read from the prepared stan_data**
   (`ls_log_lower`/`ls_log_upper`), not hardcoded — the 2nd arg is optional and
   inferred from the model dir name, but pass it explicitly if the prepared data
   is not under `prepared_data/`. If the posterior piles at a bound, widen
   `gp_length_scale` bounds in config.yaml (operator decision), re-prep, and
   re-run the pilot.
3. **Stan↔R prediction parity** (the package must reproduce the Stan field).
   MUST set `OPENBLAS_CORETYPE=Haswell` — see the OpenBLAS gotcha below; without
   it the parity gate silently reports a huge (~1e3) discrepancy that is a broken
   BLAS, not a real mismatch. Run inside the SIF with the package's R source bound
   in (the parity script sources `../leafwax-pkg/R/spatial_interpolation.R`, so
   bind the PARENT dir, not just the run dir):
   ```bash
   apptainer exec --no-home --containall \
     --bind /scratch2/fs1/$USER:/scratch2/fs1/$USER --bind /tmp:/tmp \
     --pwd /scratch2/fs1/$USER/leafwax_run \
     --env OPENBLAS_CORETYPE=Haswell \
     /scratch2/fs1/$USER/leafwax_run/leafwax-spatial.sif \
     Rscript scripts/check_stan_r_parity.R \
       model_output/baseline_sp prepared_data/stan_data_baseline_sp.rds
   ```
   PASS looks like `max|pkg - stan|` ~1e-6 on both fields.
Only when all three pass, launch the full batch.

### OpenBLAS coretype gotcha (in-SIF R linear algebra on C2)

The C2 general-cpu nodes are Intel Sapphire Rapids (cpuinfo model 143). The SIF's
R links OpenBLAS 0.3.20 (`openblas-pthread`), which PREDATES Sapphire Rapids and
mis-selects a buggy AVX-512 `dpotrf`/`dtrsm` kernel: `chol()`/`solve()` return
WRONG (and non-deterministic) results for larger matrices — e.g. `chol()` fails
"leading minor of order 33 not positive" on a matrix whose min eigenvalue is
clearly positive. This corrupts any R Gaussian-process prediction that inverts the
knot covariance (the leafwax package's `predict_one_gp_mpp`, hence the parity gate
and the re-trace prediction maps/inversions). **Stan itself is unaffected** — it
uses its own compiled Eigen, not R's OpenBLAS — so the fits are correct.

Fix: set `OPENBLAS_CORETYPE=Haswell` (correct AVX2 path) in the apptainer `--env`
for ANY in-SIF R step that does GP prediction. Must be in the environment BEFORE R
starts (OpenBLAS reads it at load time); `Sys.setenv()` inside R is too late.
This applies to the parity gate above AND to §5 re-trace when run in the SIF on C2
— a re-trace prediction number computed without it is silently wrong (data
integrity). Re-trace on Mac (Accelerate BLAS) is not affected.

### Deploy gotcha — /scratch2 purges by mtime; `touch` after rsync

C2 `/scratch2` reaps files whose **mtime** is older than ~30 days. `rsync -a`
PRESERVES the *source* mtime, so a file last edited months ago lands on C2 with
its old mtime and can be purged **mid-run** even though you just deployed it.
This bit run 2283718 (2026-07-20): `0_load_config.R` (mtime Apr 9) and
`1_extract_c4_raster.R` (Apr 12) vanished after the first 5 fits, and every later
fit died at `source("0_load_config.R")`. The SIF + frozen sediment RDS (Jun 23)
were next in line. After deploying, reset mtimes so nothing expires during the
run — but do NOT touch `prepared_data/` or `model_output/` (their mtimes drive
the resume/freshness guards):
```bash
cd /scratch2/fs1/$USER/leafwax_run
find . -type d \( -name prepared_data -o -name model_output \) -prune \
  -o -type f -exec touch {} +
```

## 2. Fresh prepared data (no stale cache)

The chordal batch needs a fresh `prepared_data/` so no stale (2D-standardized)
stan_data is reused. Two safeguards:
- `4b_stan_prep.R` now invalidates the cache on a config/prep-code signature
  change (not just mtime) — a chordal config change forces a rebuild.
- Still, prep into a clean `prepared_data/` for the run (empty it or point
  `CONFIG$output_dirs$prepared_data` at a new location). Run prep via
  `sbatch slurm/job_prep.sh` (preps all 17 model_configs; the fit array uses all 17).

## 3. Launch the fits

```bash
sbatch slurm/job_fit_chordal.sh   # array 0-16, 17 models, apptainer/1.3.6
```

Each array task fits one model and writes per-model artifacts under
`model_output/<model>/` (`fit_status.rds` always; a `DONE` marker only when the
full core set — fit.rds + diagnostics.rds + posterior_draws.rds — is on disk).
The wrapper's resume guard skips a task ONLY when that complete set is present
and newer than its stan_data, so a crash between fitting and draw-extraction
re-runs (4c resumes from fit.rds and regenerates the downstream artifacts) rather
than false-succeeding. A failed/incomplete fit exits nonzero, so SLURM marks the
task failed.

**Batch completion gate (fail-closed):**
```bash
Rscript scripts/aggregate_chordal_status.R model_output prepared_data
```
This reads the per-model status/DONE files and exits nonzero unless all 17 models
finished with a full artifact set AND every artifact is newer than its stan_data
(freshness guard against a stale DONE surviving a crashed resume). Do not assemble
the run until it passes. (There is no shared `pipeline_4c_complete.rds` for the
array — under concurrent single-model tasks it would be meaningless; the
aggregator is the batch judge.)

> **Do not use `slurm/submit.sh` for this run** — it launches the deprecated
> standardized-metric 14-model workflow and now refuses to run without
> `--legacy`. The chordal path is job_prep.sh → job_fit_chordal.sh → aggregator.

## 4. Snapshot + build the manifest

All 17 deposits come from the fit array — there is no separate assembly/copy
step (nothing is reused).

- Frozen numbers already snapshotted: `NUMBERS_c2_run_20260626.md`.
- Pull `model_output/` (and `prepared_data/`) back and snapshot to
  `results/c2_run_<date>_chordal/` (preserve `results/c2_run_20260626/` untouched).
- **Fit-time provenance** is captured during prep (`job_prep.sh` runs
  `scripts/capture_fit_provenance.R` inside the SIF → `results/fit_provenance.rds`,
  recording the SIF's real CmdStan/package versions, SIF md5, and code/input
  hashes). Snapshot it into the run root alongside `model_output/`; the manifest
  REQUIRES it.
- Build the manifest (fail-closed; hard-fails without a readable config of exactly
  17 models or the fit-provenance record; validates draw dims, the full conditional
  variable set, indexed dims, elevation content, and fit_status/DONE per model; add
  `--allow-incomplete` only for an interim snapshot):
  ```bash
  Rscript scripts/build_run_manifest.R \
    --chordal-output   results/c2_run_<date>_chordal/model_output \
    --chordal-prepared results/c2_run_<date>_chordal/_prepared_data \
    --frozen-output    results/c2_run_20260626/model_output \
    --fit-provenance   results/c2_run_<date>_chordal/fit_provenance.rds \
    --out              results/c2_run_<date>_chordal/RUN_MANIFEST_chordal
  ```
- Verify the non-spatial refits reproduce the frozen posterior on shared
  parameters (not assumed):
  ```bash
  Rscript scripts/verify_nonspatial_vs_frozen.R \
    --refit  results/c2_run_<date>_chordal/model_output \
    --frozen results/c2_run_20260626/model_output
  ```

## 5. Re-trace (spec §4) then deposits

Recompute every affected number/figure/table from the new run (slopes, GP basis
length scale, ~110×, detection thresholds via `scripts/regen_precip_space.R`,
LOOIC/Table S10, spatial figures), apply the range_factor gate to the 3
sensitivities, then build new versioned Zenodo deposits (stamp
`attr(posterior, "spatial_metric") <- "chordal"`).

## `range_factor` gate ("material" defined up front)

Escalate to more sensitivity fits only if a range_factor-off fit changes a
**manuscript-level rounded claim** or the **intercept/slope partition
interpretation**; otherwise the 3 suffice.
