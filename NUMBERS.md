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
- **δ²H_precip scale** (÷ local slope ≈0.545): single-sample SD **≈29‰**; two-sample ρ=0 threshold **≈81‰**. ⚠️ Confirm exact precip-scale values via the `leafwax` inversion re-run before quoting (currently derived via ÷slope from wax-scale).
- Fig 5 endpoints: spatial ρ=0.9 → **14.0‰**; non-spatial ρ=0 → **58.7‰** (σ_baseline 20.96‰).

## Confounding simulation (Fig 3, confounding_v2, 2026-07-02) — verified

| scenario | achieved ρ | true.std | GP post.mean | OLS | % bias absorbed |
|---|---:|---:|---:|---:|---:|
| rho00 | 0.000 | 0.474 | 0.511 | 0.734 | 85.8 |
| rho03 | 0.639 | 0.414 | 0.613 | 0.798 | 48.1 |
| empirical | 0.726 | 0.394 | 0.664 | 0.826 | 37.5 |
| rho05 | 0.752 | 0.389 | 0.682 | 0.835 | 34.2 |

Empirical headline (main text): true **0.39**, OLS **0.83**, GP **0.66**, ~40% removed / ~60% remains.

## Simulated recovery (simrecovery, 2026-07-02) — verified

- 3a uniform 0.70 → **0.68 [0.61, 0.75]** (CI covers truth).
- 3b spatially varying mean 0.65 → **0.67 [0.60, 0.75]** (CI covers truth).

## OPEN — not yet locked

1. **Prior-sensitivity (6c, SLURM 1930423)** re-running on frozen data (7 variants of baseline_veg_sp). GCA quoted posterior medians 0.614–0.638. Fold in when complete.
2. **Precip-scale inversion SD (~29‰) + two-sample threshold (~81‰)**: confirm from a `leafwax`-package inversion on the frozen fit, rather than the ÷slope derivation above.
3. **Prior tables backed up** at `results/c2_run_20260626/tables_prev_20260703/` (pre-frozen for reference).
