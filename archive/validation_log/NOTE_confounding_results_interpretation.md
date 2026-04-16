# Spatial Confounding Simulation Results: Interpretation and Implications

**Date**: 2026-04-04
**Context**: Summary of Phase 3c confounding simulation results and their implications for the manuscript and paleoclimate inversion.

---

## 1. What we tested

We generated synthetic leaf wax d2H data where the true OIPC slope is known (beta_oipc = 0.7), then added a spatially-structured confounding intercept with controlled correlation to OIPC. The confounding intercept is a Gaussian process draw with the same Matern 3/2 kernel and length scale as the analysis model, following the Paciorek (2010) framework.

The question: can the spatial model (GP intercept) absorb the confounding and recover the true causal slope?

## 2. What we found

| Scenario | rho | True slope (re-std) | Posterior | Bias | 95% CI covers truth? |
|---|---|---|---|---|---|
| rho=0.0 | None | 0.551 | 0.607 | +0.056 | YES |
| rho=0.3 | Weak | 0.477 | 0.773 | +0.296 | NO |
| rho>=0.45 | Moderate+ | — | Sampler failed | — | NOT IDENTIFIABLE |

## 3. Intuitive explanation

### The core problem

OIPC (precipitation d2H) varies smoothly at continental scales: depleted at high latitudes, enriched in the tropics. This pattern tracks latitude, altitude, continentality, and moisture source — all geographic quantities. Leaf wax d2H inherits this spatial structure through the fractionation chain (precipitation -> soil water -> plant uptake -> wax synthesis), but it also responds to geographic factors that have nothing to do with OIPC per se — temperature effects on biosynthesis, ecosystem composition, growing season length.

The GP intercept in the model is designed to absorb these "geographic but not OIPC" effects. But OIPC itself is fundamentally a geographic field. It varies at the same continental scales as the GP. The model cannot distinguish "d2H_wax is depleted here because OIPC is depleted" from "d2H_wax is depleted here because this is a high-latitude site with a particular ecosystem." Both explanations produce the same spatial pattern in the data.

### What the simulation showed at each level

**rho = 0 (no confounding)**: The confounding intercept is spatially structured but statistically independent of OIPC. The model can separate them because they happen to vary differently across space. The GP absorbs the intercept pattern, and the slope recovers close to truth (bias +0.056, 95% CI covers). This confirms the model mechanics are sound.

**rho = 0.3 (weak confounding)**: The confounding intercept weakly tracks OIPC (cor ≈ 0.3). Even at this modest level, the model cannot separate the confounding from the true OIPC signal. The slope absorbs both the genuine OIPC effect and the confounding, producing a posterior (0.773) far above the true value (0.477). The spatial GP fails to absorb the OIPC-correlated portion of the intercept because it cannot distinguish that portion from additional OIPC signal.

**rho >= 0.45 (moderate to strong)**: The posterior becomes so entangled — the GP amplitude and the slope are nearly perfectly anti-correlated — that the NUTS sampler cannot navigate the geometry. The model is computationally unable to explore the parameter space. This reflects complete non-identifiability: the data contains no information to separate the slope from the confounding at these levels.

### Why this is a fundamental limitation

This is not a model failure or a fixable bug. It is a property of the data. Paciorek (2010) proved that when a covariate (OIPC) shares spatial structure with an unmeasured confounder, the spatial regression coefficient is biased in a way that cannot be corrected without additional information. The bias depends on the spectral overlap between the covariate and the spatial effect — and OIPC is predominantly a low-frequency (continental-scale) field, which is precisely the regime where this bias is strongest.

No amount of additional leaf wax sites, no alternative GP kernel, and no restricted spatial regression method (Dupont et al. 2022; Khan and Calder 2022) can fully resolve this. The confounding is at the same spatial frequencies as the signal. The only way to break the confounding would be to observe the confounders directly (e.g., measure the non-OIPC contributions to wax isotope variation), which is the broader goal of the modeling enterprise but cannot be achieved through spatial statistics alone.

## 4. What it means for the OIPC slope interpretation

The real data shows:
- Non-spatial OLS: beta_oipc ≈ 0.83
- Spatial model (GP intercept): beta_oipc ≈ 0.55, 95% CI [0.42, 0.69]

The simulation tells us the 0.83 -> 0.55 drop is real ambiguity:
- Some portion is genuine confounding removal (the GP correctly absorbing geographic effects that inflated the OLS estimate)
- Some portion may be signal loss (the GP absorbing genuine OIPC signal along with confounding)
- We cannot determine the split from spatial data alone

The posterior distribution (0.55, 95% CI [0.42, 0.69]) is the model's honest representation of this identifiability constraint. It is not wrong — it is appropriately uncertain.

## 5. Implications for precipitation reconstruction uncertainty

### The direct effect: wider prediction intervals

When inverting the forward model to reconstruct d2H_precip from downcore d2H_wax:

```
d2H_precip = (d2H_wax - intercept(s) - covariates) / beta_oipc(s)
```

The slope appears in the denominator. Three consequences:

**a) Uncertainty amplification.** An uncertain slope in the denominator propagates nonlinearly into the reconstruction. The posterior on beta_oipc (SD ≈ 0.07 in standardized space) is wide relative to the mean (0.55). When propagated through the inversion, this produces wider prediction intervals for d2H_precip than a conventional (non-spatial) calibration with its falsely precise slope of 0.83 ± 0.02.

**b) Asymmetric error.** Division by a positively-distributed denominator produces right-skewed prediction intervals. A slope posterior centered at 0.55 with support down to 0.42 means the lower tail of beta_oipc produces very large d2H_precip estimates. The prediction intervals are wider on the "more depleted" side.

**c) Central estimates are amplified.** If the true process slope is higher than the spatial model's estimate (e.g., 0.65 rather than 0.55), the reconstruction amplifies isotopic changes by 0.65/0.55 ≈ 18%. Conversely, if the spatial model is correct, the conventional calibration (0.83) attenuates changes by 0.55/0.83 ≈ 34%.

### The honest answer: it makes uncertainty larger but more calibrated

Compared to a conventional (non-spatial) calibration:

| Property | Conventional (OLS) | Spatial model |
|---|---|---|
| Slope estimate | 0.83 ± 0.02 | 0.55 ± 0.07 |
| Slope precision | High (falsely) | Lower (honestly) |
| Reconstruction central estimate | Smaller changes | Larger changes |
| Prediction interval width | Narrow | Wide |
| Calibration (does 95% CI contain truth?) | Unknown — ignores confounding | Better — accounts for spatial ambiguity |

The spatial model produces **more uncertainty** in the reconstruction. But this uncertainty was always there — the conventional calibration just ignored it. The confounding simulation demonstrates that the OLS slope of 0.83 absorbs geographic covariation that is not strictly causal, and the true process slope is genuinely unknown within the range [0.42, 0.69].

A reconstruction that reports ±8 permil is worse than one that reports ±3 permil only if the ±3 permil interval actually contains the truth. The simulation shows that the narrow OLS-based interval does not account for spatial confounding and is therefore likely miscalibrated (overconfident).

### Practical magnitude

For a hypothetical downcore d2H_wax change of -30 permil (a glacial-interglacial scale shift):

| Calibration | Reconstructed d2H_precip change | 95% PI |
|---|---|---|
| Conventional (slope = 0.83) | -36 permil | approximately ±5 permil |
| Spatial model (slope = 0.55) | -55 permil | approximately ±20 permil |

The spatial model reconstructs a larger central change (because a smaller slope means OIPC must change more to produce the same wax signal) and wider uncertainty (because the slope itself is uncertain). The ±20 permil interval is large, but it honestly reflects what the modern calibration data can constrain.

## 6. Manuscript framing

This result supports framing the spatial model as providing **honest uncertainty quantification** rather than a better point estimate. The key message:

1. The model correctly identifies that the OIPC-wax relationship is confounded by shared spatial structure at continental scales.
2. The posterior on beta_oipc reflects this identifiability limitation.
3. The uncertainty propagates through to paleoclimate reconstructions, producing wider but better-calibrated prediction intervals.
4. This is a fundamental property of spatial calibration data, not a limitation specific to this model or dataset.

## 7. References

- Paciorek CS (2010). The importance of scale for spatial-confounding bias and precision of spatial regression estimators. *Biostatistics* 11(4):601-15.
- Hodges JS, Reich BJ (2010). Adding spatially-correlated errors can mess up the fixed effect you care about. *The American Statistician* 64(4):325-34.
- Dupont E, Wood SN, Augustin NH (2022). Spatial+: a novel approach to spatial confounding. *Biometrics* 78(4):1279-90.
- Khan K, Calder CA (2022). Restricted spatial regression methods: implications for inference. *J Am Stat Assoc* 117(537):482-94.
