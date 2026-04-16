# Spatial Confounding and the OIPC Slope: Theory, Implications, and Honest Uncertainty

**Date**: 2026-03-28
**Context**: Reference note for manuscript writing. Not for direct inclusion — draw from as needed.

---

## 1. The identifiability problem

Precipitation δ²H (OIPC) has strong spatial structure. It tracks latitude, altitude, continentality, and moisture source — all fundamentally geographic quantities. Leaf wax δ²H inherits this spatial structure through the fractionation chain: precipitation → soil water → plant uptake → leaf wax synthesis. Both variables covary across space for reasons that are partly causal (precipitation isotopes do influence wax isotopes) and partly circumstantial (high-latitude sites have both depleted precipitation and depleted wax isotopes because of the temperature-dependent fractionation environment, atmospheric distillation, and ecosystem composition).

A non-spatial regression attributes the full OIPC-d2Hwax correlation to the slope coefficient. A spatial regression (with a GP intercept) allows some of that correlation to be explained by "geographic position" rather than "OIPC per se." The slope coefficient drops because the model has an alternative explanation for part of the covariation.

## 2. Paciorek's formalization

Paciorek (2010, Biostatistics 11(4):601-15) showed that when a covariate X shares spatial structure with an unmeasured confounder Z, a spatial model systematically attenuates the estimated coefficient of X. The approximate bias is:

```
E(beta_hat) ≈ beta_true - rho × (sigma_Z / sigma_X) × S
```

where:
- rho is the correlation between X and Z at the relevant spatial scale
- sigma_Z / sigma_X is the relative variance of the confounder to the covariate
- S is a spatial smoothing factor that depends on the GP kernel and length scale

The attenuation increases when:
- The covariate is more spatially smooth (its variation is at scales the GP can absorb)
- The confounder is stronger (larger sigma_Z)
- The covariate and confounder share spatial scales (higher rho at the same frequencies)

This is not a pathology of the model. It reflects a genuine identifiability limitation: if X and Z vary at the same spatial scales, no amount of data can perfectly separate their effects. The spatial model is being honest about an ambiguity that the non-spatial model ignores.

Hodges and Reich (2010, The American Statistician 64(4):325-34) made this point bluntly: adding spatially-correlated errors can "mess up" the fixed effect you care about, not because the model is broken, but because the data genuinely cannot distinguish the fixed effect from the spatial effect at the scales where both operate.

## 3. What this means for our OIPC slope

In our models:
- Non-spatial estimate: beta_oipc ≈ 0.83 (includes both process signal and geographic covariation)
- Spatial estimate: beta_oipc ≈ 0.55 (after the intercept GP absorbs geographically-structured variance)

The ~0.28 drop reflects the portion of the OIPC-d2Hwax correlation that can be explained by geographic position alone (the intercept GP, which has a length scale of ~4,200 km — continental scale). This does not mean the "true" process slope is 0.55. It means:

**The true process slope is somewhere between 0.55 and 0.83.** The non-spatial estimate is an upper bound (it assumes no geographic confounding). The spatial estimate may be a lower bound (the GP may absorb some genuine OIPC signal along with the confounding). The posterior distribution of beta_oipc in the spatial model (0.55, 95% CI [0.42, 0.69]) represents the model's honest uncertainty about where the truth lies, given the identifiability constraint.

The staged attenuation we observed is consistent with this framework:
- GP alone drops the slope to ~0.60 (absorbs the broadest geographic covariation)
- Adding elevation drops it further to ~0.50 (elevation correlates with OIPC at r = -0.38, providing the GP with additional spatially-structured variance to absorb)
- Adding all covariates drops it to ~0.45 (each additional spatially-structured covariate gives the model another avenue to explain OIPC-correlated variance without the slope)

## 4. The two interpretations

There are two ways to read the 0.83 → 0.55 attenuation, and both may be partly correct:

**Interpretation A (confounding correction):** The non-spatial slope of 0.83 is inflated by geographic confounding. Sites at high latitudes have both depleted OIPC and depleted d2Hwax for reasons beyond the direct OIPC→wax pathway (e.g., temperature effects on biosynthesis, ecosystem composition, growing season length). The spatial model correctly identifies and removes this confound. The "true" process slope is closer to 0.55.

**Interpretation B (overcorrection):** The spatial model is too aggressive. The OIPC spatial gradient is largely *causal* — precipitation isotopes genuinely drive wax isotopes at the global scale through the fractionation chain. By allowing the intercept GP to absorb continental-scale patterns, the model removes signal along with confounding. The "true" process slope is closer to 0.83, and the attenuation represents lost signal, not removed confounding.

The reality is almost certainly a mixture. Dupont, Wood, and Augustin (2022, Biometrics 78(4):1279-90) and Zimmerman and Ver Hoef (2022, The American Statistician) have shown that the degree of overcorrection depends on the spectral overlap between the covariate and the spatial effect — confounding at low spatial frequencies (smooth, continental-scale patterns) biases the spatial model more than confounding at high frequencies (local patterns). OIPC is predominantly a low-frequency field (it varies smoothly over continental scales), which is precisely the regime where Paciorek's attenuation is strongest.

## 5. Implications for paleoclimate inversion

When inverting the model to reconstruct d2H_precip from downcore d2H_wax:

```
d2H_precip = (d2H_wax - intercept(s) - other_terms) / beta_oipc(s)
```

The slope appears in the denominator, so an attenuated slope amplifies the reconstruction. If the "true" slope is 0.70 but the model estimates 0.55, reconstructed d2H_precip changes are amplified by 0.70/0.55 = 27%.

However, the uncertainty propagation corrects for this naturally. The posterior of beta_oipc (0.55 ± 0.07) is wider than the non-spatial estimate (0.83 ± 0.02). When propagated through the inversion, this produces wider but more calibrated prediction intervals for d2H_precip.

The practical consequence: **a spatial-model-based reconstruction will show smaller central estimates of isotopic change but wider uncertainty bands compared to a conventional calibration.** Both features reflect the same underlying reality — that the modern calibration slope is confounded by geography, and we are honest about not knowing how much.

For a user of the R package who inputs a d2H_wax value at a specific location, the model returns:
1. A point estimate of d2H_precip (using the posterior mean of beta_oipc at that location)
2. A full posterior predictive distribution that propagates uncertainty from: the slope, the intercept GP, measurement error, model residual error, and (for new locations) GP prediction uncertainty

The slope uncertainty is the dominant source of inversion uncertainty, and it directly reflects the spatial confounding limitation.

## 6. What the confounding simulation will tell us

The Paciorek-framework simulation (running overnight) tests the model's behavior at four confounding levels. The expected outcomes:

- **rho = 0 (no confounding):** Model should recover beta_oipc = 0.7 exactly (like scenario 3a, which passed).
- **rho = 0.3 (weak):** Small attenuation expected. If the model still recovers ~0.7, it successfully separates weak confounding.
- **rho = 0.5 (moderate):** Some attenuation likely. The recovered slope may be 0.6-0.65. This quantifies the Paciorek bias at a standard confounding level.
- **rho = 0.45 (empirical):** This is the key result. The attenuation here directly estimates how much of the real data's 0.83 → 0.55 drop is "overcorrection" vs. genuine confounding removal.

If the model recovers beta_oipc = 0.65 at rho = 0.45 (true = 0.7), that implies a bias of -0.05. Applied to the real data: the "true" slope might be ~0.60 rather than the estimated 0.55. The 0.83 → 0.60 drop would then be genuine confounding removal, and the 0.60 → 0.55 drop would be Paciorek-type overcorrection.

## 7. Key references

- Paciorek CS (2010). The importance of scale for spatial-confounding bias and precision of spatial regression estimators. *Biostatistics* 11(4):601-15.
- Hodges JS, Reich BJ (2010). Adding spatially-correlated errors can mess up the fixed effect you care about. *The American Statistician* 64(4):325-34.
- Dupont E, Wood SN, Augustin NH (2022). Spatial+: a novel approach to spatial confounding. *Biometrics* 78(4):1279-90.
- Dupont E, Wood SN, Augustin NH (2023). Demystifying spatial confounding. arXiv:2309.16861.
- Zimmerman DL, Ver Hoef JM (2022). On deconfounding spatial confounding in linear models. *The American Statistician* 76(2):159-67.
- Khan K, Calder CA (2022). Restricted spatial regression methods: implications for inference. *J Am Stat Assoc* 117(537):482-94.
- Hanks EM, Schliep EM, Hooten MB, Hoeting JA (2015). Restricted spatial regression in practice: geostatistical models, confounding, and robustness under model misspecification. *Environmetrics* 26(4):243-54.
