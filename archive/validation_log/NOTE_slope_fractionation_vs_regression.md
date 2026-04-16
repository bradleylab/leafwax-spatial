# The Fractionation Factor Is Not the Regression Slope

**Date**: 2026-03-28
**Context**: Core conceptual note for manuscript discussion. A reviewer will ask: "if the hydrogen isotope fractionation is ~-150‰ (implying a slope of ~0.85), why does your spatial model give 0.55?"

---

The confusion comes from conflating two different quantities — the fractionation factor (a physical/biochemical constant) and the regression slope (an empirical statistical quantity). They happen to have similar values (~0.85) in the non-spatial case, which makes it tempting to equate them, but they measure different things.

## 1. The fractionation physics

The hydrogen isotope fractionation between source water and leaf wax is multiplicative:

```
(1 + d2H_wax/1000) = alpha × (1 + d2H_source/1000)
```

where alpha (the fractionation factor) ≈ 0.85, corresponding to a net fractionation epsilon = (alpha - 1) × 1000 ≈ -150‰. This follows the standard linearization (Hayes, 2001). Rearranging gives an exact linear relationship:

```
d2H_wax = alpha × d2H_source + (alpha - 1) × 1000
d2H_wax = 0.85 × d2H_source - 150
```

The slope of δ_product versus δ_reactant IS the fractionation factor alpha. This is not an approximation — it is an exact consequence of the definition of alpha. (The common shorthand δ_product ≈ δ_reactant + ε, which gives a slope of 1, IS an approximation that works for carbon isotopes where δ values are small but fails for hydrogen where δ values range over hundreds of permil.)

In a controlled experiment where you vary source water isotopes and measure wax isotopes from the same plant under the same conditions, you get a slope of ~0.85. This is the fractionation factor — a property of the biosynthetic pathway.

**The critical implication**: if the regression slope against source water is alpha, then a slope of 0.55 would imply alpha = 0.55, or epsilon = (0.55 - 1) × 1000 = **-450‰**. This is physically absurd — known net hydrogen isotope fractionation in leaf wax biosynthesis ranges from approximately -80‰ to -200‰ depending on plant type and conditions (Sachse et al., 2012). A reviewer will immediately notice this.

**But the spatial model's slope of 0.55 is not regressing against source water.** It is regressing against OIPC — a modeled estimate of precipitation δ²H, which is not the water the plant actually used. The regression slope only equals alpha when the x-variable is the actual reactant (source water). When it is a proxy for the reactant (OIPC), the slope reflects alpha convolved with everything that makes OIPC differ from source water. This is the essential point of this note.

## 2. The global regression slope is not the fractionation factor

The global regression:

```
d2H_wax_i = beta_0 + beta_oipc × OIPC_i + error_i
```

estimates the slope across 818 sites spanning the globe. This slope includes:

**(a) The fractionation factor (alpha ≈ 0.85)** — the direct biochemical relationship between source water and wax hydrogen isotopes.

**(b) Source water ≠ OIPC** — OIPC is a modeled estimate of mean annual precipitation δ²H. The actual water a plant uses differs from OIPC because of:
- Soil water evaporative enrichment (larger in arid environments)
- Seasonal bias (plants preferentially use growing-season precipitation, which differs isotopically from the annual mean)
- Soil water residence time and mixing (damping short-term isotopic variability)
- Interception and throughfall fractionation

These discrepancies between OIPC and actual source water vary geographically. They are spatially structured.

**(c) The fractionation factor itself varies geographically** — epsilon is not a constant. It depends on:
- Leaf temperature (Sachse et al., 2012)
- Relative humidity and transpiration rate (Kahmen et al., 2013)
- Plant functional type (C3 vs C4, trees vs grasses vs shrubs)
- Growth form and leaf morphology

All of these covary with OIPC across space. High-latitude sites have depleted OIPC AND different vegetation, temperature, and humidity than tropical sites.

**(d) Sedimentary integration** — lake surface sediment samples integrate lipids from an entire watershed, which may span an elevation range with heterogeneous precipitation isotopes, vegetation types, and microenvironments. The degree of integration varies by lake size, watershed topography, and vegetation density.

**The empirical slope of 0.85 in a non-spatial regression happens to approximately equal the fractionation factor, but this is partly coincidental.** It reflects the fractionation factor PLUS all the geographically-structured modifications listed above, which happen to roughly cancel in the global average. There is no reason they must equal 0.85 — they could sum to a higher or lower value depending on how the geographic modifications correlate with the OIPC gradient.

## 3. What the spatial model's slope means

The spatial model adds an intercept GP:

```
d2H_wax_i = [beta_0 + GP_intercept(s_i)] + beta_oipc × OIPC_i + error_i
```

The intercept GP absorbs variance that varies smoothly across space. This includes the geographically-structured components listed in section 2 — the source water modifications, the vegetation effects, the integration effects — insofar as they produce continental-scale spatial patterns.

The spatial slope beta_oipc estimates the **conditional** relationship: holding geographic position constant (i.e., after removing the spatial surface), how much does d2H_wax change per unit change in OIPC?

This is a fundamentally different question than the non-spatial regression asks. The non-spatial regression asks: "across all sites on Earth, what is the average relationship?" The spatial regression asks: "at a given location, if OIPC were different, what would d2H_wax be?"

The distinction is analogous to the **ecological fallacy** (Robinson 1950). A relationship observed across groups (sites) need not hold within groups (at a single location over time). The non-spatial slope is a "between-site" quantity. The spatial slope is an attempt to estimate a "within-site" quantity from cross-sectional data.

## 4. Can the biological slope be 0.85 while the field slope is lower?

Yes. Consider this concrete scenario:

Imagine two sites that happen to have the same OIPC value of -80‰, but at different elevations. Site A is a low-elevation tropical lake; Site B is a high-elevation temperate lake. Despite identical OIPC:
- Site A has stronger evapotranspiration → source water is enriched to -60‰ → d2H_wax = 0.85 × (-60) - 150 = -201‰
- Site B has less evaporation → source water ≈ precipitation → d2H_wax = 0.85 × (-80) - 150 = -218‰

These two sites have d2H_wax values that differ by 17‰ despite identical OIPC. The intercept GP captures this — it assigns different intercept values to the two locations, reflecting the different source water modifications.

Now consider the global gradient. At the global scale, OIPC ranges from ~+10‰ (tropical) to ~-200‰ (polar). But the OIPC-to-source-water offset also varies systematically: tropical sites have large evaporative enrichment (source water less negative than OIPC), while polar sites have minimal evaporation (source water ≈ OIPC). This means the effective source water range is compressed relative to the OIPC range — tropical source water is less enriched than OIPC suggests, and polar source water is about right.

If the source water range is smaller than the OIPC range, but d2H_wax faithfully tracks source water at 0.85:

```
Apparent slope = 0.85 × (source water range / OIPC range)
```

If source water range is, say, 85% of OIPC range due to tropical evaporative enrichment compressing the warm end, then:

```
Apparent slope against OIPC = 0.85 × 0.85 ≈ 0.72
```

This is a made-up number for illustration, but the point is real: **the regression slope against OIPC can be lower than the fractionation factor even without any statistical artifact, simply because OIPC is an imperfect proxy for source water, and the imperfection varies geographically.**

The spatial model detects this. The intercept GP absorbs the geographically-varying OIPC-to-source-water offset. What remains (beta_oipc ≈ 0.55) reflects the relationship after accounting for these offsets.

Whether 0.55 is the "right" number depends on how well the GP separates the genuine source water offset from the direct OIPC-wax relationship (the Paciorek overcorrection concern), but the direction is correct: the field slope should be lower than 0.85 because of source water modifications.

## 5. The three slopes

It may be useful to distinguish three quantities:

| Slope | Value | What it represents |
|---|---|---|
| **Fractionation factor** (alpha) | ~0.85 | Biochemical: d2H_wax per unit d2H_source, single plant, controlled conditions |
| **Apparent slope** (non-spatial) | ~0.83 | Empirical global: d2H_wax per unit OIPC, all geographic effects included |
| **Conditional slope** (spatial) | ~0.55 | Empirical local: d2H_wax per unit OIPC, after removing geographic patterns |

The near-equality of the fractionation factor (0.85) and the apparent slope (0.83) is not a validation of either — it reflects the coincidental near-cancellation of opposing geographic effects in the global average.

The conditional slope (0.55) is arguably the most relevant for paleoclimate application, because when inverting a downcore record at a specific site, you want the local relationship, not the global average. But it may be an underestimate of the true conditional relationship (Paciorek overcorrection, addressed in NOTE_spatial_confounding_theory.md).

## 6. What this means for the manuscript

A reviewer will reason: "the fractionation is -150‰, so alpha = 0.85 and the slope should be 0.85. Your spatial model gives 0.55, which implies alpha = 0.55, i.e., epsilon = -450‰. That's physically impossible. Something is wrong with your model."

The response:

1. The regression slope only equals alpha when regressing δ_product against δ_reactant (i.e., d2H_wax against d2H_source_water). Our model regresses d2H_wax against OIPC — a modeled estimate of precipitation δ²H, not the actual source water. The slope of d2H_wax vs. OIPC is not alpha and should not be interpreted as a fractionation factor.

2. The apparent non-spatial slope of 0.83 coincidentally approximates alpha ≈ 0.85, but this is not evidence that the regression is recovering the fractionation factor. It reflects alpha convolved with geographically-structured differences between OIPC and source water, which happen to roughly cancel in the global average.

3. The spatial model identifies that a large fraction of the OIPC-d2Hwax covariation is geographically structured (continental-scale patterns with length scale ~4,200 km). These patterns reflect the systematic, spatially-varying offset between OIPC and actual source water (evapotranspiration, seasonal bias, vegetation effects).

4. The conditional slope of 0.55 is the residual d2Hwax-OIPC relationship after the GP absorbs the geographically-structured OIPC-to-source-water offset. It does NOT imply a fractionation of -450‰. It implies that the effective slope of d2Hwax against OIPC (not against source water) is 0.55 after geographic deconfounding. The fractionation factor alpha remains ~0.85; the attenuation occurs in the OIPC-to-source-water link, not in the source-water-to-wax link.

5. For paleoclimate inversion, the model propagates the full posterior uncertainty of the slope through the inversion, producing wider but more honestly calibrated prediction intervals than conventional calibrations using the apparent slope.

## 7. What we should NOT claim

We should not claim that the "true" fractionation is -450‰ or that the slope of 0.55 is the fractionation factor. It is not. The slope of 0.55 is the conditional relationship between d2H_wax and OIPC after geographic deconfounding. It reflects the fractionation factor (alpha ≈ 0.85) attenuated by the OIPC-to-source-water offset and potentially by Paciorek-type overcorrection. These components cannot be fully separated with our data.

We should also not claim that the non-spatial slope of 0.83 is "wrong" — it is a valid estimate of the marginal (across-site) relationship, which is appropriate for applications that leverage the same global gradient (e.g., spatial reconstruction of isotope maps at a single time slice).

The contribution is identifying that these are different quantities, quantifying the difference, and providing calibrated uncertainty for both through the model framework.

## 8. References for this argument

- Hayes JM (2001). Fractionation of carbon and hydrogen isotopes in biosynthetic processes. *Rev Mineral Geochem* 43(1):225-77. [linearized fractionation equations, slope = alpha]
- Robinson WS (1950). Ecological correlations and the behavior of individuals. *Am Sociol Rev* 15(3):351-7. [ecological fallacy]
- Sachse D et al. (2012). Molecular paleohydrology: interpreting the hydrogen-isotopic composition of lipid biomarkers from photosynthesizing organisms. *Annu Rev Earth Planet Sci* 40:221-49. [variable fractionation]
- Kahmen A et al. (2013). Leaf water deuterium enrichment shapes leaf wax n-alkane δD values of angiosperm plants I: Experimental evidence and mechanistic insights. *GCA* 111:39-49. [humidity effects on fractionation]
- Paciorek CJ (2010). The importance of scale for spatial-confounding bias and precision of spatial regression estimators. *Biostatistics* 11(4):601-15. [overcorrection theory]
- Polissar PJ, D'Andrea WJ (2014). Uncertainty in paleohydrologic reconstructions from molecular δD values. *GCA* 129:146-56. [inversion uncertainty]
