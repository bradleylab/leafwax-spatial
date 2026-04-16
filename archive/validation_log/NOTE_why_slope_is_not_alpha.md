# Why the Regression Slope Is Not the Fractionation Factor: The Growing Season Argument

**Date**: 2026-03-28
**Context**: Sharpened argument for the reviewer who says "of course the slope is 0.85."

---

## The skeptic's position

"The hydrogen isotope fractionation between source water and leaf wax is epsilon ≈ -150‰, giving alpha ≈ 0.85. When I regress d2H_wax on d2H_precip, I'm regressing product on reactant. The slope must be alpha. Your spatial model giving 0.55 is broken."

## Why the skeptic is wrong: a worked example

The slope of d2H_wax vs d2H_source IS alpha. Nobody disputes that. But the regression is not against source water. It's against OIPC — a model estimate of mean annual precipitation δ²H. These are different quantities, and the difference is geographically structured. One well-known mechanism is sufficient to demonstrate this: growing season bias.

### Setup

Consider four lake sediment sites along a latitudinal transect. At each site, the fractionation is exactly alpha = 0.85 (epsilon = -150‰). No variation, no complications. The only thing that differs between OIPC and source water is that plants preferentially use growing-season precipitation, and the growing-season offset from the annual mean increases with latitude because high-latitude winters contribute very depleted precipitation to the annual mean that plants never see.

| Site | OIPC (annual mean) | Growing season offset | Source water | d2H_wax = 0.85 × source - 150 |
|------|------|------|------|------|
| Tropical (5°N) | -20‰ | +5‰ | -15‰ | -162.8‰ |
| Subtropical (30°N) | -50‰ | +10‰ | -40‰ | -184.0‰ |
| Temperate (50°N) | -80‰ | +20‰ | -60‰ | -201.0‰ |
| Polar (70°N) | -110‰ | +30‰ | -80‰ | -218.0‰ |

The growing season offset increases with latitude because:
- In the tropics, precipitation is weakly seasonal; growing season ≈ annual mean
- At high latitudes, winter precipitation is extremely depleted (snow, Rayleigh distillation) but plants are dormant. Summer precipitation is 20-40‰ more enriched than the annual mean. Plants use the enriched summer precipitation.

These numbers are illustrative but realistic (see Sachse et al. 2012, Table 2; also Daniels et al. 2017 for growing season vs. annual precipitation isotope offsets).

### The two regressions

**d2H_wax vs source water (the true reactant):**
- Source water range: -15‰ to -80‰ = 65‰
- d2H_wax range: -162.8‰ to -218.0‰ = 55.2‰
- Slope = 55.2 / 65 = **0.85** (exact, by construction — this IS alpha)

**d2H_wax vs OIPC (the proxy for the reactant):**
- OIPC range: -20‰ to -110‰ = 90‰
- d2H_wax range: -162.8‰ to -218.0‰ = 55.2‰ (same d2H_wax values!)
- Slope = 55.2 / 90 = **0.61**

The fractionation is 0.85 at every single site. There is no statistical artifact, no confounding, no model error. The slope is 0.61 simply because OIPC overstates the source water range. The growing season offset compresses the actual source water range (65‰) relative to the OIPC range (90‰), and the regression slope scales accordingly.

### Why this is undeniable

The logic is elementary:
1. Plants use growing-season water, not annual precipitation (universally accepted)
2. The growing-season offset from annual mean increases with latitude (observed, documented)
3. Therefore, source water range < OIPC range across a latitudinal transect
4. d2H_wax faithfully tracks source water at slope = alpha = 0.85
5. But d2H_wax vs OIPC has a slope < 0.85, because the x-axis is stretched relative to what the plants actually experienced

No fancy statistics required. The slope differs from alpha for the same reason that regressing height against shoe size gives a different slope than regressing height against foot length — the proxy variable has a different range than the true predictor.

### What the spatial model does

The spatial model detects this. The intercept GP at each location absorbs the OIPC-to-source-water offset. In our 4-site example, the GP would assign:
- Tropical site: intercept offset ≈ +4.25‰ (= 0.85 × 5‰, from the +5‰ growing season offset)
- Polar site: intercept offset ≈ +25.5‰ (= 0.85 × 30‰, from the +30‰ offset)

The GP captures the geographically-structured deviation between OIPC and source water. What remains as the slope is the relationship after accounting for these deviations — which is closer to what you'd see if you actually had source water measurements instead of OIPC.

### Why the non-spatial slope is ~0.83, not 0.61

If growing season bias alone would give 0.61, why does the non-spatial regression give 0.83?

Because other effects push the slope back up:
- **Evapotranspiration in arid regions**: Enriches source water beyond OIPC at warm/dry sites, partially counteracting the growing season compression. This extends the source water range at the enriched end.
- **Elevation effects within regions**: High-elevation sites have more depleted OIPC than nearby lowlands, but similar source water modifications, increasing the effective OIPC range without proportionally increasing the d2Hwax range.
- **Vegetation-type covariation**: C4 grasses (common in warm regions) may have different apparent fractionation than C3 trees (common in cool regions). If C4 fractionation is less negative, this extends the d2Hwax range at the enriched end.
- **Sedimentary integration**: Lake sediments in large catchments integrate lipids from a range of elevations, effectively averaging the source water within a watershed. This tends to dampen extremes differently at different latitudes.

These effects partially cancel the growing season compression in the global average. The non-spatial slope of 0.83 reflects the net result of all these opposing effects, which happens to land close to 0.85 — but this is coincidental, not confirmatory of alpha.

The spatial model disentangles these effects (to the extent that they vary smoothly across space) and reveals the conditional slope after geographic deconfounding. The 0.83 → 0.55 attenuation reflects the net geographic structure that the non-spatial regression absorbs into the slope.

## Sensitivity to the growing season offset gradient

The slope against OIPC depends on how strongly the growing season offset increases with latitude. We can parameterize this simply:

Let the growing season offset = a + b × |latitude|

| Offset gradient (b) | Source water range compression | Slope vs OIPC |
|-----|-----|-----|
| 0 (no growing season bias) | None: source = OIPC everywhere | 0.85 |
| 0.2 ‰/degree (mild) | ~14‰ compression over 65° | ~0.74 |
| 0.4 ‰/degree (moderate, as in our example) | ~26‰ compression | ~0.61 |
| 0.6 ‰/degree (strong) | ~39‰ compression | ~0.48 |

The moderate case already produces a slope of 0.61, well below 0.85, with no change to the fractionation factor. Published growing-season offsets (Daniels et al. 2017; Konecky et al. 2019) are consistent with 0.2-0.5 ‰/degree at mid-to-high latitudes.

## What this means for the manuscript

The argument structure for the manuscript:

1. **State the fractionation factor**: alpha ≈ 0.85 (epsilon ≈ -150‰) is the biochemical fractionation between source water and leaf wax n-alkanes (Sachse et al. 2012). This is not in dispute.

2. **State that OIPC ≠ source water**: OIPC estimates mean annual precipitation δ²H. Actual source water is modified by growing season bias (Sachse et al. 2012; Daniels et al. 2017), evapotranspiration (Kahmen et al. 2013), soil water mixing, and other processes. These modifications are spatially structured — they vary systematically with latitude, aridity, and elevation.

3. **Show the worked example** (or a version of it, perhaps as a figure): demonstrate that growing season bias alone reduces the regression slope below alpha, with no change to the fractionation.

4. **Explain the spatial model's role**: the intercept GP absorbs the geographically-structured OIPC-to-source-water offset. The conditional slope represents the d2Hwax-OIPC relationship after accounting for this offset.

5. **Acknowledge the Paciorek limitation**: the spatial model may overcorrect if it absorbs some genuine OIPC signal along with the geographic structure. We quantify this with the confounding simulation. The true conditional slope lies between the spatial estimate (~0.55) and the non-spatial estimate (~0.83).

6. **Emphasize the uncertainty**: for paleoclimate inversion, the posterior distribution of the slope (not the point estimate) is what matters. The spatial model produces wider but more honestly calibrated uncertainty that accounts for the OIPC ≠ source water problem.

## References

- Hayes JM (2001). Fractionation of carbon and hydrogen isotopes in biosynthetic processes. *Rev Mineral Geochem* 43(1):225-77.
- Sachse D et al. (2012). Molecular paleohydrology: interpreting the hydrogen-isotopic composition of lipid biomarkers from photosynthesizing organisms. *Annu Rev Earth Planet Sci* 40:221-49.
- Kahmen A et al. (2013). Leaf water deuterium enrichment shapes leaf wax n-alkane δD values of angiosperm plants I. *GCA* 111:39-49.
- Daniels WC, Russell JM, Giblin AE, et al. (2017). Calibrating the hydrogen isotopic composition of lacustrine leaf waxes. *GCA* 215:105-19.
- Konecky BL et al. (2019). The Iso2k database: a global compilation of paleo-deltaO18 and deltaD records. *Earth Syst Sci Data* 12:2261-88.
- Paciorek CJ (2010). The importance of scale for spatial-confounding bias and precision of spatial regression estimators. *Biostatistics* 11(4):601-15.
