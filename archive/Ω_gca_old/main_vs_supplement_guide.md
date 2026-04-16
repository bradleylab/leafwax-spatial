# Main Text vs. Supplement Organization Guide

## ABSTRACT (250 words max) - MAIN TEXT ONLY

### Current Abstract Assessment:
✅ Good: Already mentions spatial confounding and slope change
✅ Good: Quantifies variance partition and uncertainty reduction
❌ Needs: Stronger opening about the spatial confounding discovery
❌ Needs: More emphasis on paleoclimate implications

### Revised Abstract Structure:
```
Leaf wax hydrogen isotope ratios (δ²H_wax) are widely used to reconstruct past precipitation isotopes, but conventional calibrations overlook spatial autocorrelation, potentially biasing relationships. We compiled 818 global lake sediment n-C₂₉ alkane δ²H_wax measurements and developed Bayesian spatial models to separate geographic covariation from true process relationships. Spatial autocorrelation accounts for 54-68% of total variance, operating at continental scales (3,400-4,200 km). Critically, accounting for spatial structure reveals that the apparent global δ²H_wax-δ²H_precip slope of 0.82 drops to 0.44 in spatial models, indicating that up to half of the observed correlation in global datasets reflects spatial confounding rather than direct isotopic transfer. The δ²H_wax-δ²H_precip relationship shows near-universal slopes (SD = 0.001-0.022) but large regional intercept offsets (SD = 14.8-16.6‰), suggesting consistent fractionation processes operating on regionally-shifted baselines. Environmental predictors show surprisingly minimal influence: vegetation type contributes <2‰ per 10% cover change, and elevation effects remain highly uncertain (90% CI ~250‰) once spatial structure is included. Location-specific calibrations reduce prediction uncertainty by 30% (RMSE: 21‰ → 15‰). Models using only paleo-reconstructable variables (δ²H_precip + elevation) perform within 2% of the best model, enabling practical inversions without vegetation reconstruction. These findings demonstrate that geographic context dominates leaf wax systematics and that spatial confounding fundamentally alters calibration relationships. The R package 'leafwax' implements these spatial models for improved paleoclimate reconstructions.
```

---

## 1. INTRODUCTION

### MAIN TEXT (1,000 words):

**Paragraph 1: The proxy importance (150 words)**
- Leaf wax δ²H as molecular rainfall recorder
- Applications from seasonal to millennial timescales
- Global usage in paleoclimate studies
- *Keep from current:* Lines 43-52 (condensed)

**Paragraph 2: Known biological complexity (150 words)**
- PFT-specific fractionation patterns
- Species-level variations
- Mechanistic understanding of fractionation
- *Keep from current:* Lines 54-67 (condensed)

**Paragraph 3: Environmental influences (150 words)**
- Climate controls beyond δ²H_precip
- Temperature, humidity, soil moisture effects
- The collinearity problem
- *Keep from current:* Lines 69-86 (condensed)

**Paragraph 4: Spatial integration challenge (150 words)**
- Catchment-scale mixing
- Elevation and regional patterns
- Unknown source areas
- *Keep from current:* Lines 88-102 (condensed)

**Paragraph 5: The calibration problem (200 words)**
- OLS limitations and assumptions
- Spatial autocorrelation recognition
- **NEW: Spatial confounding concept**
- Why this matters for reconstructions

**Paragraph 6: Our approach and findings (200 words)**
- Hierarchical Bayesian spatial framework
- Global compilation approach
- **EMPHASIZE: Discovery of spatial confounding**
- Preview of key findings
- *Adapt from:* Lines 104-113

### TO SUPPLEMENT:
- Extended review of calibration literature (current lines 115-230)
- Detailed discussion of biosynthetic fractionation mechanisms
- Mathematical framework introduction
- Extended discussion of previous spatial approaches

---

## 2. METHODS

### MAIN TEXT (1,500 words):

#### 2.1 Data Compilation (300 words)
**Keep in main:**
- Sample counts and distribution (brief)
- Key variables: δ²H_precip, elevation, vegetation
- Basic filtering criteria
- *Condensed from current:* Lines 232-280

**Move to supplement:**
- Detailed filtering procedures
- Data source descriptions
- Matching procedures
- Quality control details

#### 2.2 Spatial Modeling Framework (600 words)
**Keep in main:**
- Conceptual overview (prose, minimal math)
- ONE key equation: `δ²H_wax = Xβ + GP(location) + ε`
- Explanation of Gaussian process in plain language
- Length scale interpretation
- Variance decomposition concept
- *Simplified from current:* Lines 350-450

**Move to supplement:**
- All mathematical derivations (Equations 3-13)
- Predictive process details
- Kriging equations
- Computational specifics
- Matrix formulations

#### 2.3 Model Variants (300 words)
**Keep in main:**
- Table or list of model configurations
- Basic rationale for each variant
- Prior philosophy (PC priors concept)
- Model comparison approach (LOOIC)
- *Condensed from current:* Lines 500-580

**Move to supplement:**
- Detailed prior specifications
- Stan implementation code
- Convergence diagnostics
- Sensitivity analyses

#### 2.4 Collinearity Handling (100 words)
**Keep in main:**
- One sentence on VIF screening
- Final predictor selection
- *Brief mention from current:* Lines 600-610

**Move to supplement:**
- Correlation matrices
- VIF tables
- Elastic net analysis
- Variable selection details

#### 2.5 Cross-validation (200 words)
**Keep in main:**
- Spatial block CV concept
- Performance metrics used
- *Simplified from current:* Lines 620-650

**Move to supplement:**
- Detailed CV methodology
- Regional stratification
- Leave-one-out details

---

## 3. RESULTS

### MAIN TEXT (2,000 words):

#### 3.1 Model Performance (300 words)
**Keep in main:**
- Summary table of key models (LOOIC, R², RMSE)
- Best model identification
- Brief performance comparison
- *Condensed from current:* Section 3.1

**Move to supplement:**
- Full model comparison table
- Detailed performance metrics
- Residual analyses
- Regional breakdowns

#### 3.2 Spatial Confounding Discovery (500 words)
**Keep in main:**
- **FEATURE:** Slope comparison figure (0.82 → 0.44)
- Interpretation of the change
- What this means for existing calibrations
- Geographic covariation explanation
- *Expand from current:* Lines 1235-1244

**Move to supplement:**
- Technical details of slope calculation
- Bootstrap confidence intervals
- Sensitivity to model specification

#### 3.3 Spatial Structure (400 words)
**Keep in main:**
- Variance partition results (54-68% spatial)
- Length scale estimates (3,400-4,200 km)
- Figure showing kernel weights
- Continental-scale interpretation
- *From current:* Section 3.2 key points

**Move to supplement:**
- Full posterior distributions
- Alternative parameterizations
- Sensitivity to prior choices
- Technical variance calculations

#### 3.4 Regional Patterns (400 words)
**Keep in main:**
- Asia challenge discussion
- Regional performance summary
- Map or figure of regional differences
- Universal slope, varying intercepts finding
- *From current:* Lines 1258-1306

**Move to supplement:**
- Detailed regional statistics
- Country-level analyses
- Site-specific results
- Extended regional maps

#### 3.5 Environmental Predictors (400 words)
**Keep in main:**
- Minimal vegetation effects finding
- Elevation uncertainty demonstration
- Why spatial field captures these effects
- Implications for reconstructions
- *From current:* Lines 1308-1360

**Move to supplement:**
- Full coefficient tables
- Posterior distributions
- Alternative vegetation metrics
- Detailed elevation spline analysis

---

## 4. DISCUSSION

### MAIN TEXT (2,000 words):

#### 4.1 Rethinking Global Calibrations (500 words)
**Keep in main:**
- Spatial confounding implications
- Why previous calibrations overestimate
- Geographic clustering biases
- Need for spatial approaches
- *Expand from current spatial confounding discussion*

#### 4.2 Continental Integration (400 words)
**Keep in main:**
- What 3,800 km means geochemically
- Basin-scale averaging
- Moisture transport scales
- Proxy interpretation implications
- *New geochemical interpretation*

#### 4.3 The Vegetation Paradox (400 words)
**Keep in main:**
- Lab vs field discrepancy
- Spatial confounding explanation
- Scale of integration
- Resolution limitations
- *Expand from current:* Lines 1349-1358

#### 4.4 Practical Applications (400 words)
**Keep in main:**
- When to use different calibrations
- Uncertainty propagation guidance
- Software implementation
- Example application
- *New practical guidance*

#### 4.5 Limitations and Future (300 words)
**Keep in main:**
- Sampling gaps
- Stationarity assumptions
- Future directions
- *Condensed from current discussion*

### TO SUPPLEMENT:
- Extended literature comparisons
- Technical discussion of methods
- Additional limitations
- Detailed software documentation

---

## 5. CONCLUSIONS

### MAIN TEXT (300 words):
1. Spatial confounding inflates correlations by ~50%
2. Continental-scale processes dominate (3,800 km)
3. Location > vegetation for variance explanation
4. 30% uncertainty reduction with spatial models
5. R package availability

### NO SUPPLEMENT NEEDED

---

## FIGURES ALLOCATION

### MAIN TEXT (8 figures):
1. **Figure 1:** Global sample distribution map
2. **Figure 2:** Spatial confounding demonstration (slope comparison)
3. **Figure 3:** Observed vs predicted with uncertainty
4. **Figure 4:** Variance decomposition visualization
5. **Figure 5:** Spatial kernel weights at different λ
6. **Figure 6:** Regional performance map/chart
7. **Figure 7:** Elevation uncertainty demonstration
8. **Figure 8:** Practical reconstruction example

### SUPPLEMENT (unlimited):
- S1: Correlation matrices
- S2: VIF diagnostics
- S3: Posterior distributions
- S4: Convergence diagnostics
- S5: Residual maps
- S6: Regional spatial fields
- S7: Cross-validation details
- S8: Elastic net paths
- S9: Alternative model comparisons
- S10: Site-specific predictions

---

## TABLES ALLOCATION

### MAIN TEXT (2-3 tables):
1. **Table 1:** Model comparison (top 5-6 models only)
2. **Table 2:** Key parameter estimates (slopes, intercepts, variance)
3. **Table 3:** Regional performance summary (optional)

### SUPPLEMENT:
- S1: Full model comparison (all 14 models)
- S2: Data summary statistics
- S3: Prior specifications
- S4: Detailed coefficient estimates
- S5: Cross-validation results by region

---

## EQUATIONS ALLOCATION

### MAIN TEXT (1-2 equations only):
1. Basic model: `δ²H_wax = Xβ + GP(s) + ε`
2. Maybe variance decomposition (optional)

### SUPPLEMENT (all others):
- Full hierarchical model specification
- Gaussian process formulation
- Predictive process equations
- Kriging equations
- Prior specifications
- Variance calculations
- All numbered equations from current draft

---

## SPECIFIC CONTENT TO CUT/MOVE

### Lines to DELETE entirely:
- Overly technical method descriptions
- Redundant literature review
- Excessive statistical jargon
- Detailed software implementation

### Lines to MOVE to supplement:
- 115-230: Extended calibration review
- 350-450: Mathematical framework (keep conceptual only)
- 500-580: Detailed model specifications
- 600-610: Full collinearity analysis
- Most equations (keep 1-2 max)
- Detailed tables
- Technical discussions

### Lines to CONDENSE:
- 43-102: Introduction background (tighten)
- 232-280: Data description (summarize)
- Results sections (focus on key findings)

---

## WRITING STYLE CHECKLIST

### For Main Text:
✅ Lead with geochemical insight, not statistics
✅ Explain spatial confounding clearly
✅ Use plain language for technical concepts
✅ Emphasize practical implications
✅ Connect everything to paleoclimate reconstruction
✅ Minimize equations and statistical notation
✅ Use active voice
✅ Keep paragraphs focused

### Avoid in Main Text:
❌ Statistical jargon without explanation
❌ Extensive mathematical notation
❌ Long technical discussions
❌ Implementation details
❌ Excessive model comparison
❌ Detailed diagnostics
❌ Programming/software specifics

---

## NEXT STEPS

1. **Immediate:** Rewrite abstract with spatial confounding as lead
2. **Day 1:** Restructure introduction (cut 50%, focus on problem)
3. **Day 2:** Simplify methods (move 70% to supplement)
4. **Day 3:** Reorganize results around key findings
5. **Day 4:** Rewrite discussion with geochemical focus
6. **Day 5:** Create main figures
7. **Day 6:** Organize supplement
8. **Day 7:** Polish and format
