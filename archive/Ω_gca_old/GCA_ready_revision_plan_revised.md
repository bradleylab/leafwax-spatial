# GCA Revision Plan - Updated 2025-11-03

## Executive Summary
Transform the manuscript from a statistical methods paper into a **geochemical insights paper** that uses spatial modeling to reveal fundamental truths about leaf wax isotope systematics. The headline finding is that **spatial confounding inflates apparent Î´Â²H_wax-Î´Â²H_precip relationships by ~50%**, with profound implications for paleoclimate reconstructions.

---

## 1. Core Messages & Positioning

### Primary Discovery
**Spatial confounding accounts for ~50% of the apparent global correlation** between Î´Â²H_wax and Î´Â²H_precip. When spatial structure is properly accounted for, the true process-based slope drops from 0.82 to 0.44.

### Key Geochemical Insights
1. **Continental-scale integration** (Î» = 3,400-4,200 km) dominates the Î´Â²H_wax signal
2. **Location matters more than vegetation** - spatial position explains 54-68% of variance
3. **Universal slope, regional intercepts** - the fundamental relationship is consistent globally but shifted regionally
4. **30% uncertainty reduction** possible with spatial models vs. conventional calibrations

### Target Audience
Paleoclimate reconstructors, organic geochemists, and proxy developers who need robust calibrations with proper uncertainty quantification.

---

## 2. Title Options (â‰¤20 words)

**Preferred:**
"Spatial confounding inflates leaf wax hydrogen isotope calibrations: Evidence from global Bayesian analysis"

**Alternatives:**
- "Continental-scale controls on leaf wax Î´Â²H revealed by spatial modeling of global sediments"
- "Quantifying spatial structure in the leaf wax paleohydrology proxy using hierarchical Bayesian methods"

---

## 3. Abstract Structure (250 words)

**Opening (Problem):** Leaf wax Î´Â²H_wax is widely used for paleoclimate reconstruction, but conventional calibrations ignore spatial autocorrelation, potentially biasing relationships and inflating uncertainties.

**Methods (Approach):** We compiled 818 global lake sediment n-Câ‚‚â‚‰ alkane measurements and developed Bayesian spatial models incorporating Gaussian process approximations.

**Results (Key findings):**
- Spatial autocorrelation accounts for 54-68% of total variance
- Global Î´Â²H_wax-Î´Â²H_precip slope drops from 0.82 (non-spatial) to 0.44 (spatial)
- Spatial effects operate at continental scales (3,400-4,200 km)
- Vegetation and elevation show minimal influence once spatial structure is included
- Regional offsets dominate over slope variations (intercept SD = 14.8-16.6â€°, slope SD = 0.001-0.022)

**Implications (So what?):**
- Up to 50% of apparent correlation in global datasets reflects spatial confounding
- Location-specific calibrations reduce uncertainty by ~30% (RMSE: 21â€° â†’ 15â€°)
- Models using only paleo-reconstructable variables perform within 2% of best model
- R package 'leafwax' provides tools for spatial-aware reconstructions

**Conclusion:** These findings demonstrate that geographic context dominates leaf wax isotope systematics and that accounting for spatial structure is essential for accurate paleoclimate reconstruction.

---

## 4. Main Text Structure

### 4.1 Introduction (~1,000 words)

**Paragraph 1: The proxy and its importance**
- Leaf wax Î´Â²H as molecular recorder of hydroclimate
- Applications from seasonal to millennial scales
- Current use in paleoclimate reconstructions

**Paragraph 2: Known complexities**
- Biological fractionation (PFT effects)
- Environmental influences (temperature, humidity, seasonality)
- Integration across catchments

**Paragraph 3: The calibration challenge**
- Conventional OLS approaches and assumptions
- Growing recognition of spatial autocorrelation
- Implications of ignoring spatial structure

**Paragraph 4: Spatial confounding problem**
- Geographic covariation inflates correlations
- Example: C4 plants cluster in specific climates
- Need to separate process from geography

**Paragraph 5: Our approach**
- Hierarchical Bayesian framework with spatial GP
- Global compilation of lake sediments
- Objectives: (1) quantify spatial structure, (2) reveal true process relationships, (3) improve predictions

**Paragraph 6: "Here we..." statement**
- Clear statement of what we did and found
- Preview of spatial confounding discovery

### 4.2 Methods (~1,500 words main + extensive supplement)

#### 4.2.1 Data Compilation
**Main text (300 words):**
- 818 lake sediment n-Câ‚‚â‚‰ alkane measurements
- Geographic distribution and sampling density
- Key covariates: Î´Â²H_precip (OIPC), elevation, vegetation (MODIS)
- Brief quality control

**To supplement:**
- Full filtering criteria
- Data matching procedures
- Covariate extraction methods

#### 4.2.2 Spatial Modeling Framework
**Main text (600 words):**
- Conceptual overview (prose, not equations)
- Gaussian process approximation via predictive process
- Key equation showing decomposition: observation = fixed effects + spatial field + nugget
- Interpretation of length scale parameter Î»

**To supplement:**
- Full mathematical derivation (Equations 7-13)
- Kriging details
- Computational implementation

#### 4.2.3 Model Configurations
**Main text (300 words):**
- Model variants tested (spatial vs non-spatial, with/without vegetation)
- Prior specification rationale (PC priors for variance)
- Model comparison approach (LOOIC)

**To supplement:**
- Complete prior specifications
- Stan implementation details
- Convergence diagnostics

#### 4.2.4 Handling Collinearity
**Main text (100 words):**
- Brief mention of VIF screening and elastic net
- Final predictor selection

**To supplement:**
- Correlation matrices
- VIF tables
- Elastic net coefficient paths

#### 4.2.5 Cross-validation
**Main text (200 words):**
- Spatial block CV strategy
- Performance metrics

**To supplement:**
- Detailed CV methodology
- Regional stratification

### 4.3 Results (~2,000 words)

#### 4.3.1 Model Performance (300 words)
**Focus:** Brief summary of predictive skill
- Overall RÂ² and RMSE
- Coverage of prediction intervals
- Improvement from spatial models

#### 4.3.2 Spatial Confounding Discovery (500 words)
**Focus:** The headline finding
- Slope comparison: 0.82 â†’ 0.44
- Interpretation of the change
- Implications for existing calibrations
- Figure showing spatial vs non-spatial relationships

#### 4.3.3 Spatial Structure Characterization (400 words)
**Focus:** What the spatial patterns tell us
- Variance partition: 54-68% spatial
- Length scales: 3,400-4,200 km
- Continental-scale integration
- Figure of spatial weights at different scales

#### 4.3.4 Regional Patterns (400 words)
**Focus:** Geographic heterogeneity
- Asia challenge (highest RMSE)
- Best performance in Americas/Africa
- Systematic regional offsets
- Universal slope, varying intercepts

#### 4.3.5 Environmental Predictors (400 words)
**Focus:** Why vegetation and elevation don't matter much
- Minimal vegetation effects despite expectations
- Elevation uncertainty spans ~250â€°
- Spatial field captures these effects
- Implications for paleo-reconstructions

### 4.4 Discussion (~2,000 words)

#### 4.4.1 Rethinking Global Calibrations (500 words)
**Lead with spatial confounding implications**
- Previous calibrations overestimate process strength
- Geographic clustering biases relationships
- Need for spatial-aware approaches

#### 4.4.2 Continental Integration of Signals (400 words)
**Geochemical interpretation of length scales**
- 3,800 km integration implies basin-scale averaging
- Atmospheric moisture transport scales
- Implications for proxy interpretation

#### 4.4.3 The Vegetation Paradox (400 words)
**Why lab studies don't match field observations**
- Spatial confounding explanation
- Scale of integration homogenizes signals
- Coarse resolution of global products

#### 4.4.4 Practical Applications (400 words)
**Guidance for reconstructions**
- When to use global vs regional calibrations
- Uncertainty propagation methods
- Software tools available

#### 4.4.5 Limitations and Future Directions (300 words)
- Sampling gaps (especially tropics)
- Stationarity assumptions
- Need for temporal validation

### 4.5 Conclusions (~300 words)

**Four key statements:**
1. Spatial confounding inflates apparent Î´Â²H_wax-Î´Â²H_precip relationships by ~50%
2. Continental-scale processes (Î» â‰ˆ 3,800 km) dominate leaf wax isotope patterns
3. Geographic location explains more variance than vegetation or elevation
4. Spatial models reduce reconstruction uncertainty by 30% and enable rigorous error propagation

**Final thought:** The 'leafwax' R package implements these methods for the community.

---

## 5. Figures Strategy (8 main text)

### Main Text Figures
1. **Global sample distribution** - Map showing data density and gaps
2. **Spatial confounding demonstration** - Side-by-side scatter plots showing slope change
3. **Observed vs predicted** - With uncertainty bands, comparing spatial/non-spatial
4. **Variance decomposition** - Visual showing spatial vs nugget components
5. **Length scale interpretation** - Kernel weights at different distances
6. **Regional performance** - Map or bar chart of RMSE by region
7. **Elevation uncertainty** - Spline with enormous confidence bands
8. **Practical example** - Reconstruction with proper uncertainty propagation

### Supplementary Figures
- Correlation matrices and VIF diagnostics
- Posterior distributions and traces
- Residual maps and diagnostics
- Per-site predictions
- Cross-validation details
- Regional spatial fields

---

## 6. Tables Strategy

### Main Text Tables (2-3 max)
1. **Model comparison** - LOOIC, RÂ², RMSE for key models
2. **Parameter estimates** - Slopes, intercepts, variance components
3. **Regional performance** - Summary statistics by continent (optional)

### Supplementary Tables
- Full data summary statistics
- Complete model specifications
- Cross-validation results by region
- Prior specifications

---

## 7. Supplementary Information Structure

### SI-1: Mathematical Framework
- Complete derivation of spatial model
- Gaussian process approximation details
- Kriging equations
- PC prior specifications

### SI-2: Data and Processing
- Detailed filtering criteria
- Covariate extraction methods
- Data quality assessments
- Missing data handling

### SI-3: Computational Methods
- Stan model code
- Convergence diagnostics
- Sensitivity analyses
- Cross-validation procedures

### SI-4: Extended Results
- Full posterior summaries
- Regional analyses
- Alternative model configurations
- Robustness checks

### SI-5: Additional Figures and Tables
- Diagnostic plots
- Detailed maps
- Extended statistical summaries

---

## 8. Writing Style Guidelines

### Do's
âœ“ Lead with geochemical insights, not statistical methods
âœ“ Use active voice and clear statements
âœ“ Emphasize practical implications for reconstructions
âœ“ Connect all findings back to the proxy problem
âœ“ Keep technical details minimal in main text

### Don'ts
âœ— Don't lead paragraphs with statistical jargon
âœ— Avoid extensive mathematical notation in main text
âœ— Don't assume readers are Bayesian statisticians
âœ— Minimize "we show that" constructions
âœ— No extensive model comparison in main text

---

## 9. Key Revisions from Current Draft

### Major Changes
1. **Restructure abstract** to lead with spatial confounding discovery
2. **Move most equations** to supplement (keep only 1-2 essential)
3. **Consolidate methods** to focus on concepts, not implementation
4. **Reorganize results** to feature spatial confounding prominently
5. **Reframe discussion** around geochemical implications

### Content to Move to Supplement
- Equations 3-13 and their derivations
- Detailed prior specifications
- Convergence diagnostics
- Extended model comparisons
- Technical implementation details

### New Content to Add
- Clearer explanation of spatial confounding
- More geochemical interpretation of patterns
- Practical reconstruction example
- Software implementation guide

---

## 10. Timeline and Action Items

### Immediate (Week 1)
1. Restructure abstract around spatial confounding
2. Create supplement structure and move technical content
3. Rewrite introduction with geochemical focus

### Short-term (Week 2)
1. Consolidate methods section
2. Reorganize results around key findings
3. Draft new discussion emphasizing implications

### Final (Week 3)
1. Create publication-quality figures
2. Polish prose for clarity and flow
3. Ensure all GCA formatting requirements met
4. Prepare submission package

---

## 11. Success Metrics

The revision succeeds if:
- **A geochemist can understand** the main findings without statistical background
- **The spatial confounding message** is clear and compelling
- **Practical applications** are obvious to proxy users
- **Technical readers** can access full details in supplement
- **The paper reads** as a geochemical insight paper, not a methods paper

---

## 12. Response to Reviewers Template

*For future use when addressing reviewer comments*

"We have substantially revised the manuscript to emphasize the geochemical insights rather than statistical methodology. The key discoveryâ€”that spatial confounding accounts for ~50% of apparent correlations in global calibrationsâ€”is now prominently featured. Technical details have been moved to the supplement while maintaining full reproducibility..."