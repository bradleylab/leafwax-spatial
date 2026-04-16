# Methods Section Revision Example

## BEFORE (Current Draft) - Too Technical

### 2.3 Hierarchical Spatial Model Framework

We employ a hierarchical Bayesian framework that decomposes the observed δ²H_wax signal into fixed effects, spatially-structured random effects, and unstructured error:

Level 1 (Observation model):
δ²H_wax,i ~ Normal(μᵢ, τ²)

Level 2 (Process model):
μᵢ = X₁β + Zw* + ε_spatial

where w* ~ GP(0, C*) with C*ᵢⱼ = η² × ρ(||s*ᵢ - s*ⱼ||; λ, ν)

Level 3 (Parameter model):
β ~ Normal(0, 100²)
η² ~ PC(U = 60, α = 0.01)
λ ~ PC(U = 4000, α = 0.01)
τ² ~ Exponential(1/20²)

The spatial covariance function ρ employs a Matérn formulation with smoothness parameter ν = 3/2:

ρ(d; λ, ν=3/2) = (1 + √3d/λ)exp(-√3d/λ)

We implement the predictive process approximation using m = 150 knots placed on a Fibonacci grid...

[Continues with extensive mathematical detail]

---

## AFTER (Revised for Main Text) - Geochemically Focused

### 2.2 Accounting for Spatial Structure in Leaf Wax Signals

Leaf wax δ²H values from nearby locations tend to be more similar than those from distant sites, reflecting shared regional climate patterns, moisture sources, and vegetation assemblages. To account for this spatial autocorrelation, we developed a hierarchical model that separates the δ²H_wax signal into three components:

**δ²H_wax = Environmental predictors + Spatial pattern + Local variation**

The **environmental predictors** include precipitation δ²H from OIPC, elevation, and vegetation cover—variables we can measure or estimate for both modern and paleo settings. These represent the mechanistic controls on leaf wax isotopes that operate globally.

The **spatial pattern** captures regional coherence in the isotope signal that isn't explained by our measured predictors. This component uses a Gaussian process—essentially a flexible surface that adapts to the data—to model how similarity between sites decreases with distance. The key parameter here is the "length scale" (λ), which tells us the distance over which leaf wax signals remain correlated. A length scale of 3,800 km, for example, means that sites within this distance share similar regional influences on their δ²H_wax values beyond what precipitation isotopes alone would predict.

The **local variation** represents site-specific factors we cannot measure or model explicitly—local hydrology, micro-climate effects, or sampling artifacts. This "nugget effect" sets the minimum prediction uncertainty achievable even with perfect knowledge of regional patterns.

By explicitly modeling spatial structure, we can distinguish between true process relationships (how δ²H_precip translates to δ²H_wax) and geographic covariation (sites with similar δ²H_precip happening to cluster geographically). This separation is critical because conventional calibrations conflate these two effects, potentially inflating the apparent strength of proxy relationships.

We tested 14 model configurations, varying which environmental predictors were included and whether spatial structure was considered. Model comparison used leave-one-out cross-validation (LOOIC), which evaluates predictive performance while accounting for model complexity. All models were implemented in a Bayesian framework, providing full uncertainty quantification for predictions and parameters.

---

## WHAT GOES TO SUPPLEMENT

### Supplement Section S2: Mathematical Framework

#### S2.1 Hierarchical Model Specification

The full hierarchical model is specified as follows:

**Level 1 (Observation model):**
```
δ²H_wax,i ~ Normal(μᵢ, τ²)
```

**Level 2 (Process model):**
```
μᵢ = X₁β + Zw* + ε_spatial
```

where:
- X represents the design matrix of environmental covariates
- β contains the fixed effect coefficients
- Z is the interpolation matrix from knots to observation locations
- w* represents the Gaussian process at knot locations

**Level 3 (Parameter model):**

The Gaussian process is defined as:
```
w* ~ GP(0, C*)
```

with covariance:
```
C*ᵢⱼ = η² × ρ(||s*ᵢ - s*ⱼ||; λ, ν)
```

We use penalized complexity (PC) priors for spatial parameters:
- η² ~ PC(U = 60, α = 0.01): spatial variance
- λ ~ PC(U = 4000, α = 0.01): length scale
- τ² ~ Exponential(1/20²): nugget variance

#### S2.2 Matérn Covariance Function

The spatial covariance employs a Matérn function with ν = 3/2:
```
ρ(d; λ, ν=3/2) = (1 + √3d/λ)exp(-√3d/λ)
```

This provides a flexible, once-differentiable spatial surface...

[Continue with full technical details, equations, derivations, etc.]

---

## KEY PRINCIPLES FOR REVISION

### Main Text Should:
1. **Explain concepts** without heavy mathematics
2. **Use analogies** (e.g., "flexible surface that adapts to data")
3. **Connect to geochemistry** (what it means for the proxy)
4. **Define technical terms** immediately in plain language
5. **Focus on interpretation** not implementation

### Supplement Should:
1. **Provide complete mathematics** for reproducibility
2. **Include all equations** with full notation
3. **Detail implementation** (code, convergence, etc.)
4. **Show diagnostics** and sensitivity analyses
5. **Enable technical readers** to fully understand/reproduce

### Example Transformations:

**Too Technical:** "We employ a Matérn covariance function with ν = 3/2"
**Better:** "We model how correlation between sites decreases with distance"

**Too Technical:** "The predictive process uses m = 150 knots on a Fibonacci grid"  
**Better:** "We efficiently compute spatial patterns using 150 strategically placed reference points"

**Too Technical:** "η²/(η² + τ²) represents the proportion of variance attributable to spatial structure"
**Better:** "Spatial patterns explain 54-68% of the total variation in leaf wax values"

**Too Technical:** "The PC prior on λ has U = 4000 km with probability α = 0.01"
**Better:** "We allow the spatial correlation to extend up to continental scales"
