# Model Architecture Decisions — Summary

**Source**: `12-1_choices.docx` and `12-1_new.docx` (model v12.1 development notes, Jul 2025)
**Purpose**: Records the rationale behind key design choices for the spatial model. Useful for methods writing and reviewer responses.

---

## Key Design Choices and Rationale

### 1. Predictive Process GP (not full GP)
- **Choice**: Predictive process with 120-125 knots instead of full 818×818 GP
- **Why**: Computational efficiency — O(Nk²) vs O(N³)
- **Trade-off**: Some loss of fine-scale variation for massive speedup

### 2. Fibonacci Sphere for Knot Placement
- **Choice**: Equal-area global grid using the golden angle algorithm
- **Why**: Need a generalizable model that can predict anywhere on Earth
- **Alternatives tried and discarded**:
  - `cover.design` on data locations → knots too close together due to duplicate coordinates
  - `kmeans` → poor coverage in data-sparse regions
  - `maximin` → better but still data-dependent (can't predict outside data range)
  - Simple lat/lon grid → unequal area (clustering at poles)

### 3. Exponential Decay for Multi-Scale Integration
- **Choice**: Fit a single lambda_decay parameter to weight spatial scales (1-150 km)
- **Why**: Data-driven approach to determine optimal spatial scale of environmental influence
- **Previous approach**: Fixed weights or interpolation between pre-computed scales

### 4. C4 Standardization (Fixed, not data-derived)
- **Choice**: Standardize C4 with fixed mean=20, sd=25 rather than from data
- **Why**: C4 has a known global distribution; fixed values prevent data-dependent standardization from affecting model predictions at new locations

### 5. GP Length Scale
- **Choice**: Prior centered on ~2400 km (matching ~1600 km knot spacing)
- **Why**: Ensures smooth interpolation between knots without over- or under-fitting
- **Critical lesson**: Length scale must be converted from km to standardized coordinate space using `coord_scale_km = mean(coord_scaling) * 111.0`

### 6. Elevation B-Spline
- **Choice**: Cubic B-spline with 9 interior knots (13 basis functions)
- **Why**: Capture non-linear elevation effects without over-parameterizing
- **Note**: Almost lost during a refactoring — the elevation effect is real and important

### 7. Numerical Stability
- Kernel jitter: 1e-4 (added to kernel matrix diagonal)
- Coordinate jitter: 0.0001° (~11m) to break duplicate coordinates
- adapt_delta: 0.95 (non-spatial) / 0.99 (spatial)
- All controlled via config.yaml

## Major Issues Encountered During Development

1. **Ill-conditioned kernel matrices**: Condition numbers > 400 million with cover.design knots. Root cause: duplicate locations + poor knot spacing. Fixed by Fibonacci sphere + proper length scale.

2. **Coordinate standardization warps distances**: Must account for the fact that standardized coordinates change inter-point distances. GP length scale must be in standardized units, not km.

3. **Duplicate locations cause Cholesky failures**: 818 observations but only 627 unique locations. Coordinate jitter (0.0001°) solves this.

4. **Configuration drift**: Hard-coded values scattered across scripts caused inconsistencies. Centralized in config.yaml.
