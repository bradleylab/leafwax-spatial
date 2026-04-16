# GCA-Ready Revision Plan for “geosp_leafwax_2025-10-05.pdf”

This file summarizes concrete edits to prepare the manuscript for **Geochimica et Cosmochimica Acta (GCA)**.

---

## 0) Positioning & High-level Goals
- Reframe as a **geochemical** paper that *uses* an advanced spatial–Bayesian framework, not a statistics paper with δ²H data.
- Core messages to foreground:
  - Quantified **spatial integration scale (λ)** and **variance partition** (η² vs τ²) of δ²H\_wax.
  - **Plant-type (PFT) effects** on ε\_app with uncertainty.
  - **Uncertainty propagation** from δ²H\_wax → δ²H\_precip reconstructions.
- Target main text: ≤ **8,000 words**; equations and diagnostics → **Supplementary**.

---

## 1) Title & Abstract
**Keep (main):**
- Title ≤ 20 words focused on geochemical insight.
- Abstract ≤ 250 words, single paragraph: Problem → Approach → Results (λ, variance partition, PFT offsets, prediction skill) → Geochemical implications.

**Suggested titles:**
- “Spatial structure and plant-type controls on hydrogen isotope fractionation in plant waxes from a hierarchical Bayesian model”
- “Quantifying δ²H\_wax–δ²H\_precip relationships with a spatial Bayesian framework”

**Move to SI:** none.

**Rewrite focus:**
- Avoid equations/derivations.
- Conclude with implications for paleohydrologic reconstructions.

---

## 2) Introduction (~900 words)
**Keep (main):**
- δ²H\_wax as a proxy; specific limitations of constant-fractionation, no spatial uncertainty propagation.
- Clear objectives:
  1) What fraction of δ²H\_wax variability is spatial vs local (η²/(η²+τ²))?
  2) How ε\_app varies by vegetation/PFT.
  3) Whether the framework yields process-meaningful uncertainty for reconstructions.
- “Here we…” paragraph tying goals to the model and dataset.

**Move to SI:**
- Extended review of Bayesian/geostatistical literature not specific to isotopes.

**Rewrite focus:**
- Lead with geochemical problem; keep modelling as means to an end.

---

## 3) Methods (~2,500 words in main; rest to SI)

### 3.1 Data
**Keep (main):**
- Counts, spatial distribution, PFT classes; covariates used; brief filtering logic.
- Rationale for variable selection in geochemical terms.

**Move to SI:**
- Full data cross-matching/filtering scripts; ancillary datasets.

### 3.2 Model overview
**Keep (main):**
- Concise hierarchy in prose + **1–2 key equations** (model mean & covariance kernel).
- Concept of predictive-process/knots and how λ defines integration across precomputed radii (Eqs 9–10 summarized).

**Move to SI:**
- Full math for Eqs 7–13; kriging algebra; PC prior forms; knot–site covariance details.

### 3.3 Priors & computation
**Keep (main):**
- Rationale for **PC prior** on spatial variance, brief note on priors for range and nugget; software stack (e.g., Stan/INLA).

**Move to SI:**
- Exact prior parameterizations; sampler settings; convergence diagnostics; sensitivity analyses.

### 3.4 Collinearity handling
**Keep (main):**
- One sentence: correlated predictors removed/penalized (VIF > 5; elastic net supported choices).

**Move to SI:**
- Correlation matrix, VIF table, elastic-net coefficient paths.

---

## 4) Results (~2,000 words)

### 4.1 Performance
**Keep (main):**
- Cross-validation metrics (RMSE/bias), coverage of 95% intervals.
- One compact figure: observed vs predicted with uncertainty.

**Move to SI:**
- Per-site residuals; full CV tables.

### 4.2 Spatial parameters & variance partition
**Keep (main):**
- Posterior λ (interpretation as hydrologic integration length).
- Fraction spatial = η²/(η²+τ²) with CI; brief nugget interpretation.

**Move to SI:**
- Alternative prior tests; full posterior traces.

### 4.3 Plant-type fractionation
**Keep (main):**
- Posterior ε\_app by PFT; compare to literature boxplots.
- One box/violin plot figure.

**Move to SI:**
- Per-site ε\_app posteriors.

### 4.4 δ²H\_precip predictions & uncertainty
**Keep (main):**
- Mean and SD maps; text that explains regional patterns and data-sparse behavior.

**Move to SI:**
- High-res regional panels; additional difference maps.

---

## 5) Discussion (~2,000 words)
**Keep (main):**
- Geochemical interpretation of λ (catchment/source-area scale; evapotranspiration coherence).
- Meaning of variance partition for proxy noise vs signal.
- Biological/physiological context of ε\_app offsets.
- Comparison with prior global calibrations (Sachse 2012; Tipple & Pagani 2013; Douglas 2020) and what changes for reconstructions.

**Move to SI:**
- Secondary model comparison details and derivations.

**Rewrite focus:**
- Start with key insights; put methods in background; end with practical guidance for applying the calibration and uncertainties.

---

## 6) Conclusions (~300 words)
- 3–4 concise statements:
  - Spatial process explains ~20–40% of variance (update with your result).
  - λ ≈ X km: effective environmental integration radius.
  - PFT-specific ε\_app offsets with uncertainties.
  - Framework provides formal uncertainty propagation to δ²H\_precip reconstructions.

---

## 7) Figures & Tables

**Main (≤ 8 figures):**
1. Workflow diagram (data → model → outputs).
2. Observed vs predicted δ²H\_wax (with intervals).
3. Mean map of predicted δ²H\_precip (or δ²H\_wax) and SD map.
4. Posterior λ (histogram or map of effective weights).
5. Variance partition (η² vs τ²) bar/violin with CI.
6. ε\_app by PFT (box/violin).
7. Cross-validation performance (condensed).
8. Sensitivity summary (if space allows) or move to SI.

**Supplementary:**
- Correlation/VIF matrix; elastic-net paths.
- Posterior diagnostics (R-hat, ESS, traces).
- Expanded maps, residuals, per-site posteriors.

---

## 8) Supplementary Information (SI) Structure
1. Detailed math: Eqs 7–13; kriging/predictive-process; PC priors.
2. Diagnostics: VIF tables, elastic-net, posterior traces, PPCs.
3. Expanded results: high-res maps; per-site summaries.
4. Data & code availability statement with repository DOI.

---

## 9) Style & Compliance Checklist
- δ values labeled in **‰ VSMOW** consistently.
- Keep only essential equations in main; math & proofs to SI.
- Use active voice and process framing.
- References trimmed to ~40–50; use GCA style.
- Figures: vector PDF/EPS or 600 dpi TIFF; consistent fonts; colorblind-friendly palettes.
- Prepare **Data Availability** and **Acknowledgments** statements.

---

## 10) Submission Package
- Main manuscript (Word or LaTeX in Elsevier format).
- Supplementary PDF (mathematical details, diagnostics).
- Figures as separate high-res files.
- Optional: code/data DOI (Zenodo/OSF) linked in SI.
