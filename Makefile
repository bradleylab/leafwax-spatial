# Makefile — leafwax GCA manuscript pipeline (Phase 5 W7, 2026-04-17).
#
# One entry point: `make pipeline`. From a clean checkout plus a populated
# results/c2_run_20260414/ mirror, this rebuilds every main-text figure
# and every manuscript-tables .tex, then runs the numeric audit.
#
# Rules:
#   tables    — regenerate all .tex (+ companion CSVs) from pipeline.
#   figures   — regenerate all main-figure PDFs into manuscript/figures/main_figs/.
#   audit     — check that each .tex's numbers still agree with its CSV.
#   pipeline  — tables + figures + audit.
#   clean-figures / clean-tables — remove generated artifacts.
#
# Notes:
#   - Analysis scripts are CWD-coupled. Recipes `cd` into the right
#     directory (repo root for extract_*.R / 5b_*.R, manuscript/figure_code/
#     for Figure_*.R) and copy outputs into manuscript/figures/main_figs/
#     after the run.
#   - Results are git-ignored; code + methods + TABLES.md/FIGURES.md are
#     the only tracked sources of truth. The per-artifact provenance docs
#     live next to their data under manuscript/.
#   - `R` is customizable via the `RSCRIPT` variable (`make RSCRIPT=…`).

RSCRIPT ?= Rscript

APRIL      := results/c2_run_20260414
PREP       := $(APRIL)/_prepared_data
MAIN_FIGS  := manuscript/figures/main_figs
TABLES_DIR := manuscript/tables
TABLE_CSVS := model_analysis/tables

MODELS := baseline baseline_sp baseline_env baseline_env_sp \
          baseline_veg baseline_veg_sp full full_sp \
          full_interact full_interact_sp c4_only_sp \
          elevation_only_sp elevation_c4_sp elevation_c4_interact_sp

# Per-model rds prereqs as lists.
ALL_DRAWS := $(foreach m,$(MODELS),$(APRIL)/$(m)/posterior_draws.rds)
ALL_LOO   := $(foreach m,$(MODELS),$(APRIL)/$(m)/loo.rds)
ALL_DIAG  := $(foreach m,$(MODELS),$(APRIL)/$(m)/diagnostics.rds)
ALL_STAND := $(foreach m,$(MODELS),$(PREP)/stan_data_$(m).rds)

SEDIMENT_RDS := $(APRIL)/3_sediment_ready_for_modeling.rds
TABLE_HELPERS := manuscript/table_code/table_helpers.R
POST_HELPERS_ROOT := scripts/posterior_helpers.R
POST_HELPERS_FIG  := manuscript/figure_code/posterior_helpers.R
COMMON_FIG_FUNS   := manuscript/figure_code/common_functions.R

# ─────────────────────────────────────────────────────────────────────────────
# Tables
# ─────────────────────────────────────────────────────────────────────────────

# Tables 1, 2, 3, 4 all come from extract_full_model_analysis.R in one pass.
# They share the same inputs; stamp Table 1 as the make target and include
# the rest as order-only side effects.
# Multi-output recipes are written with one canonical target carrying the
# recipe; every other output of the same script lists the canonical as a
# prereq. This is the GNU Make 3.81 idiom (macOS default). If your make
# is >= 4.3 you can collapse these pairs back to a grouped `&:` rule.

TABLE1_CANON := $(TABLES_DIR)/table1_model_performance.tex
TABLES_FULL_OTHER := $(TABLES_DIR)/table2_global_params_body.tex \
                     $(TABLES_DIR)/table3_variance_decomp.tex \
                     $(TABLES_DIR)/table4_environmental.tex
TABLES_FULL := $(TABLE1_CANON) $(TABLES_FULL_OTHER)

$(TABLE1_CANON): extract_full_model_analysis.R \
                 $(POST_HELPERS_ROOT) $(TABLE_HELPERS) \
                 $(ALL_DRAWS) $(ALL_LOO) $(ALL_DIAG) $(ALL_STAND)
	$(RSCRIPT) extract_full_model_analysis.R

$(TABLES_FULL_OTHER): $(TABLE1_CANON)

# Table S2 (body + companion CSV) from 5b_overfitting_diagnostics.R.
S2_BODY := $(TABLES_DIR)/Table_S2_regional_performance_body.tex
S2_CSV  := $(TABLES_DIR)/Table_S2_regional_performance.csv

$(S2_BODY): 5b_overfitting_diagnostics.R \
            $(POST_HELPERS_ROOT) $(TABLE_HELPERS) \
            $(ALL_DRAWS) $(ALL_LOO) $(ALL_STAND)
	$(RSCRIPT) 5b_overfitting_diagnostics.R

$(S2_CSV): $(S2_BODY)

TABLES_ALL := $(TABLES_FULL) $(S2_BODY)

.PHONY: tables
tables: $(TABLES_ALL)

# ─────────────────────────────────────────────────────────────────────────────
# Figures
# ─────────────────────────────────────────────────────────────────────────────

# Each figure script writes directly into manuscript/figures/main_figs/
# (via ../figures/main_figs/ relative to its cwd, or the full path when
# the script runs from the repo root). No intermediate scratch files —
# the target path is the output path.

# Figure 1 — OLS regression panels (canonical target = panel A).
FIG1A := $(MAIN_FIGS)/Figure1a_point_fitting.pdf
FIG1_OTHER := $(MAIN_FIGS)/Figure1b_10km_scale.pdf \
              $(MAIN_FIGS)/Figure1c_equal_weights.pdf \
              $(MAIN_FIGS)/Figure1d_bayesian_fitted.pdf
FIG1_OUT := $(FIG1A) $(FIG1_OTHER)

$(FIG1A): manuscript/figure_code/Figure_01_ols_regression.R \
          $(COMMON_FIG_FUNS) $(POST_HELPERS_FIG) \
          $(APRIL)/baseline/posterior_draws.rds \
          $(PREP)/stan_data_baseline.rds \
          $(SEDIMENT_RDS)
	mkdir -p $(MAIN_FIGS)
	cd manuscript/figure_code && $(RSCRIPT) Figure_01_ols_regression.R

$(FIG1_OTHER): $(FIG1A)

# Figure 2 — environmental variables. Runs from repo root because the
# script hardcodes repo-root-relative raster paths; emits directly into
# manuscript/figures/main_figs/.
$(MAIN_FIGS)/Figure_02_all_environmental_variables.png: \
    manuscript/figure_code/Figure_02_all_environmental_variables_noborders.R \
    $(SEDIMENT_RDS)
	mkdir -p $(MAIN_FIGS)
	$(RSCRIPT) manuscript/figure_code/Figure_02_all_environmental_variables_noborders.R

# Figure 3 — spatial confounding. Reads baseline + baseline_sp draws.
FIG3_PDF := $(MAIN_FIGS)/figure_03_spatial_confounding.pdf
FIG3_PNG := $(MAIN_FIGS)/figure_03_spatial_confounding.png
FIG3_OUT := $(FIG3_PDF) $(FIG3_PNG)

$(FIG3_PDF): manuscript/figure_code/figure_03_spatial_confounding_revised.R \
             $(POST_HELPERS_FIG) \
             $(APRIL)/baseline/posterior_draws.rds \
             $(APRIL)/baseline_sp/posterior_draws.rds \
             $(PREP)/stan_data_baseline.rds \
             $(PREP)/stan_data_baseline_sp.rds \
             $(SEDIMENT_RDS)
	mkdir -p $(MAIN_FIGS)
	cd manuscript/figure_code && $(RSCRIPT) figure_03_spatial_confounding_revised.R

$(FIG3_PNG): $(FIG3_PDF)

# Figure 4 — spatial maps (full_sp intercept / slope GP fields).
FIG4_PDF := $(MAIN_FIGS)/Figure_04.pdf
FIG4_PNG := $(MAIN_FIGS)/Figure_04.png
FIG4_OUT := $(FIG4_PDF) $(FIG4_PNG)

$(FIG4_PDF): manuscript/figure_code/Figure_04_spatial_maps.R \
             $(POST_HELPERS_FIG) \
             $(APRIL)/full_sp/posterior_draws.rds \
             $(PREP)/stan_data_full_sp.rds \
             $(SEDIMENT_RDS)
	mkdir -p $(MAIN_FIGS)
	cd manuscript/figure_code && $(RSCRIPT) Figure_04_spatial_maps.R

$(FIG4_PNG): $(FIG4_PDF)

# Figure 5 — detection thresholds. Reads full_sp and baseline draws.
$(MAIN_FIGS)/Figure_05_detection_thresholds.png: \
    manuscript/figure_code/Figure_05_detection_thresholds.R \
    $(POST_HELPERS_FIG) \
    $(APRIL)/full_sp/posterior_draws.rds \
    $(APRIL)/baseline/posterior_draws.rds \
    $(PREP)/stan_data_full_sp.rds \
    $(PREP)/stan_data_baseline.rds
	mkdir -p $(MAIN_FIGS)
	cd manuscript/figure_code && $(RSCRIPT) Figure_05_detection_thresholds.R

FIGURES_ALL := $(FIG1_OUT) \
               $(MAIN_FIGS)/Figure_02_all_environmental_variables.png \
               $(FIG3_OUT) \
               $(FIG4_OUT) \
               $(MAIN_FIGS)/Figure_05_detection_thresholds.png

# ─────────────────────────────────────────────────────────────────────────────
# Supplement figures — each script writes directly to
# manuscript/figures/supplement_figs/ (via ../figures/supplement_figs/
# relative to its cwd). No scratch files.
# ─────────────────────────────────────────────────────────────────────────────

SUP_FIGS := manuscript/figures/supplement_figs

# S1: purely analytic (no rds inputs)
$(SUP_FIGS)/Figure_S1_spatial_weighting.pdf: \
    manuscript/figure_code/Figure_S1_spatial_weighting.R
	mkdir -p $(SUP_FIGS)
	cd manuscript/figure_code && $(RSCRIPT) Figure_S1_spatial_weighting.R

# S2: pairwise correlations (reads sediment rds via posterior_helpers)
$(SUP_FIGS)/Figure_S2_pairwise_correlations.pdf: \
    manuscript/figure_code/Figure_S2_pairwise_correlations.R \
    $(POST_HELPERS_FIG) $(SEDIMENT_RDS)
	mkdir -p $(SUP_FIGS)
	cd manuscript/figure_code && $(RSCRIPT) Figure_S2_pairwise_correlations.R

# S3: global OLS residuals map (sediment rds only)
$(SUP_FIGS)/Figure_S3_residuals.pdf: \
    manuscript/figure_code/Figure_S3_residuals.R \
    $(SEDIMENT_RDS)
	mkdir -p $(SUP_FIGS)
	cd manuscript/figure_code && $(RSCRIPT) Figure_S3_residuals.R

# S4: OLS R² vs integration scale (baseline stan_data + diagnostics)
$(SUP_FIGS)/Figure_S4_ols_spatial_scale.pdf: \
    manuscript/figure_code/Figure_S4_ols_spatial_scale.R \
    $(POST_HELPERS_FIG) $(SEDIMENT_RDS) \
    $(PREP)/stan_data_baseline.rds \
    $(APRIL)/baseline/diagnostics.rds
	mkdir -p $(SUP_FIGS)
	cd manuscript/figure_code && $(RSCRIPT) Figure_S4_ols_spatial_scale.R

# S5: elevation effects across all `include_elevation == 1` models
# (reads per-model stan_data + diagnostics.rds for beta_elev_bspline)
$(SUP_FIGS)/Figure_S5_elevation.pdf: \
    manuscript/figure_code/Figure_S5_elevation.R \
    $(POST_HELPERS_FIG) $(SEDIMENT_RDS) $(ALL_DIAG) $(ALL_STAND)
	mkdir -p $(SUP_FIGS)
	cd manuscript/figure_code && $(RSCRIPT) Figure_S5_elevation.R

SUP_FIGURES_ALL := $(SUP_FIGS)/Figure_S1_spatial_weighting.pdf \
                   $(SUP_FIGS)/Figure_S2_pairwise_correlations.pdf \
                   $(SUP_FIGS)/Figure_S3_residuals.pdf \
                   $(SUP_FIGS)/Figure_S4_ols_spatial_scale.pdf \
                   $(SUP_FIGS)/Figure_S5_elevation.pdf

.PHONY: figures supplement-figures
figures: $(FIGURES_ALL)
supplement-figures: $(SUP_FIGURES_ALL)

# ─────────────────────────────────────────────────────────────────────────────
# Numeric audit + full pipeline
# ─────────────────────────────────────────────────────────────────────────────

.PHONY: audit
audit: tables
	$(RSCRIPT) manuscript/table_code/check_tables_numeric.R

.PHONY: pipeline
pipeline: tables figures supplement-figures audit

# ─────────────────────────────────────────────────────────────────────────────
# Cleanup
# ─────────────────────────────────────────────────────────────────────────────

.PHONY: clean-figures clean-tables clean
clean-figures:
	rm -rf $(MAIN_FIGS) $(SUP_FIGS)

clean-tables:
	rm -f $(TABLES_ALL) \
	      $(TABLES_DIR)/Table_S2_regional_performance.csv \
	      $(TABLE_CSVS)/table1_model_performance.csv \
	      $(TABLE_CSVS)/table2_global_params.csv \
	      $(TABLE_CSVS)/table3_variance_decomp.csv \
	      $(TABLE_CSVS)/table4_environmental.csv

clean: clean-figures clean-tables
