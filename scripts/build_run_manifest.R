#!/usr/bin/env Rscript
# build_run_manifest.R
#
# Run manifest for all 17 configured fits: 9 chordal-distance spatial models,
# 5 non-spatial models, and 3 *_rfoff sensitivity variants.
#
# The manifest is an audit record. It validates in two stages and fails closed:
#   Structural validation (ignores --allow-incomplete):
#     * config readable + EXACTLY the 17 expected model_configs;
#     * fit-time provenance record (scripts/capture_fit_provenance.R, captured
#       inside the SIF on C2) present and readable — it cannot be reconstructed.
#   Fit validation (gated by --allow-incomplete for interim snapshots):
#     * posterior readable with a valid 3-D draw structure;
#     * draw dims match the fitted config (chains, retained iterations);
#     * the COMPLETE conditional variable set 4c writes for that config is present;
#     * observation-/scale-/knot-indexed variables have the stan_data dimensions;
#     * elevation models carry beta_elev; prepared stan_data readable (+ km-reg
#       fields for spatial); fit_status.rds + DONE present.
# The written .rds records status ("authoritative"/"interim") and any problems.
#
#   Rscript scripts/build_run_manifest.R \
#     --chordal-output   results/c2_run_<date>_chordal/model_output \
#     --chordal-prepared results/c2_run_<date>_chordal/_prepared_data \
#     --fit-provenance   results/c2_run_<date>_chordal/fit_provenance.rds \
#     --out              results/c2_run_<date>_chordal/RUN_MANIFEST_chordal

suppressWarnings(suppressMessages(library(yaml)))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) { i <- match(flag, args); if (!is.na(i) && i < length(args)) args[[i + 1]] else default }
has_flag <- function(flag) flag %in% args

chordal_output   <- get_arg("--chordal-output",   "model_output")
chordal_prepared <- get_arg("--chordal-prepared", "prepared_data")
config_path      <- get_arg("--config",           "config.yaml")
fit_prov_path    <- get_arg("--fit-provenance",   file.path(dirname(chordal_output), "fit_provenance.rds"))
out_base         <- get_arg("--out",              "results/RUN_MANIFEST_chordal")
allow_incomplete <- has_flag("--allow-incomplete")

# Model taxonomy + provenance completeness/drift come from the shared module, so
# the manifest cannot disagree with capture about which fields/hashes are essential.
source("scripts/provenance_common.R")   # CHORDAL_SPATIAL / NONSPATIAL_REFIT / SENSITIVITY / EXPECTED17
is_spatial <- function(m) grepl("_sp", m)

die <- function(...) { cat("ERROR:", ..., "\n"); quit(status = 1) }
md5_or_na <- function(p) if (!is.na(p) && file.exists(p)) unname(tools::md5sum(p)) else NA_character_

# ── Structural integrity — hard fail (ignores --allow-incomplete) ─────────────
if (!file.exists(config_path)) die("config not found at '", config_path, "'.")
cfg <- tryCatch(yaml::read_yaml(config_path), error = function(e) NULL)
if (is.null(cfg) || is.null(cfg$model_configs)) die("config '", config_path, "' unreadable or has no model_configs.")
if (!setequal(names(cfg$model_configs), EXPECTED17))
  die("config model_configs != the 17 expected. Missing: ",
      paste(setdiff(EXPECTED17, names(cfg$model_configs)), collapse = ", "),
      " | Unexpected: ", paste(setdiff(names(cfg$model_configs), EXPECTED17), collapse = ", "))
warmup_ratio <- cfg$warmup_ratio
if (is.null(warmup_ratio)) die("config has no warmup_ratio — cannot compute expected retained iterations.")

fit_prov <- if (file.exists(fit_prov_path)) tryCatch(readRDS(fit_prov_path), error = function(e) NULL) else NULL
if (is.null(fit_prov))
  die("fit-time provenance not found/readable at '", fit_prov_path, "'. Run scripts/capture_fit_provenance.R ",
      "inside the SIF on C2 (job_prep.sh does this) — it cannot be reconstructed. This is required for an ",
      "authoritative manifest.")
# Validate the record's CONTENTS with the SAME complete definition capture uses
# (all versions + SIF md5 + every code/input hash + per-model prepared hashes for
# all 17), md5-format-checked. A record with any NA/invalid essential field is not
# authoritative (tier-1 hard fail).
prov_bad <- provenance_missing_fields(fit_prov, EXPECTED17)
if (length(prov_bad))
  die("fit-time provenance is incomplete/invalid (missing or bad-md5 fields): ", paste(prov_bad, collapse = ", "),
      ". Re-capture with scripts/capture_fit_provenance.R.")

# ── The conditional variable set 4c_fit_models.R writes, re-derived from config ──
# SOURCE OF TRUTH: 4c_fit_models.R (params_to_check + draws_to_save). Kept in sync
# here as the independent audit check — if 4c's saved set changes, update this too.
expected_base_vars <- function(mc) {
  v <- c("beta_0", "beta_oipc", "sigma", "lambda_decay", "effective_scale_km",
         "mu", "d2H_rep", "log_lik", "scale_weights")
  if (isTRUE(mc$include_precip))    v <- c(v, "beta_precip")
  if (isTRUE(mc$include_elevation)) v <- c(v, "beta_elev_bspline", "tau_elev_bspline")
  if (isTRUE(mc$include_c4)) { v <- c(v, "beta_c4"); if (isTRUE(mc$include_veg_interactions)) v <- c(v, "beta_oipc_x_c4") }
  if (isTRUE(mc$include_pft)) v <- c(v, "beta_tree", "beta_shrub", "beta_grass",
                                        "beta_oipc_x_tree", "beta_oipc_x_shrub", "beta_oipc_x_grass")
  if (isTRUE(mc$include_gp)) v <- c(v, "ls_intercept_km", "ls_slope_km", "sigma_intercept_spatial",
                                       "sigma_slope_spatial", "alpha_spatial", "beta_oipc_spatial",
                                       "z_intercept_spatial", "z_slope_spatial")
  v
}
var_present <- function(base, vn) any(grepl(sprintf("^%s(\\[|\\.|$)", base), vn))
var_count   <- function(base, vn) sum(grepl(sprintf("^%s(\\[|\\.)", base), vn))

# ── Provenance (manifest-builder env, plus the captured fit-time record) ─────────
git_head <- tryCatch(sub("\\s+$", "", system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)), error = function(e) NA_character_)
if (length(git_head) != 1 || !nzchar(git_head)) git_head <- NA_character_
git_dirty <- tryCatch(length(system2("git", c("status", "--porcelain"), stdout = TRUE, stderr = FALSE)) > 0, error = function(e) NA)
reg_radii <- list(density_radius_km = cfg$gp_regularization$density_radius_km,
                  oipc_range_radius_km = cfg$gp_regularization$oipc_range_radius_km)
provenance <- list(
  built_at = Sys.time(), git_head = git_head, git_dirty = git_dirty,
  manifest_builder_r_version = R.version.string,
  manifest_builder_cmdstan_version = tryCatch(cmdstanr::cmdstan_version(), error = function(e) NA_character_),
  stan_seed = cfg$stan_seed,
  fit_time = fit_prov,   # authoritative: captured inside the SIF on C2
  reg_neighborhood_radii_km = reg_radii)
# The real fitting CmdStan comes from the captured record, not a hardcode.
fitting_cmdstan <- tryCatch(fit_prov$environment$cmdstan_version, error = function(e) NA_character_)

rng <- function(x) if (is.null(x) || !length(x)) "—" else paste0("[", round(min(x), 3), ", ", round(max(x), 3), "]")

inspect_prepared <- function(model) {
  f <- file.path(chordal_prepared, paste0("stan_data_", model, ".rds"))
  if (!file.exists(f)) return(list(exists = FALSE, readable = FALSE, reg_ok = FALSE, md5 = NA_character_,
                                   N = NA, n_scales = NA, n_pp_knots = NA, reg = NULL))
  sd <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(sd)) return(list(exists = TRUE, readable = FALSE, reg_ok = FALSE, md5 = md5_or_na(f),
                               N = NA, n_scales = NA, n_pp_knots = NA, reg = NULL))
  req <- c("tau_spatial_slope", "tau_spatial_intercept", "oipc_range_at_knots", "knot_data_density", "ls_log_lower", "ls_log_upper")
  reg_ok <- all(vapply(req, function(k) !is.null(sd[[k]]), logical(1)))
  reg <- if (!is.null(sd$tau_spatial_slope)) list(
    n_pp_knots = sd$n_pp_knots, apply_range_factor = if (!is.null(sd$apply_range_factor)) as.logical(sd$apply_range_factor) else NA,
    knot_data_density = sd$knot_data_density, oipc_range_at_knots = sd$oipc_range_at_knots,
    tau_spatial_slope = sd$tau_spatial_slope, tau_spatial_intercept = sd$tau_spatial_intercept,
    range_factor_values = sd$range_factor_values) else NULL
  list(exists = TRUE, readable = TRUE, reg_ok = reg_ok, md5 = md5_or_na(f),
       N = sd$N, n_scales = sd$n_scales, n_pp_knots = sd$n_pp_knots, reg = reg)
}

# Full posterior validation: readability, draw dims vs config, complete variable
# set vs config, and indexed-variable dimensions vs stan_data. One read per model.
validate_posterior <- function(pd, mc, dims) {
  exp_chains <- mc$chains
  exp_iter   <- if (!is.null(mc$iter)) floor(mc$iter * (1 - warmup_ratio)) else NA
  base_ok <- list(exists = FALSE, readable = FALSE, struct_ok = FALSE, dims_ok = FALSE,
                  vars_ok = FALSE, indexed_ok = FALSE, has_elev = NA,
                  n_iter = NA, n_chains = NA, missing_vars = NA_character_, dim_notes = NA_character_)
  if (!file.exists(pd)) return(base_ok)
  x <- tryCatch(readRDS(pd), error = function(e) NULL)
  if (is.null(x)) return(modifyList(base_ok, list(exists = TRUE)))
  d  <- tryCatch(dim(x), error = function(e) NULL)
  vn <- tryCatch(dimnames(x)[[length(dim(x))]], error = function(e) NULL)
  struct_ok <- !is.null(d) && length(d) == 3 && all(d > 0) && !is.null(vn) && all(c("beta_0", "beta_oipc") %in% vn)
  n_iter <- if (length(d) >= 1) d[1] else NA; n_chains <- if (length(d) >= 2) d[2] else NA
  dims_ok <- isTRUE(n_chains == exp_chains) && (is.na(exp_iter) || isTRUE(n_iter == exp_iter))
  exp_vars <- expected_base_vars(mc)
  missing <- exp_vars[!vapply(exp_vars, var_present, logical(1), vn = vn)]
  vars_ok <- length(missing) == 0
  # Indexed-variable dimension checks against stan_data.
  notes <- c()
  chk <- function(base, expect, label) {
    if (!var_present(base, vn)) return(invisible())
    got <- var_count(base, vn); if (got == 0) got <- 1  # scalar written as bare name
    if (!isTRUE(got == expect)) notes <<- c(notes, sprintf("%s: %d != %s(%s)", base, got, label, expect))
  }
  # Observation-indexed [N] (4d:266-267 declare alpha_spatial/beta_oipc_spatial as
  # vector[N] — per-location, NOT knot-indexed). Knot-indexed [n_pp_knots] is only
  # z_intercept_spatial / z_slope_spatial (4d:205-206).
  n_indexed <- c("mu", "d2H_rep", "log_lik")
  if (isTRUE(mc$include_gp)) n_indexed <- c(n_indexed, "alpha_spatial", "beta_oipc_spatial")
  if (!is.na(dims$N)) for (b in n_indexed) chk(b, dims$N, "N")
  if (!is.na(dims$n_scales)) chk("scale_weights", dims$n_scales, "n_scales")
  if (isTRUE(mc$include_gp) && !is.na(dims$n_pp_knots))
    for (b in c("z_intercept_spatial", "z_slope_spatial")) chk(b, dims$n_pp_knots, "n_pp_knots")
  indexed_ok <- length(notes) == 0
  list(exists = TRUE, readable = TRUE, struct_ok = struct_ok, dims_ok = dims_ok, vars_ok = vars_ok,
       indexed_ok = indexed_ok, has_elev = if (!is.null(vn)) any(grepl("^beta_elev", vn)) else NA,
       n_iter = n_iter, n_chains = n_chains,
       missing_vars = if (length(missing)) paste(missing, collapse = ",") else "",
       dim_notes = if (length(notes)) paste(notes, collapse = "; ") else "")
}

runtime_info <- function(model) {
  f <- file.path(chordal_output, model, "runtime_info.rds")
  ri <- if (file.exists(f)) tryCatch(readRDS(f), error = function(e) NULL) else NULL
  list(elapsed_mins = if (!is.null(ri)) ri$elapsed_mins else NA)
}

# fit_status.rds must be a READABLE list for THIS model with status "completed"
# and core_complete TRUE — not merely a file that exists (a plaintext file passes
# file.exists). DONE + core artifacts must be newer than stan_data (the
# aggregator's freshness rule), so a stale DONE cannot pass here either.
validate_fit_status <- function(model) {
  f <- file.path(chordal_output, model, "fit_status.rds")
  if (!file.exists(f)) return(FALSE)
  st <- tryCatch(readRDS(f), error = function(e) NULL)
  is.list(st) && identical(st$model, model) && identical(st$status, "completed") && isTRUE(st$core_complete)
}
done_fresh <- function(model) {
  mdir <- file.path(chordal_output, model)
  sd <- file.path(chordal_prepared, paste0("stan_data_", model, ".rds"))
  paths <- file.path(mdir, c("DONE", "fit.rds", "diagnostics.rds", "posterior_draws.rds"))
  if (!file.exists(sd) || !all(file.exists(paths))) return(FALSE)
  all(file.mtime(paths) > file.mtime(sd))
}

rows <- list()
add_row <- function(model, role) {
  mc <- cfg$model_configs[[model]]
  pd <- file.path(chordal_output, model, "posterior_draws.rds")
  prep <- inspect_prepared(model)
  vp <- validate_posterior(pd, mc, prep)
  needs_elev <- isTRUE(mc$include_elevation)
  rows[[length(rows) + 1]] <<- list(
    model = model, role = role, spatial_metric = if (is_spatial(model)) "chordal" else "n/a (non-spatial)",
    fits_elevation = needs_elev,
    posterior_exists = vp$exists, posterior_readable = vp$readable, posterior_struct_ok = vp$struct_ok,
    dims_ok = vp$dims_ok, vars_ok = vp$vars_ok, indexed_ok = vp$indexed_ok,
    missing_vars = vp$missing_vars, dim_notes = vp$dim_notes,
    n_iter = vp$n_iter, n_chains = vp$n_chains, posterior_md5 = md5_or_na(pd),
    elevation_present = if (needs_elev) isTRUE(vp$has_elev) else NA,
    stan_data_md5 = prep$md5, prepared_readable = prep$readable, reg_ok = prep$reg_ok,
    fit_status_ok = validate_fit_status(model),
    done_fresh = done_fresh(model),
    runtime = runtime_info(model), reg = prep$reg)
}
for (m in CHORDAL_SPATIAL)  add_row(m, "spatial")
for (m in NONSPATIAL_REFIT) add_row(m, "nonspatial")
for (m in SENSITIVITY)      add_row(m, "sensitivity")

# ── Fit completeness — gated by --allow-incomplete ───────────────────────────
mnames <- function(pred) vapply(rows[vapply(rows, pred, logical(1))], `[[`, character(1), "model")
fully_valid <- function(r) isTRUE(r$posterior_struct_ok) && isTRUE(r$dims_ok) && isTRUE(r$vars_ok) &&
  isTRUE(r$indexed_ok) && isTRUE(r$prepared_readable) && (!is_spatial(r$model) || isTRUE(r$reg_ok)) &&
  (!isTRUE(r$fits_elevation) || isTRUE(r$elevation_present)) && isTRUE(r$fit_status_ok) && isTRUE(r$done_fresh)
n_ok <- sum(vapply(rows, fully_valid, logical(1)))

# Config/input/prepared drift: the artifacts were fit under the code/config/inputs
# the provenance record hashed. prov_is_drift is fail-CLOSED (a missing recorded
# hash OR missing current file OR mismatch all count as drift). Per-model stan_data
# is compared against its fit-time captured hash too.
drift_path <- function(now_path, rec) prov_is_drift(md5_or_na(now_path), rec)
drift <- c(
  if (drift_path(config_path, fit_prov$code_md5$config)) "config.yaml",
  if (drift_path("4a_spatial_functions.R", fit_prov$code_md5$spatial_functions)) "4a_spatial_functions.R",
  if (drift_path("4b_stan_prep.R", fit_prov$code_md5$stan_prep)) "4b_stan_prep.R",
  if (drift_path("4c_fit_models.R", fit_prov$code_md5$fit_models)) "4c_fit_models.R",
  if (drift_path("4d_leaf_wax_spatial_model.stan", fit_prov$code_md5$stan_model)) "4d_leaf_wax_spatial_model.stan",
  if (drift_path("input_data/leafwax_d2h_c29_calibration_v1.csv", fit_prov$input_md5$calibration_csv)) "calibration_csv",
  if (drift_path("results/3_sediment_ready_for_modeling.rds", fit_prov$input_md5$sediment_prepared_rds)) "sediment_prepared_rds")
# Per-model manifested stan_data md5 vs the fit-time captured stan_data md5.
stan_data_drift <- vapply(EXPECTED17, function(m) {
  cap <- fit_prov$prepared_md5[[m]]$stan_data
  now <- md5_or_na(file.path(chordal_prepared, paste0("stan_data_", m, ".rds")))
  prov_is_drift(now, cap)
}, logical(1))
if (any(stan_data_drift)) drift <- c(drift, paste0("stan_data[", paste(EXPECTED17[stan_data_drift], collapse = ","), "]"))

problems <- c(
  { b <- mnames(function(r) !isTRUE(r$posterior_exists));                              if (length(b)) paste0("missing posteriors: ", paste(b, collapse = ", ")) },
  { b <- mnames(function(r) isTRUE(r$posterior_exists) && !isTRUE(r$posterior_readable)); if (length(b)) paste0("unreadable posteriors: ", paste(b, collapse = ", ")) },
  { b <- mnames(function(r) isTRUE(r$posterior_readable) && !isTRUE(r$posterior_struct_ok)); if (length(b)) paste0("bad draw structure: ", paste(b, collapse = ", ")) },
  { b <- mnames(function(r) isTRUE(r$posterior_struct_ok) && !isTRUE(r$dims_ok));      if (length(b)) paste0("draw dims != fitted config (chains/iterations): ", paste(b, collapse = ", ")) },
  { b <- mnames(function(r) isTRUE(r$posterior_struct_ok) && !isTRUE(r$vars_ok));      if (length(b)) paste0("missing expected variables: ", paste(sprintf("%s{%s}", b, vapply(rows[match(b, vapply(rows, `[[`, "", "model"))], `[[`, "", "missing_vars")), collapse = "; ")) },
  { b <- mnames(function(r) isTRUE(r$posterior_struct_ok) && !isTRUE(r$indexed_ok));   if (length(b)) paste0("indexed-variable dim mismatch: ", paste(b, collapse = ", ")) },
  { b <- mnames(function(r) isTRUE(r$fits_elevation) && !isTRUE(r$elevation_present)); if (length(b)) paste0("elevation models NOT carrying beta_elev: ", paste(b, collapse = ", ")) },
  { b <- mnames(function(r) !isTRUE(r$prepared_readable));                             if (length(b)) paste0("unreadable/missing prepared stan_data: ", paste(b, collapse = ", ")) },
  { b <- mnames(function(r) is_spatial(r$model) && !isTRUE(r$reg_ok));                 if (length(b)) paste0("spatial models missing km-reg fields: ", paste(b, collapse = ", ")) },
  { b <- mnames(function(r) !isTRUE(r$fit_status_ok));                                 if (length(b)) paste0("fit_status.rds missing/invalid (need readable list, model match, status=completed, core_complete): ", paste(b, collapse = ", ")) },
  { b <- mnames(function(r) !isTRUE(r$done_fresh));                                    if (length(b)) paste0("DONE/core artifacts missing or not newer than stan_data: ", paste(b, collapse = ", ")) },
  if (length(drift)) paste0("working tree drifted from fit-time provenance (md5 mismatch): ", paste(drift, collapse = ", ")))

status <- if (length(problems)) "interim" else "authoritative"
cat(n_ok, "of", length(rows), "models fully valid. status =", status, "\n")
if (length(problems)) for (p in problems) cat("  -", p, "\n")
if (length(problems) && !allow_incomplete) {
  cat("\n✗ MANIFEST INVALID — not writing an authoritative record.\n")
  cat("  Fix the above, or pass --allow-incomplete for an interim (non-authoritative) snapshot.\n")
  quit(status = 1)
}

manifest <- list(
  status = status, validation_problems = problems, allow_incomplete = allow_incomplete,
  provenance = provenance,
  built_role_map = list(spatial = CHORDAL_SPATIAL, nonspatial = NONSPATIAL_REFIT, sensitivity = SENSITIVITY),
  chordal_output = chordal_output, chordal_prepared = chordal_prepared,
  fit_provenance_path = fit_prov_path, models = rows)

# ── Markdown ────────────────────────────────────────────────────────────────────
MD <- character(0); wl <- function(...) MD[[length(MD) + 1]] <<- paste0(...)
wl("# Model run manifest — status: ", toupper(status)); wl("")
if (length(problems)) wl("> **INTERIM / NON-AUTHORITATIVE** (--allow-incomplete): ", paste(problems, collapse = "; "))
wl("The run contains all 17 model configurations: 9 spatial, 5 non-spatial, and 3 range-factor sensitivity fits.")
wl(""); wl("## Provenance")
wl("- built_at: ", format(provenance$built_at))
wl("- manifest-builder: R ", provenance$manifest_builder_r_version, "; cmdstan ",
   ifelse(is.na(provenance$manifest_builder_cmdstan_version), "—", provenance$manifest_builder_cmdstan_version), " (NOT the fitting env)")
wl("- **fitting env (captured in-SIF on C2):** cmdstan ", ifelse(is.na(fitting_cmdstan), "—", fitting_cmdstan),
   "; cmdstanr ", fit_prov$environment$packages$cmdstanr, "; posterior ", fit_prov$environment$packages$posterior,
   "; R ", fit_prov$environment$r_version)
wl("- fit-time SIF md5: ", ifelse(is.null(fit_prov$sif$md5) || is.na(fit_prov$sif$md5), "—", substr(fit_prov$sif$md5, 1, 12)),
   " | fit-time code(4c) md5: ", substr(fit_prov$code_md5$fit_models, 1, 12),
   " | fit-time input md5: ", substr(fit_prov$input_md5$calibration_csv, 1, 12))
wl("- fit-time git: ", ifelse(is.null(fit_prov$git$head) || is.na(fit_prov$git$head), "—", fit_prov$git$head),
   if (isTRUE(fit_prov$git$dirty)) " (dirty)" else "")
wl("- stan_seed: ", ifelse(is.null(provenance$stan_seed), "—", provenance$stan_seed),
   " | reg radii (km): density=", reg_radii$density_radius_km, " oipc_range=", reg_radii$oipc_range_radius_km)
wl("")
wl("| model | role | struct | dims✓ | vars✓ | idx✓ | iter×chain | md5 | stan_data md5 | fits_elev | elev✓ | reg✓ | status/DONE |")
wl("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
tri <- function(x) if (is.na(x)) "—" else if (isTRUE(x)) "yes" else "NO"
for (r in rows) {
  wl("| ", r$model, " | ", r$role, " | ", tri(r$posterior_struct_ok), " | ", tri(r$dims_ok), " | ", tri(r$vars_ok),
     " | ", tri(r$indexed_ok), " | ", ifelse(is.na(r$n_iter), "—", paste0(r$n_iter, "×", r$n_chains)),
     " | ", ifelse(is.na(r$posterior_md5), "—", substr(r$posterior_md5, 1, 10)),
     " | ", ifelse(is.na(r$stan_data_md5), "—", substr(r$stan_data_md5, 1, 10)),
     " | ", if (isTRUE(r$fits_elevation)) "yes" else "no", " | ", tri(r$elevation_present),
     " | ", if (is_spatial(r$model)) tri(r$reg_ok) else "n/a",
     " | ", if (isTRUE(r$fit_status_ok) && isTRUE(r$done_fresh)) "yes" else "NO", " |")
}
wl("")
wl("Detailed km regularization diagnostics per spatial model (full per-knot vectors in the .rds; ranges here):")
for (r in rows) {
  if (is.null(r$reg)) next
  d <- r$reg; wl(""); wl("### ", r$model)
  wl("- knot_data_density range: ", rng(d$knot_data_density))
  wl("- oipc_range_at_knots range: ", rng(d$oipc_range_at_knots))
  wl("- tau_spatial_intercept range: ", rng(d$tau_spatial_intercept))
  wl("- tau_spatial_slope range: ", rng(d$tau_spatial_slope))
  wl("- range_factor values: ",
     if (!is.null(d$range_factor_values)) paste(sort(unique(round(d$range_factor_values, 2))), collapse = "/") else "—",
     " (apply_range_factor = ", d$apply_range_factor, ")")
}

rds_tmp <- paste0(out_base, ".rds.tmp"); rds_final <- paste0(out_base, ".rds")
md_tmp  <- paste0(out_base, ".md.tmp");  md_final  <- paste0(out_base, ".md")
saveRDS(manifest, rds_tmp)
if (!isTRUE(file.rename(rds_tmp, rds_final))) die("failed to write ", rds_final, " (rename failed).")
writeLines(MD, md_tmp)
if (!isTRUE(file.rename(md_tmp, md_final))) die("failed to write ", md_final, " (rename failed).")

cat("Wrote", rds_final, "and", md_final, " (status:", status, ")\n")
if (length(problems) && allow_incomplete)
  cat("NOTE: --allow-incomplete; wrote an INTERIM (non-authoritative) manifest; status recorded in the .rds.\n")
