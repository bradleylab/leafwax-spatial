# provenance_common.R
#
# Single source of truth for (a) the chordal-run model taxonomy and (b) what makes
# a fit-time provenance record COMPLETE and what counts as prep→fit DRIFT. Sourced
# by capture_fit_provenance.R, check_provenance_current.R, and build_run_manifest.R
# so the three cannot disagree about which fields and hashes are essential.

# ── Model taxonomy ────────────────────────────────────────────────────────────
CHORDAL_SPATIAL  <- c("baseline_sp", "baseline_veg_sp", "baseline_env_sp", "c4_only_sp",
                      "elevation_only_sp", "elevation_c4_sp", "elevation_c4_interact_sp",
                      "full_sp", "full_interact_sp")
NONSPATIAL_REFIT <- c("baseline", "baseline_veg", "baseline_env", "full", "full_interact")
SENSITIVITY      <- c("baseline_sp_rfoff", "baseline_env_sp_rfoff", "full_interact_sp_rfoff")
EXPECTED17       <- c(CHORDAL_SPATIAL, NONSPATIAL_REFIT, SENSITIVITY)

# ── Field / hash helpers ────────────────────────────────────────────────────────
prov_na_or_empty <- function(x) is.null(x) || length(x) == 0 || (length(x) == 1 && (is.na(x) || !nzchar(as.character(x))))
is_md5 <- function(x) is.character(x) && length(x) == 1 && grepl("^[0-9a-fA-F]{32}$", x)

# Drift is fail-closed: a missing recorded hash, a missing current file, or a
# mismatch all count as drift.
prov_is_drift <- function(now, rec) prov_na_or_empty(rec) || prov_na_or_empty(now) || !identical(now, rec)

# The COMPLETE essential-field set. Returns a character vector of missing/invalid
# field names (empty = complete). `models` is the set of models whose prepared
# hashes must be present (per-model stan_data + config). Every hash field must be
# a valid 32-hex md5; versions must be non-empty.
provenance_missing_fields <- function(p, models) {
  miss <- character(0)
  need <- function(name, val, md5 = FALSE) {
    if (prov_na_or_empty(val)) miss <<- c(miss, name)
    else if (md5 && !is_md5(val)) miss <<- c(miss, paste0(name, "(bad-md5)"))
  }
  need("environment$cmdstan_version", p$environment$cmdstan_version)
  need("packages$cmdstanr",  p$environment$packages$cmdstanr)
  need("packages$posterior", p$environment$packages$posterior)
  need("packages$loo",       p$environment$packages$loo)
  need("sif$md5", p$sif$md5, md5 = TRUE)
  for (k in c("spatial_functions", "stan_prep", "fit_models", "stan_model", "config"))
    need(paste0("code_md5$", k), p$code_md5[[k]], md5 = TRUE)
  for (k in c("calibration_csv", "sediment_prepared_rds"))
    need(paste0("input_md5$", k), p$input_md5[[k]], md5 = TRUE)
  for (mm in models) {
    need(paste0("prepared_md5$", mm, "$stan_data"), p$prepared_md5[[mm]]$stan_data, md5 = TRUE)
    need(paste0("prepared_md5$", mm, "$config"),    p$prepared_md5[[mm]]$config,    md5 = TRUE)
  }
  miss
}
