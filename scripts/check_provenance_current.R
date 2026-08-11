#!/usr/bin/env Rscript
# check_provenance_current.R
#
# Pre-fit gate: before a model is fit, confirm the captured fit-time
# provenance record is COMPLETE and still describes the CURRENT workspace — the
# files 4c actually reads. job_fit_chordal.sh runs this inside the SIF, and BEFORE
# the resume/skip guard, so no task fits (or is trusted as complete) without it.
# Guards:
#   (1) incomplete/corrupt record (shared completeness definition);
#   (2) drift in the shared code/config/modeling inputs since prep;
#   (3) drift in THIS model's prepared_data/stan_data_<m>.rds and config_<m>.rds
#       (the direct fit inputs — they can change while global hashes do not);
#   (4) the current host-side SIF md5 (passed in) vs the captured SIF md5.
# Exit 0 only if all pass; else quit(1).
#
#   Rscript scripts/check_provenance_current.R \
#     --prov results/fit_provenance.rds --model baseline_sp --current-sif-md5 <md5>

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NA_character_) {
  i <- match(flag, args); if (!is.na(i) && i < length(args)) args[[i + 1]] else default
}
prov_path <- get_arg("--prov", "results/fit_provenance.rds")
model     <- get_arg("--model")
cur_sif   <- get_arg("--current-sif-md5")

source("scripts/provenance_common.R")
die <- function(...) { cat("ERROR:", ..., "\n"); quit(status = 1) }
if (prov_na_or_empty(model)) die("--model is required.")
if (!file.exists(prov_path)) die("fit provenance not found at '", prov_path, "'. Run job_prep.sh first.")
p <- tryCatch(readRDS(prov_path), error = function(e) NULL)
if (is.null(p)) die("fit provenance at '", prov_path, "' is unreadable.")

# (1) record complete (validated against the full expected model set).
miss <- provenance_missing_fields(p, EXPECTED17)
if (length(miss)) die("provenance record incomplete/invalid — ", paste(miss, collapse = ", "))

md5_or_na <- function(f) if (file.exists(f)) unname(tools::md5sum(f)) else NA_character_

# (2) shared code/config/inputs drift.
drift <- character(0)
add_drift <- function(label, now, rec) if (prov_is_drift(now, rec))
  drift <<- c(drift, sprintf("%s (now=%s captured=%s)", label,
              ifelse(prov_na_or_empty(now), "MISSING", substr(now, 1, 12)),
              ifelse(prov_na_or_empty(rec), "MISSING", substr(rec, 1, 12))))
add_drift("4a_spatial_functions.R", md5_or_na("4a_spatial_functions.R"), p$code_md5$spatial_functions)
add_drift("4b_stan_prep.R",         md5_or_na("4b_stan_prep.R"),         p$code_md5$stan_prep)
add_drift("4c_fit_models.R",        md5_or_na("4c_fit_models.R"),        p$code_md5$fit_models)
add_drift("4d_leaf_wax_spatial_model.stan", md5_or_na("4d_leaf_wax_spatial_model.stan"), p$code_md5$stan_model)
add_drift("config.yaml",            md5_or_na("config.yaml"),            p$code_md5$config)
add_drift("calibration_csv",        md5_or_na("input_data/leafwax_d2h_c29_calibration_v1.csv"), p$input_md5$calibration_csv)
add_drift("sediment_prepared_rds",  md5_or_na("results/3_sediment_ready_for_modeling.rds"), p$input_md5$sediment_prepared_rds)

# (3) THIS model's direct fit inputs (stan_data + config).
add_drift(paste0("stan_data_", model, ".rds"),
          md5_or_na(file.path("prepared_data", paste0("stan_data_", model, ".rds"))),
          p$prepared_md5[[model]]$stan_data)
add_drift(paste0("config_", model, ".rds"),
          md5_or_na(file.path("prepared_data", paste0("config_", model, ".rds"))),
          p$prepared_md5[[model]]$config)

# (4) current SIF vs captured SIF (host computes the current md5 and passes it in).
if (prov_na_or_empty(cur_sif)) {
  drift <- c(drift, "SIF (--current-sif-md5 not provided; cannot verify the fitting image)")
} else if (prov_is_drift(cur_sif, p$sif$md5)) {
  drift <- c(drift, sprintf("SIF (now=%s captured=%s)", substr(cur_sif, 1, 12), substr(p$sif$md5, 1, 12)))
}

if (length(drift)) {
  cat("ERROR: workspace drifted from the captured fit-time provenance:\n")
  for (d in drift) cat("  ", d, "\n")
  cat("  Re-run prep to re-capture provenance, or restore the fit-time files/SIF.\n")
  quit(status = 1)
}
cat("Provenance OK for ", model, ": record complete; code/config/inputs/stan_data/config/SIF match capture (cmdstan ",
    p$environment$cmdstan_version, ").\n", sep = "")
