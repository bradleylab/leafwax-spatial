#!/usr/bin/env Rscript
# capture_fit_provenance.R
#
# Capture fit-time provenance from INSIDE the fitting SIF on Compute2, BEFORE the
# fit array runs. This is the authoritative record of the environment the fits
# actually ran in — it cannot be reconstructed as confidently after the code or
# the C2 workspace changes. build_run_manifest.R REQUIRES this record for an
# authoritative manifest.
#
# Run once via apptainer exec (so package/CmdStan versions are the SIF's, not the
# manifest-building machine's). The host passes the SIF checksum + git state as
# args (git/md5 of the SIF are host-side facts):
#   apptainer exec ... "$SIF" Rscript scripts/capture_fit_provenance.R \
#     --out results/fit_provenance.rds \
#     --sif-md5 "$SIF_MD5" --sif-path "$SIF" \
#     --git-head "$GIT_HEAD" --git-dirty "$GIT_DIRTY"
#
# Captures, from inside the SIF: R version; cmdstanr/posterior/loo/… package
# versions; the real CmdStan version (cmdstanr::cmdstan_version()); and md5 hashes
# of the actual fitting code (4a/4b/4c/4d), config, and modeling inputs
# (calibration CSV + the prepared sediment RDS). Writes results/fit_provenance.rds
# (+ .json if jsonlite is available).

suppressWarnings(suppressMessages(library(cmdstanr)))

# --containall does not reliably expose CmdStan on its default search path, so
# (like 4c's fit runner) pin it explicitly; otherwise cmdstan_version() would be
# NA and the record would silently omit the fitting CmdStan version. Honour an
# externally-provided CMDSTAN env if set, else the SIF's baked path.
cmdstan_dir <- Sys.getenv("CMDSTAN", "/root/.cmdstan/cmdstan-2.36.0")
tryCatch(cmdstanr::set_cmdstan_path(cmdstan_dir), error = function(e)
  message("set_cmdstan_path('", cmdstan_dir, "') failed: ", conditionMessage(e)))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NA_character_) {
  i <- match(flag, args); if (!is.na(i) && i < length(args)) args[[i + 1]] else default
}
out_path  <- get_arg("--out", "results/fit_provenance.rds")
sif_md5   <- get_arg("--sif-md5")
sif_path  <- get_arg("--sif-path")
git_head  <- get_arg("--git-head")
git_dirty <- get_arg("--git-dirty")
prepared_dir <- get_arg("--prepared-dir", "prepared_data")

source("scripts/provenance_common.R")   # model taxonomy + completeness/drift definitions

md5_or_na <- function(p) if (file.exists(p)) unname(tools::md5sum(p)) else NA_character_
pkg_ver   <- function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_)

# Per-model prepared-input hashes (stan_data_<m>.rds + config_<m>.rds) — the files
# 4c reads directly. Captured for every model with a prepared stan_data so drift in
# any of them (post-prep) is detectable at fit time and in the manifest.
sd_files <- list.files(prepared_dir, pattern = "^stan_data_.*\\.rds$", full.names = FALSE)
models_found <- sub("^stan_data_", "", sub("\\.rds$", "", sd_files))
prepared_md5 <- setNames(lapply(models_found, function(m) list(
  stan_data = md5_or_na(file.path(prepared_dir, paste0("stan_data_", m, ".rds"))),
  config    = md5_or_na(file.path(prepared_dir, paste0("config_", m, ".rds"))))), models_found)

record <- list(
  captured_at = as.character(Sys.time()),
  where = "inside fitting SIF on Compute2",
  environment = list(
    r_version       = R.version.string,
    cmdstan_version = tryCatch(cmdstanr::cmdstan_version(), error = function(e) NA_character_),
    packages = list(
      cmdstanr  = pkg_ver("cmdstanr"),
      posterior = pkg_ver("posterior"),
      loo       = pkg_ver("loo"),
      rstan     = pkg_ver("rstan"),
      tidyverse = pkg_ver("tidyverse"))),
  sif = list(path = sif_path, md5 = sif_md5),
  git = list(head = if (identical(git_head, "NA")) NA_character_ else git_head,
             dirty = if (identical(git_dirty, "NA")) NA else as.logical(as.integer(git_dirty))),
  code_md5 = list(
    spatial_functions = md5_or_na("4a_spatial_functions.R"),
    stan_prep         = md5_or_na("4b_stan_prep.R"),
    fit_models        = md5_or_na("4c_fit_models.R"),
    stan_model        = md5_or_na("4d_leaf_wax_spatial_model.stan"),
    config            = md5_or_na("config.yaml")),
  input_md5 = list(
    calibration_csv       = md5_or_na("input_data/leafwax_d2h_c29_calibration_v1.csv"),
    sediment_prepared_rds = md5_or_na("results/3_sediment_ready_for_modeling.rds")),
  prepared_md5 = prepared_md5)

# ── Fail closed: an authoritative provenance record must be COMPLETE ─────────────
# Uses the shared completeness definition (provenance_common.R), validated against
# the models actually prepared, so capture and the manifest cannot disagree.
missing <- provenance_missing_fields(record, models_found)
if (length(missing)) {
  cat("ERROR: fit-time provenance is INCOMPLETE — refusing to write.\n")
  cat("  missing/invalid fields:", paste(missing, collapse = ", "), "\n")
  if (prov_na_or_empty(record$environment$cmdstan_version))
    cat("  (CmdStan NA usually means the path is unset — set CMDSTAN or check the SIF.)\n")
  quit(status = 1)
}
if (!length(models_found)) { cat("ERROR: no prepared stan_data found in '", prepared_dir, "'.\n"); quit(status = 1) }

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
saveRDS(record, out_path)
cat("Wrote", out_path, "\n")
cat("  CmdStan (in-SIF):", record$environment$cmdstan_version,
    "| cmdstanr:", record$environment$packages$cmdstanr,
    "| posterior:", record$environment$packages$posterior, "\n")
cat("  SIF md5:", ifelse(is.na(sif_md5), "—", substr(sif_md5, 1, 12)),
    "| code(4c) md5:", substr(record$code_md5$fit_models, 1, 12), "\n")

if (requireNamespace("jsonlite", quietly = TRUE)) {
  json_path <- sub("\\.rds$", ".json", out_path)
  writeLines(jsonlite::toJSON(record, auto_unbox = TRUE, pretty = TRUE), json_path)
  cat("Wrote", json_path, "\n")
}
