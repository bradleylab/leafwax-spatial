#!/usr/bin/env Rscript

# Export the scientific portions of the fit-time provenance record without
# cluster-specific absolute paths or operator metadata.

input_path <- Sys.getenv(
  "LEAFWAX_FIT_PROVENANCE_RDS",
  unset = "results/c2_run_20260728_chordal/fit_provenance.rds"
)
output_path <- Sys.getenv(
  "LEAFWAX_PUBLIC_FIT_PROVENANCE",
  unset = "provenance/chordal_fit_2026-07-28/fit_environment.json"
)

if (!file.exists(input_path)) stop("Missing fit-time provenance: ", input_path)
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("The jsonlite package is required to export fit provenance")
}

record <- readRDS(input_path)
required <- c("captured_at", "environment", "sif", "code_md5", "input_md5", "prepared_md5")
missing <- setdiff(required, names(record))
if (length(missing) > 0L) {
  stop("Fit-time provenance is missing: ", paste(missing, collapse = ", "))
}

public <- list(
  fit_label = "2026-07-28 chordal-distance calibration fit",
  captured_at = record$captured_at,
  environment = record$environment,
  container_md5 = record$sif$md5,
  code_md5 = record$code_md5,
  input_md5 = record$input_md5,
  prepared_md5 = record$prepared_md5
)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(public, output_path, auto_unbox = TRUE, pretty = TRUE, na = "null")
cat("Wrote", output_path, "\n")
