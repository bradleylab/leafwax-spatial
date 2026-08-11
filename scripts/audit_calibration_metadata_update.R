#!/usr/bin/env Rscript

# Verify that the public audited calibration input differs from the exact
# chordal fit-time snapshot only in corrected provenance and location-label
# metadata. This audit does not rerun model preparation or fitting.

fit_path <- Sys.getenv(
  "LEAFWAX_FIT_INPUT_SNAPSHOT",
  unset = "provenance/chordal_fit_2026-07-28/leafwax_d2h_c29_calibration_fit.csv"
)
public_path <- Sys.getenv(
  "LEAFWAX_PUBLIC_CALIBRATION",
  unset = "input_data/leafwax_d2h_c29_calibration_v1.csv"
)
output_path <- Sys.getenv(
  "LEAFWAX_CALIBRATION_AUDIT_OUT",
  unset = "scripts/reference_outputs/calibration_metadata_update_audit.csv"
)

expected_fit_md5 <- "bb52649130a02b9b7d897325899b4f5a"
expected_public_md5 <- "25f1f279a6b487955a8adefc2f7b303a"
allowed_changed_fields <- c("source", "location", "DOI")

for (path in c(fit_path, public_path)) {
  if (!file.exists(path)) stop("Missing calibration input: ", path)
}

actual_md5 <- unname(tools::md5sum(c(fit_path, public_path)))
if (actual_md5[[1L]] != expected_fit_md5) {
  stop("Fit-time snapshot md5 mismatch: ", actual_md5[[1L]])
}
if (actual_md5[[2L]] != expected_public_md5) {
  stop("Public calibration md5 mismatch: ", actual_md5[[2L]])
}

read_as_text <- function(path) {
  read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    colClasses = "character",
    na.strings = character()
  )
}

fit <- read_as_text(fit_path)
public <- read_as_text(public_path)
if (!identical(fit$obs_id, public$obs_id)) {
  stop("Fit-time and public calibration rows do not align by obs_id")
}
if (!identical(names(fit), names(public))) {
  stop("Fit-time and public calibration columns differ")
}

audit <- data.frame(
  field = names(fit),
  n_differences = vapply(
    names(fit),
    function(field) sum(fit[[field]] != public[[field]]),
    integer(1)
  ),
  stringsAsFactors = FALSE
)

unexpected <- audit$field[
  audit$n_differences > 0L & !audit$field %in% allowed_changed_fields
]
if (length(unexpected) > 0L) {
  stop("Unexpected changed fields: ", paste(unexpected, collapse = ", "))
}
if (!identical(audit$n_differences[audit$field == "source"], 27L) ||
    !identical(audit$n_differences[audit$field == "location"], 6L) ||
    !identical(audit$n_differences[audit$field == "DOI"], 27L)) {
  stop("Expected exactly 27 corrected source/DOI values and 6 location labels")
}

audit$fit_input_md5 <- expected_fit_md5
audit$public_input_md5 <- expected_public_md5
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write.csv(audit, output_path, row.names = FALSE)
message("Wrote ", output_path)
