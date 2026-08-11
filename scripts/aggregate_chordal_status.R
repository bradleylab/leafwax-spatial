#!/usr/bin/env Rscript
# aggregate_chordal_status.R
#
# Batch-level completion gate for the chordal fit array. Each SLURM
# array task fits ONE model and writes per-model artifacts under
# model_output/<model>/ (fit_status.rds always; a DONE marker only when the full
# core artifact set — fit.rds + diagnostics.rds + posterior_draws.rds — is on
# disk). This aggregator collects those per-model files and judges the WHOLE
# batch.
#
# FAIL-CLOSED: exits nonzero unless every expected model is complete (status
# "completed" AND core artifacts present AND DONE marker present). loo.rds is
# best-effort and only warned about, not gated (loo::loo can legitimately fail).
#
# Run from the leafwax_working root after the array finishes:
#   Rscript scripts/aggregate_chordal_status.R [model_output]

args <- commandArgs(trailingOnly = TRUE)
model_output  <- if (length(args) >= 1) args[[1]] else "model_output"
prepared_data <- if (length(args) >= 2) args[[2]] else "prepared_data"

# The 17 models in the fit array (job_fit_chordal.sh MODELS), in order.
EXPECTED <- c(
  "baseline_sp", "baseline_veg_sp", "baseline_env_sp", "c4_only_sp",
  "elevation_only_sp", "elevation_c4_sp", "elevation_c4_interact_sp",
  "full_sp", "full_interact_sp",
  "baseline_sp_rfoff", "baseline_env_sp_rfoff", "full_interact_sp_rfoff",
  "baseline", "baseline_veg", "baseline_env", "full", "full_interact"
)

core_files <- c("fit.rds", "diagnostics.rds", "posterior_draws.rds")

# Freshness: every core artifact AND the DONE marker must be newer than the
# model's stan_data. Existence alone is not enough — a crash on a resumed run can
# leave a stale DONE + old diagnostics/draws (from a prior stan_data) in place.
fresh_vs_stan_data <- function(mdir, model) {
  sd <- file.path(prepared_data, paste0("stan_data_", model, ".rds"))
  if (!file.exists(sd)) return(NA)   # cannot verify freshness -> not a pass
  paths <- file.path(mdir, c(core_files, "DONE"))
  if (!all(file.exists(paths))) return(FALSE)
  all(file.mtime(paths) > file.mtime(sd))
}

rows <- lapply(EXPECTED, function(m) {
  mdir <- file.path(model_output, m)
  st_f <- file.path(mdir, "fit_status.rds")
  st   <- if (file.exists(st_f)) tryCatch(readRDS(st_f), error = function(e) NULL) else NULL
  core_ok <- all(file.exists(file.path(mdir, core_files)))
  done_ok <- file.exists(file.path(mdir, "DONE"))
  fresh   <- fresh_vs_stan_data(mdir, m)
  data.frame(
    model = m,
    status = if (!is.null(st)) st$status else "MISSING",
    core_complete = core_ok,
    done_marker = done_ok,
    fresh_vs_stan = fresh,
    loo_present = if (!is.null(st)) isTRUE(st$loo_present) else file.exists(file.path(mdir, "loo.rds")),
    divergences = if (!is.null(st) && !is.null(st$divergences)) st$divergences else NA_integer_,
    max_treedepth = if (!is.null(st) && !is.null(st$max_treedepth)) st$max_treedepth else NA_integer_,
    low_ebfmi = if (!is.null(st) && !is.null(st$low_ebfmi)) st$low_ebfmi else NA_integer_,
    max_rhat = if (!is.null(st) && !is.null(st$max_rhat)) st$max_rhat else NA_real_,
    min_ess_bulk = if (!is.null(st) && !is.null(st$min_ess_bulk)) st$min_ess_bulk else NA_real_,
    stringsAsFactors = FALSE
  )
})
tab <- do.call(rbind, rows)

# A model passes only if completed, core-complete, DONE present, AND all artifacts
# are strictly newer than its stan_data (isTRUE guards the NA "can't verify" case).
tab$ok <- tab$status == "completed" & tab$core_complete & tab$done_marker &
          vapply(tab$fresh_vs_stan, isTRUE, logical(1))

cat("== Chordal batch status (", model_output, ") ==\n", sep = "")
print(tab[, c("model", "status", "core_complete", "done_marker", "fresh_vs_stan",
              "loo_present", "divergences", "max_rhat", "min_ess_bulk")], row.names = FALSE)

n_ok <- sum(tab$ok)
n_exp <- nrow(tab)
missing_loo <- tab$model[!tab$loo_present & tab$ok]
if (length(missing_loo)) {
  cat("\nWARNING: complete models with no loo.rds (LOOIC unavailable):",
      paste(missing_loo, collapse = ", "), "\n")
}

cat(sprintf("\n%d of %d models complete.\n", n_ok, n_exp))
stale <- tab$model[tab$core_complete & tab$done_marker &
                   !vapply(tab$fresh_vs_stan, isTRUE, logical(1))]
if (length(stale)) {
  cat("STALE (artifacts/DONE not newer than stan_data, or stan_data missing):",
      paste(stale, collapse = ", "), "\n")
}
if (n_ok < n_exp) {
  incomplete <- tab$model[!tab$ok]
  cat("INCOMPLETE / FAILED / STALE:", paste(incomplete, collapse = ", "), "\n")
  cat("Re-run the affected array indices before assembling the run.\n")
  quit(status = 1)
}
cat("BATCH COMPLETE — all", n_exp, "chordal models finished with a full artifact set.\n")
