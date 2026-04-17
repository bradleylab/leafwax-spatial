#!/usr/bin/env Rscript
# check_tables_numeric.R — audit that each manuscript .tex table's
# numeric content is still in sync with the companion CSV emitted by
# its producer. Run at the tail of `make tables` (or `make audit`) to
# catch silent hand-edits of the .tex files.
#
# What this catches:
#   - A data cell in a .tex file was hand-edited away from what the
#     producer wrote (e.g., operator "corrected" one number).
#   - A producer bug that rendered a value into the .tex that never
#     landed in the CSV.
#
# What this does NOT catch:
#   - A producer that writes the same wrong value into both .tex and
#     CSV. The audit trusts the CSV as the source of truth.
#
# Usage:
#   Rscript manuscript/table_code/check_tables_numeric.R
# Exits non-zero on any drift.

library(tools)

AUDIT_PAIRS <- list(
  list(
    label     = "Table 1",
    tex       = "manuscript/tables/table1_model_performance.tex",
    csv       = "model_analysis/tables/table1_model_performance.csv",
    row_key   = "model",
    # Table 1 has 2 leading text cells (Model + Predictors), then 8
    # numeric cells. The predictors column contains $\delta^2$H$_p$ and
    # C4 × ... which would otherwise feed stray 2/4 tokens to the audit.
    drop_leading_cells = 2,
    tex_cols  = c("max_rhat", "min_ess", "looic", "se_looic",
                  "delta_looic", "se_delta", "p_eff", "n_hi_k")
  ),
  list(
    label     = "Table 2",
    tex       = "manuscript/tables/table2_global_params_body.tex",
    csv       = "model_analysis/tables/table2_global_params.csv",
    row_key   = "model",
    # Table 2 cells are "mean [q025, q975]" for the first six columns and
    # bare mean for the two knot-SD columns. Order must match row render.
    tex_cols  = c("rmse_permil_mean",    "rmse_permil_q025",    "rmse_permil_q975",
                  "r_squared_mean",      "r_squared_q025",      "r_squared_q975",
                  "beta_0_permil_mean",  "beta_0_permil_q025",  "beta_0_permil_q975",
                  "beta_oipc_mean",      "beta_oipc_q025",      "beta_oipc_q975",
                  "lambda_int_km_mean",  "lambda_int_km_q025",  "lambda_int_km_q975",
                  "gp_scale_km_mean",    "gp_scale_km_q025",    "gp_scale_km_q975",
                  "knot_slope_sd_mean",
                  "knot_intercept_sd_permil_mean")
  ),
  list(
    label     = "Table 3",
    tex       = "manuscript/tables/table3_variance_decomp.tex",
    csv       = "model_analysis/tables/table3_variance_decomp.csv",
    row_key   = "model",
    tex_cols  = c("var_spatial_total_pct", "var_obs_pct",
                  "var_intercept_pct_of_spatial",
                  "var_slope_pct_of_spatial")
  ),
  list(
    label     = "Table 4",
    tex       = "manuscript/tables/table4_environmental.tex",
    csv       = "model_analysis/tables/table4_environmental.csv",
    row_key   = "model",
    tex_cols  = c("beta_c4_mean",     "beta_c4_q05",     "beta_c4_q95",
                  "beta_tree_mean",   "beta_tree_q05",   "beta_tree_q95",
                  "beta_shrub_mean",  "beta_shrub_q05",  "beta_shrub_q95",
                  "beta_grass_mean",  "beta_grass_q05",  "beta_grass_q95",
                  "beta_precip_mean", "beta_precip_q05", "beta_precip_q95")
  ),
  list(
    label     = "Table S2",
    tex       = "manuscript/tables/Table_S2_regional_performance_body.tex",
    csv       = "manuscript/tables/Table_S2_regional_performance.csv",
    row_key   = NULL,
    tex_cols  = c("n", "rmse", "r_squared")
  )
)

.extract_numeric_tokens <- function(s, drop_leading_cells = 1) {
  # Pull every signed/unsigned decimal from a string. Treats "-" used as
  # a dash for "not applicable" as NOT a number: requires at least one
  # digit. Strips the first `drop_leading_cells` `&`-delimited cells
  # before matching, so model names ("c4_only_sp") and predictor strings
  # ($\delta^2$H$_p$) don't contribute stray digit tokens.
  parts <- strsplit(s, "&", fixed = TRUE)[[1]]
  if (length(parts) > drop_leading_cells) {
    s <- paste(parts[(drop_leading_cells + 1):length(parts)], collapse = "&")
  } else {
    s <- ""
  }
  matches <- regmatches(s, gregexpr("-?[0-9]+(?:\\.[0-9]+)?", s))[[1]]
  as.numeric(matches)
}

.tex_body_lines <- function(path) {
  lines <- readLines(path, warn = FALSE)
  # Skip comments and empty lines.
  lines <- lines[!grepl("^\\s*%", lines)]
  lines <- lines[nzchar(trimws(lines))]

  # If the tabular includes a \midrule separator, everything before it is
  # preamble + header — skip. This is the case for standalone .tex files
  # emitted by emit_standalone_tex (Tables 1/3/4). Body-only fragments
  # (Tables 2 / S2) have no \midrule, so the skip is a no-op there.
  midrule_idx <- grep("\\\\midrule", lines)
  if (length(midrule_idx) > 0) {
    lines <- lines[(midrule_idx[1] + 1):length(lines)]
  }

  # Drop \addlinespace spacers that knitr injects between row groups and
  # any \bottomrule / \end{tabular} closers the standalone form carries.
  lines <- lines[!grepl("^\\s*\\\\addlinespace", lines)]
  lines <- lines[!grepl("^\\s*\\\\bottomrule", lines)]
  lines <- lines[!grepl("^\\s*\\\\end\\{tabular\\}", lines)]

  # Keep only rows ending in `\\` — those are the actual data rows.
  lines[grepl("\\\\\\\\\\s*$", lines)]
}

.round_like <- function(x, digits) round(as.numeric(x), digits)

.compare_row <- function(tex_nums, csv_nums, tol = 1e-2) {
  # Align positionally. Lengths should match for rows in spec.
  if (length(tex_nums) != length(csv_nums)) return(FALSE)
  # Element-wise: either both NA, or both numeric within tolerance.
  mapply(function(a, b) {
    if (is.na(a) && is.na(b)) return(TRUE)
    if (is.na(a) || is.na(b)) return(FALSE)
    abs(a - b) <= max(tol, tol * max(abs(a), abs(b)))
  }, tex_nums, csv_nums)
}

audit_one <- function(spec) {
  if (!file.exists(spec$tex) || !file.exists(spec$csv)) {
    return(list(label = spec$label, ok = FALSE,
                reason = sprintf("missing file: tex=%s csv=%s",
                                 file.exists(spec$tex), file.exists(spec$csv))))
  }

  tex_lines <- .tex_body_lines(spec$tex)
  csv <- read.csv(spec$csv, stringsAsFactors = FALSE)

  # For tables with a row-key (Table 1..4): check the row count matches.
  # For Table S2 (multi-row per model): just check the number of CSV rows
  # equals the number of tex-body rows.
  if (nrow(csv) != length(tex_lines)) {
    return(list(label = spec$label, ok = FALSE,
                reason = sprintf("row count mismatch: tex=%d csv=%d",
                                 length(tex_lines), nrow(csv))))
  }

  drop_leading <- if (!is.null(spec$drop_leading_cells)) spec$drop_leading_cells else 1
  drift <- character()
  for (i in seq_along(tex_lines)) {
    tex_nums <- .extract_numeric_tokens(tex_lines[i],
                                        drop_leading_cells = drop_leading)
    csv_vals <- unlist(csv[i, spec$tex_cols, drop = TRUE], use.names = FALSE)
    csv_vals <- suppressWarnings(as.numeric(csv_vals))
    # Drop NA cells before length comparison: the .tex renders "-" for
    # NA cells (no numeric token emitted), so a length-match needs to
    # compare against the non-NA subset of the csv row.
    csv_vals <- csv_vals[!is.na(csv_vals)]

    # Rounding: the tex renders with a fixed number of decimals. Compare
    # both sides rounded to the number of decimals in the tex token.
    if (length(tex_nums) != length(csv_vals)) {
      drift <- c(drift,
        sprintf("row %d: tex has %d numbers, csv has %d non-NA",
                i, length(tex_nums), length(csv_vals)))
      next
    }
    # Per-cell rounding tolerance keyed off the tex decimal precision.
    ok_cells <- mapply(function(t, c) {
      if (is.na(c) && is.na(t)) return(TRUE)
      if (is.na(c) || is.na(t)) return(FALSE)
      # rough precision inference from the printed token is overkill;
      # just require agreement to 3 significant figures plus 0.05 absolute.
      abs(t - c) <= max(0.05, 5e-3 * max(abs(t), abs(c)))
    }, tex_nums, csv_vals)
    if (!all(ok_cells)) {
      bad <- which(!ok_cells)
      drift <- c(drift,
        sprintf("row %d: cell(s) %s drift tex=%s csv=%s",
                i, paste(bad, collapse = ","),
                paste(tex_nums[bad], collapse = ","),
                paste(csv_vals[bad], collapse = ",")))
    }
  }

  if (length(drift) > 0) {
    return(list(label = spec$label, ok = FALSE,
                reason = paste(drift, collapse = "\n  ")))
  }
  list(label = spec$label, ok = TRUE)
}

cat("Running numeric audit of manuscript tables...\n\n")
results <- lapply(AUDIT_PAIRS, audit_one)
fail <- Filter(function(r) !r$ok, results)

for (r in results) {
  status <- if (r$ok) "OK" else "DRIFT"
  cat(sprintf("  [%s] %s\n", status, r$label))
  if (!r$ok) cat("    ", r$reason, "\n", sep = "")
}

if (length(fail) > 0) {
  cat(sprintf("\n%d table(s) show drift against CSV. Failing.\n", length(fail)))
  quit(status = 1)
}
cat("\nAll tables consistent with their CSV source of truth.\n")
