#!/usr/bin/env Rscript
# check_r_deps.R — Lint that fails CI when an analysis script's library() call
# isn't covered by the Dockerfile's install.packages list (or base R).
#
# Usage: Rscript .github/scripts/check_r_deps.R
# Exits 1 with a diff if drift is found, 0 otherwise. No external R deps.
#
# Why this exists: build-container.yml only rebuilds on Dockerfile changes.
# Adding library(X) to a script silently drifts the script's dependency
# surface away from the container's frozen install list. This lint closes
# that gap on every PR + push to master.

# ── 1. Collect .R files (only git-tracked — matches what CI checks out) ──
# `git ls-files "*.R"` via system2 doesn't recurse into subdirs the way the
# shell-quoted pathspec does. Pull all tracked paths and filter in R.
tracked <- system2("git", "ls-files", stdout = TRUE)
tracked <- tracked[grepl("\\.R$", tracked)]
all_r <- file.path(".", tracked)
all_r <- all_r[!grepl("^\\./\\.github/", all_r)]

# ── 2. Extract library()/require() calls from each .R file ────────────────
extract_libs <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^\\s*#", lines)]              # drop comment-only lines
  pat <- "(library|require|requireNamespace)\\(['\"]?([A-Za-z][A-Za-z0-9._]*)['\"]?"
  m <- regmatches(lines, regexec(pat, lines))
  vapply(m, function(x) if (length(x) >= 3) x[3] else NA_character_,
         character(1))
}

used <- sort(unique(unlist(lapply(all_r, extract_libs))))
used <- used[!is.na(used)]

# ── 3. Parse Dockerfile install.packages list ─────────────────────────────
docker <- readLines("Dockerfile", warn = FALSE)
quoted <- regmatches(docker, gregexpr('"[a-zA-Z][a-zA-Z0-9._]*"', docker))
installed <- sort(unique(gsub('"', "", unlist(quoted))))

# Strip Dockerfile string literals that aren't package names (CMDSTAN_VERSION
# is one — anything that's clearly a setting). Also drop the LABEL values
# we know about; everything else stays.
not_packages <- c("CMDSTAN_VERSION")
installed <- setdiff(installed, not_packages)

# ── 4. Base + recommended R packages (always available, no install needed) ──
base_pkgs <- c(
  # base
  "base", "utils", "stats", "tools", "methods", "datasets", "parallel",
  "splines", "grid", "graphics", "grDevices", "compiler",
  # recommended (ship with most R distributions including rocker/r-ver)
  "nlme", "KernSmooth", "Matrix", "lattice", "mgcv", "cluster",
  "codetools", "foreign", "survival", "boot", "rpart", "nnet",
  "class", "MASS", "spatial"
)

# Optional output backends may be checked with requireNamespace() without being
# required for the tracked pipeline. The calibration freezer always writes CSV
# and RDS; arrow adds a local Parquet copy only when already installed.
optional_pkgs <- c("arrow")

# ── 4b. Meta-package expansion ────────────────────────────────────────────
# When the Dockerfile installs `tidyverse` it also installs its component
# packages (ggplot2, dplyr, etc.) plus common transitive deps that scripts
# library() directly. Account for these so the lint doesn't false-positive.
meta_expand <- list(
  tidyverse = c("ggplot2", "dplyr", "tidyr", "readr", "purrr", "tibble",
                "stringr", "forcats", "lubridate",
                "scales", "RColorBrewer"),
  ggplot2   = c("scales", "RColorBrewer"),
  ggpubr    = c("ggplot2", "scales", "RColorBrewer"),
  bayesplot = c("ggplot2", "scales"),
  cowplot   = c("ggplot2", "scales"),
  patchwork = c("ggplot2"),
  rnaturalearth = c("rnaturalearthdata", "sp"),
  fields    = c("spam", "viridisLite"),
  geoR      = c("splancs"),
  sf        = c("DBI"),
  terra     = c("Rcpp")
)
for (m in intersect(names(meta_expand), installed)) {
  installed <- c(installed, meta_expand[[m]])
}
installed <- sort(unique(installed))

# ── 5. Diff and report ────────────────────────────────────────────────────
covered <- c(installed, base_pkgs, optional_pkgs)
drift   <- setdiff(used, covered)

cat(sprintf("Scanned %d .R files; found %d unique library()/require() targets.\n",
            length(all_r), length(used)))
cat(sprintf("Dockerfile install.packages: %d packages.\n", length(installed)))

if (length(drift) > 0) {
  cat("\nERROR: the following packages are loaded by .R scripts but are NOT\n")
  cat("in the Dockerfile install.packages list and NOT base/recommended R:\n")
  for (p in drift) cat("  - ", p, "\n", sep = "")
  cat("\nFix: add them to the install.packages(c(...)) list in Dockerfile,\n")
  cat("then commit. The build-container workflow will rebuild + push to GHCR.\n")
  quit(status = 1)
}

cat("\nOK: all library()/require() calls are covered.\n")
quit(status = 0)
