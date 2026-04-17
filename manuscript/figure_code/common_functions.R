# common_functions.R
# Shared setup for manuscript/figure_code/ scripts.
# Ported from Ω_archive_leafwax_model/spatial_model_current/figures/common_functions.R
# and re-rooted for the Apr-2026 layout: repo root is two levels up
# (manuscript/figure_code/ → manuscript/ → repo root).

library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(loo)
library(patchwork)
library(viridis)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
# `gt` and `fields` were in the original archive copy but are not required
# by any current figure script; load lazily if a script needs them.
suppressWarnings({
  if (requireNamespace("gt", quietly = TRUE)) library(gt)
  if (requireNamespace("fields", quietly = TRUE)) library(fields)
})

# PROJECT_ROOT = manuscript/ (parent of figure_code/)
PROJECT_ROOT <- normalizePath("..", mustWork = FALSE)
# REPO_ROOT holds config.yaml, input_data/, results/, pipeline scripts.
REPO_ROOT <- normalizePath(file.path("..", ".."), mustWork = FALSE)

cat("=== Directory Information ===\n")
cat("Current directory:", getwd(), "\n")
cat("Project root (manuscript/):", PROJECT_ROOT, "\n")
cat("Repo root:", REPO_ROOT, "\n")

load_project_config <- function() {
  possible_paths <- c(
    file.path(REPO_ROOT, "0_load_config.R"),
    file.path(PROJECT_ROOT, "0_load_config.R"),
    "../../0_load_config.R",
    "../0_load_config.R"
  )
  config_path <- NULL
  for (path in possible_paths) {
    if (file.exists(path)) { config_path <- path; break }
  }
  if (is.null(config_path)) stop("0_load_config.R not found")
  cat("Loading config from:", config_path, "\n")
  current_dir <- getwd()
  tryCatch({
    setwd(dirname(config_path))
    source(basename(config_path))
  }, finally = {
    setwd(current_dir)
  })
}

create_directories <- function() {
  dir.create("main", showWarnings = FALSE, recursive = TRUE)
  dir.create("supplement", showWarnings = FALSE, recursive = TRUE)
  dir.create("tables", showWarnings = FALSE, recursive = TRUE)
}

get_world_map <- function() {
  ne_countries(scale = "medium", returnclass = "sf")
}

print_figure_message <- function(figure_name) {
  cat("\n", strrep("=", 60), "\n")
  cat("Creating", figure_name, "\n")
  cat(strrep("=", 60), "\n")
}

cat("=== Common functions loaded ===\n\n")
