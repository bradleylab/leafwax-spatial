# Dockerfile for leafwax spatial modeling pipeline
# Builds an image with R 4.4.x, CmdStan, terra, sf, and all dependencies
# needed for the Bayesian spatial d2H calibration pipeline.
#
# Target: GHCR (ghcr.io/bradleylab/leafwax-spatial)
# Runtime: Apptainer on WashU Compute2 HPC
#
# Build: docker build -t leafwax-spatial .
# Or via GitHub Actions (see .github/workflows/build-container.yml)

FROM rocker/r-ver:4.4.1

LABEL org.opencontainers.image.source="https://github.com/bradleylab/leafwax-spatial"
LABEL org.opencontainers.image.description="Bayesian spatial calibration of leaf wax hydrogen isotopes (d2H) against precipitation and environmental covariates. R 4.4.1, CmdStan 2.36.0, terra, sf."
LABEL org.opencontainers.image.authors="Alexander S. Bradley"
LABEL org.opencontainers.image.licenses="MIT"
LABEL org.opencontainers.image.title="leafwax-calibration"

# System dependencies for terra, sf, and CmdStan
# terra needs: GDAL, PROJ, GEOS
# sf needs: GDAL, PROJ, GEOS, libudunits2
# CmdStan needs: g++, make
RUN apt-get update && apt-get install -y --no-install-recommends \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libudunits2-dev \
    libnetcdf-dev \
    libsqlite3-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libuv1-dev \
    libtbb-dev \
    tcl8.6-dev \
    tk8.6-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    cmake \
    g++ \
    make \
    && rm -rf /var/lib/apt/lists/*

# Pin the CRAN releases used by the frozen R 4.4.1 container. Deriv 4.3.0 uses
# a newer R API, while terra 1.9-34 requires a newer GDAL multidimensional API
# than Ubuntu 22.04 provides. Install terra separately so its compiled
# geospatial dependency layer is cached and verified independently.
RUN R -e ' \
    install.packages( \
        "https://cloud.r-project.org/src/contrib/Archive/Deriv/Deriv_4.2.0.tar.gz", \
        repos = NULL, \
        type = "source" \
    ); \
    if (!requireNamespace("Deriv", quietly = TRUE) || \
        as.character(packageVersion("Deriv")) != "4.2.0") { \
        stop("Deriv 4.2.0 installation failed") \
    }; \
    install.packages("Rcpp", repos = "https://cloud.r-project.org", Ncpus = 1); \
    if (!requireNamespace("Rcpp", quietly = TRUE)) { \
        stop("Rcpp installation failed") \
    }; \
    install.packages( \
        "https://cloud.r-project.org/src/contrib/Archive/terra/terra_1.9-27.tar.gz", \
        repos = NULL, \
        type = "source", \
        Ncpus = 1 \
    ); \
    if (!requireNamespace("terra", quietly = TRUE) || \
        as.character(packageVersion("terra")) != "1.9.27") { \
        stop("terra 1.9-27 installation failed") \
    } \
'

# Install the R packages used across preparation, fitting, diagnostics, and
# manuscript-number regeneration. Compilation is deliberately serial: parallel
# source builds of terra and its dependency tree can exceed GitHub runner memory.
# Arrow is optional in 2i_freeze_calibration.R and is not required by the tracked
# pipeline, so it is omitted from the image.
RUN R -e ' \
    packages <- c( \
        "yaml", \
        "digest", \
        "jsonlite", \
        "tidyverse", \
        "sf", \
        "fields", \
        "geosphere", \
        "ncdf4", \
        "posterior", \
        "loo", \
        "bayesplot", \
        "geoR", \
        "ape", \
        "ggpubr", \
        "corrplot", \
        "cowplot", \
        "patchwork", \
        "viridis", \
        "rnaturalearth", \
        "rnaturalearthdata", \
        "broom", \
        "car", \
        "factoextra", \
        "glmnet", \
        "gstat", \
        "knitr", \
        "ncf", \
        "openxlsx", \
        "spdep" \
    ); \
    install.packages(packages, repos = "https://cloud.r-project.org", Ncpus = 1); \
    missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]; \
    if (length(missing)) stop("R package installation failed: ", paste(missing, collapse = ", ")) \
'

# Install CmdStanR from R-universe (official distribution channel)
RUN R -e ' \
    install.packages("cmdstanr", repos = c( \
        "https://stan-dev.r-universe.dev", \
        "https://cloud.r-project.org" \
    )) \
'

# Install CmdStan (the C++ toolchain that cmdstanr calls)
# This compiles Stan's math library — takes ~10 min
ENV CMDSTAN_VERSION=2.36.0
RUN R -e 'cmdstanr::install_cmdstan(version = "2.36.0", cores = 2)'

# Verify the installation
RUN R -e ' \
    library(terra); cat("terra", as.character(packageVersion("terra")), "\n"); \
    library(sf); cat("sf", as.character(packageVersion("sf")), "\n"); \
    library(cmdstanr); cat("cmdstanr", as.character(packageVersion("cmdstanr")), "\n"); \
    cmdstanr::cmdstan_path(); \
    cat("CmdStan OK\n"); \
    cat("GDAL:", terra::gdal(), "\n"); \
    cat("GEOS:", sf::sf_extSoftVersion()["GEOS"], "\n"); \
    cat("PROJ:", sf::sf_extSoftVersion()["PROJ"], "\n") \
'

# Default working directory — will be bind-mounted at runtime
WORKDIR /pipeline

CMD ["R"]
