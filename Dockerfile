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

# Install core R packages needed by the pipeline
# Order matters: dependencies first, then packages that depend on them
# Post-fit scripts (5c, 5e, 7) added geoR/ape/ggpubr/corrplot/cowplot/
# patchwork/rnaturalearth/viridis over time; backfilled here so the full
# pipeline runs in-container (not just the fit step).
RUN R -e ' \
    install.packages(c( \
        "yaml", \
        "tidyverse", \
        "terra", \
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
        "rnaturalearthdata" \
    ), repos = "https://cloud.r-project.org", Ncpus = 4) \
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
RUN R -e 'cmdstanr::install_cmdstan(version = "2.36.0", cores = 4)'

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
