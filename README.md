# arcpt: AR Changepoint Detection

<!-- badges: start -->
[![R-CMD-check](https://github.com/atharv-1511/arcpt/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/atharv-1511/arcpt/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/atharv-1511/arcpt/branch/main/graph/badge.svg)](https://app.codecov.io/gh/atharv-1511/arcpt?branch=main)
<!-- badges: end -->

An R package for detecting changepoints in autoregressive (AR) time series.

## Overview

`arcpt` provides tools for detecting changepoints in AR(1) and AR(2) time series, with optional trend components. It uses the regression-based changepoint framework from the `EnvCpt` package.

This package was developed as part of the **GSoC 2026** application for the "Generalized Changepoint Regression" project under the R Project for Statistical Computing.

## Installation

```r
# Install from GitHub
# install.packages("devtools")
devtools::install_github("atharv-1511/arcpt")
```

## Usage

### Basic AR(1) Changepoint Detection

```r
library(arcpt)

# Generate AR(1) data with a changepoint
set.seed(42)
x1 <- arima.sim(model = list(ar = 0.3), n = 100)
x2 <- arima.sim(model = list(ar = 0.9), n = 100)
data <- c(x1, x2)

# Detect changepoints
result <- ar_cpt(data, order = 1)
print(result)
```

### AR(2) with Trend

```r
# AR(2) model with trend component
result <- ar_cpt(data, order = 2, trend = TRUE)
summary(result)
```

### Visualization

```r
# Plot results
plot(result)
```

### Extract Results

```r
# Get changepoint locations
cpts.arcpt(result)

# Get number of changepoints
ncpts.arcpt(result)

# Get coefficients per segment
coef(result)
```

## Features

- **AR(1) and AR(2) models**: Detect changes in autoregressive structure
- **Optional trend**: Include linear trend component
- **Multiple methods**: PELT (default) or AMOC
- **Various penalties**: MBIC, BIC, SIC, AIC, Hannan-Quinn, Manual
- **Comprehensive output**: Changepoint locations, segment coefficients, fitted model
- **S3 methods**: `print`, `summary`, `plot`, `coef`, `cpts.arcpt`, `ncpts.arcpt`
- **Comprehensive validation**: Extensive input checking with informative error messages
- **52 unit tests**: Thorough test coverage for reliability

## Test Coverage

The package includes 52 comprehensive tests covering:

- Input validation (data types, NA values, invalid parameters)
- AR(1) and AR(2) model functionality
- PELT and AMOC methods
- All penalty types
- S3 methods (print, summary, plot, coef)
- Edge cases (constant data, white noise, near-unit-root)
- Numerical stability (small/large variance)
- Multiple changepoint detection
- Trend + AR combinations

## Dependencies

- `changepoint`: Core changepoint detection algorithms
- `EnvCpt`: Changepoint regression framework

## Author

**Atharv Raskar** - GSoC 2026 Contributor

## License

GPL-3

## Acknowledgments

- Rebecca Killick (GSoC Mentor, author of `changepoint` and `EnvCpt`)
- Colin Gallagher (GSoC Mentor)
