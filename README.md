# arcpt: Changepoint Detection in Regression Models

<!-- badges: start -->
[![R-CMD-check](https://github.com/atharv-1511/arcpt/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/atharv-1511/arcpt/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/atharv-1511/arcpt/branch/main/graph/badge.svg)](https://app.codecov.io/gh/atharv-1511/arcpt?branch=main)
<!-- badges: end -->

A comprehensive R package for detecting changepoints in autoregressive time series, mixed effects models, and non-parametric regression.

## Overview

`arcpt` provides tools for detecting changepoints in various regression settings:

1. **AR Changepoint Detection** (`ar_cpt`): AR(1) and AR(2) time series with optional trend
2. **Mixed Effects Changepoint Detection** (`me_cpt`): Hierarchical data with random effects
3. **Non-parametric Changepoint Detection** (`np_cpt`): Flexible spline-based regression

This package was developed as part of the **GSoC 2026** application for the "Generalized Changepoint Regression" project under the R Project for Statistical Computing.

## Installation

```r
# Install from GitHub
# install.packages("devtools")
devtools::install_github("atharv-1511/arcpt")
```

## Usage

### 1. AR Changepoint Detection

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
plot(result)

# AR(2) model with trend component
result <- ar_cpt(data, order = 2, trend = TRUE)
summary(result)
```

### 2. Mixed Effects Changepoint Detection

```r
# Generate hierarchical data
set.seed(42)
n_groups <- 10
n_per_group <- 50

data <- data.frame(
  group = rep(1:n_groups, each = n_per_group),
  time = rep(1:n_per_group, n_groups)
)

# Random intercepts per group
random_int <- rnorm(n_groups, 0, 2)

# Fixed effect changes at time 25
data$y <- ifelse(data$time <= 25,
  0.5 * data$time + random_int[data$group],
  1.5 * data$time - 25 + random_int[data$group]
) + rnorm(nrow(data), 0, 1)

# Detect changepoints with mixed effects
result <- me_cpt(y ~ time + (1|group), data = data)
print(result)
summary(result)
plot(result, type = "data")
plot(result, type = "random")
```

### 3. Non-parametric Changepoint Detection

```r
# Generate non-linear data with changepoint
set.seed(42)
n <- 200
x <- seq(0, 4*pi, length.out = n)

# First half: sine wave; Second half: different pattern
y <- c(
  sin(x[1:100]) + rnorm(100, 0, 0.2),
  2 * sin(2 * x[101:200]) + rnorm(100, 0, 0.2)
)

# Detect changepoints with spline regression
result <- np_cpt(x, y, spline_type = "ns", df = 6)
print(result)
plot(result, type = "fit")
plot(result, type = "segments")

# Predictions
predict(result, newx = seq(0, 4*pi, length.out = 50))
```

## Features

### Core Functions

| Function | Description | Key Arguments |
|----------|-------------|---------------|
| `ar_cpt()` | AR changepoint detection | `order`, `trend`, `method`, `penalty` |
| `me_cpt()` | Mixed effects changepoints | `formula`, `data`, `method`, `penalty` |
| `np_cpt()` | Non-parametric changepoints | `spline_type`, `df`, `method`, `penalty` |

### AR Changepoint Detection (`ar_cpt`)

- **AR(1) and AR(2) models**: Detect changes in autoregressive structure
- **Optional trend**: Include linear trend component
- **Multiple methods**: PELT (default) or AMOC
- **Various penalties**: MBIC, BIC, SIC, AIC, Hannan-Quinn, Manual
- **S3 methods**: `print`, `summary`, `plot`, `coef`, `cpts.arcpt`, `ncpts.arcpt`

### Mixed Effects Changepoint Detection (`me_cpt`)

- **Hierarchical data**: Support for grouped/nested data structures
- **Random intercepts**: Group-specific baseline effects
- **Fixed effects that change**: Coefficients that vary at changepoints
- **lme4-style syntax**: `y ~ x + (1|group)` formula notation
- **Variance components**: Estimates of random effect and residual variance
- **S3 methods**: `print`, `summary`, `plot`, `coef`

### Non-parametric Changepoint Detection (`np_cpt`)

- **B-splines**: Flexible basis functions with local support
- **Natural splines**: Splines with linear extrapolation at boundaries
- **Automatic df selection**: Cross-validation for degrees of freedom
- **Segment-wise fits**: Different smooth functions per segment
- **Prediction**: Interpolate/extrapolate to new x values
- **S3 methods**: `print`, `summary`, `plot`, `fitted`, `residuals`, `predict`

## Test Coverage

The package includes **97 comprehensive tests** covering:

### AR Model Tests (52 tests)
- Input validation (data types, NA values, invalid parameters)
- AR(1) and AR(2) model functionality
- PELT and AMOC methods
- All penalty types
- S3 methods (print, summary, plot, coef)
- Edge cases (constant data, white noise, near-unit-root)
- Numerical stability (small/large variance)
- Multiple changepoint detection
- Trend + AR combinations

### Advanced Methods Tests (45 tests)
- Mixed effects input validation
- Hierarchical data handling
- Formula parsing
- Variance component estimation
- Non-parametric spline fitting
- B-spline and natural spline support
- Prediction functionality
- Edge cases and robustness

## Dependencies

- `changepoint`: Core changepoint detection algorithms
- `EnvCpt`: Changepoint regression framework
- `splines`: B-spline and natural spline basis functions
- `stats`: Statistical functions

## Why This Package?

Current changepoint packages force all regression coefficients to change at each changepoint. This package demonstrates:

1. **Mixed effects**: Random effects stay constant while fixed effects change
2. **Non-parametric**: Flexible relationships that can differ by segment
3. **Extensible framework**: Foundation for the full "fixed vs changing covariates" approach

This addresses a genuine gap that practitioners have requested for years.

## Author

**Atharv Raskar** - GSoC 2025 Selectee | CRAN Package Author (Grepreaper)

## License

GPL-3

## Acknowledgments

- Rebecca Killick (GSoC Mentor, author of `changepoint` and `EnvCpt`)
- Colin Gallagher (GSoC Mentor)

## Related Links

- [GSoC 2026 Project Page](https://github.com/rstats-gsoc/gsoc2026/wiki/Generalized-changepoint-regression)
- [EnvCpt Package](https://github.com/rkillick/EnvCpt)
- [changepoint Package](https://github.com/rkillick/changepoint)
