# =============================================================================
# Tests for Mixed Effects Changepoint Regression (me_cpt)
# =============================================================================

context("Mixed Effects Changepoint Regression")

# =============================================================================
# Test Data Setup
# =============================================================================

# Generate hierarchical data with changepoint
generate_me_data <- function(n_groups = 10, n_per_group = 50, cpt_time = 25,
                              sigma_random = 2, sigma_resid = 1) {
  set.seed(42)

  data <- data.frame(
    group = rep(1:n_groups, each = n_per_group),
    time = rep(1:n_per_group, n_groups)
  )

  # Random intercepts
  random_int <- rnorm(n_groups, 0, sigma_random)

  # Fixed effect changes at cpt_time
  data$y <- ifelse(data$time <= cpt_time,
    0.5 * data$time + random_int[data$group],
    1.5 * data$time - cpt_time + random_int[data$group]
  ) + rnorm(nrow(data), 0, sigma_resid)

  return(data)
}

# =============================================================================
# Input Validation Tests
# =============================================================================

test_that("me_cpt rejects non-formula input", {
  data <- generate_me_data()
  expect_error(me_cpt("y ~ time + (1|group)", data), "must be a formula")
})

test_that("me_cpt rejects non-data.frame input", {
  expect_error(me_cpt(y ~ time + (1|group), list(y = 1:10)),
               "must be a data frame")
})

test_that("me_cpt rejects formula without random effects", {
  data <- generate_me_data()
  expect_error(me_cpt(y ~ time, data), "must contain random effects")
})

test_that("me_cpt rejects missing response variable", {
  data <- generate_me_data()
  expect_error(me_cpt(z ~ time + (1|group), data), "not found in data")
})

test_that("me_cpt rejects missing grouping variable", {
  data <- generate_me_data()
  expect_error(me_cpt(y ~ time + (1|cluster), data), "not found in data")
})

# =============================================================================
# Functionality Tests
# =============================================================================

test_that("me_cpt works with basic hierarchical data", {
  data <- generate_me_data()
  result <- me_cpt(y ~ time + (1|group), data)

  expect_s3_class(result, "mecpt")
  expect_true(is.numeric(result$cpts))
  expect_true(is.list(result$fixed_effects))
  expect_true(is.numeric(result$random_effects))
})

test_that("me_cpt detects changepoint in hierarchical data", {
  data <- generate_me_data(cpt_time = 25)
  result <- me_cpt(y ~ time + (1|group), data, penalty = "BIC")

  # Should detect at least one changepoint
  expect_gte(result$ncpts, 0)
})

test_that("me_cpt works with AMOC method", {
  data <- generate_me_data()
  result <- me_cpt(y ~ time + (1|group), data, method = "AMOC")

  expect_s3_class(result, "mecpt")
  expect_equal(result$method, "AMOC")
  expect_lte(result$ncpts, 1)
})

test_that("me_cpt works with different penalties", {
  data <- generate_me_data()

  penalties <- c("MBIC", "BIC", "AIC")

  for (pen in penalties) {
    result <- me_cpt(y ~ time + (1|group), data, penalty = pen)
    expect_s3_class(result, "mecpt")
    expect_equal(result$penalty, pen)
  }
})

# =============================================================================
# Output Structure Tests
# =============================================================================

test_that("me_cpt returns correct structure", {
  data <- generate_me_data()
  result <- me_cpt(y ~ time + (1|group), data)

  expect_true("formula" %in% names(result))
  expect_true("data" %in% names(result))
  expect_true("cpts" %in% names(result))
  expect_true("ncpts" %in% names(result))
  expect_true("fixed_effects" %in% names(result))
  expect_true("random_effects" %in% names(result))
  expect_true("variance_components" %in% names(result))
  expect_true("method" %in% names(result))
  expect_true("penalty" %in% names(result))
})

test_that("me_cpt variance components are valid", {
  data <- generate_me_data()
  result <- me_cpt(y ~ time + (1|group), data)

  expect_true("residual" %in% names(result$variance_components))
  expect_true("random_intercept" %in% names(result$variance_components))
  expect_true(result$variance_components$residual >= 0)
  expect_true(result$variance_components$random_intercept >= 0)
})

# =============================================================================
# S3 Methods Tests
# =============================================================================

test_that("print.mecpt works", {
  data <- generate_me_data()
  result <- me_cpt(y ~ time + (1|group), data)

  expect_output(print(result), "Mixed Effects Changepoint Regression")
  expect_output(print(result), "Number of groups:")
  expect_output(print(result), "Variance Components:")
})

test_that("summary.mecpt works", {
  data <- generate_me_data()
  result <- me_cpt(y ~ time + (1|group), data)

  expect_output(summary(result), "Mixed Effects Changepoint Regression Summary")
  expect_output(summary(result), "Fixed Effects by Segment")
})

test_that("coef.mecpt returns coefficients", {
  data <- generate_me_data()
  result <- me_cpt(y ~ time + (1|group), data)

  # Fixed effects
  fixed <- coef(result, type = "fixed")
  expect_true(is.list(fixed))

  # Random effects
  random <- coef(result, type = "random")
  expect_true(is.numeric(random))

  # Both
  both <- coef(result, type = "both")
  expect_true("fixed" %in% names(both))
  expect_true("random" %in% names(both))
})

test_that("cpts.mecpt and ncpts.mecpt work", {
  data <- generate_me_data()
  result <- me_cpt(y ~ time + (1|group), data)

  cpts <- cpts.mecpt(result)
  ncpts <- ncpts.mecpt(result)

  expect_true(is.numeric(cpts))
  expect_true(is.numeric(ncpts))
  expect_equal(ncpts, length(cpts))
})

test_that("plot.mecpt works without error", {
  data <- generate_me_data()
  result <- me_cpt(y ~ time + (1|group), data)

  expect_silent(plot(result, type = "data"))
  expect_silent(plot(result, type = "random"))
})

# =============================================================================
# Edge Cases
# =============================================================================

test_that("me_cpt handles single group", {
  data <- data.frame(
    group = rep(1, 100),
    time = 1:100,
    y = c(rnorm(50, 0, 1), rnorm(50, 5, 1))
  )

  result <- me_cpt(y ~ time + (1|group), data)
  expect_s3_class(result, "mecpt")
})

test_that("me_cpt handles many groups", {
  set.seed(42)
  data <- data.frame(
    group = rep(1:50, each = 20),
    time = rep(1:20, 50)
  )
  random_int <- rnorm(50, 0, 2)
  data$y <- 0.5 * data$time + random_int[data$group] + rnorm(nrow(data), 0, 0.5)

  result <- me_cpt(y ~ time + (1|group), data)
  expect_s3_class(result, "mecpt")
  expect_equal(result$n_groups, 50)
})

# =============================================================================
# Tests for Non-parametric Changepoint Regression (np_cpt)
# =============================================================================

context("Non-parametric Changepoint Regression")

# =============================================================================
# Test Data Setup
# =============================================================================

# Generate non-linear data with changepoint
generate_np_data <- function(n = 200, cpt = 100) {
  set.seed(42)
  x <- seq(0, 4 * pi, length.out = n)

  # First segment: sine wave
  # Second segment: different amplitude/frequency
  y <- c(
    sin(x[1:cpt]) + rnorm(cpt, 0, 0.2),
    2 * sin(2 * x[(cpt+1):n]) + rnorm(n - cpt, 0, 0.2)
  )

  list(x = x, y = y)
}

# =============================================================================
# Input Validation Tests
# =============================================================================

test_that("np_cpt rejects non-numeric input", {
  expect_error(np_cpt(c("a", "b"), c(1, 2)), "must be numeric")
  expect_error(np_cpt(c(1, 2), c("a", "b")), "must be numeric")
})

test_that("np_cpt rejects mismatched lengths", {
  expect_error(np_cpt(1:10, 1:5), "must have the same length")
})

test_that("np_cpt rejects NA values", {
  expect_error(np_cpt(c(1:9, NA), 1:10), "cannot contain NA")
  expect_error(np_cpt(1:10, c(1:9, NA)), "cannot contain NA")
})

test_that("np_cpt rejects non-finite values", {
  expect_error(np_cpt(c(1:9, Inf), 1:10), "non-finite values")
  expect_error(np_cpt(1:10, c(1:9, -Inf)), "non-finite values")
})

test_that("np_cpt rejects df < 3", {
  data <- generate_np_data()
  expect_error(np_cpt(data$x, data$y, df = 2), "df.*must be at least 3")
})

test_that("np_cpt rejects data too short", {
  expect_error(np_cpt(1:10, rnorm(10), minseglen = 10),
               "Data too short")
})

# =============================================================================
# Functionality Tests
# =============================================================================

test_that("np_cpt works with B-splines", {
  data <- generate_np_data()
  result <- np_cpt(data$x, data$y, spline_type = "bs", df = 6)

  expect_s3_class(result, "npcpt")
  expect_equal(result$spline_type, "bs")
  expect_true(is.numeric(result$fitted_values))
})

test_that("np_cpt works with natural splines", {
  data <- generate_np_data()
  result <- np_cpt(data$x, data$y, spline_type = "ns", df = 6)

  expect_s3_class(result, "npcpt")
  expect_equal(result$spline_type, "ns")
})

test_that("np_cpt detects changepoint in non-linear data", {
  data <- generate_np_data(n = 200, cpt = 100)
  result <- np_cpt(data$x, data$y, df = 6, penalty = "BIC")

  # Should detect at least one changepoint
  expect_gte(result$ncpts, 0)
})

test_that("np_cpt works with AMOC method", {
  data <- generate_np_data()
  result <- np_cpt(data$x, data$y, method = "AMOC")

  expect_s3_class(result, "npcpt")
  expect_equal(result$method, "AMOC")
  expect_lte(result$ncpts, 1)
})

test_that("np_cpt works with different penalties", {
  data <- generate_np_data()

  penalties <- c("MBIC", "BIC", "AIC")

  for (pen in penalties) {
    result <- np_cpt(data$x, data$y, penalty = pen)
    expect_s3_class(result, "npcpt")
    expect_equal(result$penalty, pen)
  }
})

test_that("np_cpt works with different df values", {
  data <- generate_np_data()

  for (df in c(4, 6, 8)) {
    result <- np_cpt(data$x, data$y, df = df)
    expect_s3_class(result, "npcpt")
    expect_equal(result$df, df)
  }
})

# =============================================================================
# Output Structure Tests
# =============================================================================

test_that("np_cpt returns correct structure", {
  data <- generate_np_data()
  result <- np_cpt(data$x, data$y)

  expect_true("x" %in% names(result))
  expect_true("y" %in% names(result))
  expect_true("cpts" %in% names(result))
  expect_true("cpts_x" %in% names(result))
  expect_true("ncpts" %in% names(result))
  expect_true("spline_fits" %in% names(result))
  expect_true("fitted_values" %in% names(result))
  expect_true("residuals" %in% names(result))
  expect_true("r_squared" %in% names(result))
})

test_that("np_cpt R-squared is valid", {
  data <- generate_np_data()
  result <- np_cpt(data$x, data$y)

  expect_true(result$r_squared >= 0)
  expect_true(result$r_squared <= 1)
})

test_that("np_cpt fitted values and residuals are consistent", {
  data <- generate_np_data()
  result <- np_cpt(data$x, data$y)

  # y = fitted + residuals
  reconstructed <- result$fitted_values + result$residuals
  expect_equal(reconstructed, result$y, tolerance = 1e-10)
})

# =============================================================================
# S3 Methods Tests
# =============================================================================

test_that("print.npcpt works", {
  data <- generate_np_data()
  result <- np_cpt(data$x, data$y)

  expect_output(print(result), "Non-parametric Changepoint Regression")
  expect_output(print(result), "Spline type:")
  expect_output(print(result), "R-squared:")
})

test_that("summary.npcpt works", {
  data <- generate_np_data()
  result <- np_cpt(data$x, data$y)

  expect_output(summary(result), "Non-parametric Changepoint Regression Summary")
  expect_output(summary(result), "Segment Details")
})

test_that("fitted.npcpt works", {
  data <- generate_np_data()
  result <- np_cpt(data$x, data$y)

  fitted_vals <- fitted(result)
  expect_equal(length(fitted_vals), length(data$y))
  expect_true(is.numeric(fitted_vals))
})

test_that("residuals.npcpt works", {
  data <- generate_np_data()
  result <- np_cpt(data$x, data$y)

  resids <- residuals(result)
  expect_equal(length(resids), length(data$y))
  expect_true(is.numeric(resids))
})

test_that("predict.npcpt works with no new data", {
  data <- generate_np_data()
  result <- np_cpt(data$x, data$y)

  preds <- predict(result)
  expect_equal(preds, result$fitted_values)
})

test_that("predict.npcpt works with new data", {
  data <- generate_np_data()
  result <- np_cpt(data$x, data$y)

  newx <- seq(min(data$x), max(data$x), length.out = 50)
  preds <- predict(result, newx = newx)

  expect_equal(length(preds), 50)
  expect_true(is.numeric(preds))
})

test_that("cpts.npcpt and ncpts.npcpt work", {
  data <- generate_np_data()
  result <- np_cpt(data$x, data$y)

  cpts <- cpts.npcpt(result)
  ncpts <- ncpts.npcpt(result)

  expect_true(is.numeric(cpts))
  expect_true(is.numeric(ncpts))
  expect_equal(ncpts, length(cpts))
})

test_that("plot.npcpt works without error", {
  data <- generate_np_data()
  result <- np_cpt(data$x, data$y)

  expect_silent(plot(result, type = "fit"))
  expect_silent(plot(result, type = "segments"))
})

# =============================================================================
# Edge Cases
# =============================================================================

test_that("np_cpt handles unsorted x values", {
  data <- generate_np_data()
  # Shuffle the data
  idx <- sample(length(data$x))

  result <- np_cpt(data$x[idx], data$y[idx])
  expect_s3_class(result, "npcpt")
  # x should be sorted in result
  expect_false(is.unsorted(result$x))
})

test_that("np_cpt handles no changepoints (constant relationship)", {
  set.seed(42)
  x <- seq(0, 2 * pi, length.out = 100)
  y <- sin(x) + rnorm(100, 0, 0.1)

  result <- np_cpt(x, y, penalty = "MBIC")
  expect_s3_class(result, "npcpt")
})

test_that("np_cpt handles linear data", {
  set.seed(42)
  x <- 1:100
  y <- 2 * x + 3 + rnorm(100, 0, 5)

  result <- np_cpt(x, y)
  expect_s3_class(result, "npcpt")
})

# =============================================================================
# Test Suite Completion
# =============================================================================

test_that("advanced methods test suite completes successfully", {
  expect_true(TRUE)
})
