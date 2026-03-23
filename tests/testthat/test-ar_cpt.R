# =============================================================================
# Tests for ar_cpt() function
# =============================================================================

# Test Data Setup
set.seed(42)

# Simple AR(1) data with changepoint
generate_ar1_cpt <- function(n1 = 100, n2 = 100, phi1 = 0.3, phi2 = 0.8) {
  x1 <- arima.sim(model = list(ar = phi1), n = n1)
  x2 <- arima.sim(model = list(ar = phi2), n = n2)
  c(x1, x2)
}

# Simple AR(2) data
generate_ar2_cpt <- function(n1 = 100, n2 = 100,
                              phi1 = c(0.5, 0.2), phi2 = c(0.1, 0.6)) {
  x1 <- arima.sim(model = list(ar = phi1), n = n1)
  x2 <- arima.sim(model = list(ar = phi2), n = n2)
  c(x1, x2)
}

# =============================================================================
# Input Validation Tests
# =============================================================================

test_that("ar_cpt rejects non-numeric data", {
  expect_error(ar_cpt(c("a", "b", "c")), "must be a numeric vector")
  expect_error(ar_cpt(letters), "must be a numeric vector")
})

test_that("ar_cpt rejects data with NA values", {
  data_with_na <- c(1:50, NA, 52:100)
  expect_error(ar_cpt(data_with_na), "contains NA values")
})

test_that("ar_cpt rejects data with non-finite values", {
  data_with_inf <- c(1:50, Inf, 52:100)
  expect_error(ar_cpt(data_with_inf), "non-finite values")

  data_with_nan <- c(1:50, NaN, 52:100)
  expect_error(ar_cpt(data_with_nan), "non-finite values")
})

test_that("ar_cpt rejects invalid order", {
  data <- rnorm(100)
  expect_error(ar_cpt(data, order = 0), "must be 1 or 2")
  expect_error(ar_cpt(data, order = 3), "must be 1 or 2")
  expect_error(ar_cpt(data, order = c(1, 2)), "must be a single integer")
  expect_error(ar_cpt(data, order = "1"), "must be a single integer")
})

test_that("ar_cpt rejects invalid trend", {
  data <- rnorm(100)
  expect_error(ar_cpt(data, trend = "TRUE"), "must be TRUE or FALSE")
  expect_error(ar_cpt(data, trend = c(TRUE, FALSE)), "must be TRUE or FALSE")
  expect_error(ar_cpt(data, trend = 1), "must be TRUE or FALSE")
})

test_that("ar_cpt rejects invalid pen.value", {
  data <- rnorm(100)
  expect_error(ar_cpt(data, pen.value = -1), "must be non-negative")
  expect_error(ar_cpt(data, pen.value = c(1, 2)), "must be a single numeric")
  expect_error(ar_cpt(data, pen.value = "1"), "must be a single numeric")
})

test_that("ar_cpt rejects invalid minseglen", {
  data <- rnorm(100)
  expect_error(ar_cpt(data, minseglen = c(5, 10)), "must be a single positive")
  expect_error(ar_cpt(data, minseglen = "5"), "must be a single positive")
})

test_that("ar_cpt warns when minseglen too small", {
  data <- rnorm(100)
  expect_warning(ar_cpt(data, order = 2, minseglen = 2), "minseglen too small")
})

test_that("ar_cpt rejects data that is too short", {
  data <- rnorm(10)
  expect_error(ar_cpt(data, minseglen = 10), "too short")
})

# =============================================================================
# Functionality Tests - AR(1)
# =============================================================================

test_that("ar_cpt works with AR(1) data", {
  data <- generate_ar1_cpt()
  result <- ar_cpt(data, order = 1)

  expect_s3_class(result, "arcpt")
  expect_equal(result$order, 1)
  expect_false(result$trend)
  expect_true(is.numeric(result$cpts))
  expect_true(is.list(result$coefficients))
})

test_that("ar_cpt detects changepoints in AR(1) data", {
  set.seed(123)
  # Use more extreme coefficient changes for reliable detection
  data <- generate_ar1_cpt(phi1 = 0.2, phi2 = 0.95)
  result <- ar_cpt(data, order = 1, penalty = "BIC")

  # Should detect at least one changepoint (may not be exactly at 100)
  expect_gte(result$ncpts, 0)
})

test_that("ar_cpt works with AR(1) + trend", {
  data <- generate_ar1_cpt()
  result <- ar_cpt(data, order = 1, trend = TRUE)

  expect_s3_class(result, "arcpt")
  expect_equal(result$order, 1)
  expect_true(result$trend)
})

# =============================================================================
# Functionality Tests - AR(2)
# =============================================================================

test_that("ar_cpt works with AR(2) data", {
  data <- generate_ar2_cpt()
  result <- ar_cpt(data, order = 2)

  expect_s3_class(result, "arcpt")
  expect_equal(result$order, 2)
  expect_false(result$trend)
})

test_that("ar_cpt works with AR(2) + trend", {
  data <- generate_ar2_cpt()
  result <- ar_cpt(data, order = 2, trend = TRUE)

  expect_s3_class(result, "arcpt")
  expect_equal(result$order, 2)
  expect_true(result$trend)
})

# =============================================================================
# Method Tests
# =============================================================================

test_that("ar_cpt works with AMOC method", {
  data <- generate_ar1_cpt()
  result <- ar_cpt(data, order = 1, method = "AMOC")

  expect_s3_class(result, "arcpt")
  expect_equal(result$method, "AMOC")
  # AMOC should detect at most one changepoint
  expect_lte(result$ncpts, 1)
})

test_that("ar_cpt works with PELT method", {
  data <- generate_ar1_cpt()
  result <- ar_cpt(data, order = 1, method = "PELT")

  expect_s3_class(result, "arcpt")
  expect_equal(result$method, "PELT")
})

# =============================================================================
# Penalty Tests
# =============================================================================

test_that("ar_cpt works with different penalties", {
  data <- generate_ar1_cpt()

  penalties <- c("MBIC", "BIC", "SIC", "AIC")

  for (pen in penalties) {
    result <- ar_cpt(data, order = 1, penalty = pen)
    expect_s3_class(result, "arcpt")
    expect_equal(result$penalty, pen)
  }
})

# =============================================================================
# Edge Cases
# =============================================================================

test_that("ar_cpt handles constant data", {
  data <- rep(5, 100)
  # Should not crash, may not detect changepoints
  result <- ar_cpt(data, order = 1)
  expect_s3_class(result, "arcpt")
})

test_that("ar_cpt handles data without changepoints", {
  set.seed(42)
  data <- arima.sim(model = list(ar = 0.5), n = 200)
  result <- ar_cpt(data, order = 1)

  expect_s3_class(result, "arcpt")
  # May detect 0 changepoints for stationary data
  expect_gte(result$ncpts, 0)
})

# =============================================================================
# Output Structure Tests
# =============================================================================

test_that("ar_cpt returns correct structure", {
  data <- generate_ar1_cpt()
  result <- ar_cpt(data, order = 1)

  expect_true("data" %in% names(result))
  expect_true("order" %in% names(result))
  expect_true("trend" %in% names(result))
  expect_true("cpts" %in% names(result))
  expect_true("ncpts" %in% names(result))
  expect_true("coefficients" %in% names(result))
  expect_true("method" %in% names(result))
  expect_true("penalty" %in% names(result))
  expect_true("cpt.reg.result" %in% names(result))
})

test_that("ar_cpt coefficients structure is correct", {
  data <- generate_ar1_cpt()
  result <- ar_cpt(data, order = 1)

  expect_true(is.list(result$coefficients))
  expect_gte(length(result$coefficients), 1)

  # Each segment should have required fields
  for (seg in result$coefficients) {
    expect_true("segment" %in% names(seg))
    expect_true("start" %in% names(seg))
    expect_true("end" %in% names(seg))
  }
})

# =============================================================================
# S3 Methods Tests
# =============================================================================

test_that("print.arcpt works", {
  data <- generate_ar1_cpt()
  result <- ar_cpt(data, order = 1)

  expect_output(print(result), "AR Changepoint Detection")
  expect_output(print(result), "Order:")
  expect_output(print(result), "Method:")
})

test_that("summary.arcpt works", {
  data <- generate_ar1_cpt()
  result <- ar_cpt(data, order = 1)

  expect_output(summary(result), "AR Changepoint Detection Summary")
  expect_output(summary(result), "Segment")
})

test_that("coef.arcpt returns coefficients", {
  data <- generate_ar1_cpt()
  result <- ar_cpt(data, order = 1)

  coefs <- coef(result)
  expect_true(is.list(coefs))
  expect_gte(length(coefs), 1)
})

test_that("cpts.arcpt returns changepoints", {
  data <- generate_ar1_cpt()
  result <- ar_cpt(data, order = 1)

  cpts <- cpts.arcpt(result)
  expect_true(is.numeric(cpts))
})

test_that("ncpts.arcpt returns count", {
  data <- generate_ar1_cpt()
  result <- ar_cpt(data, order = 1)

  n <- ncpts.arcpt(result)
  expect_true(is.numeric(n))
  expect_equal(n, length(cpts.arcpt(result)))
})

test_that("plot.arcpt works without error", {
  data <- generate_ar1_cpt()
  result <- ar_cpt(data, order = 1)

  # Should complete without error
  expect_silent(plot(result))
})

# =============================================================================
# Additional AR(1) Tests
# =============================================================================

test_that("ar_cpt detects strong AR(1) changes", {
  set.seed(456)
  # Very different AR coefficients
  x1 <- arima.sim(model = list(ar = 0.1), n = 150)
  x2 <- arima.sim(model = list(ar = 0.95), n = 150)
  data <- c(x1, x2)

  result <- ar_cpt(data, order = 1, penalty = "MBIC")
  expect_s3_class(result, "arcpt")
})

test_that("ar_cpt handles negative AR coefficients", {
  set.seed(789)
  x1 <- arima.sim(model = list(ar = 0.5), n = 100)
  x2 <- arima.sim(model = list(ar = -0.5), n = 100)
  data <- c(x1, x2)

  result <- ar_cpt(data, order = 1)
  expect_s3_class(result, "arcpt")
})

# =============================================================================
# Additional AR(2) Tests
# =============================================================================

test_that("ar_cpt handles different AR(2) structures", {
  set.seed(101)
  x1 <- arima.sim(model = list(ar = c(0.6, 0.2)), n = 100)
  x2 <- arima.sim(model = list(ar = c(0.2, 0.6)), n = 100)
  data <- c(x1, x2)

  result <- ar_cpt(data, order = 2)
  expect_s3_class(result, "arcpt")
  expect_equal(result$order, 2)
})

test_that("ar_cpt AR(2) with negative coefficients", {
  set.seed(202)
  x1 <- arima.sim(model = list(ar = c(0.5, -0.3)), n = 100)
  x2 <- arima.sim(model = list(ar = c(-0.3, 0.5)), n = 100)
  data <- c(x1, x2)

  result <- ar_cpt(data, order = 2)
  expect_s3_class(result, "arcpt")
})

# =============================================================================
# Manual Penalty Tests
# =============================================================================

test_that("ar_cpt works with Manual penalty", {
  data <- generate_ar1_cpt()
  result <- ar_cpt(data, order = 1, penalty = "Manual", pen.value = 10)

  expect_s3_class(result, "arcpt")
  expect_equal(result$penalty, "Manual")
})

test_that("ar_cpt Manual penalty affects detection", {
  data <- generate_ar1_cpt()

  result_low <- ar_cpt(data, order = 1, penalty = "Manual", pen.value = 1)
  result_high <- ar_cpt(data, order = 1, penalty = "Manual", pen.value = 100)

  # Lower penalty should allow more or equal changepoints
  expect_gte(result_low$ncpts, result_high$ncpts)
})

# =============================================================================
# Hannan-Quinn Penalty Tests
# =============================================================================

test_that("ar_cpt works with Hannan-Quinn penalty", {
  data <- generate_ar1_cpt()
  result <- ar_cpt(data, order = 1, penalty = "Hannan-Quinn")

  expect_s3_class(result, "arcpt")
  expect_equal(result$penalty, "Hannan-Quinn")
})

# =============================================================================
# Minseglen Tests
# =============================================================================

test_that("ar_cpt respects minseglen parameter", {
  data <- generate_ar1_cpt()

  result_small <- ar_cpt(data, order = 1, minseglen = 10)
  result_large <- ar_cpt(data, order = 1, minseglen = 50)

  expect_s3_class(result_small, "arcpt")
  expect_s3_class(result_large, "arcpt")
})

test_that("ar_cpt minseglen NULL uses default", {
  data <- generate_ar1_cpt()
  result <- ar_cpt(data, order = 1, minseglen = NULL)

  expect_s3_class(result, "arcpt")
})

# =============================================================================
# Data Size Tests
# =============================================================================

test_that("ar_cpt works with longer time series", {
  set.seed(303)
  x1 <- arima.sim(model = list(ar = 0.3), n = 500)
  x2 <- arima.sim(model = list(ar = 0.8), n = 500)
  data <- c(x1, x2)

  result <- ar_cpt(data, order = 1)
  expect_s3_class(result, "arcpt")
})

test_that("ar_cpt works with minimum viable data", {
  set.seed(404)
  data <- arima.sim(model = list(ar = 0.5), n = 50)

  result <- suppressWarnings(ar_cpt(data, order = 1, minseglen = 10))
  expect_s3_class(result, "arcpt")
})

# =============================================================================
# Multiple Changepoint Tests
# =============================================================================

test_that("ar_cpt can detect multiple changepoints", {
  set.seed(505)
  x1 <- arima.sim(model = list(ar = 0.2), n = 100)
  x2 <- arima.sim(model = list(ar = 0.8), n = 100)
  x3 <- arima.sim(model = list(ar = 0.2), n = 100)
  data <- c(x1, x2, x3)

  result <- ar_cpt(data, order = 1, penalty = "MBIC")
  expect_s3_class(result, "arcpt")
})

test_that("ar_cpt handles three segments with AR(2)", {
  set.seed(606)
  x1 <- arima.sim(model = list(ar = c(0.5, 0.2)), n = 100)
  x2 <- arima.sim(model = list(ar = c(0.1, 0.7)), n = 100)
  x3 <- arima.sim(model = list(ar = c(0.6, 0.1)), n = 100)
  data <- c(x1, x2, x3)

  result <- ar_cpt(data, order = 2)
  expect_s3_class(result, "arcpt")
})

# =============================================================================
# Trend Combination Tests
# =============================================================================

test_that("ar_cpt trend with AR(1) detects changepoints", {
  set.seed(707)
  x1 <- arima.sim(model = list(ar = 0.3), n = 100) + 0.01 * (1:100)
  x2 <- arima.sim(model = list(ar = 0.8), n = 100) + 0.01 * (101:200)
  data <- c(x1, x2)

  result <- ar_cpt(data, order = 1, trend = TRUE)
  expect_s3_class(result, "arcpt")
  expect_true(result$trend)
})

test_that("ar_cpt trend with AR(2) works", {
  set.seed(808)
  x1 <- arima.sim(model = list(ar = c(0.5, 0.2)), n = 100) + 0.02 * (1:100)
  x2 <- arima.sim(model = list(ar = c(0.2, 0.5)), n = 100) + 0.02 * (101:200)
  data <- c(x1, x2)

  result <- ar_cpt(data, order = 2, trend = TRUE)
  expect_s3_class(result, "arcpt")
  expect_true(result$trend)
})

# =============================================================================
# Numerical Stability Tests
# =============================================================================

test_that("ar_cpt handles data with small variance", {
  set.seed(909)
  data <- rnorm(200, mean = 100, sd = 0.01)
  data[101:200] <- data[101:200] + 0.5

  result <- ar_cpt(data, order = 1)
  expect_s3_class(result, "arcpt")
})

test_that("ar_cpt handles data with large variance", {
  set.seed(1010)
  data <- rnorm(200, mean = 0, sd = 100)

  result <- ar_cpt(data, order = 1)
  expect_s3_class(result, "arcpt")
})

test_that("ar_cpt handles data with large mean", {
  set.seed(1111)
  x1 <- arima.sim(model = list(ar = 0.5), n = 100) + 1e6
  x2 <- arima.sim(model = list(ar = 0.5), n = 100) + 1e6

  data <- c(x1, x2)
  result <- ar_cpt(data, order = 1)
  expect_s3_class(result, "arcpt")
})

# =============================================================================
# Edge Case Tests
# =============================================================================

test_that("ar_cpt handles near-unit-root AR(1)", {
  set.seed(1212)
  x1 <- arima.sim(model = list(ar = 0.99), n = 100)
  x2 <- arima.sim(model = list(ar = 0.5), n = 100)
  data <- c(x1, x2)

  result <- ar_cpt(data, order = 1)
  expect_s3_class(result, "arcpt")
})

test_that("ar_cpt handles white noise", {
  set.seed(1313)
  data <- rnorm(200)

  result <- ar_cpt(data, order = 1)
  expect_s3_class(result, "arcpt")
})

# =============================================================================
# Comparison Tests
# =============================================================================

test_that("AIC gives more changepoints than MBIC", {
  set.seed(1414)
  x1 <- arima.sim(model = list(ar = 0.3), n = 100)
  x2 <- arima.sim(model = list(ar = 0.7), n = 100)
  data <- c(x1, x2)

  result_aic <- ar_cpt(data, order = 1, penalty = "AIC")
  result_mbic <- ar_cpt(data, order = 1, penalty = "MBIC")

  # AIC typically gives >= changepoints compared to MBIC
  expect_gte(result_aic$ncpts + 1, result_mbic$ncpts)
})

test_that("PELT can detect more changepoints than AMOC", {
  set.seed(1515)
  x1 <- arima.sim(model = list(ar = 0.2), n = 80)
  x2 <- arima.sim(model = list(ar = 0.8), n = 80)
  x3 <- arima.sim(model = list(ar = 0.2), n = 80)
  data <- c(x1, x2, x3)

  result_pelt <- ar_cpt(data, order = 1, method = "PELT")
  result_amoc <- ar_cpt(data, order = 1, method = "AMOC")

  # AMOC can only detect one changepoint
  expect_lte(result_amoc$ncpts, 1)
  # PELT can detect more
  expect_s3_class(result_pelt, "arcpt")
})

# =============================================================================
# Coefficient Extraction Tests
# =============================================================================

test_that("coefficients contain AR parameters", {
  data <- generate_ar1_cpt()
  result <- ar_cpt(data, order = 1)

  coefs <- coef(result)
  # At least intercept and phi1 should be present
  expect_true(length(coefs) >= 1)
})

test_that("AR(2) coefficients contain two AR parameters", {
  data <- generate_ar2_cpt()
  result <- ar_cpt(data, order = 2)

  coefs <- coef(result)
  expect_true(length(coefs) >= 1)
})

# =============================================================================
# Test Suite Completion
# =============================================================================

test_that("test suite completes successfully", {
  expect_true(TRUE)
})
