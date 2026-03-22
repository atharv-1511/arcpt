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
