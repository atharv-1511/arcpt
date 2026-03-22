# =============================================================================
# Tests for arcpt S3 Methods
# =============================================================================

# Test Data Setup
set.seed(42)
test_data <- c(arima.sim(model = list(ar = 0.3), n = 100),
               arima.sim(model = list(ar = 0.8), n = 100))

# =============================================================================
# Print Method Tests
# =============================================================================

test_that("print.arcpt produces output", {
  result <- ar_cpt(test_data, order = 1)
  expect_output(print(result), "AR Changepoint Detection Results")
  expect_output(print(result), "AR\\(1\\)")
  expect_output(print(result), "Method:")
})

test_that("print.arcpt returns invisibly", {
  result <- ar_cpt(test_data, order = 1)
  expect_invisible(print(result))
})

test_that("print.arcpt shows trend when used", {
  result <- ar_cpt(test_data, order = 1, trend = TRUE)
  expect_output(print(result), "with trend")
})

# =============================================================================
# Summary Method Tests
# =============================================================================

test_that("summary.arcpt produces output", {
  result <- ar_cpt(test_data, order = 1)
  expect_output(summary(result), "AR Changepoint Detection Summary")
  expect_output(summary(result), "Model Configuration")
  expect_output(summary(result), "Segment Coefficients")
})

test_that("summary.arcpt returns invisibly", {
  result <- ar_cpt(test_data, order = 1)
  output <- capture.output(summ <- summary(result))
  expect_true(is.list(summ))
  expect_true("order" %in% names(summ))
  expect_true("ncpts" %in% names(summ))
})

# =============================================================================
# Plot Method Tests
# =============================================================================

test_that("plot.arcpt creates a plot", {
  result <- ar_cpt(test_data, order = 1)
  expect_silent(plot(result))
})
test_that("plot.arcpt returns invisibly", {
  result <- ar_cpt(test_data, order = 1)
  expect_invisible(plot(result))
})

test_that("plot.arcpt works with show_segments FALSE", {
  result <- ar_cpt(test_data, order = 1)
  expect_silent(plot(result, show_segments = FALSE))
})

# =============================================================================
# Coef Method Tests
# =============================================================================

test_that("coef.arcpt returns a list", {
  result <- ar_cpt(test_data, order = 1)
  coefs <- coef(result)

  expect_true(is.list(coefs))
  expect_true(length(coefs) >= 1)
  expect_true(all(grepl("segment_", names(coefs))))
})

# =============================================================================
# cpts.arcpt Tests
# =============================================================================

test_that("cpts.arcpt returns changepoint locations", {
  result <- ar_cpt(test_data, order = 1)
  cpts_result <- cpts.arcpt(result)

  expect_true(is.numeric(cpts_result))
  expect_equal(cpts_result, result$cpts)
})

# =============================================================================
# ncpts.arcpt Tests
# =============================================================================

test_that("ncpts.arcpt returns number of changepoints", {
  result <- ar_cpt(test_data, order = 1)
  ncpts_result <- ncpts.arcpt(result)

  expect_true(is.numeric(ncpts_result))
  expect_equal(ncpts_result, result$ncpts)
  expect_equal(ncpts_result, length(result$cpts))
})
