#' Changepoint Detection in AR Time Series
#'
#' Detects changepoints in autoregressive (AR) time series. Supports AR(1) and
#' AR(2) models with optional trend components.
#'
#' @param data A numeric vector containing the time series data.
#' @param order Integer specifying the AR order. Must be 1 or 2.
#' @param trend Logical indicating whether to include a linear trend component.
#'   Default is FALSE.
#' @param method Character string specifying the changepoint detection method.
#'   Either "PELT" (default) or "AMOC" (At Most One Change).
#' @param penalty Character string specifying the penalty type. Options include
#'   "MBIC" (default), "BIC", "SIC", "AIC", "Hannan-Quinn", "Manual".
#' @param pen.value Numeric value for manual penalty. Only used when
#'   penalty = "Manual".
#' @param minseglen Integer specifying the minimum segment length. Default is
#'   calculated based on AR order and trend.
#'
#' @return An object of class "arcpt" containing:
#' \itemize{
#'   \item \code{data} - The original time series
#'   \item \code{order} - The AR order used
#'   \item \code{trend} - Whether trend was included
#'   \item \code{cpts} - Vector of detected changepoint locations
#'   \item \code{ncpts} - Number of changepoints detected
#'   \item \code{coefficients} - List of AR coefficients per segment
#'   \item \code{method} - The method used
#'   \item \code{penalty} - The penalty type used
#'   \item \code{cpt.reg.result} - The underlying cpt.reg object
#' }
#'
#' @details
#' This function uses the changepoint regression framework from the EnvCpt
#' package to detect changes in AR structure. The data is transformed into
#' a regression problem where the current observation is regressed on lagged
#' values (and optionally time for trend).
#'
#' For AR(1): \eqn{y_t = \mu + \phi_1 y_{t-1} + \epsilon_t}
#'
#' For AR(2): \eqn{y_t = \mu + \phi_1 y_{t-1} + \phi_2 y_{t-2} + \epsilon_t}
#'
#' With trend: \eqn{y_t = \mu + \beta t + \phi_1 y_{t-1} + ... + \epsilon_t}
#'
#' @examples
#' # Generate AR(1) data with a changepoint
#' set.seed(42)
#' x1 <- arima.sim(model = list(ar = 0.3), n = 100)
#' x2 <- arima.sim(model = list(ar = 0.9), n = 100)
#' data <- c(x1, x2)
#'
#' # Detect changepoints
#' result <- ar_cpt(data, order = 1)
#' print(result)
#'
#' # With trend
#' result_trend <- ar_cpt(data, order = 1, trend = TRUE)
#'
#' # AR(2) model
#' result_ar2 <- ar_cpt(data, order = 2)
#'
#' @export
ar_cpt <- function(data,
                   order = 1,
                   trend = FALSE,
                   method = c("PELT", "AMOC"),
                   penalty = c("MBIC", "BIC", "SIC", "AIC", "Hannan-Quinn", "Manual"),
                   pen.value = 0,
                   minseglen = NULL) {


  # ==========================================================================
  # Input Validation
  # ==========================================================================


  # Validate data

if (!is.numeric(data)) {
    stop("Argument 'data' must be a numeric vector.")
  }
  if (!is.vector(data) && !is.ts(data)) {
    stop("Argument 'data' must be a vector or time series object.")
  }
  if (any(is.na(data))) {
    stop("Argument 'data' contains NA values. Please remove or impute missing values.")
  }
  if (any(!is.finite(data))) {
    stop("Argument 'data' contains non-finite values (Inf or NaN).")
  }

  n <- length(data)

  # Validate order
  if (!is.numeric(order) || length(order) != 1) {
    stop("Argument 'order' must be a single integer (1 or 2).")
  }
  order <- as.integer(order)
  if (order < 1 || order > 2) {
    stop("Argument 'order' must be 1 or 2.")
  }

  # Validate trend
  if (!is.logical(trend) || length(trend) != 1) {
    stop("Argument 'trend' must be TRUE or FALSE.")
  }

  # Validate method
  method <- match.arg(method)

  # Validate penalty
  penalty <- match.arg(penalty)

  # Validate pen.value
  if (!is.numeric(pen.value) || length(pen.value) != 1) {
    stop("Argument 'pen.value' must be a single numeric value.")
  }
  if (pen.value < 0) {
    stop("Argument 'pen.value' must be non-negative.")
  }

  # Validate and set minseglen
  # Minimum segment length should be at least number of parameters + 1
  min_required <- order + 1 + as.integer(trend) + 1
  if (is.null(minseglen)) {
    minseglen <- max(min_required + 2, 5)  # Default: a bit more than minimum
  } else {
    if (!is.numeric(minseglen) || length(minseglen) != 1) {
      stop("Argument 'minseglen' must be a single positive integer.")
    }
    minseglen <- as.integer(minseglen)
    if (minseglen < min_required) {
      warning(sprintf("minseglen too small for AR(%d)%s model. Setting to %d.",
                      order, if(trend) " with trend" else "", min_required))
      minseglen <- min_required
    }
  }

  # Check data length
  if (n < 2 * minseglen) {
    stop(sprintf("Data length (%d) is too short. Need at least %d observations for minseglen=%d.",
                 n, 2 * minseglen, minseglen))
  }

  # ==========================================================================
  # Build Design Matrix
  # ==========================================================================

  data_matrix <- build_ar_design_matrix(data, order = order, trend = trend)

  # ==========================================================================
  # Run Changepoint Detection
  # ==========================================================================

  # Use EnvCpt's cpt.reg function
  cpt_result <- EnvCpt:::cpt.reg(
    data = data_matrix,
    method = method,
    penalty = penalty,
    pen.value = pen.value,
    minseglen = minseglen
  )

  # Extract changepoints
  detected_cpts <- changepoint::cpts(cpt_result)
  # Adjust for removed observations due to lags
  detected_cpts_adjusted <- detected_cpts + order
  # Remove the final point (always equals n)
  detected_cpts_adjusted <- detected_cpts_adjusted[detected_cpts_adjusted < n]

  # ==========================================================================
  # Estimate Coefficients per Segment
  # ==========================================================================

  coefficients <- estimate_segment_coefficients(data, detected_cpts_adjusted, order, trend)

  # ==========================================================================
  # Build Result Object
  # ==========================================================================

  result <- list(
    data = data,
    order = order,
    trend = trend,
    cpts = detected_cpts_adjusted,
    ncpts = length(detected_cpts_adjusted),
    coefficients = coefficients,
    method = method,
    penalty = penalty,
    pen.value = pen.value,
    minseglen = minseglen,
    cpt.reg.result = cpt_result
  )

  class(result) <- "arcpt"
  return(result)
}


#' Build AR Design Matrix
#'
#' Internal function to construct the design matrix for AR changepoint detection.
#'
#' @param data Numeric vector of time series data.
#' @param order AR order (1 or 2).
#' @param trend Logical indicating whether to include trend.
#'
#' @return A matrix suitable for cpt.reg().
#'
#' @keywords internal
build_ar_design_matrix <- function(data, order, trend) {
  n <- length(data)

  if (order == 1) {
    # AR(1): y[t] ~ intercept (+ trend) + y[t-1]
    response <- data[-1]  # y[2:n]
    intercept <- rep(1, n - 1)
    lag1 <- data[-n]  # y[1:(n-1)]

    if (trend) {
      time_index <- 2:n
      data_matrix <- cbind(response, intercept, time_index, lag1)
    } else {
      data_matrix <- cbind(response, intercept, lag1)
    }

  } else if (order == 2) {
    # AR(2): y[t] ~ intercept (+ trend) + y[t-1] + y[t-2]
    response <- data[-(1:2)]  # y[3:n]
    intercept <- rep(1, n - 2)
    lag1 <- data[-c(1, n)]  # y[2:(n-1)]
    lag2 <- data[-c(n - 1, n)]  # y[1:(n-2)]

    if (trend) {
      time_index <- 3:n
      data_matrix <- cbind(response, intercept, time_index, lag1, lag2)
    } else {
      data_matrix <- cbind(response, intercept, lag1, lag2)
    }
  }

  return(data_matrix)
}


#' Estimate Segment Coefficients
#'
#' Internal function to estimate AR coefficients for each segment.
#'
#' @param data Original time series data.
#' @param cpts Vector of changepoint locations.
#' @param order AR order.
#' @param trend Whether trend was included.
#'
#' @return List of coefficient vectors per segment.
#'
#' @keywords internal
estimate_segment_coefficients <- function(data, cpts, order, trend) {
  n <- length(data)
  all_cpts <- c(0, cpts, n)
  n_segments <- length(all_cpts) - 1

  coefficients <- vector("list", n_segments)

  for (i in seq_len(n_segments)) {
    start_idx <- all_cpts[i] + 1
    end_idx <- all_cpts[i + 1]

    segment_data <- data[start_idx:end_idx]
    segment_length <- length(segment_data)

    # Need enough data to fit AR model
    min_required <- order + 3

    if (segment_length >= min_required) {
      # Fit AR model using stats::arima
      ar_fit <- tryCatch(
        stats::arima(segment_data, order = c(order, 0, 0), method = "CSS-ML"),
        error = function(e) NULL
      )

      if (!is.null(ar_fit)) {
        coefs <- ar_fit$coef
        names(coefs) <- gsub("ar", "phi", names(coefs))
        coefficients[[i]] <- list(
          segment = i,
          start = start_idx,
          end = end_idx,
          coefficients = coefs,
          sigma2 = ar_fit$sigma2
        )
      } else {
        coefficients[[i]] <- list(
          segment = i,
          start = start_idx,
          end = end_idx,
          coefficients = NULL,
          sigma2 = NULL,
          note = "Failed to estimate coefficients"
        )
      }
    } else {
      coefficients[[i]] <- list(
        segment = i,
        start = start_idx,
        end = end_idx,
        coefficients = NULL,
        sigma2 = NULL,
        note = "Segment too short for coefficient estimation"
      )
    }
  }

  return(coefficients)
}
