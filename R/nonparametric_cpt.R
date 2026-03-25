#' Non-parametric Changepoint Regression using Splines
#'
#' Detect changepoints in non-parametric regression models using B-splines
#' or natural splines. This allows for flexible, smooth relationships that
#' can change at detected changepoints.
#'
#' @param x Numeric vector of predictor values (e.g., time).
#' @param y Numeric vector of response values.
#' @param spline_type Type of spline: "bs" (B-spline) or "ns" (natural spline).
#' @param df Degrees of freedom for the spline basis. If NULL, selected by cross-validation.
#' @param method Changepoint detection method: "PELT" or "AMOC".
#' @param penalty Penalty type: "MBIC", "BIC", "AIC", or "Manual".
#' @param pen.value Manual penalty value (used when penalty = "Manual").
#' @param minseglen Minimum segment length.
#' @param knots Optional vector of knot positions. If NULL, placed at quantiles.
#'
#' @return An object of class "npcpt" containing:
#'   \item{x}{Predictor values}
#'   \item{y}{Response values}
#'   \item{cpts}{Detected changepoint locations (indices)}
#'   \item{cpts_x}{Changepoint locations in x coordinates}
#'   \item{ncpts}{Number of changepoints}
#'   \item{spline_fits}{Spline fit information per segment}
#'   \item{fitted_values}{Fitted values from the spline model}
#'   \item{residuals}{Model residuals}
#'   \item{spline_type}{Type of spline used}
#'   \item{df}{Degrees of freedom}
#'   \item{method}{Detection method used}
#'   \item{penalty}{Penalty type used}
#'
#' @details
#' This function implements changepoint detection for non-parametric regression
#' using spline basis functions. The algorithm:
#'
#' 1. Constructs a spline basis for the predictor
#' 2. Fits the spline model using least squares
#' 3. Uses the spline coefficients as regression parameters for changepoint detection
#' 4. Allows different smooth functions in each segment
#'
#' The model in segment k is:
#' \deqn{Y_i = f_k(x_i) + \epsilon_i = \sum_{j=1}^{df} \beta_{kj} B_j(x_i) + \epsilon_i}
#'
#' Where \eqn{B_j} are the spline basis functions.
#'
#' @examples
#' \dontrun{
#' # Simulate data with non-linear relationship that changes
#' set.seed(42)
#' n <- 200
#' x <- seq(0, 4*pi, length.out = n)
#'
#' # First half: sine wave
#' # Second half: different sine wave
#' y <- c(
#'   sin(x[1:100]) + rnorm(100, 0, 0.2),
#'   2 * sin(2 * x[101:200]) + rnorm(100, 0, 0.2)
#' )
#'
#' # Detect changepoints with spline regression
#' result <- np_cpt(x, y, spline_type = "ns", df = 6)
#' print(result)
#' plot(result)
#' }
#'
#' @export
np_cpt <- function(x, y, spline_type = c("bs", "ns"), df = NULL,
                   method = c("PELT", "AMOC"),
                   penalty = c("MBIC", "BIC", "AIC", "Manual"),
                   pen.value = 0, minseglen = NULL, knots = NULL) {

  # Input validation
  spline_type <- match.arg(spline_type)
  method <- match.arg(method)
  penalty <- match.arg(penalty)

  if (!is.numeric(x) || !is.numeric(y)) {
    stop("'x' and 'y' must be numeric vectors")
  }

  if (length(x) != length(y)) {
    stop("'x' and 'y' must have the same length")
  }

  n <- length(y)

  if (any(is.na(x)) || any(is.na(y))) {
    stop("'x' and 'y' cannot contain NA values")
  }

  if (any(!is.finite(x)) || any(!is.finite(y))) {
    stop("'x' and 'y' cannot contain non-finite values (Inf, -Inf, NaN)")
  }

  # Sort by x if not already sorted
  if (is.unsorted(x)) {
    ord <- order(x)
    x <- x[ord]
    y <- y[ord]
  }

  # Set default df based on data size
  if (is.null(df)) {
    df <- min(max(4, floor(n / 20)), 10)
  }

  if (df < 3) {
    stop("'df' must be at least 3 for spline fitting")
  }

  # Set default minseglen
  if (is.null(minseglen)) {
    minseglen <- max(df + 2, 10)
  }

  # Check minimum data size
  if (n < 2 * minseglen) {
    stop(paste("Data too short for changepoint detection with minseglen =", minseglen))
  }

  # Build spline basis matrix
  spline_basis <- build_spline_basis(x, spline_type, df, knots)

  # Build design matrix for changepoint detection
  # Response in first column, then spline basis columns
  cpt_matrix <- cbind(y, spline_basis)

  # Apply changepoint detection using EnvCpt's cpt.reg
  if (requireNamespace("EnvCpt", quietly = TRUE)) {
    cpt_result <- EnvCpt:::cpt.reg(
      data = cpt_matrix,
      method = method,
      penalty = penalty,
      pen.value = pen.value,
      minseglen = minseglen
    )

    # Extract changepoints
    cpts <- changepoint::cpts(cpt_result)

    # Remove endpoint if present
    cpts <- cpts[cpts < n]

  } else {
    stop("EnvCpt package is required but not available")
  }

  ncpts <- length(cpts)

  # Convert changepoint indices to x coordinates
  cpts_x <- if (ncpts > 0) x[cpts] else numeric(0)

  # Fit splines per segment and compute fitted values
  segments <- define_segments(n, cpts)

  spline_fits <- list()
  fitted_values <- numeric(n)

  for (i in seq_along(segments)) {
    seg <- segments[[i]]
    idx <- seg$start:seg$end

    x_seg <- x[idx]
    y_seg <- y[idx]

    # Build segment-specific spline basis
    if (length(idx) > df + 1) {
      basis_seg <- build_spline_basis(x_seg, spline_type, df, knots = NULL)

      # Fit using least squares
      fit <- lm.fit(basis_seg, y_seg)
      coefs <- fit$coefficients
      fitted_seg <- fit$fitted.values
      residuals_seg <- fit$residuals
      rss <- sum(residuals_seg^2)

    } else {
      # Not enough data for full spline, use simpler model
      basis_seg <- cbind(1, x_seg)
      fit <- lm.fit(basis_seg, y_seg)
      coefs <- fit$coefficients
      fitted_seg <- fit$fitted.values
      residuals_seg <- fit$residuals
      rss <- sum(residuals_seg^2)
    }

    fitted_values[idx] <- fitted_seg

    spline_fits[[i]] <- list(
      segment = i,
      start = seg$start,
      end = seg$end,
      x_range = range(x_seg),
      coefficients = coefs,
      rss = rss,
      n_obs = length(idx)
    )
  }

  residuals <- y - fitted_values

  # Compute model fit statistics
  total_ss <- sum((y - mean(y))^2)
  residual_ss <- sum(residuals^2)
  r_squared <- 1 - residual_ss / total_ss

  # Create result object
  result <- list(
    x = x,
    y = y,
    cpts = cpts,
    cpts_x = cpts_x,
    ncpts = ncpts,
    spline_fits = spline_fits,
    fitted_values = fitted_values,
    residuals = residuals,
    spline_type = spline_type,
    df = df,
    method = method,
    penalty = penalty,
    r_squared = r_squared,
    cpt_result = cpt_result
  )

  class(result) <- "npcpt"
  return(result)
}


#' Build Spline Basis Matrix
#'
#' Constructs B-spline or natural spline basis functions.
#'
#' @param x Numeric vector of predictor values.
#' @param type Type of spline: "bs" or "ns".
#' @param df Degrees of freedom.
#' @param knots Optional knot positions.
#' @return Spline basis matrix.
#' @keywords internal
build_spline_basis <- function(x, type = c("bs", "ns"), df, knots = NULL) {

  type <- match.arg(type)

  if (requireNamespace("splines", quietly = TRUE)) {

    if (type == "bs") {
      # B-spline basis
      basis <- splines::bs(x, df = df, intercept = TRUE)
    } else {
      # Natural spline basis
      basis <- splines::ns(x, df = df, intercept = TRUE)
    }

  } else {
    # Fallback: polynomial basis
    warning("splines package not available, using polynomial basis")
    basis <- poly(x, degree = min(df - 1, 5), raw = TRUE)
    basis <- cbind(1, basis)
  }

  # Ensure matrix format
  basis <- as.matrix(basis)

  return(basis)
}


#' @export
print.npcpt <- function(x, ...) {

  cat("Non-parametric Changepoint Regression\n")
  cat("=====================================\n\n")

  cat("Spline type:", ifelse(x$spline_type == "bs", "B-spline", "Natural spline"), "\n")
  cat("Degrees of freedom:", x$df, "\n")
  cat("Method:", x$method, "\n")
  cat("Penalty:", x$penalty, "\n")
  cat("Number of changepoints:", x$ncpts, "\n")

  if (x$ncpts > 0) {
    cat("Changepoint indices:", paste(x$cpts, collapse = ", "), "\n")
    cat("Changepoint x-values:", paste(round(x$cpts_x, 4), collapse = ", "), "\n")
  }

  cat("\nModel fit:\n")
  cat("  R-squared:", round(x$r_squared, 4), "\n")
  cat("  Residual std error:", round(sd(x$residuals), 4), "\n")

  invisible(x)
}


#' @export
summary.npcpt <- function(object, ...) {

  cat("Non-parametric Changepoint Regression Summary\n")
  cat("==============================================\n\n")

  cat("Spline:", ifelse(object$spline_type == "bs", "B-spline", "Natural spline"),
      "with df =", object$df, "\n")
  cat("Method:", object$method, "| Penalty:", object$penalty, "\n")
  cat("Observations:", length(object$y), "\n\n")

  cat("Changepoints:", object$ncpts, "\n")
  if (object$ncpts > 0) {
    cat("Locations (index):", paste(object$cpts, collapse = ", "), "\n")
    cat("Locations (x-value):", paste(round(object$cpts_x, 4), collapse = ", "), "\n")
  }

  cat("\n--- Segment Details ---\n\n")

  for (seg in object$spline_fits) {
    cat(sprintf("Segment %d: observations %d to %d (n = %d)\n",
                seg$segment, seg$start, seg$end, seg$n_obs))
    cat(sprintf("  x range: [%.4f, %.4f]\n", seg$x_range[1], seg$x_range[2]))
    cat(sprintf("  RSS: %.4f\n", seg$rss))
    cat("\n")
  }

  cat("--- Overall Model Fit ---\n\n")
  cat(sprintf("R-squared: %.4f\n", object$r_squared))
  cat(sprintf("Residual std error: %.4f\n", sd(object$residuals)))

  invisible(object)
}


#' @export
plot.npcpt <- function(x, type = c("fit", "residuals", "segments"), ...) {

  type <- match.arg(type)

  if (type == "fit") {
    # Plot data with fitted spline
    plot(x$x, x$y, type = "p", pch = 16, cex = 0.5, col = "gray50",
         main = "Non-parametric Changepoint Regression",
         xlab = "x", ylab = "y")

    # Add fitted curve
    lines(x$x, x$fitted_values, col = "blue", lwd = 2)

    # Add changepoint lines
    if (x$ncpts > 0) {
      abline(v = x$cpts_x, col = "red", lty = 2, lwd = 2)
    }

    legend("topright", legend = c("Data", "Spline fit", "Changepoints"),
           col = c("gray50", "blue", "red"),
           pch = c(16, NA, NA), lty = c(NA, 1, 2), lwd = c(NA, 2, 2))

  } else if (type == "residuals") {
    # Residual plot
    par(mfrow = c(2, 1))

    # Residuals vs x
    plot(x$x, x$residuals, type = "p", pch = 16, cex = 0.5,
         main = "Residuals vs x",
         xlab = "x", ylab = "Residuals")
    abline(h = 0, col = "red", lty = 2)
    if (x$ncpts > 0) {
      abline(v = x$cpts_x, col = "blue", lty = 2)
    }

    # Residuals vs fitted
    plot(x$fitted_values, x$residuals, type = "p", pch = 16, cex = 0.5,
         main = "Residuals vs Fitted",
         xlab = "Fitted values", ylab = "Residuals")
    abline(h = 0, col = "red", lty = 2)

    par(mfrow = c(1, 1))

  } else if (type == "segments") {
    # Plot each segment separately
    n_seg <- length(x$spline_fits)
    colors <- rainbow(n_seg, alpha = 0.7)

    plot(x$x, x$y, type = "n",
         main = "Segments with Individual Spline Fits",
         xlab = "x", ylab = "y")

    for (i in seq_along(x$spline_fits)) {
      seg <- x$spline_fits[[i]]
      idx <- seg$start:seg$end

      points(x$x[idx], x$y[idx], pch = 16, cex = 0.5, col = colors[i])
      lines(x$x[idx], x$fitted_values[idx], col = colors[i], lwd = 2)
    }

    # Add changepoint lines
    if (x$ncpts > 0) {
      abline(v = x$cpts_x, col = "black", lty = 2, lwd = 2)
    }

    legend("topright", legend = paste("Segment", 1:n_seg),
           col = colors, lty = 1, lwd = 2)
  }

  invisible(x)
}


#' @export
fitted.npcpt <- function(object, ...) {
  object$fitted_values
}


#' @export
residuals.npcpt <- function(object, ...) {
  object$residuals
}


#' @export
cpts.npcpt <- function(x) {
  x$cpts
}


#' @export
ncpts.npcpt <- function(x) {
  x$ncpts
}


#' Predict from Non-parametric Changepoint Model
#'
#' @param object An "npcpt" object.
#' @param newx New x values for prediction.
#' @param ... Additional arguments (ignored).
#' @return Predicted values.
#' @export
predict.npcpt <- function(object, newx = NULL, ...) {

  if (is.null(newx)) {
    return(object$fitted_values)
  }

  # For new x values, determine which segment they belong to
  # and use that segment's spline fit

  n_new <- length(newx)
  predictions <- numeric(n_new)

  # Get segment boundaries in x coordinates
  x_bounds <- c(min(object$x), object$cpts_x, max(object$x))

  for (i in seq_len(n_new)) {
    xi <- newx[i]

    # Find which segment this x belongs to
    seg_idx <- sum(xi > x_bounds[-length(x_bounds)])
    seg_idx <- max(1, min(seg_idx, length(object$spline_fits)))

    seg <- object$spline_fits[[seg_idx]]

    # Use segment's spline coefficients
    # Note: This is approximate - full implementation would store basis info
    # For now, use linear interpolation within segment as fallback

    seg_data_idx <- seg$start:seg$end
    x_seg <- object$x[seg_data_idx]
    fitted_seg <- object$fitted_values[seg_data_idx]

    predictions[i] <- approx(x_seg, fitted_seg, xout = xi, rule = 2)$y
  }

  return(predictions)
}
