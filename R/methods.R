#' Print Method for arcpt Objects
#'
#' @param x An object of class "arcpt".
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.arcpt <- function(x, ...) {
  cat("AR Changepoint Detection Results\n")
  cat("================================\n\n")

  cat("Model: AR(", x$order, ")", sep = "")
  if (x$trend) cat(" with trend")
  cat("\n")

  cat("Method:", x$method, "\n")
  cat("Penalty:", x$penalty, "\n")
  cat("Data length:", length(x$data), "\n\n")

  if (x$ncpts == 0) {
    cat("No changepoints detected.\n")
  } else {
    cat("Number of changepoints:", x$ncpts, "\n")
    cat("Changepoint locations:", paste(x$cpts, collapse = ", "), "\n")
  }

  invisible(x)
}


#' Summary Method for arcpt Objects
#'
#' @param object An object of class "arcpt".
#' @param ... Additional arguments (ignored).
#'
#' @return A list containing summary information (invisibly).
#'
#' @export
summary.arcpt <- function(object, ...) {
  cat("AR Changepoint Detection Summary\n")
  cat("================================\n\n")

  cat("Model Configuration:\n")
  cat("  - AR order:", object$order, "\n")
  cat("  - Trend:", if(object$trend) "Yes" else "No", "\n")
  cat("  - Method:", object$method, "\n")
  cat("  - Penalty:", object$penalty, "\n")
  cat("  - Min segment length:", object$minseglen, "\n\n")

  cat("Data:\n")
  cat("  - Length:", length(object$data), "\n")
  cat("  - Range: [", round(min(object$data), 3), ", ",
      round(max(object$data), 3), "]\n\n", sep = "")

  cat("Results:\n")
  cat("  - Changepoints detected:", object$ncpts, "\n")

  if (object$ncpts > 0) {
    cat("  - Changepoint locations:", paste(object$cpts, collapse = ", "), "\n")
  }

  cat("\nSegment Coefficients:\n")
  cat("---------------------\n")

  for (seg in object$coefficients) {
    cat(sprintf("Segment %d (t=%d to %d):\n", seg$segment, seg$start, seg$end))

    if (!is.null(seg$coefficients)) {
      for (j in seq_along(seg$coefficients)) {
        cat(sprintf("  %s = %.4f\n", names(seg$coefficients)[j], seg$coefficients[j]))
      }
      if (!is.null(seg$sigma2)) {
        cat(sprintf("  sigma^2 = %.4f\n", seg$sigma2))
      }
    } else if (!is.null(seg$note)) {
      cat(sprintf("  Note: %s\n", seg$note))
    }
    cat("\n")
  }

  invisible(list(
    order = object$order,
    trend = object$trend,
    ncpts = object$ncpts,
    cpts = object$cpts,
    coefficients = object$coefficients
  ))
}


#' Plot Method for arcpt Objects
#'
#' Creates a visualization of the time series with detected changepoints.
#'
#' @param x An object of class "arcpt".
#' @param show_segments Logical indicating whether to color-code segments.
#'   Default is TRUE.
#' @param show_fit Logical indicating whether to show fitted values per segment.
#'   Default is FALSE.
#' @param ... Additional arguments passed to plot().
#'
#' @return Invisibly returns the input object.
#'
#' @export
plot.arcpt <- function(x, show_segments = TRUE, show_fit = FALSE, ...) {
  # Extract data
  data <- x$data
  n <- length(data)
  cpts <- x$cpts
  all_cpts <- c(0, cpts, n)

  # Set up colors for segments
  n_segments <- length(all_cpts) - 1
  segment_colors <- grDevices::rainbow(n_segments, alpha = 0.7)

  # Create main plot
  plot(data, type = "n",
       main = paste0("AR(", x$order, ") Changepoint Detection",
                     if(x$trend) " with Trend" else ""),
       xlab = "Time", ylab = "Value", ...)

  # Plot segments with different colors
  if (show_segments && n_segments > 1) {
    for (i in seq_len(n_segments)) {
      start_idx <- all_cpts[i] + 1
      end_idx <- all_cpts[i + 1]
      graphics::lines(start_idx:end_idx, data[start_idx:end_idx],
            col = segment_colors[i], lwd = 1.5)
    }
  } else {
    graphics::lines(data, col = "steelblue", lwd = 1.5)
  }

  # Add changepoint lines
  if (length(cpts) > 0) {
    graphics::abline(v = cpts, col = "red", lwd = 2, lty = 2)
  }

  # Add legend
  if (length(cpts) > 0) {
    graphics::legend("topright",
           legend = c("Data", "Changepoints"),
           col = c("steelblue", "red"),
           lty = c(1, 2),
           lwd = c(1.5, 2),
           bg = "white",
           cex = 0.8)
  }

  invisible(x)
}


#' Extract Changepoint Locations
#'
#' @param x An object of class "arcpt".
#'
#' @return Integer vector of changepoint locations.
#'
#' @export
cpts.arcpt <- function(x) {
  return(x$cpts)
}


#' Extract Number of Changepoints
#'
#' @param x An object of class "arcpt".
#'
#' @return Integer indicating the number of changepoints.
#'
#' @export
ncpts.arcpt <- function(x) {
  return(x$ncpts)
}


#' Extract Coefficients from arcpt Objects
#'
#' @param object An object of class "arcpt".
#' @param ... Additional arguments (ignored).
#'
#' @return A list of coefficient vectors per segment.
#'
#' @export
coef.arcpt <- function(object, ...) {
  coefs <- lapply(object$coefficients, function(seg) {
    if (!is.null(seg$coefficients)) {
      return(seg$coefficients)
    } else {
      return(NA)
    }
  })
  names(coefs) <- paste0("segment_", seq_along(coefs))
  return(coefs)
}
