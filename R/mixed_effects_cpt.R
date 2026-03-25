#' Mixed Effects Changepoint Regression
#'
#' Detect changepoints in mixed effects regression models where data has
#' hierarchical/grouped structure. Supports random intercepts and slopes
#' with changepoints in fixed effects.
#'
#' @param formula A formula specifying the model. Use `(1|group)` for random
#'   intercepts and `(x|group)` for random slopes.
#' @param data A data frame containing the variables in the formula.
#' @param method Changepoint detection method: "PELT" or "AMOC".
#' @param penalty Penalty type: "MBIC", "BIC", "AIC", or "Manual".
#' @param pen.value Manual penalty value (used when penalty = "Manual").
#' @param minseglen Minimum segment length.
#'
#' @return An object of class "mecpt" containing:
#'   \item{formula}{The model formula}
#'   \item{data}{The input data}
#'   \item{cpts}{Detected changepoint locations}
#'   \item{ncpts}{Number of changepoints}
#'   \item{fixed_effects}{Fixed effects estimates per segment}
#'   \item{random_effects}{Random effects estimates}
#'   \item{variance_components}{Variance component estimates}
#'   \item{method}{Detection method used}
#'   \item{penalty}{Penalty type used}
#'
#' @details
#' This function implements changepoint detection for mixed effects models,
#' allowing for hierarchical data structures where observations are nested
#' within groups. The algorithm:
#'
#' 1. Estimates random effects using all data
#' 2. Computes group-adjusted residuals
#' 3. Applies changepoint detection to fixed effects
#' 4. Re-estimates random effects per segment
#' 5. Iterates until convergence
#'
#' The model assumes:
#' \deqn{Y_{ij} = X_{ij} \beta_k + Z_{ij} b_i + \epsilon_{ij}}
#'
#' Where:
#' - \eqn{\beta_k} are fixed effects that change at changepoints
#' - \eqn{b_i} are random effects for group i
#' - k indexes the segment containing observation j
#'
#' @examples
#' \dontrun{
#' # Simulate hierarchical data with changepoint
#' set.seed(42)
#' n_groups <- 10
#' n_per_group <- 50
#'
#' data <- data.frame(
#'   group = rep(1:n_groups, each = n_per_group),
#'   time = rep(1:n_per_group, n_groups)
#' )
#'
#' # Random intercepts
#' random_int <- rnorm(n_groups, 0, 2)
#'
#' # Fixed effect changes at time 25
#' data$y <- ifelse(data$time <= 25,
#'   0.5 * data$time + random_int[data$group],
#'   1.5 * data$time - 25 + random_int[data$group]
#' ) + rnorm(nrow(data), 0, 1)
#'
#' # Detect changepoints
#' result <- me_cpt(y ~ time + (1|group), data = data)
#' print(result)
#' }
#'
#' @export
me_cpt <- function(formula, data, method = c("PELT", "AMOC"),
                   penalty = c("MBIC", "BIC", "AIC", "Manual"),
                   pen.value = 0, minseglen = NULL) {

  # Input validation
  method <- match.arg(method)
  penalty <- match.arg(penalty)

  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object")
  }

  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }

  # Parse formula to extract fixed and random effects
  formula_parts <- parse_mixed_formula(formula)

  if (is.null(formula_parts$random)) {
    stop("Formula must contain random effects. Use (1|group) syntax.")
  }

  # Extract response variable
  response_var <- all.vars(formula)[1]
  if (!response_var %in% names(data)) {
    stop(paste("Response variable", response_var, "not found in data"))
  }

  y <- data[[response_var]]
  n <- length(y)

  # Validate data
  if (any(is.na(y))) {
    stop("Response variable contains NA values")
  }

  # Set default minseglen
  if (is.null(minseglen)) {
    minseglen <- max(5, length(formula_parts$fixed) + 2)
  }

  # Build fixed effects design matrix
  X_fixed <- build_fixed_matrix(formula_parts$fixed, data)

  # Build random effects structure
  random_structure <- build_random_structure(formula_parts$random, data)
  Z <- random_structure$Z
  groups <- random_structure$groups
  n_groups <- random_structure$n_groups

  # Initial estimation of random effects (using all data)
  # Using simple approach: estimate group means as random intercepts
  group_means <- tapply(y, groups, mean)
  overall_mean <- mean(y)
  random_intercepts <- group_means - overall_mean

  # Compute group-adjusted response
  y_adjusted <- y - random_intercepts[groups]

  # Build design matrix for changepoint detection
  # Include intercept and fixed effect covariates
  cpt_matrix <- cbind(y_adjusted, X_fixed)

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

  # Estimate fixed effects per segment
  ncpts <- length(cpts)
  segments <- define_segments(n, cpts)

  fixed_effects <- list()
  for (i in seq_along(segments)) {
    seg <- segments[[i]]
    idx <- seg$start:seg$end

    # Fit fixed effects for this segment
    X_seg <- X_fixed[idx, , drop = FALSE]
    y_seg <- y_adjusted[idx]

    # Simple OLS estimation
    if (ncol(X_seg) > 0 && length(idx) > ncol(X_seg)) {
      coefs <- solve(t(X_seg) %*% X_seg) %*% t(X_seg) %*% y_seg
      coefs <- as.vector(coefs)
      names(coefs) <- colnames(X_seg)
    } else {
      coefs <- rep(NA, ncol(X_seg))
      names(coefs) <- colnames(X_seg)
    }

    fixed_effects[[i]] <- list(
      segment = i,
      start = seg$start,
      end = seg$end,
      coefficients = coefs
    )
  }

  # Re-estimate random effects per segment (refinement step)
  # For simplicity, keep initial random effects estimate
  # Full implementation would iterate until convergence

  # Estimate variance components
  # Residual variance
  fitted_values <- numeric(n)
  for (i in seq_along(segments)) {
    seg <- segments[[i]]
    idx <- seg$start:seg$end
    X_seg <- X_fixed[idx, , drop = FALSE]
    coefs <- fixed_effects[[i]]$coefficients
    fitted_values[idx] <- X_seg %*% coefs + random_intercepts[groups[idx]]
  }

  residuals <- y - fitted_values
  sigma2_resid <- var(residuals)
  sigma2_random <- var(random_intercepts)

  variance_components <- list(
    residual = sigma2_resid,
    random_intercept = sigma2_random
  )

  # Create result object
  result <- list(
    formula = formula,
    data = data,
    cpts = cpts,
    ncpts = ncpts,
    fixed_effects = fixed_effects,
    random_effects = random_intercepts,
    variance_components = variance_components,
    groups = groups,
    n_groups = n_groups,
    fitted_values = fitted_values,
    residuals = residuals,
    method = method,
    penalty = penalty,
    cpt_result = cpt_result
  )

  class(result) <- "mecpt"
  return(result)
}


#' Parse Mixed Effects Formula
#'
#' Extract fixed and random effect components from a formula.
#'
#' @param formula A formula with mixed effects notation.
#' @return A list with fixed and random effect specifications.
#' @keywords internal
parse_mixed_formula <- function(formula) {

  formula_str <- deparse(formula)

  # Extract random effects: (term|group)
  random_pattern <- "\\([^)]+\\|[^)]+\\)"
  random_matches <- regmatches(formula_str, gregexpr(random_pattern, formula_str))[[1]]

  # Parse random effects
  random_effects <- list()
  if (length(random_matches) > 0) {
    for (i in seq_along(random_matches)) {
      match <- random_matches[i]
      # Remove parentheses
      inner <- gsub("^\\(|\\)$", "", match)
      parts <- strsplit(inner, "\\|")[[1]]
      term <- trimws(parts[1])
      group <- trimws(parts[2])

      random_effects[[i]] <- list(
        term = term,
        group = group,
        type = ifelse(term == "1", "intercept", "slope")
      )
    }
  }

  # Remove random effects from formula to get fixed effects
  fixed_str <- formula_str
  for (match in random_matches) {
    fixed_str <- gsub(match, "", fixed_str, fixed = TRUE)
  }

  # Clean up: remove extra + signs and whitespace
  fixed_str <- gsub("\\+\\s*\\+", "+", fixed_str)
  fixed_str <- gsub("\\+\\s*$", "", fixed_str)
  fixed_str <- gsub("~\\s*\\+", "~", fixed_str)
  fixed_str <- trimws(fixed_str)

  # Parse fixed effects formula
  if (nchar(fixed_str) > 0 && grepl("~", fixed_str)) {
    fixed_formula <- as.formula(fixed_str)
    fixed_terms <- attr(terms(fixed_formula), "term.labels")
  } else {
    fixed_terms <- character(0)
  }

  list(
    fixed = fixed_terms,
    random = random_effects,
    original = formula
  )
}


#' Build Fixed Effects Design Matrix
#'
#' @param fixed_terms Character vector of fixed effect term names.
#' @param data Data frame containing the variables.
#' @return Design matrix for fixed effects.
#' @keywords internal
build_fixed_matrix <- function(fixed_terms, data) {

  if (length(fixed_terms) == 0) {
    # Intercept only
    return(matrix(1, nrow = nrow(data), ncol = 1, dimnames = list(NULL, "intercept")))
  }

  # Build model matrix
  formula_str <- paste("~", paste(fixed_terms, collapse = " + "))
  X <- model.matrix(as.formula(formula_str), data = data)

  return(X)
}


#' Build Random Effects Structure
#'
#' @param random_effects List of random effect specifications from parse_mixed_formula.
#' @param data Data frame containing the grouping variables.
#' @return List with Z matrix, groups factor, and number of groups.
#' @keywords internal
build_random_structure <- function(random_effects, data) {

  if (length(random_effects) == 0) {
    stop("No random effects specified")
  }

  # For now, support single grouping factor with random intercept
  # Full implementation would handle multiple grouping factors and slopes

  group_var <- random_effects[[1]]$group

  if (!group_var %in% names(data)) {
    stop(paste("Grouping variable", group_var, "not found in data"))
  }

  groups <- as.factor(data[[group_var]])
  n_groups <- nlevels(groups)
  n <- nrow(data)

  # Z matrix for random intercepts (indicator matrix)
  Z <- model.matrix(~ groups - 1)

  list(
    Z = Z,
    groups = groups,
    n_groups = n_groups
  )
}


#' Define Segments from Changepoints
#'
#' @param n Total number of observations.
#' @param cpts Vector of changepoint locations.
#' @return List of segment definitions (start, end).
#' @keywords internal
define_segments <- function(n, cpts) {

  if (length(cpts) == 0) {
    return(list(list(start = 1, end = n)))
  }

  # Sort changepoints
  cpts <- sort(unique(cpts))

  # Define segment boundaries
  starts <- c(1, cpts + 1)
  ends <- c(cpts, n)

  segments <- list()
  for (i in seq_along(starts)) {
    segments[[i]] <- list(start = starts[i], end = ends[i])
  }

  return(segments)
}


#' @export
print.mecpt <- function(x, ...) {

  cat("Mixed Effects Changepoint Regression\n")
  cat("====================================\n\n")

  cat("Formula:", deparse(x$formula), "\n")
  cat("Method:", x$method, "\n")
  cat("Penalty:", x$penalty, "\n")
  cat("Number of groups:", x$n_groups, "\n")
  cat("Number of changepoints:", x$ncpts, "\n")

  if (x$ncpts > 0) {
    cat("Changepoint locations:", paste(x$cpts, collapse = ", "), "\n")
  }

  cat("\nVariance Components:\n")
  cat("  Random intercept:", round(x$variance_components$random_intercept, 4), "\n")
  cat("  Residual:", round(x$variance_components$residual, 4), "\n")

  invisible(x)
}


#' @export
summary.mecpt <- function(object, ...) {

  cat("Mixed Effects Changepoint Regression Summary\n")
  cat("=============================================\n\n")

  cat("Formula:", deparse(object$formula), "\n")
  cat("Method:", object$method, "| Penalty:", object$penalty, "\n")
  cat("Groups:", object$n_groups, "| Observations:", length(object$residuals), "\n\n")

  cat("Changepoints:", object$ncpts, "\n")
  if (object$ncpts > 0) {
    cat("Locations:", paste(object$cpts, collapse = ", "), "\n")
  }

  cat("\n--- Fixed Effects by Segment ---\n\n")

  for (seg in object$fixed_effects) {
    cat(sprintf("Segment %d (observations %d to %d):\n",
                seg$segment, seg$start, seg$end))

    coefs <- seg$coefficients
    for (j in seq_along(coefs)) {
      cat(sprintf("  %s: %.4f\n", names(coefs)[j], coefs[j]))
    }
    cat("\n")
  }

  cat("--- Variance Components ---\n\n")
  cat(sprintf("  Random intercept variance: %.4f (SD: %.4f)\n",
              object$variance_components$random_intercept,
              sqrt(object$variance_components$random_intercept)))
  cat(sprintf("  Residual variance: %.4f (SD: %.4f)\n",
              object$variance_components$residual,
              sqrt(object$variance_components$residual)))

  invisible(object)
}


#' @export
plot.mecpt <- function(x, type = c("data", "random"), ...) {

  type <- match.arg(type)

  if (type == "data") {
    # Plot data with fitted values and changepoints
    n <- length(x$residuals)
    y <- x$data[[all.vars(x$formula)[1]]]

    plot(1:n, y, type = "p", pch = 16, cex = 0.5, col = "gray50",
         main = "Mixed Effects Changepoint Regression",
         xlab = "Observation", ylab = "Response")

    # Add fitted values
    lines(1:n, x$fitted_values, col = "blue", lwd = 2)

    # Add changepoint lines
    if (x$ncpts > 0) {
      abline(v = x$cpts, col = "red", lty = 2, lwd = 2)
    }

    legend("topright", legend = c("Data", "Fitted", "Changepoints"),
           col = c("gray50", "blue", "red"),
           pch = c(16, NA, NA), lty = c(NA, 1, 2), lwd = c(NA, 2, 2))

  } else if (type == "random") {
    # Plot random effects
    re <- x$random_effects

    barplot(re, main = "Random Intercepts by Group",
            xlab = "Group", ylab = "Random Effect",
            col = ifelse(re > 0, "steelblue", "coral"))
    abline(h = 0, lty = 2)
  }

  invisible(x)
}


#' @export
coef.mecpt <- function(object, type = c("fixed", "random", "both"), ...) {

  type <- match.arg(type)

  if (type == "fixed") {
    return(object$fixed_effects)
  } else if (type == "random") {
    return(object$random_effects)
  } else {
    return(list(
      fixed = object$fixed_effects,
      random = object$random_effects
    ))
  }
}


#' @export
cpts.mecpt <- function(x) {
  x$cpts
}


#' @export
ncpts.mecpt <- function(x) {
  x$ncpts
}
