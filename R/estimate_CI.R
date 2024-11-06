#' Estimate Confidence Interval
#'
#' This function estimates the confidence interval using the Fortran subroutine.
#' It serves as the user-facing interface for the package.
#'
#' @param time A numeric vector containing the time points.
#' @param x.series A numeric vector containing the x-series data.
#' @param y.series A numeric vector containing the y-series data.
#' @param alpha The significance level for the confidence interval (default is 0.05).
#' @return A list containing the estimated correlation, confidence interval values, and persistence times.
#' @export
estimate_CI <- function(time, x.series, y.series, alpha = 0.05) {
  # Validate inputs
  if (!is.numeric(time) || length(time) < 2) {
    rlang::abort("Time must be a numeric vector with at least two elements.")
  }
  if (!is.numeric(x.series) || length(x.series) != length(time)) {
    rlang::abort("X-series must be a numeric vector of the same length as time.")
  }
  if (!is.numeric(y.series) || length(y.series) != length(time)) {
    rlang::abort("Y-series must be a numeric vector of the same length as time.")
  }
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    rlang::abort("Alpha must be a value between 0 and 1.")
  }

  # Call the Fortran subroutine
  result <- .Fortran("_pearsont3sub",
                     as.numeric(time),
                     as.numeric(x.series),
                     as.numeric(y.series),
                     as.numeric(alpha),
                     r = double(1),
                     ci = double(2),
                     taux = double(1),
                     tauy = double(1))

  # Pack results into value output
  list(
    correlation = result$r,
    confidence_interval = result$ci,
    persistence_time_x = result$taux,
    persistence_time_y = result$tauy
  )
}
