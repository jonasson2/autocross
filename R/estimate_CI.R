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
  n = length(time)
  # Validate inputs
  if (!is.numeric(time) || n < 2) {
    rlang::abort("Time must be a numeric vector with at least two elements.")
  }
  if (!is.numeric(x.series) || length(x.series) != n) {
    rlang::abort("X-series must be a numeric vector of the same length as time.")
  }
  if (!is.numeric(y.series) || length(y.series) != n) {
    rlang::abort("Y-series must be a numeric vector of the same length as time.")
  }
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    rlang::abort("Alpha must be a value between 0 and 1.")
  }

  # Call the Fortran subroutine
  print(paste('n=', n))
  print(paste('time=', time))
  result <- .Fortran("pearsont3sub", n, time, x.series, y.series, alpha,
                     r = double(1), ci = double(2), taux = double(1),
                     tauy = double(1))

  # Pack results into value output
  # result = list(
  #   r = 0.5,
  #   ci = c(0.4, 0.6),
  #   taux = 3.0,
  #   tauy = 4.0
  # )
  list(
    correlation = result$r,
    confidence_interval = result$ci,
    persistence_time_x = result$taux,
    persistence_time_y = result$tauy
  )
}
