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
  print(paste('n=', n))
  # Call the Fortran subroutine
  print('Calling Fortran')
  print(typeof(n))
  print(typeof(time))
  print(typeof(x.series))
  print(typeof(y.series))
  print(typeof(alpha))
  flush.console()
  #result = .Fortran('ss', x = x.series, s = double(1))
  result <- .Fortran("pt3", n = n, time = time, x = x.series,
                     y = y.series, alpha = alpha,
                     r = double(1), ci = double(2), taux = double(1),
                     tauy = double(1))
  str(result)
  result$s
  # list(
  #   correlation = result$r,
  #   confidence_interval = result$ci,
  #   persistence_time_x = result$taux,
  #   persistence_time_y = result$tauy
  # )
}
