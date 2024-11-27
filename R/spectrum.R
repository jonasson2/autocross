#' Estimate spectrum...
#'
#' This function ...
#'
#' @param ...
#' @export
spectrum <- function(x, y, config) {
  result <- .Fortran("rx_subroutine",
                     nx = nrow(x),
                     ny = nrow(y),
                     tx = x$time,
                     x = x$series,
                     ty = y$time,
                     y = y$series,
                     nsim = config$nsim,
                     ofac = config$ofac,
                     hifac = config$hifac,
                     n50 = config$n50,
                     alpha = config$alpha,
                     iwin = config$iwin)
  str(result)
}
