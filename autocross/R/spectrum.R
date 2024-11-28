#' Estimate spectrum...
#'
#' This function ...
#'
#' @param x  Data frame with...
#' @param y  Data frame...
#' @param config  Configuration...
#' @return named list...
#' @export
spectrum <- function(x, y, config) {
  print('config:')
  str(config)
  dim_result <- .Fortran("rx_setdim",
                     nx = nrow(x),
                     ny = nrow(y),
                     tx = x$time,
                     ty = y$time,
                     ofac = config$ofac,
                     hifac = config$hifac,
                     n50 = as.integer(config$n50),
                     nout = integer(1)
                     )
  nout = dim_result$nout
  print('dim_result:')
  str(dim_result)
  print(paste('n50=', config$n50))
  str(y$time)
  print('nout:')
  str(nout)
  result <- .Fortran("rx_subroutine",
                     nx = as.integer(nrow(x)),
                     ny = as.integer(nrow(y)),
                     nout = as.integer(nout),
                     tx = x$time,
                     x = x$series,
                     ty = y$time,
                     y = y$series,
                     nsim = config$nsim,
                     ofac = config$ofac,
                     hifac = config$hifac,
                     n50 = config$n50,
                     alpha = config$alpha,
                     iwin = config$iwin,
                     rhox = double(1),
                     rhoy = double(1),
                     taux = double(1),
                     tauy = double(1),
                     dof = double(1),
                     db6 = double(1),
                     false_alarm = double(1),
                     faccritx = double(1),
                     faccrity = double(1),
                     alphacritx = double(1),
                     alphacrity = double(1),
                     data_x = matrix(0, nout, 12),
                     data_y = matrix(0, nout, 12),
                     data_xy = matrix(0, nout, 2),
                     data_cxy = matrix(0, nout, 7),
                     data_phxy = matrix(0, nout, 6)
                     )
  return(result)
}
