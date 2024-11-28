#' Estimate spectrum...
#'
#' This function ...
#'
#' @param ...
#' @export
spectrum <- function(x, y, config) {
  nout = 141
  str(x$time)
  str(nout)
  str(config)
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
  str(result)
}
