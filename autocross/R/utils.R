# triangle_average.R

#' Compute triangular weighted average of vector x with apex at index k
#'
#' @param x numeric vector of length n
#' @param k integer in [1, n] indicating the index of apex
#' @return scalar weighted average
#’ @keywords internal
triangle_average <- function(x, k) {
  n <- length(x)
  if (k < 1 || k > n) stop("k must be between 1 and length(x)")
  positions <- seq_len(n)
  weights_raw <- ifelse(
    positions <= k,
    positions / k,
    (n - positions + 1) / (n - k + 1)
  )
  weights <- weights_raw / sum(weights_raw)
  sum(weights * x)
}

#' Resample a variable by triangular averaging over moving windows defined by edges
#'
#' @param yr integer vector of all years (from start to end)
#' @param df data frame containing columns 'yr' and the target variable
#' @param edges integer vector of breakpoints (must lie within yr range)
#' @param variable string, name of the variable in df to resample
#' @return data frame with 'yr' (midpoints) and the resampled variable
#' @examples
#' df <- data.frame(yr = c(2000, 2005, 2010), value = c(1, 3, 2))
#' triangle_resampling(2000:2010, df, c(2000, 2005, 2010), "value")
triangle_resampling <- function(yr, df, edges, variable) {
  # a) interpolate the target variable over full year sequence
  x <- approx(x = df$yr, y = df[[variable]], xout = yr, rule=2)$y
  
  # b) compute edge indices within yr
  edgeidx <- edges - edges[1] + 1
  
  # c) compute midpoints between successive edges and prepare output
  midpoints <- (head(edges, -1) + tail(edges, -1)) / 2
  n_mid     <- length(midpoints)
  newdf     <- data.frame(yr = midpoints)
  newdf[[variable]] <- numeric(n_mid)
  n_edges   <- length(edges)
  
  # d) for each midpoint, define window and apply triangle_average
  for (i in seq_len(n_mid)) {
    # determine the inclusive index range for this window
    idx <- if (i == 1L)          seq(edgeidx[1],           edgeidx[2])
    else if (i == n_mid)  seq(edgeidx[n_edges - 1], edgeidx[n_edges])
    else                  seq(edgeidx[i - 1],       edgeidx[i + 2])
    
    # apex position within the window: midpoint of the two edge indices
    mid_idx <- (edgeidx[i] + edgeidx[i + 1]) %/% 2
    k       <- mid_idx - idx[1] + 1
    
    # compute triangular average
    newdf[i, variable] <- triangle_average(x[idx], k)
  }
  
  newdf
}