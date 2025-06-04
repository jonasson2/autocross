# FUNCTION TO COMPUTE CORRELATION BETWEEN A PAIR OF TIME SERIES
#' @export
find_CI = function(df1, var1, df2, var2, beta, n) {
  require(pracma)
  ref_yr = min(df1$yr[1], df2$yr[1])
  min_yr = max(tail(df1$yr, 1), tail(df2$yr, 1))
  max_age = ref_yr - min_yr
  edges = generate_resampling_edges(ref_yr, max_age, beta, n)
  resamp1 = apply_resampling(df1, var1, edges)
  resamp2 = apply_resampling(df2, var2, edges)
  common.res = common_time(resamp1, var1, resamp2, var2)
  time = 2000 - common.res$yr
  x = drop(detrend(common.res$x))
  y = drop(detrend(common.res$y))
  CI = estimate_CI(time, x, y, alpha=0.01)
  res = data.frame(time=time, x=x, y=y)
  info = list(avg_int = max_age/n, 
              max_int = edges[n-1] - edges[n],
              min_int = edges[1] - edges[2])
  result = list(CI=CI, xy=res, max_age=max_age, info=info)
}

#' @export
simulate_shifted <- function(n, rho, max_shift) {
  result <- simulate_bivar(n, rho = rho)
  k = n - max_shift
  x = seq(0, k-1)*((n-1)/(k-1))
  xp = seq(0, n-1)
  yp = approx(x, result$x2[1:k], xout = xp)$y
  list(yr = 2000 - as.numeric(seq(1, n)), x = result$x1, y = yp)
}

simulate_bivar <- function(n, 
                           phi_vec = c(phi1 = 0.8, phi2 = 0.5), 
                           rho = 0.3, 
                           burnin = 100) {
  # n       : number of post–burn‐in observations you want
  # phi_vec : length-2 vector of AR(1) coefficients, named c(phi1, phi2)
  # rho     : desired correlation between X1_t and X2_t at lag 0
  # burnin  : how many initial draws to discard for stationarity
  
  stopifnot(length(phi_vec) == 2, abs(rho) < 1)
  # Stationary var of univariate AR(1): Var(X_i) = sigma_i^2/(1 - phi_i^2)
  # We set Var(X_i)=1 so that each marginal has unit variance ⇒
  # Var(ε_i) = 1 - phi_i^2
  sigma_vec <- sqrt(1 - phi_vec^2)
  
  # Build the innovation covariance 
  Sigma <- matrix(c(
    (1 - phi_vec[1]^2),  rho * (1 - phi_vec[1] * phi_vec[2]),
    rho * (1 - phi_vec[1] * phi_vec[2]), (1 - phi_vec[2]^2)
  ), 2, 2, byrow = TRUE)
  
  # coefficient matrix A = diag(phi1, phi2)
  A <- diag(phi_vec)
  
  # simulate
  require(MASS)  # for mvrnorm
  set.seed(42)
  total <- n + burnin
  eps  <- mvrnorm(total, mu = c(0,0), Sigma = Sigma)
  X    <- matrix(0, nrow = total, ncol = 2)
  
  for (t in 2:total) {
    X[t, ] <- A %*% X[t - 1, ] + eps[t, ]
  }
  
  # drop burn-in
  X <- X[(burnin + 1):total, , drop = FALSE]
  colnames(X) <- c("x1","x2")
  as.data.frame(X)
}

#' @export
generate_resampling_edges = function(ref_yr, max_age, beta, num_edges) {
  # Creates unevenly spaced edges with increasing spacing back in time as
  # controlled by an exponent beta, that should be >= 1. 
  #
  # Parameters
  #   ref_yr:    Reference year for age
  #   max_age:   Age used for the last edge
  #   beta:      The spacing is proportional to age^(1/beta)
  time_pts = seq(0, max_age^(1/beta), length.out = num_edges)^beta
  edges = ref_yr - round(time_pts)
  edges
}

#' @export
apply_resampling <- function(df, variables, edges) {
  # Resamples one or several time series at time points specified by edges.
  # The resampling is obtained by taking a triangle weighted average of the
  # original series, interpolating it linearly if necessary. The triangles
  # extend two edged-intervals except at the end points, where a single interval
  # is used. Works both for series going forward in time and back in time.
  #
  # Parameters
  #   df:         Dataframe with one or more time series. The column yr gives the
  #               time points.
  #   variables:  List of column names in df.
  #   edges:      Time points where the resampled values are computed, as
  #               returned by generate_resampling_edges.
  #
  # Value:
  #   A new data frame with yr = edges and columns specified by variables
  #   containing the resampled time series.
  
  newdf <- data.frame(yr = edges)
  n = length(edges)
  w = c(3, 4, 5, 6, 6, 6, 6)
  tri = c(head(w, 6), rev(w))
  tri1 = rev(w)
  trin = w
  tri = tri/sum(tri)
  tri1 = tri1/sum(tri1)
  trin = trin/sum(trin)
  for (v in variables) {
    newdf[[v]] <- NA
    for (i in 1:n) {
      if (i == 1) {
        e0 = edges[1]
        e1 = edges[2]
        nout = 5
        wgt = tri1
      }
      else if (i == n) {
        e0 = edges[n-1]
        e1 = edges[n]
        nout = 5
        wgt = trin
      }
      else {
        e0 = edges[i-1]
        e1 = edges[i+1]
        nout = 9
        wgt = tri
      }
      xout = seq(e0, e1, length.out = nout)
      yout = approx(x = df$yr, y = df[[v]], xout = xout)$y
      newdf[i, v] <- sum(wgt*yout)  # triangle weighted average
    }
  }
  newdf
}

#' @export
common_time <- function(df1, var1, df2, var2) {
  # df1 and df2 should have the same yr columns (from apply_resampling).
  # Returns the values of yr where both var1 and var2 are specified
  # (don't have NaN value), along with the corresponding var1 and var2
  stopifnot(identical(df1$yr, df2$yr))
  I <- !is.nan(df1[[var1]]) & !is.nan(df2[[var2]])
  x = df1[I, var1]
  y = df2[I, var2]
  list(yr = df1$yr[I], x = x, y = y)
}
