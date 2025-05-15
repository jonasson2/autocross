generate_resampling_edges = function(ref_yr, max_age, alpha, num_edges) {
  # Creates unevenly spaced edges with increasing spacing back in time as
  # controlled by an exponent alpha
  #
  # Parameters
  #   ref_yr:    Reference year for age
  #   max_age:   Age used for the last edge
  #   alpha:     The spacing is proportional to age^(1/alpha)
  time_pts = max_age - seq(max_age^(1/alpha), 0, length.out = num_edges)^alpha
  edges = ref_yr - time_pts
  edges
}

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
  tri = c(seq(1,4), 5, seq(4, 1))
  tri = tri/sum(tri)
  tri1 = seq(5, 1)/sum(seq(1, 5))
  trin = seq(1, 5)/sum(seq(5, 1))
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
