library(pracma)

resample = function(begyr, endyr, alpha, n) {
  d = endyr - begyr
  edges = endyr - seq(0^alpha, d^alpha, length.out=n+1)^(1/alpha)
  edges[n+1] = begyr
  edges = rev(edges)
  # OUTPUT THE RESAMPLING INTERVALS:
  print('Resampling intervals:')
  reportyr = c(-8000,-4000,-2000,0,1000,1500,1700,1800,1900,2000)
  reportyr <- sort(unique(c(reportyr, begyr, endyr)))
  reportyr <- reportyr[begyr <= reportyr & reportyr <= endyr]
  for (yr in reportyr) {
    i <- findInterval(yr, edges, rightmost.closed = TRUE)
    len = edges[i+1] - edges[i]
    cat(sprintf("%5d  %5.1f\n", yr, len))
  }
  edges
}

apply_resampling <- function(df, edges, variable) {
  # Resamples the specified variable to the intervals specified by edges
  # from the function "resample". Includes interpolated values at the
  # edges so that cases where no points fall in the intervals are also
  # handled. Returns a resampled dataframe.
  midpoints <- (head(edges, -1) + tail(edges, -1)) / 2
  newdf <- data.frame(yr = midpoints)
  newdf[[variable]] <- NA
  for (i in 1:(length(edges) - 1)) {
    left_val <- approx(x = df$yr, y = df[[variable]], xout = edges[i])$y
    right_val <- approx(x = df$yr, y = df[[variable]], xout = edges[i+1])$y
    inside <- (edges[i] < df$yr) & (df$yr < edges[i+1])
    vals <- c(left_val, df[[variable]][inside], right_val)
    newdf[i, variable] <- mean(vals, na.rm = TRUE)
  }
  newdf
}

common_time <- function(df1, var1, df2, var2) {
  # df1 and df2 should have the same yr columns (from apply_resampling).
  # Returns the values of yr where both var1 and var2 are specified
  # (don't have NaN value), along with the corresponding var1 and var2
  stopifnot(identical(df1$yr, df2$yr))
  I <- !is.nan(df1[[var1]]) & !is.nan(df2[[var2]])
  #x = detrend(df1[I,var1])
  #y = detrend(df2[I,var2])
  x = df1[I, var1]
  y = df2[I, var2]
  list(yr = df1$yr[I], x = x, y = y)
}
