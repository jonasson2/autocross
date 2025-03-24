resample = function(beginyear, endyear, alpha) {
  # Returns edges of sampling intervals, where the length of each interval is
  # proportional to estimated uncertainty in the ice core timescale. If alpha =
  # 2 the interval lengths are max(1, 2*age/100) corresponding 2 standard
  # deviations and an uncertainty increase of 1 year per century. The first
  # interval may be cut short as the first edge is beginyear. The edges go
  # back in time.
  yr = endyear
  edges = endyear
  age = 0
  while (yr > beginyear) {
    len = max(1, round(age/100*alpha))
    yr = max(beginyear, yr - len)
    edges = c(edges, yr)
    age = age + len
  }
  edges
}

apply_resampling = function(df, edges, variable) {
  # Resamples the specified variable to the intervals specified by edges
  # from the function "resample". Returns a resampled dataframe.
  midpoints <- (head(edges, -1) + tail(edges, -1))/2
  newdf = data.frame(yr = midpoints)
  newdf[,variable] = NA
  for (i in 1:length(midpoints)) {  # nrow(rsint)
    len = edges[i] - edges[i+1]
    I = edges[i+1] <= df$yr & df$yr <= edges[i+1]
    newdf[i, variable] = mean(df[[variable]][I])
  }
  newdf
}

common_time <- function(df1, var1, df2, var2) {
  # df1 and df2 should have the same yr columns (from apply_resampling).
  # Returns the values of yr where both var1 and var2 are specified
  # (don't have NaN value), along with the corresponding var1 and var2
  stopifnot(identical(df1$yr, df2$yr))
  I <- !is.nan(df1[[var1]]) & !is.nan(df2[[var2]])
  list(yr = df1$yr[I], x = df1[I,var1], y = df2[I,var2])
}
