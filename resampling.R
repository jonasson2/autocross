resample = function(beginyear, endyear, alpha) {
  # Returns midpoints and lengths of sampling intervals, where the length of
  # each interval is proportional to estimated uncertainty in the ice core
  # timescale. If alpha = 2 the interval lengths are max(1, 2*age/100)
  # corresponding 2 standard deviations and an uncertainty increase of 1 year
  # per century.
  age = 0
  yr = endyear
  lengths = numeric(0)
  midyr = numeric(0)
  while (yr > beginyear) {
    len = max(1, round(age/100*alpha))
    lengths = c(lengths, len)
    mp = age - 0.5 + len/2
    midyr = c(midyr, endyear - mp)
    age = age + len
    yr = yr - len
  }
  data.frame(length=lengths, midyr=midyr)
}

apply_resampling = function(df, rsint, variable, reference_yr, n_points) {
  newdf = data.frame(yr = rsint$midyr[1:n_points])
  newdf[,variable] = NA
  for (i in 1:n_points) {  # nrow(rsint)
    row = rsint[i,]
    len = row$length
    I = row$midyr - len/2 <= df$yr & df$yr <= row$midyr + len/2
    newdf[i, variable] = mean(df[[variable]][I])
  }
  newdf
}
