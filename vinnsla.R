setwd("~/autocross")
library(ggplot2)
library(autocross)
library(pracma)

# All dataframes go backwards time and have a column yr

# READ TEMPERATUE DATA
debugSource('read_data.R')
t.res = get_temperature()
tasi = t.res$tasi
tasim = t.res$tasim
sty = t.res$sty
græn = t.res$græn

# READ ICE CORE DATA
dye3.res = get_dye3()
dye3.monthly = dye3.res$monthly
dye3.annual = dye3.res$annual
renland = get_renland()

# READ MARINE CORE DATA
sst = get_sst()

# READ HAK-HVT TEMPERATURE PROXY
lakes = get_lakes()
#lakes$tproxy = detrend(lakes$tproxy)

# PLOT MONTHLY AND ANNUAL MEANS
source('plot_d18_t.R')  # To do: búa til tvö plott.
plot_monthly(tasim, dye3.monthly)
plot_annual(dye3.annual, tasi, græn, sty)

# DETERMINE EDGES, I.E. END POINTS OF RESAMPLING INTERVALS
#debugSource('resampling.R')
# debug_all()
#obs.edges = resample(old_yr, new_yr, 0.9, 1

# FUNCTION TO COMPUTE CORRELATION BETWEEN A PAIR OF TIME SERIES
find_CI = function(df1, var1, df2, var2, alpha, n) {
  new_yr = min(df1$yr[1], df2$yr[1])
  #old_yr = max(tail(df1$yr, 1), tail(df2$yr, 1))
  ref_yr = 1979
  max_age = 5779
  edges = generate_resampling_edges(ref_yr, max_age, alpha, n)
  resamp1 = apply_resampling(df1, var1, edges)
  resamp2 = apply_resampling(df2, var2, edges)
  common.res = common_time(resamp1, var1, resamp2, var2)
  time = 2000 - common.res$yr
  x = drop(detrend(common.res$x))
  y = drop(detrend(common.res$y))
  #x = common.res$x
  #y = common.res$y
  CI = estimate_CI(time, x, y)
  res = data.frame(time=time, x=x, y=y)
  result = list(CI=CI, xy=res, max_age = max_age)
}

pretty.print <- function(interval, alpha, M) {
  F <- formatC(M, format = "f", digits = 3)
  F <- matrix(F, nrow = nrow(M), ncol = ncol(M))
  alpha_labels <- formatC(alpha, format = "f", digits = 3)
  interval_labels <- formatC(interval, format = "f", digits = 3)
  col_width <- max(nchar(c("interval", alpha_labels, interval_labels, F))) + 2
  pad <- function(x) sprintf(paste0("%-", col_width, "s"), x)
  cat(pad("interval"), paste(sapply(alpha_labels, pad), collapse = ""), "\n")
  total_cols <- 1 + length(alpha_labels)
  cat(rep("-", col_width * total_cols), sep = "", "\n")
  for (i in seq_along(interval_labels)) {
    cat(pad(interval_labels[i]), paste(sapply(F[i, ], pad), collapse = ""), "\n")
  }
}

#s = find_CI(sty, 't', dye3.annual, 'd18', alpha=1, n=10)
#s = find_CI(sty, 't', dye3.annual, 'd18', alpha=1, n=10)

compute.table <- function(series1, var1, series2, var2) {
  cat("Computing tables for correlations")
  alpha.values = seq(0.60, 1.00, 0.1)
  n.values = c(10, 20, 50, 100, 200, 400, 800)
  M = length(n.values)
  N = length(alpha.values)
  R = matrix(0, M, N)
  ci.width = matrix(0, M, N)
  interval = rep(0, M)
  for (i in 1:M) {
    n = n.values[i]
    for (j in 1:N) {
      alpha = alpha.values[j]
      CI.result = find_CI(series1, var1, series2, var2, alpha, n)
      sxy = CI.result$CI
      R[i,j] = sxy$r
      interval[i] = CI.result$max_age/n
      ci.width[i,j] = sxy$ci[2] - sxy$ci[1]
    }
  }
  signif.indicator = ci.width/2/R
  cat('\nCORRELATIONS:\n')
  pretty.print(interval, alpha.values, R)
  cat('\nCONFIDENCE INTERVAL WIDHTS\n')
  pretty.print(interval, alpha.values, ci.width)
  cat('\nSIGNIFICANCE INDICATORS\n')
  pretty.print(interval, alpha.values, signif.indicator)
}
compute.table(dye3.annual, 'd18', lakes, 'tproxy')
# compute.table(lakes, 'tproxy', sst, 'SST')
# CIxy = find_CI(dye3.annual, 'd18', sst, 'SST', 0.8, 100)
# CIxy = find_CI(renland, 'd18', dye3.annual, 'd18', 0.9, 10)
# CI = CIxy$CI
# xy = CIxy$xy
# str(CI)
