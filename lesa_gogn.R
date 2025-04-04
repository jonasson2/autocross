#setwd("C:/Users/kbo/Documents/autocross")
setwd("~/autocross")
library(ggplot2)
library(autocross)

# READ TEMPERATURE DATA
source('read_data.R')
t.res = get_temperature()
tasi = t.res$tasi
tasim = t.res$tasim
sty = t.res$sty
græn = t.res$græn

# READ DYE3 DATA
dye3.res = get_dye3()
dye3.monthly = dye3.res$monthly
dye3.annual = dye3.res$annual

# READ HAK-HVT TEMPERATURE PROXY
lakes = get_lakes()

# PLOT MONTHLY AND ANNUAL MEANS
source('plot_d18_t.R')  # To do: búa til tvö plott.
plot_monthly(tasim, dye3.monthly)
plot_annual(dye3.annual, tasi, græn, sty)

# DETERMINE EDGES, I.E. END POINTS OF RESAMPLING INTERVALS
debugSource('resampling.R')
# debug_all()
end_yr = tail(dye3.annual$yr, 1)
beg_yr = -3800
#obs.edges = resample(beg_yr, end_yr, 0.9, 100)

# FUNCTION TO COMPUTE CORRELATION BETWEEN A PAIR OF TIME SERIES
find_CI = function(edges, df1, var1, df2, var2) {
  resamp1 = apply_resampling(df1, edges, var1)
  resamp2 = apply_resampling(df2, edges, var2)
  common.res = common_time(resamp1, var1, resamp2, var2)
  time = rev(2000 - common.res$yr)
  x = common.res$x
  y = common.res$y
  result = list(CI=estimate_CI(time, x, y, alpha=0.01), xy=common.res)
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

#s = find_CI(edges, sty, 't', dye3.annual, 'd18')
#s = find_CI(edges, sty, 't', dye3.annual, 'd18')

compute.table <- function() {
  alpha.values = seq(0.60, 1.00, 0.05)
  n.values = c(10, 20, 30, 50, 100, 150, 200, 300, 400, 600, 800)
  M = length(n.values)
  N = length(alpha.values)
  R = matrix(0, M, N)
  ci.width = matrix(0, M, N)
  interval = rep(0, M)
  for (i in 1:M) {
    n = n.values[i]
    for (j in 1:N) {
      alpha = alpha.values[j]
      proxy.edges = resample(beg_yr, end_yr, alpha, n)
      CI.result = find_CI(proxy.edges, dye3.annual, 'd18', lakes, 'tproxy')
      sxy = CI.result$CI
      R[i,j] = sxy$r
      interval[i] = (end_yr - beg_yr)/n
      ci.width[i,j] = sxy$ci[2] - sxy$ci[1]
    }
  }
  signif.indicator = ci.width/2/R
  cat('\nCORRELATIONS:')
  pretty.print(interval, alpha.values, R)
  cat('\nCONFIDENCE INTERVAL WIDHTS')
  pretty.print(interval, alpha.values, ci.width)
  pretty.print(interval, alpha.values, signif.indicator)
}
# str(s)
# print(s$ci[2] - s$ci[1])
compute.table()

proxy.edges = resample(beg_yr, end_yr, 0.6, 10)
CIxy = find_CI(proxy.edges, dye3.annual, 'd18', lakes, 'tproxy')
CI = CIxy$CI
xy = CIxy$xy
str(CI)
