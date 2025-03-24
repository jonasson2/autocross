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

# PLOT MONTHLY AND ANNUAL MEANS
source('plot_d18_t.R')  # To do: búa til tvö plott.
plot_monthly(tasim, dye3.monthly)
plot_annual(dye3.annual, tasi, græn, sty)

# DETERMINE EDGES, I.E. END POINTS OF RESAMPLING INTERVALS
source('resampling.R')
end_yr = tail(dye3.annual$yr, 1)
beg_yr = 1780
edges = resample(beg_yr, end_yr, 2)

# FUNCTION TO COMPUTE CORRELATION BETWEEN A PAIR OF TIME SERIES
find_CI = function(edges, df1, var1, df2, var2) {
  resamp1 = apply_resampling(df1, edges, var1)
  resamp2 = apply_resampling(df2, edges, var2)
  common.res = common_time(resamp1, var1, resamp2, var2)
  time = 2000 - common.res$yr
  x = common.res$x
  y = common.res$y
  result = estimate_CI(time, x, y)
  str(result)
}

#s = find_CI(edges, sty, 't', dye3.annual, 'd18')
#s = find_CI(edges, sty, 't', dye3.annual, 'd18')
s = find_CI(edges, græn, 't', tasi, 'd18')
str(s)
