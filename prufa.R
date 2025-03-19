if (Sys.info()['sysname'] == 'Windows') {
  setwd("C:/Users/kbo/Documents/autocross")
  if (is.loaded('pearsont3sub')) dyn.unload('pearsont3.dll')
  dyn.load('pearsont3.dll')
} else {
  setwd('~/autocross')
  if (is.loaded('pearsont3sub')) dyn.unload('pearsont3.so')
  dyn.load('pearsont3.so')
}

data = read.table('test_data.txt', col.names=c('time', 'x.series', 'y.series'))

print(data$time)
#source('R/estimate_CI.R')
result = estimate_CI(data$time, data$x.series, data$y.series)
str(result)

