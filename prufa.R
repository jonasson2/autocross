data = read.table('test_data.txt', col.names=c('time', 'x.series', 'y.series'))
source('R/estimate_CI.R')
dyn.load('pearsont3sub.so')
result = estimate_CI(data$time, data$x.series, data$y.series)
