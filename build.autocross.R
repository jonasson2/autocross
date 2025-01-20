#setwd("~/autocross/autocross")
setwd("C:/Users/kbo/Documents/autocross/autocross")
library(devtools)
devtools::build()
devtools::document()
install.packages("~/autocross/autocross_0.0.0.9000.tar.gz", repos=NULL, type="source")
library(autocross)
data = read.table('../test_data.txt', col.names=c('time', 'x.series', 'y.series'))
result = estimate_CI(data$time, data$x.series, data$y.series)
str(result)


x = read.table('../x.dat', col.names=c('time', 'series'))
y = read.table('../y.dat', col.names=c('time', 'series'))
config = list(
  nsim = as.integer(1000),
  ofac = 4.0,
  hifac = 1.0,
  n50 = as.integer(7),
  alpha = 0.05,
  iwin = as.integer(1))
result = spectrum(x,y,config)
str(result)

