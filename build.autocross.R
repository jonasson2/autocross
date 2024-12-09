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
