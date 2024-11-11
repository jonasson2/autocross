estimate_CI <- function(time, x.series, y.series, alpha = 0.05) {
  n = length(time)
  print(paste('n=', n))
  # Call the Fortran subroutine
  print('Calling Fortran')
  print(typeof(n))
  print(typeof(time))
  print(typeof(x.series))
  print(typeof(y.series))
  print(typeof(alpha))
  flush.console()
  #result = .Fortran('ss', x = x.series, s = double(1))
  print(paste('n=', n))
  result <- .Fortran("pearsont3sub", n = n, time = time, x = x.series,
                     y = y.series, alpha = alpha,
                     r = double(1), ci = double(2), taux = double(1),
                     tauy = double(1))
  str(result)
  result$s
  # list(
  #   correlation = result$r,
  #   confidence_interval = result$ci,
  #   persistence_time_x = result$taux,
  #   persistence_time_y = result$tauy
  # )
}

setwd('~/drive/autocross')
data = read.table('test_data.txt', col.names=c('time', 'x.series', 'y.series'))
if (is.loaded('pearsont3sub')) dyn.unload('pearsont3.so')
dyn.load('pearsont3.so')
print(data$time)
#source('R/estimate_CI.R')
dyn.load('pearsont3.so')
result = estimate_CI(data$time, data$x.series, data$y.series)
str(result)
if (is.loaded('prufa')) dyn.unload('prufa.so')
dyn.load('prufa.so')
n = length(data$x.series)
summa = double(1)
.Fortran('estimate_CI', n, data$x.series, data$y.series, summa)
print(paste('summa=', summa))

