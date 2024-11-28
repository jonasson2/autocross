if (Sys.info()['sysname'] == 'Windows') {
  setwd("C:/Users/kbo/Documents/autocross")
  if (is.loaded('rx_subroutine')) dyn.unload('libredfitx.dll')
  dyn.load('libredfitx.dll')
} else {
  setwd('~/autocross')
  if (is.loaded('rx_subroutine')) dyn.unload('libredfitx.so')
  dyn.load('libredfitx.so')
}

x = read.table('x.dat', col.names=c('time', 'series'))
y = read.table('y.dat', col.names=c('time', 'series'))
config = list(
  nsim = 1000,
  ofac = 4.0,
  hifac = 1.0,
  n50 = 7,
  alpha = 0.05,
  iwin = 1)

source('R/spectrum.R')
str(x)
str(y)
str(config)
result = spectrum(x,y,config)
str(result)
