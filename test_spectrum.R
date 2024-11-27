

if (Sys.info()['sysname'] == 'Windows') {
  setwd("C:/Users/kbo/Documents/autocross")
  if (is.loaded('pearsont3sub')) dyn.unload('pearsont3.dll')
  dyn.load('pearsont3.dll')
} else {
  setwd('~/drive/autocross')
  if (is.loaded('pearsont3sub')) dyn.unload('pearsont3.so')
  dyn.load('pearsont3.so')
}

x = read.table('x.dat', col.names=c('time', 'series'))
y = read.table('y.dat', col.names=c('time', 'series'))


source('R/Specrum.R')
result = Spectrum(x,y,config)
str(result)
