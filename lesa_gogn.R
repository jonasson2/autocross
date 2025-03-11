#setwd("C:/Users/kbo/Documents/autocross")
library(ggplot2)

source('read_data.R')
result = get_tasi()
tasi = result$tasi
tasim = result$tasim
sty = result$sty
græn = result$græn

result = get_dye3()
dye3 = result$dye3
annual = result$annual

core = subset(dye3, year > 1890, select=c('d18', 'year'))
core$d18s = (core$d18 - mean(core$d18))/sd(core$d18)

dev.off()
yrmin = 1920
yrmax = 1980
plot(tasim$yrmonth, tasim$t, 'l', xlim=c(yrmin,yrmax), ylim=c(-15, 8))
lines(core$year, core$d18s, 'l', col='red')
lines(tasi$yr, tasi$t, 'l', col='blue')
abline(v = seq(yrmin, yrmax, by = 1), col = "gray", lty = 1)

dev.off()
yrmin = 1798
yrmax = 1980
d18shifted = annual$d18 + 27
plot(tasi$yr, tasi$t, 'l', xlim=c(yrmin,yrmax), ylim=c(-5, 5), lwd=2)
lines(annual$yr, d18shifted, 'l', col='red', lwd=2)
lines(sty$yr, sty$t, 'l', col='dodgerblue', lwd=2)
lines(græn$yr, græn$t, 'l', col='mediumseagreen', lwd=2)
abline(v = seq(yrmin, yrmax, by = 5), col = "gray", lty = 1)
abline(h = seq(-5,5), col = "gray", lty = 1)
