plot_monthly = function(tasim, dye3) {
  if (!is.null(dev.list())) {
    dev.off()
  }
  core = subset(dye3, year > 1890, select=c('d18', 'year'))
  core$d18scaled = (core$d18 - mean(core$d18))/sd(core$d18)*5
  yrmin = 1920
  yrmax = 1980
  plot(tasim$yrmonth, tasim$t, 'l', xlim=c(yrmin,yrmax), ylim=c(-15, 15), 
       xlab = 'Year', ylab = 'Temp. °C / Delta-O18 (scaled)')
  lines(core$year, core$d18scaled, 'l', col='red')
  
  # lines(tasi$yr, tasi$t, 'l', col='blue')
  abline(v = seq(yrmin, yrmax, by = 1), col = "gray", lty = 1)
}

plot_annual = function(annual, tasi, græn, sty) {
  dev.off()
  yrmin = 1790
  yrmax = 1980
  d18shifted = annual$d18 + 27
  plot(tasi$yr, tasi$t, 'l', xlim=c(yrmin,yrmax), ylim=c(-5, 5), lwd=2, 
       xlab = 'Year', ylab = 'Temp. °C')
  lines(annual$yr, d18shifted, 'l', col='red', lwd=2)
  lines(sty$yr, sty$t, 'l', col='dodgerblue', lwd=2)
  lines(græn$yr, græn$t, 'l', col='mediumseagreen', lwd=2)
  abline(v = seq(yrmin, yrmax, by = 10), col = "gray", lty = 1)
  abline(h = seq(-5,5), col = "gray", lty = 1)
}