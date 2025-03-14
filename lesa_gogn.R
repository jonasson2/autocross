#setwd("C:/Users/kbo/Documents/autocross")
setwd("~/autocross")
library(ggplot2)

# READ TEMPERATURE DATA
source('read_data.R')
result = get_temperature()
tasi = result$tasi
tasim = result$tasim
sty = result$sty
græn = result$græn

# READ DYE3 DATA
result = get_dye3()
dye3 = result$dye3
annual = result$annual

# PLOT MONTHLY AND ANNUAL MEANS
source('plot_d18_t.R')  # To do: búa til tvö plott.
plot_monthly(tasim, dye3)
plot_annual(annual, tasi, græn, sty)


