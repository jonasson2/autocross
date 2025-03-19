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

# TEST VARIABLE TIME STEP CORRELATION
source('read_data.R')
#source('resampling.R')
library(autocross)
hiti = get_temperature()
sty = hiti$sty

dye3 = get_dye3()
end_yr = tail(dye3$annual$yr, 1)
beg_yr = 1780
n_points = 119
rsint = resample(beg_yr, end_yr, 2)
sty_resamp = apply_resampling(sty, rsint, "t", reference_yr, n_points)
dye3_resamp = apply_resampling(dye3$annual, rsint, "d18", reference_yr, n_points)

result = estimate_CI(time = dye3_resamp$yr,
                     x.series = dye3_resamp$d18,
                     y.series = sty_resamp$t)
