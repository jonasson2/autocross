library(dplyr)
get_dye3 = function() {

  tscale=read.table('GICC05_time_scale.tab',header=TRUE,row.names=NULL, sep='\t')

  tscale <- read.table("GICC05_time_scale.tab",
                       header = TRUE,        # The first row contains column names
                       sep = "\t",           # Tab-separated values
                       row.names = NULL,
                       check.names = FALSE,
                       stringsAsFactors = FALSE,
                       comment.char = "")

  dye3 <- read.table("DYE-3maincore_water_isotopes.tab",
                     header = TRUE,        # The first row contains column names
                     sep = "\t",           # Tab-separated values
                     row.names = NULL)

  names(dye3)[names(dye3) == "d18O"] = "d18"
  # types <- read.table("NGRIP2_layer_types.txt",
  #                     header = TRUE,        # The first row contains column names
  #                     sep = "\t",           # Tab-separated values
  #                     row.names = NULL)

  dye3$time <- approx(x=tscale$DYE3, y=tscale$Age, xout=dye3$depth)$y - 1/12
  dye3$year <- approx(x=tscale$DYE3, y=tscale$year, xout=dye3$depth)$y + 1/12
  dye3 = subset(dye3, !is.na(time) & !is.na(d18))
  annual <- dye3 %>%
    group_by(yr = as.integer(year)) %>%
    summarise(d18 = mean(d18, na.rm = TRUE))
  list(monthly=dye3, annual=annual)
}

get_lakes = function() {
  lakes <- read.table("lakes.txt",
              header = TRUE,        # The first row contains column names
              sep = "\t",           # Tab-separated values
              row.names = NULL,
              col.names = c("age", "tproxy"))
  lakes$yr = 2000 - lakes$age
  lakes
}

get_temperature = function() {
  # Reads monthly temperature data from Tasiilaq
  # Returns a dataframe of monthly means with time in fractional years
  # (Jan = 1/24, Feb = 3/24 etc., corresponding to middle of month),
  # and a dataframe of annual means with time in whole years.

  #setwd('~/autocross')
  tasifile=read.table('tasiilaq_1894.txt',header=FALSE)[,3:15]
  yr = tasifile[,1]
  T = as.matrix(tasifile[,2:13])/10
  tasimonth = as.vector(t(T))
  yrmonth <- as.vector(sapply(yr, function(y) y + ((2 * (1:12) - 1) / 24)))
  tasidf = data.frame(yrmonth = yrmonth, t = tasimonth)
  tasidf = subset(tasidf, t > -999)
  tasi = tasidf %>%
    group_by(yr = as.integer(yrmonth)) %>%
    summarise(t = mean(t, na.rm = TRUE))
  sty=read.table('t_sty_1798.txt',header=TRUE)[c('AR','T')]
  græn=read.table('SW_Greenland1840.txt',header=FALSE)
  names(sty)=c('yr','t')
  names(græn)=c('yr','t')
  list(tasim=tasidf, tasi=tasi, sty=sty, græn=græn)
}
