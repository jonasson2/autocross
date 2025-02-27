setwd("C:/Users/kbo/Documents/autocross")
library(ggplot2)

sty=read.table('t_sty_1798.txt',header=TRUE)
names(sty)=c('stod','year','t')

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

types <- read.table("NGRIP2_layer_types.txt",
                   header = TRUE,        # The first row contains column names
                   sep = "\t",           # Tab-separated values
                   row.names = NULL)


plot(d18O$depth,d18O$d18O,'l')

dye3$time <- approx(x=tscale$DYE3, y=tscale$Age, xout=dye3$depth)$y
dye3$year <- as.integer(approx(x=tscale$DYE3, y=tscale$year, xout=dye3$depth)$y)

library(dplyr)

annual <- dye3 %>%
  filter(!is.na(year)) %>%  # Exclude rows with missing year values
  group_by(year) %>%
  summarise(d18O = mean(d18O, na.rm = TRUE))

print(annual)

ig=merge(annual,sty, by='year')
print(ig)

ggplot(tail(ig,80),aes(year)) + geom_line(aes(year,d18O))+
  geom_line(aes(year,t),color='blue')+
  scale_x_continuous(breaks=seq(1900,1990,4))+
  scale_y_continuous(name="d18O", sec.axis=sec_axis(~., name="t"))



ggplot(tail(annual,80),aes(year,d18O)) + geom_line(size=1)+
     scale_x_continuous(breaks=seq(1900,1990,4))



ggplot(tail(head(d18O,400),200),aes(time,d18O)) + geom_line(size=1)+
  +   scale_x_continuous(breaks=seq(40,90,1))



ggplot(sty,aes(AR,T)) + geom_line(size=1)







plot(types$Depth, cumsum(types$Type-1))

sum(types$Type==1)
sum(types$Type==2)
sum(types$Type==3)


d18O$time <- approx(x=tscale$DYE3, y=tscale$Age, xout=d18O$depth)$y
