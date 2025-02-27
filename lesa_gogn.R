#setwd("C:/Users/kbo/Documents/autocross")
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

names(dye3)[names(dye3) == "d18O"] = "d18"
types <- read.table("NGRIP2_layer_types.txt",
                   header = TRUE,        # The first row contains column names
                   sep = "\t",           # Tab-separated values
                   row.names = NULL)


dye3$time <- approx(x=tscale$DYE3, y=tscale$Age, xout=dye3$depth)$y
dye3$year <- as.integer(approx(x=tscale$DYE3, y=tscale$year, xout=dye3$depth)$y)

library(dplyr)

annual <- dye3 %>%
  filter(!is.na(year)) %>%  # Exclude rows with missing year values
  group_by(year) %>%
  summarise(d18 = mean(d18, na.rm = TRUE))

print(annual)

ig=merge(annual,sty, by='year')
print(ig)

ig$sty_std = (ig$t - mean(ig$t))/sd(ig$t)
ig$d18_std = (ig$d18 - mean(ig$d18))/sd(ig$d18)

ggplot(ig, aes(year)) +
  geom_line(aes(year,d18_std),color='darkred')+
  geom_line(aes(year,sty_std),color='blue')+
  scale_x_continuous(breaks=seq(1800,1990,4))+
  scale_y_continuous(name="d18", sec.axis=sec_axis(~., name="t"))

ggplot(tail(annual,80),aes(year,d18)) + geom_line(size=1)+
     scale_x_continuous(breaks=seq(1900,1990,4))

ggplot(sty,aes(AR,T)) + geom_line(size=1)

plot(types$Depth, cumsum(types$Type-1))
sum(types$Type==1)
sum(types$Type==2)
sum(types$Type==3)
str(types)
types$age <- as.integer(approx(x=tscale$NGRIP2, y=tscale$Age, xout=types$Depth)$y)
t
