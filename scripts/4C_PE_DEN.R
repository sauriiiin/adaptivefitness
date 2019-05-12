##### PLATE EFFECT DENSITY PLOTS (4 Control)
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 05/12/2019

##### INITIALIZE
#install.packages("RMariaDB")
#install.packages("ggplot2")
library(RMariaDB)
library(ggplot2)
library(gridExtra)
source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

##### SET VARIABLES
expt_name = '4C3_GA1'
expt = 'FS1-1'
out_path = 'figs/illustrations/';
density = 6144;

tablename_fit = sprintf('%s_%d_FITNESS',expt_name,density);
tablename_p2o = '4C3_pos2orf_name1';
tablename_bpos = '4C3_borderpos';

p2c_info = NULL
p2c_info[1] = '4C3_pos2coor6144'
p2c_info[2] = '6144plate'
p2c_info[3] = '6144col'
p2c_info[4] = '6144row'

##### GET DATA

hours = dbGetQuery(conn, sprintf('select distinct hours from %s order by hours asc', tablename_fit))
p2c = dbGetQuery(conn, sprintf('select * from %s a order by a.%s, a.%s, a.%s',
                               p2c_info[1],
                               p2c_info[2],
                               p2c_info[3],
                               p2c_info[4]))
n_plates = dbGetQuery(conn, sprintf('select distinct %s from %s a order by %s asc',
                                    p2c_info[2],
                                    p2c_info[1],
                                    p2c_info[2]))

hr = hours[[1]][9]
pl = n_plates[[1]][1]

fitdat = dbGetQuery(conn, sprintf('select * from %s a, %s b where a.hours = %d and a.pos = b.pos and b.%s = %d order by b.%s, b.%s',
                                  tablename_fit,p2c_info[1],hr,p2c_info[2],
                                  pl,p2c_info[3],p2c_info[4]))
fitdat$bg[is.na(fitdat$average)] = NA
min = min(fitdat$average, na.rm=T)
max = max(fitdat$average, na.rm=T)

fitdat$source[fitdat$`6144row`%%2==1 & fitdat$`6144col`%%2==1] = 'TL'
fitdat$source[fitdat$`6144row`%%2==0 & fitdat$`6144col`%%2==1] = 'BL'
fitdat$source[fitdat$`6144row`%%2==1 & fitdat$`6144col`%%2==0] = 'TR'
fitdat$source[fitdat$`6144row`%%2==0 & fitdat$`6144col`%%2==0] = 'BR'


a <- ggplot(data = fitdat, aes(x=average, col = source)) +
  geom_density(lwd = 1.2) + 
  scale_colour_manual(name="Source",
                      values=c("TL"="#D32F2F","TR"="#536DFE","BL"="#388E3C","BR"="#795548"),
                      breaks=c("TL","TR","BL","BR"),
                      labels=c("Top Left","Top Right","Bottom Left","Bottom Right"),guide = F) +
  scale_x_continuous(breaks = seq(0,1000,100),
                     minor_breaks = seq(0,1000,20),
                     limits = c(min,max)) +
  scale_y_continuous(breaks = seq(0,0.02,0.005),
                     minor_breaks = seq(0,0.02,0.001),
                     limits = c(0,0.02)) +
  labs(x = 'Observed Pixel Count', y = 'Density') +
  theme_linedraw()

f <- ggplot(data = fitdat, aes(x=fitness, col = source)) +
  geom_density(lwd = 1.2) + 
  scale_colour_manual(name="Source",
                      values=c("TL"="#D32F2F","TR"="#536DFE","BL"="#388E3C","BR"="#795548"),
                      breaks=c("TL","TR","BL","BR"),
                      labels=c("Top Left","Top Right","Bottom Left","Bottom Right"),guide = F) +
  scale_x_continuous(breaks = seq(0,2,0.05),
                     minor_breaks = seq(0,2,0.025),
                     limits = c(0.8,1.2)) +
  scale_y_continuous(breaks = seq(0,15,1),
                     minor_breaks = seq(0,15,0.5),
                     limits = c(0,12)) +
  labs(x = 'Fitness', y = 'Density') +
  theme_linedraw()

#add NIL fit plot here, but firt modify fetch query to get everything at once, or just add another one, whichever simpler.

grid.arrange(a,f,nrow=1)





