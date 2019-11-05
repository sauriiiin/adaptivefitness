##### MPP1_EZ
##### Multipin Pilot #1
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 11/05/2019

##### INSTALL PACKAGES
# comment this section out once you have installed them.
# you do not need to run this section everytime.
install.packages("RMariaDB")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("gridExtra")
install.packages("ggExtra")
install.packages("grid")
install.packages("tidyverse")
install.packages("ggpubr")
install.packages("stringr")

##### INITIALIZE
library(RMariaDB)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggExtra)
library(grid)
library(tidyverse)
library(ggpubr)
library(stringr)

sql.usr <- ''
sql.pwd <- ''
sql.db <- ''

conn <- dbConnect(MariaDB(), dbname = sql.db,
                  usr = sql.usr,
                  password = sql.pwd,
                  host = 'paris.csb.pitt.edu')

##### GETTING DATA FROM SQL
tblname_JPEG <- '_JPEG'
tblname_p2o <- '_pos2orf_name'
tblname_p2c <- '_pos2coor'

data <- dbGetQuery(conn, sprintf('select c.*, b.orf_name, a.hours, a.average
          from %s a, %s b, %s c
          where a.pos = b.pos and a.pos = c.pos
          order by a.hours, c.384plate, c.384col, c.384row',
                                 tblname_JPEG, tblname_p2o, tblname_p2c))

##### PLOTTING DATA
ggplot(data) +
  geom_point(aes(x = hours, y = average, col = orf_name)) +
  geom_smooth(aes(x = hours, y = average, col = orf_name),
              method = 'loess') +
  scale_x_continuous(breaks = seq(0,80,4)) +
  theme_linedraw() 
ggsave("growth.jpg",
       width = 7, height = 7,
       dpi = 300)

##### REMOVING OUTLIERS USING MEDIAN ADJUSTED DEVIATIONS
for (hr in unique(data$hours)) {
  for (orf in unique(data$orf_name[data$hours == hr])) {
    m <- median(data$average[data$hours == hr & data$orf_name == orf])
    madev <- mad(data$average[data$hours == hr & data$orf_name == orf])
    ul <- m + 2*madev
    ll <- m -2*madev
    
    data$average[data$hours == hr & data$orf_name == orf & data$average > ul] <- NA
    data$average[data$hours == hr & data$orf_name == orf & data$average < ll] <- NA
  } 
}

ggplot(data) +
  geom_point(aes(x = hours, y = average, col = orf_name)) +
  geom_smooth(aes(x = hours, y = average, col = orf_name),
              method = 'loess') +
  scale_x_continuous(breaks = seq(0,80,4)) +
  theme_linedraw() 
ggsave("growth_clean.jpg",
       width = 7, height = 7,
       dpi = 300)

             