##### MPP1
##### Multipin Pilot #1
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 10/25/2019

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
source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

# conn <- dbConnect(MariaDB(), dbname = '',
#                   usr = '',
#                   password = '',
#                   host = 'paris.csb.pitt.edu')

out_path = 'figs/multipin_pilot/';

##### GETTING DATA FROM SQL

data <- dbGetQuery(conn, 'select c.*, b.orf_name, a.hours, a.average
          from Branden.pilot_1_rep1_1_SC_Com_Glu_384_RAW a,
          Branden.pilot_1_rep1_1_SC_Com_Glu_pos2orf_name b, Branden.pilot_1_pos2coor384 c
          where a.pos = b.pos and a.pos = c.pos
          order by a.hours, c.384plate, c.384col, c.384row')

##### CORRECTING FOR ZOOM LEVEL
# Don't need to do this once precorrected data is uploaded to SQL
data$average[data$hours %in% c(46,48,56,59,61)] <- data$average[data$hours %in% c(46,48,56,59,61)] * 1.543034963

##### PLOTTING DATA
ggplot(data[data$hours > 6 & data$hours != 64,]) +
  geom_point(aes(x = hours, y = average, col = orf_name)) +
  geom_smooth(aes(x = hours, y = average, col = orf_name),
              method = 'loess') +
  scale_x_continuous(breaks = seq(0,80,4)) +
  theme_linedraw() 
ggsave(sprintf("%sgrowth.jpg",out_path),
       width = 7, height = 7,
       dpi = 300)

ggplot(data[data$hours == 48,]) +
  geom_line(aes(x = average, col = orf_name), stat = 'density', trim = T)


##### REMOVING OUTLIERS USING 2MAD
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

ggplot(data[data$hours == 48,]) +
  geom_line(aes(x = average, col = orf_name), stat = 'density', trim = T)

ggplot(data[data$hours > 6 & data$hours != 64,]) +
  geom_point(aes(x = hours, y = average, col = orf_name)) +
  geom_smooth(aes(x = hours, y = average, col = orf_name),
              method = 'loess') +
  scale_x_continuous(breaks = seq(0,80,4)) +
  theme_linedraw() 
ggsave(sprintf("%sgrowth_clean.jpg",out_path),
       width = 7, height = 7,
       dpi = 300)

             