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

out_path = 'figs/multipin_pilot/';

##### ZOOM2PIXEL

z2p = readxl::read_xlsx('figs/multipin_pilot/z2p.xlsx')
z2p$ratio <- z2p$width/z2p$height
z2p$height <- z2p$width/median(z2p$ratio)
z2p$pixcount <- z2p$width * z2p$height
# z2p$pixcount <- z2p$pixcount - z2p$pixcount[39]
z2p$pixcount <- z2p$pixcount[z2p$focalLength == 55]/z2p$pixcount

pix_cor = NULL
i = 1
for (fl in unique(z2p$focalLength)) {
  pix_cor$focalLength[i] <- fl
  pix_cor$correction[i] <- mean(z2p$pixcount[z2p$focalLength == fl])
  i <- i + 1
}
pix_cor <- data.frame(pix_cor)
write.table(pix_cor, file = sprintf('%spix_cor.csv',out_path),
            sep = ",", row.names = F)

z2p.model <- lm(formula = pixcount ~ focalLength + I(focalLength^2) + I(focalLength^3), data = z2p)
z2p.model
ggplot(z2p) +
  geom_point(aes(x = focalLength, y = pixcount)) +
  geom_point(data = data.frame(seq(18,55,1)),
             aes(x = seq(18,55,1), y = predict(z2p.model, data.frame(focalLength = seq(18,55,1)))))

##### PLOTTING DATA

data <- dbGetQuery(conn, 'select c.*, b.orf_name, a.hours, a.average
          from Branden.pilot_1_rep1_1_SC_Com_Glu_384_RAW a,
          Branden.pilot_1_rep1_1_SC_Com_Glu_pos2orf_name b, Branden.pilot_1_pos2coor384 c
          where a.pos = b.pos and a.pos = c.pos
          order by a.hours, c.384plate, c.384col, c.384row')

data$average[data$hours %in% c(46,48,56,59,61)] <- data$average[data$hours %in% c(46,48,56,59,61)] * 1.543034963
# data$hours[data$hours %in% c(31,32,34,36,38,41)] <- data$hours[data$hours %in% c(31,32,34,36,38,41)] + 12

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

             