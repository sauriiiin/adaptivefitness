##### SDPG FITNESS RESULTS
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 10/05/2020

##### INITIALIZE
library(ggplot2)
library(ggridges)
library(gridExtra)
library(grid)
library(tidyverse)
library(ggpubr)
library(stringr)
library(scales)
library(egg)
library(zoo)
library(ggrepel)
library(reshape2)
library(RMariaDB)

source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")
out_path = 'figs/SDPG/';

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 7
txt <- 7
lbls <- 9

##### GATHER DATA
data1536 <- dbGetQuery(conn, 'select * from SDPG_GLU_FS_R1_1536_FITNESS')

for (hr in unique(data1536$hours)) {
  for (orf in unique(data1536$orf_name[data1536$hours == hr])) {
    m <- median(data1536$fitnessit[data1536$hours == hr & data1536$orf_name == orf])
    madev <- mad(data1536$fitness[data1536$hours == hr & data1536$orf_name == orf])
    ul <- m + 2*madev
    ll <- m -2*madev
    
    data1536$fitness[data1536$hours == hr & data1536$orf_name == orf & data1536$fitness > ul] <- NA
    data1536$fitness[data1536$hours == hr & data1536$orf_name == orf & data1536$fitness < ll] <- NA
  } 
}

ggplot(data1536[data1536$hours == 8,],
       aes(x = orf_name, y = fitness)) +
  geom_violin() +
  facet_wrap(~hours)


##### COMPARE 1536 and 6144 FITNESS
data <- dbGetQuery(conn, 'select a.orf_name, a.cs_mean fitness_6144, b.cs_mean fitness_1536
                   from
                   SDPG_GLU_FS_R1_6144_FITNESS_STATS a, SDPG_GLU_FS_R1_1536_FITNESS_STATS b
                   where a.orf_name = b.orf_name
                   and a.hours = 8 and b.hours = 24
                   and a.orf_name not in ("BOR")')

ggplot(data, aes(x = fitness_6144, y = fitness_1536)) +
  geom_abline() +
  geom_point() +
  geom_text_repel(aes(label = orf_name),
            nudge_y      = 0.05,
            direction    = "y",
            angle        = 0,
            vjust        = 0,
            segment.size = 0.05) +
  # stat_cor(method =  "pearson", label.x = 1.05, label.y = 0.95) +
  coord_cartesian(xlim = c(0.9, 1.1),
                  ylim = c(0.9, 1.1)) +
  labs(x = 'Mean Fitness @ 6144-density',
       y = 'Mean Fitness @ 1536-density')

ggsave(sprintf("%sFITNESS_CORR.jpg",out_path),
       height = two.c, width = two.c, units = 'mm',
       dpi = 300)
