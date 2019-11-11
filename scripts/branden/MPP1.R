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

##### GETTING DATA FROM SQL
dat.r1p1 <- dbGetQuery(conn, 'select c.*, b.orf_name, a.hours, a.average
                      from MPP1_R1P1_384_JPEG a,
                      MPP1_pos2orf_name b, MPP1_pos2coor c
                      where a.pos = b.pos and a.pos = c.pos
                      order by a.hours, c.plate, c.col, c.row')
dat.r1p1$rep <- 'Replicate 1'
dat.r1p1$stage <- '1.SC_COM+GLU' 

dat.r1p2 <- dbGetQuery(conn, 'select c.*, b.orf_name, a.hours, a.average
                      from MPP1_R1P2_384_JPEG a,
                      MPP1_pos2orf_name b, MPP1_pos2coor c
                      where a.pos = b.pos and a.pos = c.pos
                      order by a.hours, c.plate, c.col, c.row')
dat.r1p2$rep <- 'Replicate 1'
dat.r1p2$stage <- '2.SD+MET+GLU' 

dat.r1p3 <- dbGetQuery(conn, 'select c.*, b.orf_name, a.hours, a.average
                      from MPP1_R1P3_1536_JPEG a,
                      MPP1_pos2orf_name b, MPP1_pos2coor c
                      where a.pos = b.pos and a.pos = c.pos
                      order by a.hours, c.plate, c.col, c.row')
dat.r1p3$rep <- 'Replicate 1'
dat.r1p3$stage <- '3.SD+MET+GAL' 

dat.r1p4 <- dbGetQuery(conn, 'select c.*, b.orf_name, a.hours, a.average
                      from MPP1_R1P4_1536_JPEG a,
                      MPP1_pos2orf_name b, MPP1_pos2coor c
                      where a.pos = b.pos and a.pos = c.pos
                      order by a.hours, c.plate, c.col, c.row')
dat.r1p4$rep <- 'Replicate 1'
dat.r1p4$stage <- '4.SD-MET+GAL' 

dat.r2p1 <- dbGetQuery(conn, 'select c.*, b.orf_name, a.hours, a.average
                      from MPP1_R2P1_384_JPEG a,
                       MPP1_pos2orf_name b, MPP1_pos2coor c
                       where a.pos = b.pos and a.pos = c.pos
                       order by a.hours, c.plate, c.col, c.row')
dat.r2p1$rep <- 'Replicate 2'
dat.r2p1$stage <- '1.SC_COM+GLU' 

dat.r2p2 <- dbGetQuery(conn, 'select c.*, b.orf_name, a.hours, a.average
                       from MPP1_R2P2_384_JPEG a,
                       MPP1_pos2orf_name b, MPP1_pos2coor c
                       where a.pos = b.pos and a.pos = c.pos
                       order by a.hours, c.plate, c.col, c.row')
dat.r2p2$rep <- 'Replicate 2'
dat.r2p2$stage <- '2.SD+MET+GLU' 

dat.r2p3 <- dbGetQuery(conn, 'select c.*, b.orf_name, a.hours, a.average
                       from MPP1_R2P3_1536_JPEG a,
                       MPP1_pos2orf_name b, MPP1_pos2coor c
                       where a.pos = b.pos and a.pos = c.pos
                       order by a.hours, c.plate, c.col, c.row')
dat.r2p3$rep <- 'Replicate 2'
dat.r2p3$stage <- '3.SD+MET+GAL' 

dat.r2p4 <- dbGetQuery(conn, 'select c.*, b.orf_name, a.hours, a.average
                       from MPP1_R2P4_1536_JPEG a,
                       MPP1_pos2orf_name b, MPP1_pos2coor c
                       where a.pos = b.pos and a.pos = c.pos
                       order by a.hours, c.plate, c.col, c.row')
dat.r2p4$rep <- 'Replicate 2'
dat.r2p4$stage <- '4.SD-MET+GAL' 

data <- rbind(dat.r1p1, dat.r1p2, dat.r1p3, dat.r1p4,
              dat.r2p1, dat.r2p2, dat.r2p3, dat.r2p4)

# dbWriteTable(conn, "MPP1_ALL_JPEG", data, overwrite = T)

##### PLOTTING DATA
growth <- ggplot(data) +
  geom_point(aes(x = hours, y = average, col = orf_name)) +
  geom_smooth(aes(x = hours, y = average, col = orf_name),
              method = 'loess') +
  labs(x = 'Time (hrs.)',
       y = 'Colony Size (pix.)') +
  theme_linedraw() +
  scale_color_manual(name = 'Strain',
                         breaks = c('BY4741', 'BY4742'),
                         values = c('BY4741' = '#3F51B5', 'BY4742' = '#FFC107')) +
  facet_wrap(.~rep * stage, scales = 'free',
             ncol = 4)
ggsave(sprintf("%sgrowth.jpg",out_path),
       growth,
       width = 15, height = 8,
       dpi = 300)

##### REMOVING OUTLIERS USING 2MAD
data$cln_avg <- data$average
for (r in unique(data$rep)) {
  for (s in unique(data$stage[data$rep == r])) {
    for (hr in unique(data$hours[data$rep == r & data$stage == s])) {
      for (orf in unique(data$orf_name[data$rep == r & data$stage == s & data$hours == hr])) {
        m <- median(data$cln_avg[data$rep == r & data$stage == s & data$hours == hr & data$orf_name == orf], na.rm = T)
        madev <- mad(data$cln_avg[data$rep == r & data$stage == s & data$hours == hr & data$orf_name == orf], na.rm = T)
        ul <- m + 2*madev
        ll <- m - 2*madev
        
        data$cln_avg[data$rep == r & data$stage == s & data$hours == hr & data$orf_name == orf & data$cln_avg > ul] <- NA
        data$cln_avg[data$rep == r & data$stage == s & data$hours == hr & data$orf_name == orf & data$cln_avg < ll] <- NA
      } 
    }
  }
}

growth.cln <- ggplot(data) +
  geom_point(aes(x = hours, y = cln_avg, col = orf_name), alpha = 0.7) +
  geom_smooth(aes(x = hours, y = cln_avg, col = orf_name),
              method = 'loess') +
  labs(x = 'Time (hrs.)',
       y = 'Colony Size (pix.)') +
  theme_linedraw() +
  scale_color_manual(name = 'Strain',
                     breaks = c('BY4741', 'BY4742'),
                     values = c('BY4741' = '#536DFE', 'BY4742' = '#009688')) +
  facet_wrap(.~rep * stage, scales = 'free',
             ncol = 4)
ggsave(sprintf("%sgrowth_cln.jpg",out_path),
       growth.cln,
       width = 15, height = 8,
       dpi = 300)

ggplot(data) +
  geom_boxplot(aes(x = stage, y = cln_avg,
                   fill = orf_name, col = rep))


##### SATURATION POINT
sat.dat <- NULL
for (r in unique(data$rep)) {
  for (s in unique(data$stage[data$rep == r])) {
    sat.hr <- max(data$hours[data$rep == r & data$stage == s])
    sat.dat <- rbind(sat.dat,
                     data[data$rep == r & data$stage == s & data$hours == sat.hr,])
  }
}

ggplot(data) +
  geom_boxplot(aes(x = orf_name, y = cln_avg, fill = orf_name)) +
  stat_compare_means(aes(x = orf_name, y = cln_avg, group = orf_name),
                     geom = 'text', paired  = F) +
  facet_wrap(.~stage)



             