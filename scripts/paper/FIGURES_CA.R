##### LID COVER ART
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 12/17/2020

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
library(wesanderson)
library(ggrepel)
library(reshape2)
library(qvalue)
library(RMariaDB)

source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")
out_path = 'figs/paper/coverart/';

##### GET PLATE MAPS
p2c_384 <- dbGetQuery(conn, 'select * from 4C4_pos2coor where density = 384')
p2o_all <- dbGetQuery(conn, 'select a.*, b.orf_name from 4C4_pos2coor a, 4C4_pos2orf_name b
                      where a.pos = b.pos')

p2o_all$pos_source <- as.numeric(str_trunc(as.character(p2o_all$pos), 4, side = 'left', ellipsis = ''))
p2o_all <- merge(p2o_all, p2c_384[c('pos','plate')], by.x = 'pos_source', by.y = 'pos', suffixes = c('','_source'))
p2o_all$title1 <- sprintf('%d-Density',p2o_all$density)
p2o_all$title1 <- factor(p2o_all$title1, sprintf('%d-Density',sort(unique(p2o_all$density))))
p2o_all$title2 <- sprintf('Plate Number %d',p2o_all$plate)
p2o_all$title2 <- factor(p2o_all$title2, sprintf('Plate Number %d',sort(unique(p2o_all$plate))))

p2o_all$plate_source[is.na(p2o_all$orf_name)] <- NA
##### PLOT
plot.pm <- ggplot(p2o_all,
       aes(x = col, y = row)) +
  geom_tile(aes(fill = as.factor(plate_source)), col = 'black') +
  scale_y_reverse() +
  scale_fill_manual(guide = F,
                    values = c('1' = "#E45129",
                               '2' = "#009688",
                               '3' = "#FFC107",
                               '4' = "#3F51B5")) +
  # scale_fill_brewer(guide = F, palette="Accent") +
  facet_wrap(.~title1*title2, scale = 'free') +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        # panel.background = element_blank(),
        rect = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        # strip.text = element_text(size = 7,
        #                           margin = margin(0.1,0,0.1,0, "mm")))
        strip.background = element_blank(),
        strip.text = element_blank())
ggsave(sprintf("%sPLATEMAPS.jpg",out_path), plot.pm,
       height = 6, width = 11, dpi = 300)


rgb(228,81,41, max=255)

