##### F28FU - INIT
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 02/09/2021

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
out_path = 'figs/f28fu/'

##### FIGURE SIZE
one.c <- 90 #single column 
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 7
txt <- 7
lbls <- 9

##### LOAD DATA
expt.name <- 'F28FU_FOA'
map.data <- dbGetQuery(conn, sprintf('select a.*, c.*
                       from %s_pos2coor a, %s_pos2strainid b, %s_strainid2orf_name c
                       where a.pos = b.pos and b.strain_id = c.strain_id
                       order by a.density, a.plate, a.col, a.row', expt.name, expt.name, expt.name))

##### PLOT PLATEMAPS
map.plot <- ggplot(map.data[!(map.data$density == 6144 & map.data$plate > 3),]) +
  geom_tile(aes(x = col, y = row, fill = orf_name), col = 'black') +
  labs(x = 'Column', y = 'Row', fill = 'Control\nStrain') +
  # scale_fill_discrete(breaks = c('O1','O2','O3',
  #                              'C19','C20','C21',
  #                              'C4','C5','C6',
  #                              'REF'),
  #                     labels = c('O1' = 'Original 1','O2' = 'Original 2','O3' = 'Original 3',
  #                                'C19' = 'Copy 1','C20' = 'Copy 2','C21' = 'Copy 3',
  #                                'C4' = 'New 1','C5' = 'New 2','C6' = 'New 3',
  #                                'REF')) +
  # scale_fill_discrete(na.translate = F) +
  scale_y_reverse() +
  facet_wrap(.~density*plate, scales = 'free', ncol = 6) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%s%s_PLATEMAPS.jpg",out_path,expt.name), map.plot,
       height = 160, width = 400, units = 'mm',
       dpi = 300)





