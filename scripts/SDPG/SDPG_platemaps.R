##### SDPG PLATEMAPS
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
data <- dbGetQuery(conn, 'select a.*, b.orf_name from
                   SDPG_pos2coor a,
                   SDPG_pos2orf_name b
                   where a.pos = b.pos
                   order by a.density, a.plate, a.col, a.row')

##### PLATEMAP
platemaps <- ggplot(data, aes(x = col, y = row, fill = orf_name)) +
  geom_tile(col = 'black') +
  scale_y_reverse() +
  scale_fill_discrete(na.value="transparent") +
  facet_wrap(~density*plate, scales = 'free',
             ncol = 2) +
  theme_linedraw() +
  labs(title = 'SDPG Platemaps',
       x = 'Columns',
       y = 'Rows') +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")),
        legend.key.size = unit(3, "mm")) +
  guides(fill = guide_legend(nrow = 5))

ggsave(sprintf("%sPLATEMAP.jpg",out_path), platemaps,
       height = two.c, width = two.c, units = 'mm',
       dpi = 300)




