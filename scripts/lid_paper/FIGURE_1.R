##### FIGURE 1
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 09/21/2020

##### INITIALIZE
library(ggplot2)
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
out_path = 'figs/paper/';

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 7
txt <- 5
lbls <- 9

##### PLATE ILLUSTRATIONS
# 1.0, 1.4, 2.9, 4.0, 4.9, 6.1, 6.9, 7.8, 9.0, 10.0, 11.0 hours
# name = 'Phenotype',
# breaks = c('Beneficial','Neutral','Deleterious'),
# values = c('Deleterious'='#3F51B5',
#            'Neutral'='#212121',
#            'Beneficial'='#FFC107')

plates <- dbGetQuery(conn, 'select * from 4C4_FS_RND2_6144_DATA
                           where plate = 2 and hours = 11.04')
plates$colony[plates$orf_name == 'BF_control'] <- 'Reference'
plates$colony[plates$orf_name != 'BF_control'] <- 'Mutant'
plates$colony[is.na(plates$colony)] <- 'Gap'

ggplot(plates[plates$colony != 'Gap',]) +
  geom_point(aes(x = col, y = row,
                 shape = colony, col = colony)) +
  scale_y_continuous(trans = 'reverse') +
  scale_shape_manual(name = 'Colony Type',
                     breaks = c("Reference","Mutant"),
                     values = c("Reference" = 16,
                                "Mutant" = 17)) +
  scale_color_manual(name = 'Time',
                    breaks = c("Reference","Mutant"),
                    labels = c("9.0 hr","10.0 hr"),
                    values = c("Reference" = "#212121",
                               "Mutant" = "#FFC107")) +
  labs(x = 'Column',
       y = 'Row') +
  theme_linedraw() +
  coord_cartesian(xlim = c(10,34),
                  ylim = c(10,26))





