##### 4C4 Analysis
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 12/20/2019
##### Analyzing the repeat LID validation

##### INITIALIZE
library(RMariaDB)
library(ggplot2)
library(ggExtra)
library(gridExtra)
library(ggpubr)
source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

out_path <- 'figs/4C4/'

expt <- '4C4'
pattern <- 'UP1'

##### PS2 DATA ANALYSIS
stage <- 'PS2'
data <- dbGetQuery(conn, 'select * from 4C4_PS2_1536_FITNESS a, 4C4_pos2coor b
                   where a.pos = b.pos
                   order by hours, b.plate, b.col, b.row')

# TECHNICAL SOURCE
data$tech_src[data$row%%2==1 & data$col%%2==1] = '1TL'
data$tech_src[data$row%%2==0 & data$col%%2==1] = '3BL'
data$tech_src[data$row%%2==1 & data$col%%2==0] = '2TR'
data$tech_src[data$row%%2==0 & data$col%%2==0] = '4BR'
# BIOLOGICAL SOURCE
data$bio_src[data$plate == 1 & data$tech_src == '1TL'] <- 'one'
data$bio_src[data$plate == 1 & data$tech_src == '2TR'] <- 'two'
data$bio_src[data$plate == 1 & data$tech_src == '3BL'] <- 'three'
data$bio_src[data$plate == 1 & data$tech_src == '4BR'] <- 'four'
data$bio_src[data$plate == 2 & data$tech_src == '1TL'] <- 'four'
data$bio_src[data$plate == 2 & data$tech_src == '2TR'] <- 'one'
data$bio_src[data$plate == 2 & data$tech_src == '3BL'] <- 'two'
data$bio_src[data$plate == 2 & data$tech_src == '4BR'] <- 'three'
data$bio_src[data$plate == 3 & data$tech_src == '1TL'] <- 'three'
data$bio_src[data$plate == 3 & data$tech_src == '2TR'] <- 'four'
data$bio_src[data$plate == 3 & data$tech_src == '3BL'] <- 'one'
data$bio_src[data$plate == 3 & data$tech_src == '4BR'] <- 'two'
data$bio_src[data$plate == 4 & data$tech_src == '1TL'] <- 'two'
data$bio_src[data$plate == 4 & data$tech_src == '2TR'] <- 'three'
data$bio_src[data$plate == 4 & data$tech_src == '3BL'] <- 'four'
data$bio_src[data$plate == 4 & data$tech_src == '4BR'] <- 'one'
# COLONY TYPE
data$col_type[data$orf_name == 'BF_control'] <- 'control'
data$col_type[is.na(data$col_type)] <- 'query'

# PLOTTING PS2 DATA
ggplot(data[data$hours == 21,]) +
  geom_violin(aes(x = bio_src, y = average), fill = 'grey') +
  geom_boxplot(aes(x = bio_src, y = average),
               width = 0.35) +
  labs(title = sprintf('%s_%s - %s', expt, pattern, stage),
       subtitle = 'Comparing Colony Size of Biological Sources',
       x = 'Biological Source',
       y = 'Colony Size (pix.)') +
  scale_x_discrete(breaks = c('one','two','three','four'),
                   limits = c('one','two','three','four'),
                   labels = c('X','Y','Z','C')) +
  theme_linedraw()
ggsave(sprintf("%s%s%s_%s_biosrc_pix.jpg",out_path,expt,pattern,stage),
       width = 5, height = 5,
       dpi = 300)

ggplot(data[data$hours == 21,]) +
  geom_violin(aes(x = bio_src, y = fitness), fill = 'grey') +
  geom_boxplot(aes(x = bio_src, y = fitness),
               width = 0.35) +
  labs(title = sprintf('%s_%s - %s', expt, pattern, stage),
       subtitle = 'Comparing Fitness of Biological Sources',
       x = 'Biological Source',
       y = 'Fitness') +
  scale_x_discrete(breaks = c('one','two','three','four'),
                   limits = c('one','two','three','four'),
                   labels = c('X','Y','Z','C')) +
  theme_linedraw()
ggsave(sprintf("%s%s%s_%s_biosrc_fit.jpg",out_path,expt,pattern,stage),
       width = 5, height = 5,
       dpi = 300)


ggplot(data[data$hours == 21,]) +
  geom_violin(aes(x = tech_src, y = average), fill = 'grey') +
  geom_boxplot(aes(x = tech_src, y = average),
               width = 0.35) +
  labs(title = sprintf('%s_%s - %s', expt, pattern, stage),
       subtitle = 'Comparing Colony Size of Technical Sources',
       x = 'Technical Source',
       y = 'Colony Size (pix.)') +
  scale_x_discrete(breaks = c('1TL','2TR','3BL','4BR'),
                   limits = c('1TL','2TR','3BL','4BR'),
                   labels = c('Top Left','Top Right','Bottom Left','Bottom Right')) +
  theme_linedraw()
ggsave(sprintf("%s%s%s_%s_techsrc_pix.jpg",out_path,expt,pattern,stage),
       width = 5, height = 5,
       dpi = 300)

ggplot(data[data$hours == 21,]) +
  geom_violin(aes(x = tech_src, y = fitness), fill = 'grey') +
  geom_boxplot(aes(x = tech_src, y = fitness),
               width = 0.35) +
  labs(title = sprintf('%s_%s - %s', expt, pattern, stage),
       subtitle = 'Comparing Fitness of Technical Sources',
       x = 'Technical Source',
       y = 'Fitness') +
  scale_x_discrete(breaks = c('1TL','2TR','3BL','4BR'),
                   limits = c('1TL','2TR','3BL','4BR'),
                   labels = c('Top Left','Top Right','Bottom Left','Bottom Right')) +
  theme_linedraw()
ggsave(sprintf("%s%s%s_%s_techsrc_fit.jpg",out_path,expt,pattern,stage),
       width = 5, height = 5,
       dpi = 300)

