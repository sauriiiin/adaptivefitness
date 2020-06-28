

##### INITIALIZE
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyverse)
library(ggpubr)
library(stringr)
library(egg)
library(zoo)
library(ggrepel)
library(plotly)
library(scales)
library(reshape2)
library(RMariaDB)

source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")
out_path = 'figs/PAYM/';

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 7
txt <- 5
lbls <- 9

ul <- 635
ll <- 288.0675
m <- 463

##### RND PLATE SURFACE
rnd.data <- dbGetQuery(conn, 'select * from 4C4_FS_RND_6144_DATA')
rnd.data$rnd_avg[rnd.data$rnd_avg == 0] <- NA
rnd.data$average[rnd.data$average == 0] <- NA

rnd.data$source[rnd.data$row%%2==1 & rnd.data$col%%2==1] = 'Top Left'
rnd.data$source[rnd.data$row%%2==0 & rnd.data$col%%2==1] = 'Top Right'
rnd.data$source[rnd.data$row%%2==1 & rnd.data$col%%2==0] = 'Bottom Left'
rnd.data$source[rnd.data$row%%2==0 & rnd.data$col%%2==0] = 'Bottom Right'

rnd.data$colony[rnd.data$orf_name == 'BF_control'] <- 'Reference'
rnd.data$colony[rnd.data$orf_name != 'BF_control'] <- 'Query'
rnd.data$colony[is.na(rnd.data$colony)] <- 'Gap'

rnd.data$source <- factor(rnd.data$source, levels = c('Top Left','Top Right','Bottom Left','Bottom Right'))

rnd.data$truth[rnd.data$hours < rnd.data$rnd_hrs] <- 'Beneficial'
rnd.data$truth[rnd.data$hours > rnd.data$rnd_hrs] <- 'Deleterious'
rnd.data$truth[rnd.data$hours == rnd.data$rnd_hrs] <- 'Reference'
rnd.data$truth[rnd.data$colony == 'Gap'] <- 'Gap'

dat.raw <- rnd.data[rnd.data$hours == 6.14 & rnd.data$plate == 3,]

plt.avg <- ggplot(dat.raw[dat.raw$plate == 3,],
                  aes(x = col, y = row)) +
  geom_tile(aes(fill = average)) +
  scale_size_continuous(range = c(0,2.5), guide = F) +
  scale_x_continuous(breaks = seq(0,96,4)) +
  scale_y_continuous(breaks = seq(0,64,4),trans = 'reverse') +
  scale_fill_gradient(name = 'Colony Size (pix.)',
                      low = "#FF9800", high = "black", na.value = "white",
                      limits = c(ll,ul), oob = squish) +
  labs(x = 'Column', y = 'Row') +
  coord_cartesian(xlim = c(1,96), ylim = c(64,1)) +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")),
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"))

plt.avg.ref <- ggplot(dat.raw[dat.raw$plate == 3 & dat.raw$colony == 'Reference',],
                  aes(x = col, y = row)) +
  geom_tile(aes(fill = average)) +
  scale_size_continuous(range = c(0,2.5), guide = F) +
  scale_x_continuous(breaks = seq(0,96,4)) +
  scale_y_continuous(breaks = seq(0,64,4),trans = 'reverse') +
  scale_fill_gradient(name = 'Colony Size (pix.)',
                      low = "#FF9800", high = "black", na.value = "white",
                      limits = c(ll,ul), oob = squish) +
  labs(x = 'Column', y = 'Row') +
  coord_cartesian(xlim = c(1,96), ylim = c(64,1)) +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")),
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"))

plt.rndavg <- ggplot(dat.raw[dat.raw$plate == 3,],
                      aes(x = col, y = row)) +
  geom_tile(aes(fill = rnd_avg)) +
  scale_size_continuous(range = c(0,2.5), guide = F) +
  scale_x_continuous(breaks = seq(0,96,4)) +
  scale_y_continuous(breaks = seq(0,64,4),trans = 'reverse') +
  scale_fill_gradient(name = 'Colony Size (pix.)',
                      low = "#FF9800", high = "black", na.value = "white",
                      limits = c(ll,ul), oob = squish) +
  labs(x = 'Column', y = 'Row') +
  coord_cartesian(xlim = c(1,96), ylim = c(64,1)) +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")),
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"))

plt.truth <- ggplot(dat.raw[dat.raw$plate == 3,],
                     aes(x = col, y = row)) +
  geom_tile(aes(fill = truth), col = 'black') +
  scale_size_continuous(range = c(0,2.5), guide = F) +
  scale_x_continuous(breaks = seq(0,96,4)) +
  scale_y_continuous(breaks = seq(0,64,4),trans = 'reverse') +
  scale_fill_manual(name = '',
                    breaks = c('Reference','Beneficial','Deleterious','Gap'),
                    values = c('Reference' = '#303F9F',
                               'Beneficial' = '#388E3C',
                               'Deleterious' = '#F44336',
                               'Gap' = 'white')) +
  labs(x = 'Column', y = 'Row') +
  coord_cartesian(xlim = c(1,96), ylim = c(64,1)) +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")),
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"))

plt.hours <- ggplot(dat.raw[dat.raw$plate == 3,],
                    aes(x = col, y = row)) +
  geom_tile(aes(fill = as.factor(rnd_hrs)), col = 'black') +
  scale_size_continuous(range = c(0,2.5), guide = F) +
  scale_x_continuous(breaks = seq(0,96,4)) +
  scale_y_continuous(breaks = seq(0,64,4),trans = 'reverse') +
  # scale_fill_manual(name = '',
  #                   breaks = c('Reference','Beneficial','Deleterious','Gap'),
  #                   values = c('Reference' = '#303F9F',
  #                              'Beneficial' = '#388E3C',
  #                              'Deleterious' = '#F44336',
  #                              'Gap' = 'white')) +
  labs(x = 'Column', y = 'Row') +
  coord_cartesian(xlim = c(1,96), ylim = c(64,1)) +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")),
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"))


plt.srf <- ggarrange(plt.avg, plt.avg.ref, plt.rndavg, plt.truth,
                                 ncol = 2, nrow = 2)

ggsave(sprintf("%sVP2_6_3.jpg",out_path), plt.srf,
       height = 155, width = two.c, units = 'mm',
       dpi = 300)


##### VP2 Hourly and PIX ES wise results
load('figs/paper/VP2PIXES.RData')

rnd.es$lid_result2[rnd.es$lid_result == 'Ben/Ben'] <- 'True Beneficial' 
rnd.es$lid_result2[rnd.es$lid_result == 'Del/Del'] <- 'True Deleterious' 
rnd.es$lid_result2[rnd.es$lid_result == 'Neu/Ben' | rnd.es$lid_result == 'Neu/Del'] <- 'False Negative' 
rnd.es$lid_result2[rnd.es$lid_result == 'Del/Ben' | rnd.es$lid_result == 'Ben/Del'] <- 'False Positve' 

rnd.es$mcat_result2[rnd.es$mcat_result == 'Ben/Ben'] <- 'True Beneficial' 
rnd.es$mcat_result2[rnd.es$mcat_result == 'Del/Del'] <- 'True Deleterious' 
rnd.es$mcat_result2[rnd.es$mcat_result == 'Neu/Ben' | rnd.es$mcat_result == 'Neu/Del'] <- 'False Negative' 
rnd.es$mcat_result2[rnd.es$mcat_result == 'Del/Ben' | rnd.es$mcat_result == 'Ben/Del'] <- 'False Positve' 


rnd.es$lid_result2 <- factor(rnd.es$lid_result2, levels = c('True Beneficial',
                                                           'True Deleterious',
                                                           'False Negative',
                                                           'False Positve'))

rnd.es$mcat_result2 <- factor(rnd.es$mcat_result2, levels = c('True Beneficial',
                                                            'True Deleterious',
                                                            'False Negative',
                                                            'False Positve'))

mcatbar <- ggplot(rnd.es) +
  geom_bar(aes(x = as.factor(hours), fill = mcat_result2)) +
  labs(title = 'MCAT',
       x = 'Reference Population Time Point (hours)',
       y = '') +
  scale_fill_manual(name = 'Effects',
                    breaks = c('True Beneficial',
                               'True Deleterious',
                               'False Negative',
                               'False Positve'),
                    values = c('True Beneficial' = '#388E3C',
                               'True Deleterious' = '#F44336',
                               'False Negative' = 'black',
                               'False Positve' = 'grey80'),
                    guide = F) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm"))

lidbar <- ggplot(rnd.es) +
  geom_bar(aes(x = as.factor(hours), fill = lid_result2)) +
  labs(title = 'LID',
       x = 'Reference Population Time Point (hours)',
       y = '') +
  scale_fill_manual(name = 'Effects',
                    breaks = c('True Beneficial',
                               'True Deleterious',
                               'False Negative',
                               'False Positve'),
                    values = c('True Beneficial' = '#388E3C',
                               'True Deleterious' = '#F44336',
                               'False Negative' = 'black',
                               'False Positve' = 'grey80'),
                    guide = F) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm"))

truthbar <- ggplot(rnd.es,
                   aes(x = as.factor(hours))) +
  geom_bar(aes(fill = truth)) +
  labs(title = 'TRUTH',
       x = 'Reference Population Time Point (hours)',
       y = 'Number of Mock Mutants') +
  scale_fill_manual(name = 'Result',
                    breaks = c('Beneficial',
                               'Deleterious'),
                    values = c('Deleterious'='#F44336',
                               'Beneficial'='#388E3C'),
                    labels = c('Deleterious'='True Deleterious',
                               'Beneficial'='True Beneficial'),
                    guide = F) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm"))


lidpixes <- ggplot(rnd.es[!(rnd.es$lid_pix_es > 0 & rnd.es$truth == 'Deleterious') &
                            !(rnd.es$lid_pix_es < 0 & rnd.es$truth == 'Beneficial'),],
                   aes(x = lid_pix_es * 100, fill = lid_result2)) +
  geom_histogram(binwidth = 5) +
  labs(title = 'LID',
       x = 'Effect Size',
       y = 'Number of Mock Mutants') +
  scale_x_continuous(breaks = seq(-1000,1000,100),
                     minor_breaks = seq(-1000,1000,10),
                     labels = paste(seq(-1000,1000,100),'%', sep = '')) +
  scale_fill_manual(name = 'Result',
                    breaks = c('True Beneficial',
                               'True Deleterious',
                               'False Negative',
                               'False Positve'),
                    values = c('True Beneficial' = '#388E3C',
                               'True Deleterious' = '#F44336',
                               'False Negative' = 'black',
                               'False Positve' = 'grey80'),
                    drop = F) +
  theme_linedraw() + 
  # facet_wrap(.~hours) +
  theme(plot.title = element_text(size = titles),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm"))

mcatpixes <- ggplot(rnd.es[!(rnd.es$lid_pix_es > 0 & rnd.es$truth == 'Deleterious') &
                             !(rnd.es$lid_pix_es < 0 & rnd.es$truth == 'Beneficial'),],
                    aes(x = lid_pix_es*100, fill = mcat_result2)) +
  geom_histogram(binwidth = 5) +
  labs(title = 'MCAT',
       x = 'Effect Size',
       y = '') +
  scale_x_continuous(breaks = seq(-1000,1000,100),
                     minor_breaks = seq(-1000,1000,10),
                     labels = paste(seq(-1000,1000,100),'%', sep = '')) +
  scale_fill_manual(name = 'Result',
                    breaks = c('True Beneficial',
                               'True Deleterious',
                               'False Negative',
                               'False Positve'),
                    values = c('True Beneficial' = '#388E3C',
                               'True Deleterious' = '#F44336',
                               'False Negative' = 'black',
                               'False Positve' = 'grey80'),
                    drop = F) +
  theme_linedraw() +
  # facet_wrap(.~hours) +
  theme(plot.title = element_text(size = titles),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm"))

hrlybars <- ggpubr::ggarrange(truthbar, lidbar, mcatbar,
                              nrow = 1)
pixes <- ggpubr::ggarrange(lidpixes, mcatpixes,
                           common.legend = T,
                           legend = 'bottom')
pixesbars <- ggarrange(hrlybars, pixes,
                       nrow = 2,
                       heights = c(1,1.4),
                       labels = c('a','b'),
                       label.args = list(gp=gpar(font = 2, fontsize = lbls),
                                         hjust=-1))
ggsave(sprintf("%sVP2_RESULTS.jpg",out_path), pixesbars,
       height = two.c, width = two.c, units = 'mm',
       dpi = 300)


######
rnd.es[!(rnd.es$lid_pix_es < 0 & rnd.es$truth == 'Deleterious') &
         !(rnd.es$lid_pix_es > 0 & rnd.es$truth == 'Beneficial'),]
