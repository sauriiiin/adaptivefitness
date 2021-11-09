##### F28FUR  - PLOTS
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 03/19/2021

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
out_path = 'figs/f28fu/repeat/'

##### FIGURE SIZE
one.c <- 90 #single column 
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 7
txt <- 7
lbls <- 9


#####
d <- 6144
rr <- 'R1'
h <- 36

temp.boxplot <- foa.data %>%
  filter(density == d & replicate == rr & hours == h &
           !(orf_name %in% c('BOR','REF','YHR021W-A')) & !is.na(orf_name)) %>%
  ggplot(aes(y = fitness_clean, x = orf_name, fill = orf_name)) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(method = 'wilcox.test', ref.group = "BF_control",
                     label = "p.signif", label.y = 1.2,
                     size = 1.5, angle = 90, vjust = 0.5) +
  scale_y_continuous(breaks = seq(0,2,0.05),
                     minor_breaks = seq(0,2,0.01)) +
  coord_cartesian(ylim = c(0.8,1.2)) +
  scale_fill_discrete(guide = F) +
  labs(x = 'Strains', y = 'Fitness') +
  facet_wrap(.~density*hours*condition,
             ncol = 2) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))

ggsave(sprintf("%sBOXPLOT_COLONYSIZES_%d_%s_%d.jpg",out_path, d, rr, h), temp.boxplot,
       height = two.c, width = two.c, units = 'mm',
       dpi = 300)


e <- c(997, 998)
temp.sd.boxplot <- foa.sd.data %>%
  filter(exp_id %in% e & orf_name != 'YNR015W') %>%
  ggplot(aes(y = fitness, x = orf_name, fill = orf_name)) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(method = 'wilcox.test', ref.group = "BF_control",
                     label = "p.signif", label.y = 1.2,
                     size = 1.5, angle = 90, vjust = 0.5) +
  scale_y_continuous(breaks = seq(0,2,0.05),
                     minor_breaks = seq(0,2,0.01)) +
  coord_cartesian(ylim = c(0.8,1.2)) +
  scale_fill_discrete(guide = F) +
  labs(title = 'SanDiego SC+GAL+G418',
       x = 'Strains', y = 'Fitness') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
# ggsave(sprintf("%sBOXPLOT_COLONYSIZES_SD_6144.jpg",out_path), temp.sd.boxplot,
#        height = two.c/2, width = two.c/2, units = 'mm',
#        dpi = 300)

temp.gal <- foa.data %>%
  filter(density == 6144 & hours == 36 & condition == 'GAL' &
           !(orf_name %in% c('BOR','REF','YHR021W-A')) & !is.na(orf_name)) %>%
  ggplot(aes(y = fitness_clean, x = orf_name, fill = orf_name)) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(method = 'wilcox.test', ref.group = "BF_control",
                     label = "p.signif", label.y = 1.2,
                     size = 1.5, angle = 90, vjust = 0.5) +
  scale_y_continuous(breaks = seq(0,2,0.05),
                     minor_breaks = seq(0,2,0.01)) +
  coord_cartesian(ylim = c(0.8,1.2)) +
  scale_fill_discrete(guide = F) +
  labs(title = 'Pittsburgh SC+GAL+G418',
       x = 'Strains', y = 'Fitness') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        # axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
gal.gal <- ggarrange(temp.gal, temp.sd.boxplot, nrow = 2)
ggsave(sprintf("%sBOXPLOT_COLONYSIZES_6144_GAL.jpg",out_path), gal.gal,
       height = two.c, width = two.c, units = 'mm',
       dpi = 300)

gal.gal.fit <- merge(foa.es %>%
  filter(density == 6144 & hours == 36 & condition == 'GAL' &
           !(group2 %in% c('BOR','REF','YHR021W-A'))) %>%
  group_by(group2) %>%
  summarise(mean = mean(orf_fit)) %>%
  data.frame(),
  foa.sd.es[!(foa.sd.es$orf_name %in% c('YHR021W-A','BF_control','YNR015W')),] %>%
  group_by(orf_name) %>%
  summarise(mean = mean(mean)) %>%
  data.frame(), 
  by.x = "group2", by.y = "orf_name")

##### GAL GAL CORR
ggplot(gal.gal.fit, aes(x = mean.x, y = mean.y)) +
  geom_abline() +
  geom_point() +
  # geom_text_repel(aes(label = group2), max.overlaps = 30) +
  geom_smooth(method = 'lm') +
  stat_cor(method = 'spearman') +
  coord_cartesian(xlim = c(0.9,1.12),
                  ylim = c(0.9,1.12)) +
  labs(x = 'Pittsburgh', y = 'SanDiego') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        # axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        # axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%sSCATTER_PITT_SD.jpg",out_path),
       height = one.c, width = one.c, units = 'mm',
       dpi = 300)
