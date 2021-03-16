##### F28FU - CONT EXP
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 02/19/2021

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
library(effsize)
library(qvalue)
library(RMariaDB)

source("R/functions/isoutlier.R")
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

##### CONTROL STRAIN EXPERIMENT
##### GATHER DATA
expt.name <- 'CSE'

tblname_dat <- sprintf('F28FU_%s_FS_6144_FITNESS', expt.name)
tblname_p2o <- sprintf('F28FU_%s_pos2orf_name', expt.name)
tblname_p2c <- sprintf('F28FU_%s_pos2coor', expt.name)

cs.data <- dbGetQuery(conn, sprintf('select c.*, b.orf_name, a.hours, a.average, a.fitness
                                    from %s a, %s b, %s c
                                    where a.pos = b.pos and a.pos = c.pos
                                    order by a.hours, c.plate, c.col, c.row',
                                    tblname_dat, tblname_p2o, tblname_p2c))

cs.data$condition[cs.data$plate %in% c(1,2,3)] <- 'GLU'
cs.data$condition[cs.data$plate %in% c(4,5,6)] <- 'CAS'
cs.data$condition[cs.data$plate %in% c(7,8,9)] <- 'SDA'
cs.data$condition[cs.data$plate %in% c(10,11,12)] <- 'GAL'
cs.data$condition[cs.data$plate %in% c(13,14,15)] <- 'GLY'

cs.data$source[cs.data$row%%2==1 & cs.data$col%%2==1] = 'TL'
cs.data$source[cs.data$row%%2==0 & cs.data$col%%2==1] = 'BL'
cs.data$source[cs.data$row%%2==1 & cs.data$col%%2==0] = 'TR'
cs.data$source[cs.data$row%%2==0 & cs.data$col%%2==0] = 'BR'
cs.data$source <- factor(unique(cs.data$source), levels = c('TL','TR','BL','BR'))

cs.data$rep <- as.numeric(str_trunc(as.character(cs.data$pos), 4, side = 'left', ellipsis = ''))

cs.data$orf_name[cs.data$orf_name == 'O1'] <- 'Original 1'
cs.data$orf_name[cs.data$orf_name == 'O2'] <- 'Original 2'
cs.data$orf_name[cs.data$orf_name == 'O3'] <- 'Original 3'
cs.data$orf_name[cs.data$orf_name == 'C4'] <- 'New 1'
cs.data$orf_name[cs.data$orf_name == 'C5'] <- 'New 2'
cs.data$orf_name[cs.data$orf_name == 'C6'] <- 'New 3'
cs.data$orf_name[cs.data$orf_name == 'C19'] <- 'Copy 1'
cs.data$orf_name[cs.data$orf_name == 'C20'] <- 'Copy 2'
cs.data$orf_name[cs.data$orf_name == 'C21'] <- 'Copy 3'

cs.data$orf_name <- factor(unique(cs.data$orf_name),
                           levels = c('Original 1','Original 2','Original 3',
                                      'Copy 1','Copy 2','Copy 3',
                                      'New 1','New 2','New 3','REF'))

##### CLEAN DATA
head(cs.data)

for (h in unique(cs.data$hours)) {
  for (p in unique(cs.data$plate[cs.data$hours == h])) {
    temp <- cs.data[cs.data$hours == h & cs.data$plate == p,]
    temp$average[isoutlier(temp$average)] <- NA
    cs.data$average_clean[cs.data$hours == h & cs.data$plate == p] <- temp$average
    
    temp$fitness[isoutlier(temp$fitness)] <- NA
    cs.data$fitness_clean[cs.data$hours == h & cs.data$plate == p] <- temp$fitness
    
    for (o in unique(cs.data$orf_name[cs.data$hours == h & cs.data$plate == p])) {
      temp <- cs.data[cs.data$hours == h & cs.data$plate == p & cs.data$orf_name == o,]
      temp$average_clean[isoutlier(temp$average_clean)] <- NA
      cs.data$average_clean[cs.data$hours == h & cs.data$plate == p & cs.data$orf_name == o] <- temp$average_clean
      
      temp$fitness_clean[isoutlier(temp$fitness_clean)] <- NA
      cs.data$fitness_clean[cs.data$hours == h & cs.data$plate == p & cs.data$orf_name == o] <- temp$fitness_clean
      
      for (r in unique(cs.data$rep[cs.data$hours == h & cs.data$plate == p & cs.data$orf_name == o])) {
        temp <- cs.data[cs.data$hours == h & cs.data$plate == p & cs.data$orf_name == o & cs.data$rep == r,]
        temp$average_clean[isoutlier(temp$average_clean)] <- NA
        cs.data$average_clean[cs.data$hours == h & cs.data$plate == p & cs.data$orf_name == o & cs.data$rep == r] <- temp$average_clean
        
        temp$fitness_clean[isoutlier(temp$fitness_clean)] <- NA
        cs.data$fitness_clean[cs.data$hours == h & cs.data$plate == p & cs.data$orf_name == o & cs.data$rep == r] <- temp$fitness_clean
      }
    }
  }
}


##### PLOT COLONY SIZE
plot.den.cs <- ggplot(cs.data[cs.data$condition != 'GLY' & cs.data$orf_name != 'REF' & cs.data$hours > 0,],
       aes(x = average_clean, y = orf_name, fill = orf_name)) +
  geom_density_ridges(quantile_lines = TRUE,
                      quantiles = c(0.25, 0.5, 0.75),
                      scale = 3, alpha = 0.8, size = 0.3,
                      vline_size = 0.2, vline_color = "black",
                      na.rm = T) +
  # geom_vline(xintercept = c(0.99,1.01), linetype = 'dashed', col = 'red', lwd = 0.5) +
  # scale_x_continuous(breaks = seq(0,2,0.1),
  #                    minor_breaks = seq(0,2,0.01)) +
  labs(y = 'Strains', x = 'Colony Size (pix.)') +
  scale_fill_discrete(guide = F) +
  facet_wrap(.~hours*condition,
             ncol = 4) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.y = element_text(angle = -45, vjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))

plot.den.fit <- ggplot(cs.data[cs.data$condition != 'GLY' & cs.data$orf_name != 'REF' & cs.data$hours > 0,],
                      aes(x = fitness_clean, y = orf_name, fill = orf_name)) +
  geom_density_ridges(quantile_lines = TRUE,
                      quantiles = c(0.25, 0.5, 0.75),
                      scale = 3, alpha = 0.8, size = 0.3,
                      vline_size = 0.2, vline_color = "black",
                      na.rm = T) +
  # geom_vline(xintercept = c(0.99,1.01), linetype = 'dashed', col = 'red', lwd = 0.5) +
  scale_x_continuous(breaks = seq(0,2,0.1),
                     minor_breaks = seq(0,2,0.01)) +
  labs(y = 'Strains', x = 'Fitness') +
  scale_fill_discrete(guide = F) +
  facet_wrap(.~hours*condition,
             ncol = 4) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.y = element_text(angle = -45, vjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))

plot.box.cs <- ggplot(cs.data[cs.data$condition != 'GLY' & cs.data$orf_name != 'REF' & cs.data$hours > 0,],
       aes(y = average_clean, x = orf_name, fill = orf_name)) +
  geom_boxplot(outlier.shape = NA) +
  # scale_y_continuous(breaks = seq(0,2,0.05),
  #                    minor_breaks = seq(0,2,0.01)) +
  # coord_cartesian(ylim = c(0.8,1.2)) +
  scale_fill_discrete(guide = F) +
  labs(x = 'Strains', y = 'Colony Size (pix.)') +
  facet_wrap(.~hours*condition,
             ncol = 4) +
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

ggsave(sprintf("%s%s_DENSITY_COLONYSIZES.jpg",out_path,expt.name), plot.den.cs,
       height = one.c*1.5, width = two.c*1.5, units = 'mm',
       dpi = 300)
ggsave(sprintf("%s%s_DENSITY_FITNESS.jpg",out_path,expt.name), plot.den.fit,
       height = one.c*1.5, width = two.c*1.5, units = 'mm',
       dpi = 300)
ggsave(sprintf("%s%s_BOXPLOT_COLONYSIZES.jpg",out_path,expt.name), plot.box.cs,
       height = one.c*1.5, width = two.c*1.5, units = 'mm',
       dpi = 300)


##### CALCULATING EFFECT SIZE (CLIFF's DELTA)
cs.es <- NULL
i <- 1
for (d in unique(cs.data$density)) {
  for (h in unique(cs.data$hours[cs.data$density == d])){
    for (c in unique(cs.data$condition[cs.data$density == d & cs.data$hours == h])) {
      for (o1 in unique(cs.data$orf_name[cs.data$density == d & cs.data$hours == h & cs.data$condition == c])) {
        cont_fit <- cs.data$fitness_clean[cs.data$density == d & cs.data$hours == h &
                                            cs.data$condition == c & cs.data$orf_name == o1]
        cont_mean <- median(cs.data$fitness_clean[cs.data$density == d & cs.data$hours == h &
                                                    cs.data$condition == c & cs.data$orf_name == o1], na.rm = T)
        for (o2 in unique(cs.data$orf_name[cs.data$density == d & cs.data$hours == h & cs.data$condition == c &
                                           cs.data$orf_name != o1])){
          orf_fit <- cs.data$fitness_clean[cs.data$density == d & cs.data$hours == h &
                                             cs.data$condition == c & cs.data$orf_name == o2]
          orf_mean <- median(cs.data$fitness_clean[cs.data$density == d & cs.data$hours == h &
                                                     cs.data$condition == c & cs.data$orf_name == o2], na.rm = T)
          temp_cd <- cliff.delta(orf_fit, cont_fit, conf.level=.95,
                                 use.unbiased=TRUE, use.normal=FALSE, return.dm=FALSE)
          temp_p <- wilcox.test(cont_fit, orf_fit, alternative = "two.sided")
          cs.es$density[i] <- d
          cs.es$hours[i] <- h
          cs.es$condition[i] <- c
          cs.es$group1[i] <- o1
          cs.es$group2[i] <- o2
          cs.es$cliff.delta[i] <- temp_cd$estimate
          cs.es$magnitude[i] <- as.character(temp_cd$magnitude)
          cs.es$wilcox[i] <- temp_p$p.value
          cs.es$effect_size[i] <- (orf_mean - cont_mean)/cont_mean * 100
          i <- i + 1
        }
      }
    }
  } 
}
cs.es <- data.frame(cs.es)
cs.es$wilcox.adj <- p.adjust(cs.es$wilcox, method = 'BH')
head(cs.es)

cs.es$group1 <- factor(unique(cs.es$group1),
                       levels = c('Original 1','Original 2','Original 3',
                                  'Copy 1','Copy 2','Copy 3',
                                  'New 1','New 2','New 3','REF'))
cs.es$group2 <- factor(unique(cs.es$group2),
                       levels = c('Original 1','Original 2','Original 3',
                                  'Copy 1','Copy 2','Copy 3',
                                  'New 1','New 2','New 3','REF'))


plot.es.mag <- ggplot(cs.es[cs.es$condition != 'GLY' & cs.es$group1 != 'REF' & cs.es$group2 != 'REF' & cs.es$hours > 0,]) +
  geom_tile(aes(x = group1, y = group2, fill = magnitude), col = 'black') +
  facet_wrap(.~density*hours*condition, ncol = 4) +
  theme_linedraw() +
  scale_fill_discrete(name = 'Magnitude') +
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

plot.es.val <- ggplot(cs.es[cs.es$condition != 'GLY' & cs.es$group1 != 'REF' & cs.es$group2 != 'REF' & cs.es$hours > 0,]) +
  geom_tile(aes(x = group1, y = group2, fill = cliff.delta), col = 'black') +
  facet_wrap(.~density*hours*condition, ncol = 4) +
  theme_linedraw() +
  scale_fill_continuous(name = "Cliff's Delta") +
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

ggsave(sprintf("%s%s_ES_MAGNITUDE.jpg",out_path,expt.name), plot.es.mag,
       height = one.c*1.5, width = two.c*1.5, units = 'mm',
       dpi = 300)
ggsave(sprintf("%s%s_ES_VALUE.jpg",out_path,expt.name), plot.es.val,
       height = one.c*1.5, width = two.c*1.5, units = 'mm',
       dpi = 300)


##### 5FOA EXPERIMENT
##### GATHER DATA
expt.name <- 'FOA'

tblname_dat1 <- sprintf('F28FU_%s_FS_1536_FITNESS', expt.name)
tblname_dat2 <- sprintf('F28FU_%s_FS_6144_FITNESS', expt.name)
tblname_p2o <- sprintf('F28FU_%s_pos2orf_name', expt.name)
tblname_p2c <- sprintf('F28FU_%s_pos2coor', expt.name)

foa.data <- dbGetQuery(conn, sprintf('(select c.*, b.orf_name, a.hours, a.average, a.fitness
                                    from %s a, %s b, %s c
                                    where a.pos = b.pos and a.pos = c.pos)
                                    union
                                    (select c.*, b.orf_name, a.hours, a.average, a.fitness
                                    from %s a, %s b, %s c
                                    where a.pos = b.pos and a.pos = c.pos)
                                    order by density, hours, plate, col, row',
                                    tblname_dat1, tblname_p2o, tblname_p2c,
                                    tblname_dat2, tblname_p2o, tblname_p2c))

foa.data$condition[foa.data$plate == 1] <- 'GLU'
foa.data$condition[foa.data$plate == 2] <- 'CAS'
foa.data$condition[foa.data$plate == 3] <- 'SDA'
foa.data$condition[foa.data$plate == 4] <- 'GAL'
foa.data$condition[foa.data$plate == 5] <- 'GLY'

foa.data$source[foa.data$row%%2==1 & foa.data$col%%2==1] = 'TL'
foa.data$source[foa.data$row%%2==0 & foa.data$col%%2==1] = 'BL'
foa.data$source[foa.data$row%%2==1 & foa.data$col%%2==0] = 'TR'
foa.data$source[foa.data$row%%2==0 & foa.data$col%%2==0] = 'BR'
foa.data$source <- factor(unique(foa.data$source), levels = c('TL','TR','BL','BR'))

foa.data$rep <- as.numeric(str_trunc(as.character(foa.data$pos), 4, side = 'left', ellipsis = ''))

##### CLEAN DATA
head(foa.data)

for (d in unique(foa.data$density)) {
  for (h in unique(foa.data$hours[foa.data$density == d])) {
    for (p in unique(foa.data$plate[foa.data$hours == h & foa.data$density == d])) {
      temp <- foa.data[foa.data$hours == h & foa.data$density == d & foa.data$plate == p,]
      temp$average[isoutlier(temp$average)] <- NA
      foa.data$average_clean[foa.data$hours == h & foa.data$density == d & foa.data$plate == p] <- temp$average
      
      temp$fitness[isoutlier(temp$fitness)] <- NA
      foa.data$fitness_clean[foa.data$hours == h & foa.data$density == d & foa.data$plate == p] <- temp$fitness
      
      for (o in unique(foa.data$orf_name[foa.data$hours == h & foa.data$density == d & foa.data$plate == p & foa.data$orf_name != 'BOR' & !is.na(foa.data$orf_name)])) {
        temp <- foa.data[foa.data$hours == h & foa.data$density == d & foa.data$plate == p & foa.data$orf_name == o & !is.na(foa.data$orf_name),]
        temp$average_clean[isoutlier(temp$average_clean)] <- NA
        foa.data$average_clean[foa.data$hours == h & foa.data$density == d & foa.data$plate == p & foa.data$orf_name == o & !is.na(foa.data$orf_name)] <- temp$average_clean
        
        temp$fitness_clean[isoutlier(temp$fitness_clean)] <- NA
        foa.data$fitness_clean[foa.data$hours == h & foa.data$density == d & foa.data$plate == p & foa.data$orf_name == o & !is.na(foa.data$orf_name)] <- temp$fitness_clean
        
        for (r in unique(foa.data$rep[foa.data$hours == h & foa.data$density == d & foa.data$plate == p & foa.data$orf_name == o & !is.na(foa.data$orf_name)])) {
          temp <- foa.data[foa.data$hours == h & foa.data$density == d & foa.data$plate == p & foa.data$orf_name == o & !is.na(foa.data$orf_name) & foa.data$rep == r,]
          temp$average_clean[isoutlier(temp$average_clean)] <- NA
          foa.data$average_clean[foa.data$hours == h & foa.data$density == d & foa.data$plate == p & foa.data$orf_name == o & !is.na(foa.data$orf_name) & foa.data$rep == r] <- temp$average_clean
          
          temp$fitness_clean[isoutlier(temp$fitness_clean)] <- NA
          foa.data$fitness_clean[foa.data$hours == h & foa.data$density == d & foa.data$plate == p & foa.data$orf_name == o & !is.na(foa.data$orf_name) & foa.data$rep == r] <- temp$fitness_clean
        }
      }
    }
  }
}

##### PLOT COLONY SIZE
plot.den.cs <- ggplot(foa.data[foa.data$condition != 'GLY' & !(foa.data$orf_name %in% c('BOR','REF')) & !is.na(foa.data$orf_name) & foa.data$hours > 0,],
                      aes(x = average_clean, y = orf_name, fill = orf_name)) +
  geom_density_ridges(quantile_lines = TRUE,
                      quantiles = c(0.25, 0.5, 0.75),
                      scale = 3, alpha = 0.8, size = 0.3,
                      vline_size = 0.2, vline_color = "black",
                      na.rm = T) +
  # geom_vline(xintercept = c(0.99,1.01), linetype = 'dashed', col = 'red', lwd = 0.5) +
  # scale_x_continuous(breaks = seq(0,2,0.1),
  #                    minor_breaks = seq(0,2,0.01)) +
  labs(y = 'Strains', x = 'Colony Size (pix.)') +
  scale_fill_discrete(guide = F) +
  facet_wrap(.~density*hours*condition,
             ncol = 4, scales = 'free_x') +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.y = element_text(angle = -45, vjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))

plot.den.fit <- ggplot(foa.data[foa.data$condition != 'GLY' & !(foa.data$orf_name %in% c('BOR','REF')) & !is.na(foa.data$orf_name) & foa.data$hours > 0,],
                       aes(x = fitness_clean, y = orf_name, fill = orf_name)) +
  geom_density_ridges(quantile_lines = TRUE,
                      quantiles = c(0.25, 0.5, 0.75),
                      scale = 3, alpha = 0.8, size = 0.3,
                      vline_size = 0.2, vline_color = "black",
                      na.rm = T) +
  # geom_vline(xintercept = c(0.99,1.01), linetype = 'dashed', col = 'red', lwd = 0.5) +
  scale_x_continuous(breaks = seq(0,2,0.1),
                     minor_breaks = seq(0,2,0.01)) +
  labs(y = 'Strains', x = 'Fitness') +
  scale_fill_discrete(guide = F) +
  facet_wrap(.~density*hours*condition,
             ncol = 4) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.y = element_text(angle = -45, vjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))

plot.box.fit <- ggplot(foa.data[foa.data$condition != 'GLY' & !(foa.data$orf_name %in% c('BOR','REF')) & !is.na(foa.data$orf_name) & foa.data$hours > 0,],
                      aes(y = fitness_clean, x = orf_name, fill = orf_name)) +
  geom_boxplot(outlier.shape = NA) +
  # scale_y_continuous(breaks = seq(0,2,0.05),
  #                    minor_breaks = seq(0,2,0.01)) +
  # coord_cartesian(ylim = c(0.8,1.2)) +
  scale_fill_discrete(guide = F) +
  labs(x = 'Strains', y = 'Fitness') +
  facet_wrap(.~density*hours*condition,
             ncol = 4) +
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

ggsave(sprintf("%s%s_DENSITY_COLONYSIZES.jpg",out_path,expt.name), plot.den.cs,
       height = one.c*1.5*3, width = two.c*1.5, units = 'mm',
       dpi = 300)
ggsave(sprintf("%s%s_DENSITY_FITNESS.jpg",out_path,expt.name), plot.den.fit,
       height = one.c*1.5*3, width = two.c*1.5, units = 'mm',
       dpi = 300)
ggsave(sprintf("%s%s_BOXPLOT_FITNESS.jpg",out_path,expt.name), plot.box.fit,
       height = one.c*1.5*3, width = two.c*1.5, units = 'mm',
       dpi = 300)

# plot.hm.cs <- ggplot(foa.data[foa.data$hours > 0,]) +
#   geom_tile(aes(x = col, y = row, fill = average), col = 'black') +
#   scale_y_reverse() +
#   labs(y = 'Row', x = 'Column', fill = 'Pixel Counts (pix.)') +
#   facet_wrap(.~hours*condition, ncol = 5) +
#   theme_linedraw() +
#   theme(plot.title = element_blank(),
#         axis.title = element_text(size = titles),
#         axis.text = element_text(size = txt),
#         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#         legend.title = element_text(size = titles),
#         legend.text = element_text(size = txt),
#         legend.key.size = unit(3, "mm"),
#         legend.position = "bottom",
#         legend.box.spacing = unit(0.5,"mm"),
#         strip.text = element_text(size = txt,
#                                   margin = margin(0.1,0,0.1,0, "mm")))
# ggsave(sprintf("%s%s_HEATMAP_1536COLONYSIZES.jpg",out_path,expt.name), plot.hm.cs,
#        height = one.c*1.5*2, width = two.c*1.5, units = 'mm',
#        dpi = 300)

##### CALCULATING EFFECT SIZE (CLIFF'S DELTA)
foa.es <- NULL
i <- 1
for (d in unique(foa.data$density)) {
  for (h in unique(foa.data$hours[foa.data$density == d])){
    for (c in unique(foa.data$condition[foa.data$density == d & foa.data$hours == h])) {
      cont_fit <- foa.data$fitness_clean[foa.data$density == d & foa.data$hours == h &
                                           foa.data$condition == c & foa.data$orf_name == 'BF_control']
      cont_mean <- median(foa.data$fitness_clean[foa.data$density == d & foa.data$hours == h &
                                                   foa.data$condition == c & foa.data$orf_name == 'BF_control'], na.rm = T)
      for (o in unique(foa.data$orf_name[foa.data$density == d & foa.data$hours == h & foa.data$condition == c &
                                          foa.data$orf_name != 'BF_control' & foa.data$orf_name != 'BOR' & !is.na(foa.data$orf_name)])){
        orf_fit <- foa.data$fitness_clean[foa.data$density == d & foa.data$hours == h &
                                            foa.data$condition == c & foa.data$orf_name == o]
        orf_mean <- median(foa.data$fitness_clean[foa.data$density == d & foa.data$hours == h &
                                                    foa.data$condition == c & foa.data$orf_name == o], na.rm = T)
        if (!is.na(orf_mean)) {
          temp_cd <- cliff.delta(orf_fit, cont_fit, conf.level=.95,
                                 use.unbiased=TRUE, use.normal=FALSE, return.dm=FALSE)
          foa.es$density[i] <- d
          foa.es$hours[i] <- h
          foa.es$condition[i] <- c
          foa.es$group1[i] <- 'BF_control'
          foa.es$group2[i] <- o
          foa.es$cliff.delta[i] <- temp_cd$estimate
          foa.es$magnitude[i] <- as.character(temp_cd$magnitude)
          foa.es$effect_size[i] <- (orf_mean - cont_mean)/cont_mean * 100
          i <- i + 1
        }
      }
    }
  } 
}
foa.es <- data.frame(foa.es)
head(foa.es)

plot.es.mag <- ggplot(foa.es[foa.es$condition != 'GLY' & foa.es$hours > 0,]) +
  geom_tile(aes(x = group1, y = group2, fill = magnitude), col = 'black') +
  facet_wrap(.~density*hours*condition, ncol = 4) +
  theme_linedraw() +
  scale_fill_discrete(name = 'Magnitude') +
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

plot.es.val <- ggplot(foa.es[foa.es$condition != 'GLY' & foa.es$hours > 0,]) +
  geom_tile(aes(x = group1, y = group2, fill = cliff.delta), col = 'black') +
  facet_wrap(.~density*hours*condition, ncol = 4) +
  theme_linedraw() +
  scale_fill_continuous(name = "Cliff's Delta") +
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

ggsave(sprintf("%s%s_ES_MAGNITUDE.jpg",out_path,expt.name), plot.es.mag,
       height = one.c*3, width = two.c, units = 'mm',
       dpi = 300)
ggsave(sprintf("%s%s_ES_VALUE.jpg",out_path,expt.name), plot.es.val,
       height = one.c*3, width = two.c, units = 'mm',
       dpi = 300)


