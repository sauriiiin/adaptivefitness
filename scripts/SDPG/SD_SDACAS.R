##### SAN DIEGO CAS SD
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 11/18/2020

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

out_path = 'figs/SDPG/SD/';

##### FIGURE TEXT SIZE
titles <- 7
txt <- 7
lbls <- 9

##### GATHER DATA
sdacas <- dbGetQuery(conn, 'select a.*, b.orf_name, b.fitness
                    from brian_031918.JPEG_RESULTS_v2_BF5_CASEGAL_screen a, brian_031918.FITNESS_v2_BF5_CASEGAL_screen b
                    where a.exp_id = 91 and a.hours = 35
                    and a.exp_id = b.exp_id and a.hours = b.hours and a.pid = b.pid
                    UNION
                    select e.*, f.orf_name, f.fitness
                    from brian_031918.JPEG_RESULTS_v2_BF5_CASEGAL_screen e, brian_031918.FITNESS_v2_BF5_CASEGAL_screen f
                    where e.exp_id = 93 and e.hours = 58
                    and e.exp_id = f.exp_id and e.hours = f.hours and e.pid = f.pid')

##### FITNESS PIX SCATTER PLOT
ggplot(sdacas,
       aes(x = fitness, y = average)) +
  geom_point() +
  stat_cor(method = 'spearman') +
  facet_grid(.~exp_id)

#####
head(sdacas)


##### LID VS MCAT at PITT
lid_mcat <- dbGetQuery(conn, 'select a.orf_name, a.hours, a.cs_median lid_fitness, b.cs_median mcat_fitness
                      from SDPG_CAS_FS_R1_6144_FITNESS_STATS a,
                       SDPG_CAS_FS_R1_6144_MCAT_FITNESS_STATS b
                       where a.orf_name = b.orf_name and a.hours = b.hours
                       and a.hours = 24')

ggplot(lid_mcat,
       aes(x = lid_fitness, y = mcat_fitness)) +
  geom_abline() +
  geom_point() + 
  stat_cor(method = 'spearman')

load('figs/SDPG/results.RDATA')

mcatres <- mcatres[mcatres$orf_name != 'YHR021W-A' &
                   mcatres$orf_name != 'BOR' &
                   mcatres$orf_name != 'REF' &
                   !is.na(mcatres$orf_name),]
# Removing plates with bad colony grids
mcatres <- mcatres[!(mcatres$density == 6144 & mcatres$arm == 'GLU' & mcatres$rep != 'R3' & mcatres$hours >= 20),]
mcatres <- mcatres[!(mcatres$density == 6144 & mcatres$arm == 'CAS' & mcatres$rep == 'R2' & mcatres$hours >= 20),]

# Removing outliers
for (d in unique(mcatres$density)) {
  for (a in unique(mcatres$arm[mcatres$density == d])) {
    for (r in unique(mcatres$rep[mcatres$density == d & 
                                mcatres$arm == a])) {
      for (h in unique(mcatres$hours[mcatres$density == d & 
                                    mcatres$arm == a & 
                                    mcatres$rep == r])) {
        for (o in unique(mcatres$orf_name[mcatres$density == d & 
                                         mcatres$arm == a & 
                                         mcatres$rep == r &
                                         mcatres$hours == h])) {
          temp <- mcatres$fitness[mcatres$density == d & 
                                   mcatres$arm == a & 
                                   mcatres$rep == r &
                                   mcatres$hours == h &
                                   mcatres$orf_name == o]
          m <- median(temp, na.rm = T)
          madev <- mad(temp)
          ul <- m + 3*madev
          ll <- m - 3*madev
          
          mcatres$fitness[mcatres$density == d & 
                           mcatres$arm == a & 
                           mcatres$rep == r &
                           mcatres$hours == h &
                           mcatres$orf_name == o &
                           mcatres$fitness > ul] <- NA
          mcatres$fitness[mcatres$density == d & 
                           mcatres$arm == a & 
                           mcatres$rep == r &
                           mcatres$hours == h &
                           mcatres$orf_name == o &
                           mcatres$fitness < ll] <- NA
          
          
        }
        cont_median <- median(mcatres$average[mcatres$density == d & 
                                               mcatres$arm == a & 
                                               mcatres$rep == r &
                                               mcatres$hours == h &
                                               mcatres$orf_name == 'BF_control'], na.rm = T)
        mcatres$ccs[mcatres$density == d & 
                     mcatres$arm == a & 
                     mcatres$rep == r &
                     mcatres$hours == h] <- mcatres$cs_median[mcatres$density == d & 
                                                              mcatres$arm == a & 
                                                              mcatres$rep == r &
                                                              mcatres$hours == h] * cont_median
      }
    }
  }
}

# Saturation information
mcatres$time <- NA
mcatres$time[mcatres$density == 1536 & mcatres$hours == 48 & mcatres$arm != 'SDA'] <- "Saturated"
mcatres$time[mcatres$density == 1536 & mcatres$hours %in% c(72, 78) & mcatres$arm == 'SDA'] <- "Saturated"
mcatres$time[mcatres$density == 6144 & mcatres$hours == 24 & mcatres$arm != 'SDA'] <- "Saturated"
mcatres$time[mcatres$density == 6144 & mcatres$hours == 36 & mcatres$arm == 'SDA'] <- "Saturated"
mcatres$time[is.na(mcatres$time)] <- "Other"

# Factorize the Experimental Arms
mcatres$arm <- factor(mcatres$arm, levels = c('GLU','CAS','SDA'))
