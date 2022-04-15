##### DELETION SCREEN ANALYSIS
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 12/09/2021 

##### INITIALIZE
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyverse)
library(ggpubr)
library(ggridges)
library(stringr)
library(egg)
library(zoo)
library(ggrepel)
library(ggforce)
library(plotly)
library(scales)
library(reshape2)
library(locfit)
library(growthcurver)
library(rstatix)
library(gtools)
library(growthrates)
library(RMariaDB)
library(genefilter)
library(apeglm)
library(clusterProfiler)
library(org.Sc.sgd.db)

out_path <- "~/R/Projects/adaptivefitness/output/translatome/"
fig_path <- "~/R/Projects/adaptivefitness/figs/translatome/"

expt.name <- "tr_del_all"

source("~/R/Projects/adaptivefitness/R/functions/isoutlier.R")
source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

`%notin%` <- Negate(`%in%`)

# load('output/translatome/AllData.RData')

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 8
txt <- 8
lbls <- 9

##### INITIALIZE
tr.conds <- data.frame(arms = c('R1','R2','R1','R2','R1','R2','R1','R2','R3','R1','R2','R3','R1','R2','R1','R2','R3'),
                       conds = c('YPDA','YPDA','DM','DM','HU','HU','HO','HO','HO','TN','TN','TN','FL','FL','SA','SA','SA'))

borders <- dbGetQuery(conn, 'select * from TR_DEL_borderpos')

##### GATHER DATA
data.fit.all <- NULL
for (c in unique(tr.conds$conds)) {
  for (a in tr.conds$arms[tr.conds$conds == c]) {
    temp.fit <- dbGetQuery(conn, sprintf('select a.*, b.density, b.plate, b.col, b.row
                                         from TR_DEL_FS_%s_%s_1536_FITNESS a, TR_DEL_pos2coor b
                                         where a.pos = b.pos
                                         order by a.hours, b.plate, b.col, b.row', c, a))
    temp.fit$condition <- c
    temp.fit$arm <- a
    temp.fit$saturation <- max(temp.fit$hours)
    data.fit.all <- rbind(data.fit.all, temp.fit)
  }
}
data.fit.all$average[data.fit.all$pos %in% borders$pos] <- NA
data.fit.all <- data.fit.all[!(data.fit.all$condition == 'FL' & data.fit.all$arm == 'R1'),]
data.fit.all <- data.fit.all[!(data.fit.all$condition == 'FL' & data.fit.all$arm == 'R3'),]
data.fit.all <- data.fit.all[!(data.fit.all$condition == 'HO' & data.fit.all$arm == 'R3'),]

data.fit.all$rep <- as.numeric(str_trunc(as.character(data.fit.all$pos), 4, side = 'left', ellipsis = ''))

data.fit.all %>%
  filter(orf_name != 'BOR') %>%
  ggplot(aes(x = hours, y = average)) +
  geom_line(aes(group = pos)) +
  facet_grid(arm ~ condition)

data.fit.sum <- data.fit.all %>%
  filter(orf_name != 'BOR') %>%
  group_by(arm, condition, hours, rep, strain_id, orf_name) %>%
  summarize(avg.median = median(average, na.rm = T), avg.mad = mad(average, na.rm = T),
            fitness.median = median(fitness, na.rm = T), fitness.mad = mad(fitness, na.rm = T),
            .groups = 'keep') %>%
  data.frame()

##### REMOVE OUTLIERS
data.fit.all <- merge(data.fit.all, data.fit.sum,
                      by = c('arm','condition','hours','rep','strain_id','orf_name'), all = T)

data.fit.all$average[data.fit.all$average < (data.fit.all$avg.median - 2*data.fit.all$avg.mad) |
                       data.fit.all$average > (data.fit.all$avg.median + 2*data.fit.all$avg.mad)] <- NA
data.fit.all$fitness[data.fit.all$fitness < (data.fit.all$fitness.median - 2*data.fit.all$fitness.mad) |
                       data.fit.all$fitness > (data.fit.all$fitness.median + 2*data.fit.all$fitness.mad)] <- NA

data.fit.sum <- data.fit.all %>%
  group_by(condition, arm, hours, rep, strain_id, orf_name) %>%
  summarize(avg.median = median(average, na.rm = T), avg.mad = mad(average, na.rm = T),
            fitness.median = median(fitness, na.rm = T), fitness.mad = mad(fitness, na.rm = T),
            .groups = 'keep') %>%
  data.frame()
data.fit.sum <- data.fit.sum[!is.na(data.fit.sum$orf_name) & data.fit.sum$orf_name != 'BOR',]

data.fit.all <- merge(data.fit.all[,c(1:15)], data.fit.sum,
                      by = c('arm','condition','hours','rep','strain_id','orf_name'), all = T)


##### REFERENCE LIMITS
data.fit.lim <- data.fit.all %>%
  filter(orf_name == 'HO') %>%
  group_by(arm, condition, hours, orf_name, rep) %>%
  summarize(average = median(average, na.rm = T),
            fitness = median(fitness, na.rm = T),
            .groups = 'keep') %>%
  group_by(arm, condition, hours, orf_name) %>%
  summarize(avg_ll = quantile(average, 0.025, na.rm = T),
            avg_m = median(average, na.rm = T),
            avg_ul = quantile(average, 0.975, na.rm = T),
            fitness_ll = quantile(fitness, 0.025, na.rm = T),
            fitness_m = median(fitness, na.rm = T),
            fitness_ul = quantile(fitness, 0.975, na.rm = T),
            .groups = 'keep') %>%
  data.frame()

data.fit.all <- merge(data.fit.all, data.fit.lim[,-4], by = c('arm','condition','hours'))

##### SAVE PRELIM DATA
write.csv(data.fit.all[,c(1:16)], file = 'output/translatome/tr_del_fitandcs_highres_all.csv')
write.csv(data.fit.sum, file = 'output/translatome/tr_del_fitandcs_highres_summary_armwise.csv')

##### GROWTH CURVES FOR AARON'S HITS
del.hits <- read_delim(file = 'output/translatome/deletion/deletion_mutant_highly_deleterious.txt', delim = ' ')
match.ids <- read_delim(file = 'rawdata/translatome/match_orf_ids_2012', delim = ' ')
match.ids$full_orf_id <- str_replace(match.ids$full_orf_id, 'chr_','chr')

data.fit.hit.all <- merge(data.fit.all, del.hits[,c(2,3,5:7)], by = c('condition','arm','rep','strain_id','orf_name'))

data.fit.all %>%
  ggplot(aes(x = hours, y = average)) +
  geom_line(aes(group = pos), alpha = 0.7) +
  geom_line(data = data.fit.hit.all, aes(x = hours, y = average, group = pos),
            col = 'blue') +
  geom_line(data = data.fit.lim, aes(x = hours, y = avg_ul),
            col = 'red', linetype = 'dashed') +
  geom_line(data = data.fit.lim, aes(x = hours, y = avg_ll),
            col = 'red', linetype = 'dashed') +
  facet_grid(arm ~ condition)


##### PLATE MAPS
data.fit.all %>%
  filter(condition == 'SA', hours == 32) %>%
  ggplot(aes(x = col, y = row)) +
  geom_tile(aes(fill = fitness), col = 'black') +
  scale_y_reverse() +
  facet_wrap(.~arm)

#####
s2o.tbl <- readxl::read_excel('/home/sbp29/MATLAB/TR_DEL_S2O.xlsx') %>% data.frame()

write.csv(merge(s2o.tbl, match.ids, by.x = 'orf_name', by.y = 'full_orf_id', all.x = T),
          file = 'output/translatome/deletion/del_mutants.csv')
write.csv(merge(merge(s2o.tbl, match.ids, by.x = 'orf_name', by.y = 'full_orf_id', all.x = T), del.hits[,c(2,3,7)], by = 'orf_name'),
          file = 'output/translatome/deletion/del_mutants_hits.csv')

write.csv(oe.hits[,c(3,7)], file = 'output/translatome/oe_mutant_fitness_hits.csv')
