##### OVEREXPRESSION SCREEN ANALYSIS
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 10/28/2021 

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

expt.name <- "tr_oe_all"

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
controls <- read.csv(file = '/home/sbp29/R/Projects/adaptivefitness/rawdata/translatome/condition_controls.csv', stringsAsFactors = F)

orf_types <- read.csv(file = '/home/sbp29/R/Projects/adaptivefitness/rawdata/translatome/oe_transient')
orf_types$category[orf_types$is_transient + orf_types$translated + orf_types$is_candidate == 3] <- 'Transient'
orf_types$category[orf_types$is_conserved + orf_types$translated + orf_types$is_candidate == 3] <- 'Conserved'
orf_types$category[is.na(orf_types$category)] <- 'Others'

tr.conds <- data.frame(arms = c('ONE','ONE','ONE','ONE','TWO','TWO','TWO','TWO'),
                       conds = c('GA','SA','HO','HU','GA','DM','FL','TN'))

borders <- dbGetQuery(conn, 'select * from TRANS_OE_borderpos')
smudge <- dbGetQuery(conn, 'select * from TR_OE_ONE_FS_GA_smudgebox')


##### GATHER DATA
data.fit.all <- NULL
for (a in unique(tr.conds$arms)) {
  for (c in tr.conds$conds[tr.conds$arms == a]) {
    temp.fit <- dbGetQuery(conn, sprintf('select a.*, b.density, b.plate_no, b.plate_col, b.plate_row
                                         from TR_OE_FS_%s_%s_6144_FITNESS a, TRANS_OE_pos2coor b
                                         where a.pos = b.pos
                                         order by a.hours, b.plate_no, b.plate_col, b.plate_row', a, c))
    temp.fit$arm <- a
    temp.fit$condition <- c
    temp.fit$saturation <- max(temp.fit$hours)
    data.fit.all <- rbind(data.fit.all, temp.fit)
  }
}
data.fit.all$average[data.fit.all$pos %in% borders$pos] <- NA
data.fit.all$average[data.fit.all$condition == 'HO' & data.fit.all$plate %in% c(2,3,10,18,19,23,26,32)] <- NA
data.fit.all$average[data.fit.all$arm == 'ONE' & data.fit.all$plate == 14] <- NA
data.fit.all$average[data.fit.all$arm == 'ONE' & data.fit.all$pos %in% smudge$pos] <- NA

data.fit.all$fitness[data.fit.all$pos %in% borders$pos] <- NA
data.fit.all$fitness[data.fit.all$condition == 'HO' & data.fit.all$plate %in% c(2,3,10,18,19,23,26,32)] <- NA
data.fit.all$fitness[data.fit.all$arm == 'ONE' & data.fit.all$plate == 14] <- NA
data.fit.all$fitness[data.fit.all$arm == 'ONE' & data.fit.all$pos %in% smudge$pos] <- NA

data.fit.all <- merge(data.fit.all, orf_types, by = 'orf_name', all.x = T)
data.fit.all$category[data.fit.all$orf_name == 'BF_control'] <- 'Reference'
# data.fit.all[is.na(data.fit.all$category),] %>%
#   group_by(orf_name) %>%
#   count() %>% data.frame()
data.fit.all$rep <- as.numeric(str_trunc(as.character(data.fit.all$pos), 6, side = 'left', ellipsis = ''))

data.fit.all$source[data.fit.all$plate_row%%2==1 & data.fit.all$plate_col%%2==1] = 'TL'
data.fit.all$source[data.fit.all$plate_row%%2==0 & data.fit.all$plate_col%%2==1] = 'BL'
data.fit.all$source[data.fit.all$plate_row%%2==1 & data.fit.all$plate_col%%2==0] = 'TR'
data.fit.all$source[data.fit.all$plate_row%%2==0 & data.fit.all$plate_col%%2==0] = 'BR'

data.fit.all <- data.fit.all[!(data.fit.all$condition == 'FL' & data.fit.all$hours == 192) &
                               !(data.fit.all$condition == 'GA' & data.fit.all$hours == 67) &
                               !(data.fit.all$condition == 'SA' & data.fit.all$hours %in% c(51,63,72)),]

data.fit.sum <- data.fit.all %>%
  group_by(arm, condition, hours, rep, strain_id, orf_name, category) %>%
  summarize(avg.median = median(average, na.rm = T), avg.mad = mad(average, na.rm = T),
            fitness.median = median(fitness, na.rm = T), fitness.mad = mad(fitness, na.rm = T),
            .groups = 'keep') %>%
  data.frame()

##### REMOVE OUTLIERS
data.fit.all <- merge(data.fit.all, data.fit.sum,
                      by = c('arm','condition','hours','rep','strain_id','orf_name','category'), all = T)

data.fit.all$fitness[data.fit.all$fitness < (data.fit.all$fitness.median - 2*data.fit.all$fitness.mad) |
                       data.fit.all$fitness > (data.fit.all$fitness.median + 2*data.fit.all$fitness.mad)] <- NA
data.fit.all$average[is.na(data.fit.all$fitness)] <- NA

data.fit.sum <- data.fit.all %>%
  group_by(arm, condition, hours, rep, strain_id, orf_name, category) %>%
  summarize(avg.median = median(average, na.rm = T),
            fitness.median = median(fitness, na.rm = T),
            .groups = 'keep') %>%
  data.frame()
head(data.fit.sum)

data.fit.all <- merge(data.fit.all[,c(1:22)], data.fit.sum,
                      by = c('arm','condition','hours','rep','strain_id','orf_name','category'), all = T)

##### REFERENCE LIMITS
data.fit.lim <- data.fit.all %>%
  filter(orf_name == 'BF_control') %>%
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

data.fit.sum <- merge(data.fit.sum, data.fit.lim[,-4], by = c('arm','condition','hours'))

##### GROWTH CURVE ANALYSIS
data.fit.lim.src <- data.fit.all %>%
  filter(orf_name == 'BF_control') %>%
  group_by(arm, condition, hours, source) %>%
  summarize(avg_ll = quantile(average, 0.025, na.rm = T),
            avg_m = median(average, na.rm = T),
            avg_ul = quantile(average, 0.975, na.rm = T),
            fitness_ll = quantile(fitness, 0.025, na.rm = T),
            fitness_m = median(fitness, na.rm = T),
            fitness_ul = quantile(fitness, 0.975, na.rm = T),
            .groups = 'keep') %>%
  data.frame()

data.fit.all <- merge(data.fit.all[,c(1:22)], data.fit.lim.src, by = c('arm','condition','hours','source'))
data.fit.all$norm_cs <- data.fit.all$fitness * data.fit.all$avg_m
data.fit.all$pos <- as.numeric(data.fit.all$pos)

data.strains <- data.fit.all %>%
  group_by(arm, condition, rep, strain_id, orf_name, category) %>%
  count() %>% data.frame()

data.gr.res <- NULL
data.gc.res <- NULL
for (a in unique(data.fit.all$arm)) {
  for (c in unique(data.fit.all$condition[data.fit.all$arm == a])) {
    sprintf('Now analyzing %s %s.', a, c)
    data.pred <- NULL
    col.names <- NULL
    data.temp <- data.fit.all[data.fit.all$arm == a & data.fit.all$condition == c,]
    for (p in unique(data.temp$pos)) {
      r <- unique(data.temp$rep[data.temp$pos == p])
      s <- unique(data.temp$strain_id[data.temp$pos == p])
      temp <- data.temp[data.temp$pos == p,]
      if (sum(!is.na(temp$average)) >= length(temp$average) * 0.7) {
        lo <- loess.smooth(temp$hours, log(temp$average),
                           span = 0.6, evaluation = 100, degree = 2,
                           family = 'gaussian')
        data.pred <- cbind(data.pred,exp(lo$y))
        # data.pred <- cbind(data.pred, temp$average)
        col.names <- cbind(col.names,paste(a,c,p,r,s,sep = ','))
      }
    }
    data.pred <- cbind(lo$x, data.pred)
    # data.pred <- cbind(temp$hours, data.pred)
    data.pred <- data.frame(data.pred)
    colnames(data.pred) <- c('Time',col.names)
    # save(data.pred, file = sprintf('/home/sbp29/R/Projects/adaptivefitness/output/translatome/data.pred.pos.%s.%s.RData',a,c))

    ## GROWTH CURVE ANALYSIS
    temp.gr.res <- NULL
    for (i in 2:dim(data.pred)[2]) {
      fit0 <- fit_easylinear(data.pred$Time[2:dim(data.pred)[1]], data.pred[2:dim(data.pred)[1],i], h = 8, quota = 1);

      temp_res <- data.frame(colnames(data.pred[i]), maxgr = coef(fit0)[[3]],
                             dtime = log(2)/coef(fit0)[[3]], ltime = coef(fit0)[[4]])
      temp.gr.res <- rbind(temp.gr.res, temp_res)
    }
    temp.gr.res <- data.frame(temp.gr.res)
    colnames(temp.gr.res) <- c('sample','gr','dtime','lag')
    temp <- str_split(temp.gr.res$sample, ',', simplify = T)
    colnames(temp) <- c('arm','condition','pos','rep','strain_id')
    temp.gr.res <- cbind(temp, temp.gr.res)

    temp.gc.res <- SummarizeGrowthByPlate(data.pred)
    temp <- str_split(temp.gc.res$sample, ',', simplify = T)
    colnames(temp) <- c('arm','condition','pos','rep','strain_id')
    temp.gc.res <- cbind(temp, temp.gc.res)

    data.gr.res <- rbind(data.gr.res, temp.gr.res)
    data.gc.res <- rbind(data.gc.res, temp.gc.res)
  }
}
# data.gcgr.res <- merge(data.gc.res, data.gr.res, by = c('arm','condition','pos','rep','strain_id','sample'))
# data.gcgr.res <- merge(data.strains[,-7], data.gcgr.res, by = c('arm','condition','pos','rep','strain_id'))

save.image(file = '/home/sbp29/R/Projects/adaptivefitness/output/translatome/pinning_artifact_correction/INIT_CS_GROWTH_CURVE_ALL.RData')




