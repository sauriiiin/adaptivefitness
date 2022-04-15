##### OVEREXPRESSION SCREEN HITS
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 01/06/2022 

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
library(reshape2)
library(RMariaDB)

source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

load("~/R/Projects/adaptivefitness/output/translatome/220106.RData")
tr.oe.fit.hits <- read.table(file = '/home/acwach/overexp/oe_hits_05fdr.txt')
tr.oe.fit.diff.hits <- read.table(file = '/home/acwach/overexp/oe_hits_fit_diff_fdr05.txt')

`%notin%` <- Negate(`%in%`)

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 8
txt <- 8
lbls <- 9

##### GATHER DATA FROM MySQL
platemap <- dbGetQuery(conn, 'select * from BARFLEX_SPACE_AGAR_180313
                       union
                       select * from PROTOGENE_COLLECTION')

##### COUNTS
tr.oe.fit.hits %>%
  group_by(arm, condition, hours) %>%
  count()
# write.csv(tr.oe.fit.hits, file = 'output/translatome/tr_oe_fit_hits.csv')

tr.oe.fit.diff.hits %>%
  group_by(arm, condition, hours) %>%
  count()
# write.csv(tr.oe.fit.diff.hits, file = 'output/translatome/tr_oe_fit_diff_hits.csv')

merge(tr.oe.fit.hits[,c('orf_id','strain_id','arm','condition')],
      tr.oe.fit.diff.hits[,c('orf_id','strain_id','arm','condition')],
      by = c('orf_id','strain_id','arm','condition')) %>%
  group_by(arm, condition) %>%
  count()

head(tr.oe.fit.diff.hits)
head(platemap)

merge(tr.oe.fit.diff.hits, platemap, by.x = c('orf_id','strain_id'), by.y = c('orf_name','strain_id')) %>%
  group_by(arm, condition, `384plate`) %>%
  count() %>%
  data.frame()

merge(tr.oe.fit.diff.hits, platemap, by.x = c('orf_id','strain_id'), by.y = c('orf_name','strain_id')) %>%
  group_by(arm, condition, `384plate`) %>%
  filter(`384plate` == 27, condition == 'TN') %>%
  data.frame() %>%
  ggplot(aes(x = X384col, y = X384row)) +
  geom_tile(col = 'black') +
  coord_cartesian(xlim = c(0,24),
                  ylim = c(0,16))

##### FILTERING THE TOP HITS
rbind(merge(tr.oe.fit.hits[,c('orf_id','strain_id','arm','condition')],
      tr.oe.fit.diff.hits[,c('orf_id','strain_id','arm','condition')],
      by = c('orf_id','strain_id','arm','condition')) %>%
  group_by(arm, condition),
  tr.oe.fit.hits[,c('orf_id','strain_id','arm','condition')] %>%
    filter(condition == 'HU')) %>%
  group_by(arm, condition) %>%
  count()

temp.sa.u <- tr.oe.fit.diff.hits[str_detect(tr.oe.fit.diff.hits$orf_id, 'smor') &
                                   tr.oe.fit.diff.hits$orf_id %notin% tr.oe.fit.hits[tr.oe.fit.hits$condition == 'SA','orf_id'],
                                 c('orf_id','strain_id','arm','condition','mean_fitness_diff')] %>% filter(condition == 'SA')
temp.sa.u <- head(temp.sa.u[order(-temp.sa.u$mean_fitness_diff),c('orf_id','strain_id','arm','condition')],2)
temp.sa.a <- tr.oe.fit.diff.hits[!str_detect(tr.oe.fit.diff.hits$orf_id, 'smor') &
                                   tr.oe.fit.diff.hits$orf_id %notin% tr.oe.fit.hits[tr.oe.fit.hits$condition == 'SA','orf_id'],
                                 c('orf_id','strain_id','arm','condition','mean_fitness_diff')] %>% filter(condition == 'SA')
temp.sa.a <- head(temp.sa.a[order(-temp.sa.a$mean_fitness_diff),c('orf_id','strain_id','arm','condition')],2)

temp.fl.u <- tr.oe.fit.diff.hits[str_detect(tr.oe.fit.diff.hits$orf_id, 'smor') &
                                   tr.oe.fit.diff.hits$orf_id %notin% tr.oe.fit.hits[tr.oe.fit.hits$condition == 'FL','orf_id'],
                                 c('orf_id','strain_id','arm','condition','mean_fitness_diff')] %>% filter(condition == 'FL')
temp.fl.u <- head(temp.fl.u[order(-temp.fl.u$mean_fitness_diff),c('orf_id','strain_id','arm','condition')],3)
temp.fl.a <- tr.oe.fit.diff.hits[!str_detect(tr.oe.fit.diff.hits$orf_id, 'smor') &
                                   tr.oe.fit.diff.hits$orf_id %notin% tr.oe.fit.hits[tr.oe.fit.hits$condition == 'FL','orf_id'],
                                 c('orf_id','strain_id','arm','condition','mean_fitness_diff')] %>% filter(condition == 'FL')
temp.fl.a <- head(temp.fl.a[order(-temp.fl.a$mean_fitness_diff),c('orf_id','strain_id','arm','condition')],3)

temp.tn.u <- tr.oe.fit.diff.hits[str_detect(tr.oe.fit.diff.hits$orf_id, 'smor') &
                                   tr.oe.fit.diff.hits$orf_id %notin% tr.oe.fit.hits[tr.oe.fit.hits$condition == 'TN','orf_id'],
                                 c('orf_id','strain_id','arm','condition','mean_fitness_diff')] %>% filter(condition == 'TN')
temp.tn.u <- head(temp.tn.u[order(-temp.tn.u$mean_fitness_diff),c('orf_id','strain_id','arm','condition')],5)
temp.tn.a <- tr.oe.fit.diff.hits[!str_detect(tr.oe.fit.diff.hits$orf_id, 'smor') &
                                   tr.oe.fit.diff.hits$orf_id %notin% tr.oe.fit.hits[tr.oe.fit.hits$condition == 'TN','orf_id'],
                                 c('orf_id','strain_id','arm','condition','mean_fitness_diff')] %>% filter(condition == 'TN')
temp.tn.a <- head(temp.tn.a[order(-temp.tn.a$mean_fitness_diff),c('orf_id','strain_id','arm','condition')],4)


val.hits <- rbind(merge(tr.oe.fit.hits[,c('orf_id','strain_id','arm','condition')],
            tr.oe.fit.diff.hits[,c('orf_id','strain_id','arm','condition')],
            by = c('orf_id','strain_id','arm','condition')) %>%
        group_by(arm, condition),
      tr.oe.fit.hits[,c('orf_id','strain_id','arm','condition')] %>%
        filter(condition == 'HU'),
      temp.sa.u, temp.sa.a,
      temp.fl.u, temp.fl.a,
      temp.tn.u, temp.tn.a) %>%
  group_by(arm, condition) %>%
  data.frame()
val.hits <- val.hits[order(val.hits$condition),]
# val.hits <- val.hits[order(val.hits$orf_id),]
# write.csv(val.hits, file = 'output/translatome/tr_oe_hits2validate.csv')

platemap <- platemap[platemap$strain_id %in% val.hits$strain_id,]
platemap <- platemap[order(platemap$`384plate`, platemap$`384col`, platemap$`384row`),]
# write.csv(platemap, file = 'output/translatome/tr_oe_hits2validate_location.csv')
