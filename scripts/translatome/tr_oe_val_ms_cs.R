##### OVEREXPRESSION VALIDATION MINI SCREEN ANALYSIS
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 03/17/2021 

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

expt.name <- "tr_oe_val_ms"

source("~/R/Projects/adaptivefitness/R/functions/isoutlier.R")
source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

`%notin%` <- Negate(`%in%`)

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 8
txt <- 8
lbls <- 9

##### INITIALIZE
tr.conds <- readxl::read_xlsx(path = '/home/sbp29/RAW_Data/TranslatomeOE_VAL/TR_OE_VAL_MS_INFO2.xlsx') %>% data.frame()
borders <- dbGetQuery(conn, 'select * from TR_OE_VAL_MS_borderpos')

orf_types <- read.csv(file = '/home/sbp29/R/Projects/adaptivefitness/rawdata/translatome/oe_transient')
orf_types$category[orf_types$is_transient + orf_types$translated + orf_types$is_candidate == 3] <- 'Transient'
orf_types$category[orf_types$is_conserved + orf_types$translated + orf_types$is_candidate == 3] <- 'Conserved'
orf_types$category[is.na(orf_types$category)] <- 'Others'

##### GATHER DATA
data.cs.all <- NULL
for (a in unique(tr.conds$arm)) {
  temp.cs <- dbGetQuery(conn, sprintf('select a.*, b.density, b.plate_no, b.plate_col, b.plate_row,
                                         d.strain_id, c.orf_name
                                         from TR_OE_VAL_MS_FS_%s_6144_RAW a, TR_OE_VAL_MS_pos2coor b,
                                         TR_OE_VAL_MS_pos2orf_name c, TR_OE_VAL_MS_pos2strainid d
                                         where a.pos = b.pos and b.pos = c.pos and c.pos = d.pos
                                         order by a.hours, b.plate_no, b.plate_col, b.plate_row', a))
  temp.cs$condition <- a
  temp.cs$saturation <- max(temp.cs$hours)
  data.cs.all <- rbind(data.cs.all, temp.cs)
}
data.cs.all$average[data.cs.all$pos %in% borders$pos] <- NA
data.cs.all$average[data.cs.all$average <= 300] <- NA

# GENERATING CLEAN TABLES
# for (a in unique(data.cs.all$condition)) {
#   dbWriteTable(conn, sprintf('TR_OE_VAL_MS_FS_%s_6144_CLEAN',a), data.cs.all[data.cs.all$condition == a,
#                                                                               c('pos','hours','average')],
#                overwrite = T)
# }

##### ADDING FEATURES
data.cs.all <- merge(data.cs.all, orf_types, by = 'orf_name', all.x = T)
data.cs.all$category[data.cs.all$orf_name == 'BF_control'] <- 'Reference'

data.cs.all$rep <- as.numeric(str_trunc(as.character(data.cs.all$pos), 5, side = 'left', ellipsis = ''))

# ##### REMOVING OUTLIER DATA
# data.cs.sum <- data.cs.all %>%
#   group_by(condition, hours, rep, strain_id, orf_name, category) %>%
#   summarize(avg.median = median(average, na.rm = T), avg.mad = mad(average, na.rm = T),
#             .groups = 'keep') %>%
#   data.frame()
# 
# data.cs.all <- merge(data.cs.all, data.cs.sum,
#                       by = c('condition','hours','rep','strain_id','orf_name','category'), all = T)
# 
# data.cs.all$average[data.cs.all$average < (data.cs.all$avg.median - 2*data.cs.all$avg.mad) |
#                        data.cs.all$average > (data.cs.all$avg.median + 2*data.cs.all$avg.mad)] <- NA
# 
## REPLICATE COUNT PER REP
data.cs.cnt <- data.cs.all[!is.na(data.cs.all$average),] %>%
  group_by(condition, hours, rep) %>%
  count() %>%
  data.frame()
data.cs.all <- merge(data.cs.all, data.cs.cnt, by = c('condition','hours','rep'), all.x = T)
# data.cs.all$n[is.na(data.cs.all$n)] <- 0

# ## REMOVE DATA WITH LESS THAN 4 REPLICATES
# data.cs.all$average[data.cs.all$n < 4] <- NA
# 
# ##### COLONY SIZE SUMMARY
# data.cs.sum<- data.cs.all %>%
#   group_by(condition, hours, rep, strain_id, orf_name, n, category) %>%
#   summarize(avg.median = median(average, na.rm = T),
#             .groups = 'keep') %>%
#   data.frame()
# head(data.cs.sum)

##### SAVING COLONY SIZE AND SUMMARY DATA
# write.csv(data.cs.sum, file = '/home/sbp29/R/Projects/adaptivefitness/output/translatome/pinning_artifact_correction/tr_oe_val_ms_cs_summary.csv')
write.csv(data.cs.all[,c(1:12,21)] %>% filter(orf_name != 'NULL'),
          file = '/home/sbp29/R/Projects/adaptivefitness/output/translatome/pinning_artifact_correction/tr_oe_val_ms_cs_all.csv')

##### PLOTTING CS DATA
data.cs.all %>%
  filter(orf_name != 'NULL') %>%
  ggplot(aes(x = hours, y = average, group = rep)) +
  stat_summary(geom = 'line', fun.x = 'mean') +
  facet_wrap(.~condition)

