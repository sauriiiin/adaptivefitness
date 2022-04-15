##### OVEREXPRESSION SCREEN ANALYSIS
##### for data used in the paper
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 02/07/2022 

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
library(ggforce)
library(plotly)
library(scales)
library(reshape2)
library(RMariaDB)

out_path <- "~/R/Projects/adaptivefitness/output/translatome/"
fig_path <- "~/R/Projects/adaptivefitness/figs/translatome/"

expt.name <- "tr_oe_paper"

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
controls <- read.csv(file = '/home/sbp29/R/Projects/adaptivefitness/rawdata/translatome/condition_controls.csv', stringsAsFactors = F)
pgs <- dbGetQuery(conn, 'select orf_name from PROTOGENES where pg_2012 = 1')

orf_types <- read.csv(file = '/home/sbp29/R/Projects/adaptivefitness/rawdata/translatome/oe_transient')
orf_types$category[orf_types$is_transient + orf_types$translated + orf_types$is_candidate == 3] <- 'Transient'
orf_types$category[orf_types$is_conserved + orf_types$translated + orf_types$is_candidate == 3] <- 'Conserved'
orf_types$category[is.na(orf_types$category)] <- 'Others'

tr.conds <- data.frame(arms = c('ONE','ONE','ONE','ONE','TWO','TWO','TWO'),
                       conds = c('GA','SA','HO','HU','DM','FL','TN'))

borders <- dbGetQuery(conn, 'select * from TRANS_OE_borderpos')
smudge <- dbGetQuery(conn, 'select * from TR_OE_ONE_FS_GA_smudgebox')

##### GATHER DATA
data.fit.all <- NULL
for (a in unique(tr.conds$arms)) {
  for (c in tr.conds$conds[tr.conds$arms == a]) {
    temp.fit <- dbGetQuery(conn, sprintf('select a.*, b.density, b.plate_no, b.plate_col, b.plate_row
                                         from TR_OE_PAPER_FS_%s_%s_6144_FITNESS a, TRANS_OE_pos2coor b
                                         where a.pos = b.pos
                                         order by a.hours, b.plate_no, b.plate_col, b.plate_row', a, c))
    temp.fit$arm <- a
    temp.fit$condition <- c
    temp.fit$saturation <- max(temp.fit$hours)
    data.fit.all <- rbind(data.fit.all, temp.fit)
  }
}
data.fit.all$average[data.fit.all$condition == 'HO' & data.fit.all$plate_no %in% c(2,3,10,18,19,23,26,32)] <- NA
data.fit.all$average[data.fit.all$arm == 'ONE' & data.fit.all$plate_no == 14] <- NA
data.fit.all$average[data.fit.all$arm == 'ONE' & data.fit.all$pos %in% smudge$pos] <- NA

data.fit.all$fitness[data.fit.all$condition == 'HO' & data.fit.all$plate_no %in% c(2,3,10,18,19,23,26,32)] <- NA
data.fit.all$fitness[data.fit.all$arm == 'ONE' & data.fit.all$plate_no == 14] <- NA
data.fit.all$fitness[data.fit.all$arm == 'ONE' & data.fit.all$pos %in% smudge$pos] <- NA

data.fit.all <- merge(data.fit.all, orf_types, by = 'orf_name', all.x = T)
data.fit.all$category[data.fit.all$orf_name == 'BF_control'] <- 'Reference'
data.fit.all$category[data.fit.all$strain_id >= 1000000] <- 'ConditionControls'

data.fit.all$rep <- as.numeric(str_trunc(as.character(data.fit.all$pos), 5, side = 'left', ellipsis = ''))

data.fit.all$source[data.fit.all$plate_row%%2==1 & data.fit.all$plate_col%%2==1] = 'TL'
data.fit.all$source[data.fit.all$plate_row%%2==0 & data.fit.all$plate_col%%2==1] = 'BL'
data.fit.all$source[data.fit.all$plate_row%%2==1 & data.fit.all$plate_col%%2==0] = 'TR'
data.fit.all$source[data.fit.all$plate_row%%2==0 & data.fit.all$plate_col%%2==0] = 'BR'

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

##### REPLICATE COUNT PER REP
data.fit.cnt <- data.fit.all[!is.na(data.fit.all$fitness),] %>%
  group_by(condition, rep) %>% 
  count() %>%
  data.frame()
data.fit.all <- merge(data.fit.all, data.fit.cnt, by = c('condition','rep'), all.x = T)
data.fit.all$n[is.na(data.fit.all$n)] <- 0

##### REMOVE DATA WITH LESS THAN 4 REPLICATES
data.fit.all$fitness[data.fit.all$n < 4] <- NA

##### FITNESS SUMMARY
data.fit.sum <- data.fit.all %>%
  group_by(arm, condition, hours, rep, strain_id, orf_name, n, category) %>%
  summarize(avg.median = median(average, na.rm = T),
            fitness.median = median(fitness, na.rm = T),
            .groups = 'keep') %>%
  data.frame()
head(data.fit.sum)

data.fit.all <- merge(data.fit.all[,c(1:22,27)], data.fit.sum,
                      by = c('arm','condition','hours','rep','strain_id','orf_name','n','category'), all = T)

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

# write.csv(data.fit.sum[,c(1:10)], file = '/home/sbp29/R/Projects/adaptivefitness/output/translatome/pinning_artifact_correction/tr_oe_fitandcs_paper_hours_summary.csv')
# write.csv(data.fit.all[,c(1:16)], file = '/home/sbp29/R/Projects/adaptivefitness/output/translatome/pinning_artifact_correction/tr_oe_fitandcs_paper_hours_all.csv')

##### PLOTTING REP-WISE RESULTS
data.fit.sum %>%
  filter(category %in% c('Reference','Transient')) %>%
  ggplot(aes(x = rep, y = fitness.median)) +
  geom_point(aes(col = category), size = 0.5) +
  labs(x = 'Index',
       y = 'Fitness') +
  facet_wrap(.~condition) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_blank(),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  guides(color = guide_legend(nrow=1, byrow=TRUE, order = 1,
                              override.aes = list(size = 2)))


##### FITNESS DIFFERENTIAL
data.fit.diff <- merge(data.fit.sum[data.fit.sum$condition %notin% c('GA','DM'),
             c('arm','condition','hours','rep','strain_id','orf_name','n','category','avg.median','fitness.median')] %>%
               filter(n > 3),
      data.fit.sum[data.fit.sum$condition %in% c('GA','DM'),
                   c('arm','condition','hours','rep','n','avg.median','fitness.median')] %>%
        filter(n > 3),
      by = c('arm','rep'), suffixes = c('_stress','_nonstress'))
head(data.fit.diff)

data.fit.diff$avg.diff <- data.fit.diff$avg.median_stress - data.fit.diff$avg.median_nonstress
data.fit.diff$fitness.diff <- data.fit.diff$fitness.median_stress - data.fit.diff$fitness.median_nonstress

data.fit.diff %>%
  filter(category %in% c('Reference','Transient')) %>%
  ggplot(aes(x = rep, y = fitness.diff)) +
  geom_point(aes(col = category), size = 0.5) +
  labs(x = 'Index',
       y = 'Fitness') +
  facet_wrap(.~condition_stress) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_blank(),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  guides(color = guide_legend(nrow=1, byrow=TRUE, order = 1,
                              override.aes = list(size = 2)))

data.fit.diff %>%
  filter(fitness.diff > 0.05, category != 'ConditionControls') %>%
  group_by(condition_stress, category) %>%
  count() %>% data.frame()


##### DISTRIBUTION OF FITNESS EFFECTS OVER PLATEMAPS
p2s <- dbGetQuery(conn, 'select a.*, b.strain_id
                from TRANS_OE_pos2coor a, TRANS_OE_pos2strainid b
                where a.pos = b.pos')
p2s$rep <- as.numeric(str_trunc(as.character(p2s$pos), 5, side = 'left', ellipsis = ''))
p2s$pos <- as.numeric(p2s$pos)
# write.csv(p2s, file = '/home/sbp29/R/Projects/adaptivefitness/output/translatome/positionid2reps.csv')


data.fit.diff.p2s <- merge(p2s, data.fit.diff, by = c('rep','strain_id'))
head(data.fit.diff.p2s)

data.fit.diff.p2s %>%
  filter(density == 384, condition_stress == 'FL', plate_no < 25,
         fitness.diff > 0.05) %>%
  ggplot(aes(x = plate_col, y = plate_row, fill = fitness.diff)) +
  geom_tile(col = 'black') +
  labs(x = 'Columns', y = 'Rows') +
  scale_y_reverse() +
  facet_wrap(.~plate_no) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_blank(),
        legend.text = element_text(size = txt),
        legend.position = 'right',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm")))


##### HITS ACCORDING TO AARON
# tr.oe.fit.diff.hits <- read.table(file = '/home/acwach/overexp/oe_hits_fit_diff_fdr05.txt')
tr.oe.fit.diff.hits <- read.table(file = '/home/acwach/overexp/oe_fit_diff_fdr05_all')
tr.oe.fit.diff.hits$hit <- 'Yes'
data.fit.diff.p2s <- merge(data.fit.diff.p2s, tr.oe.fit.diff.hits[,c('strain_id','condition','hit')],
                      by.x = c('strain_id','condition_stress'),
                      by.y = c('strain_id','condition'),
                      all.x = T)
data.fit.diff.p2s$hit[is.na(data.fit.diff.p2s$hit)] <- 'No'

platemap.hits <- data.fit.diff.p2s %>%
  filter(density == 384, plate_no < 25) %>%
  ggplot(aes(x = plate_col, y = plate_row, fill = hit)) +
  geom_tile(col = 'black') +
  labs(x = 'Columns', y = 'Rows') +
  scale_y_reverse() +
  scale_fill_manual(name = 'Hits',
                    values = c('Yes' = 'red',
                                  'No' = 'white')) +
  facet_grid(condition_stress~plate_no) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        panel.grid = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'right',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%s/PLATEMAP_HITS.jpg",fig_path), platemap.hits,
       height = 100, width = 700, units = 'mm',
       dpi = 300)


##### ODD DM PATTERN
data.fit.all$pos <- as.numeric(data.fit.all$pos)
data.fit.odd <- data.fit.all %>% filter(condition == 'DM',
  pos %in% c(4454008166,19282008166,2482008166,33404008166,24704008166,12554008167,21254008167,33404008167,29354008167,6532008167,4454008167,16604008167,14632008167,27382008167,19282008167,8504008167,11182008167,31432008167,23332008167,2482008167,24704008167,4454008180,14632008180,8504008180,23332008180,29354008180,31432008180,11182008180,19282008180,6532008180,33404008180,12554008180,21254008180,24704008180,27382008180,2482008180,16604008180,16604008197,12554008197,11182008197,23332008197,4454008197,24704008197,14632008197,33404008197,8504008197,27382008197,31432008197,6532008197,29354008197,19282008197,2482008197,21254008197,4454008199,19282008199,33404008199,27382008199,11182008199,12554008199,2482008199,6532008199,8504008199,24704008199,16604008199,23332008199,14632008199,21254008199,29354008199,31432008199,8504008201,14632008201,12554008201,2482008201,11182008201,33404008201,23332008201,16604008201,29354008201,24704008201,6532008201,4454008201,31432008201,19282008201,21254008201,27382008201,4454008214,2482008214,12554008214,31432008214,16604008214,27382008214,33404008214,14632008214,19282008214,8504008214,6532008214,21254008214,24704008214,11182008214,23332008214,29354008214,4454008219,14632008219,11182008219,33404008219,29354008219,16604008219,31432008219,6532008219,27382008219,12554008219,24704008219,23332008219,19282008219,21254008219,8504008219,2482008219,24704008223,14632008223,6532008223,12554008223,11182008223,19282008223,21254008223,2482008223,8504008223,33404008223,4454008223,16604008223,31432008223,23332008223,29354008223,27382008223,21254008230,29354008230,27382008230,12554008230,19282008230,8504008230,33404008230,31432008230,4454008230,6532008230,2482008230,14632008230,11182008230,24704008230,23332008230,16604008230,2482008232,8504008232,6532008232,31432008232,12554008232,14632008232,16604008232,23332008232,29354008232,11182008232,24704008232,33404008232,19282008232,21254008232,4454008232,27382008232,2482008233,8504008233,11182008233,12554008233,21254008233,31432008233,4454008233,19282008233,29354008233,6532008233,24704008233,23332008233,14632008233,33404008233,27382008233,16604008233,29354008254,21254008254,2482008254,24704008254,27382008254,12554008254,19282008254,16604008254,4454008254,14632008254,8504008254,23332008254,33404008254,11182008254,31432008254,6532008254,21254008258,19282008258,27382008258,24704008258))
data.fit.odd

data.fit.odd %>%
  ggplot(aes(x = seq(1:201))) +
  geom_point(aes(y = fitness, col = as.character(rep))) +
  scale_color_discrete(guide = 'none')


data.fit.diff.p2s %>%
  filter(density == 384, condition_nonstress == 'DM', rep %in% unique(data.fit.odd$rep)) %>%
  ggplot(aes(x = plate_col, y = plate_row, fill = as.character(rep))) +
  geom_tile(col = 'black') +
  labs(x = 'Columns', y = 'Rows') +
  scale_fill_discrete(guide = 'none') +
  scale_y_reverse() +
  facet_wrap(.~plate_no) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_blank(),
        legend.text = element_text(size = txt),
        legend.position = 'right',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm")))


data.fit.odd2 <- data.fit.all %>% filter(condition == 'DM',
                                         rep %in% p2s$rep[p2s$density == 384 & p2s$plate_no == 22])
data.fit.odd2$odd[data.fit.odd2$rep %in% unique(data.fit.odd$rep)] <- 'Yes'
data.fit.odd2$odd[is.na(data.fit.odd2$odd)] <- 'No'

data.fit.odd2 %>%
  ggplot(aes(x = seq(1:4640), y = fitness)) +
  geom_point(aes(col = odd)) +
  scale_color_discrete(guide = 'none')

data.fit.odd2$pospos <- sample(4640, 4640)

data.fit.odd3 <- data.fit.all %>% filter(condition == 'DM')
data.fit.odd3$pospos <- sample(160104, 160104)

data.fit.odd3 %>%
  ggplot(aes(x = as.factor(pos), y = fitness)) +
  geom_point() +
  scale_color_discrete(guide = 'none') +
  theme(axis.text = element_blank())

unique(data.fit.odd2$orf_name)

##### AARON'S GROWTH ANALYSIS
data.growth <- read.table(file = '/home/acwach/overexp/colony_growth_oe_full', sep = ' ', header = T)

temp.growth <- NULL
for (c in unique(data.growth$condition)) {
  temp <- data.growth[data.growth$condition == c,]
  temp <- merge(p2s[p2s$density == 6144,], temp[,c('condition','position','growth')], by.x = 'pos', by.y = 'position', all.x = T)
  temp$condition[is.na(temp$condition)] <- c
  temp.growth <- rbind(temp, temp.growth)
}
data.growth <- temp.growth

data.growth$hours[data.growth$condition == 'DM'] <- 1
data.growth$hours[data.growth$condition == 'FL'] <- 2
data.growth$hours[data.growth$condition == 'GA'] <- 3
data.growth$hours[data.growth$condition == 'HU'] <- 4
data.growth$hours[data.growth$condition == 'SA'] <- 5
data.growth$hours[data.growth$condition == 'TN'] <- 6
data.growth$hours[data.growth$condition == 'HO'] <- 7

data.growth <- data.growth[,c('pos','hours','growth')]
colnames(data.growth) <- c('pos','hours','average')

dbWriteTable(conn, 'TR_OE_GROWTH_FS_ALL_6144_CLEAN', data.growth, overwrite = T)

##### POST LID ON GROWTH DATA