##### OVEREXPRESSION SCREEN HITS VALIDATION
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 01/30/2022 

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
library(locfit)
library(growthcurver)
library(rstatix)
library(gtools)
library(growthrates)
library(readxl)

`%notin%` <- Negate(`%in%`)


##### OVEREXPRESSION HITS
platemap <- read_xlsx('/home/sbp29/R/Projects/adaptivefitness/rawdata/translatome/validation/TR_HIT_VAL_PR_PLATEMAP.xlsx') %>% data.frame()
res.biotek <- read_xlsx('/home/sbp29/R/Projects/adaptivefitness/rawdata/translatome/validation/TR_HIT_VAL_PR_RESULTS.xlsx') %>% data.frame()
data.val1 <- read_xlsx('/home/sbp29/R/Projects/adaptivefitness/rawdata/translatome/validation/TR_HIT_VAL_PR_VALUES_1.xlsx') %>% data.frame()
data.val2 <- read_xlsx('/home/sbp29/R/Projects/adaptivefitness/rawdata/translatome/validation/TR_HIT_VAL_PR_VALUES_2.xlsx') %>% data.frame()
data.val3 <- read_xlsx('/home/sbp29/R/Projects/adaptivefitness/rawdata/translatome/validation/TR_HIT_VAL_PR_VALUES_3.xlsx') %>% data.frame()

head(data.val1)
data.val <- merge(data.frame(rbind(cbind(data.val1 %>% melt(id.vars = 'Time', variable.name = 'Sample', value.name = 'OD'), val_id = 'PR1'),
                 cbind(data.val2 %>% melt(id.vars = 'Time', variable.name = 'Sample', value.name = 'OD'), val_id = 'PR2'),
                 cbind(data.val3 %>% melt(id.vars = 'Time', variable.name = 'Sample', value.name = 'OD'), val_id = 'PR3'))),
      platemap,
      by.x = c('val_id','Sample'), by.y = c('val_id','pos_id')) %>%
  filter(orf_name %notin% c('BLANK','REAL BLANK'))
write.csv(data.val, file = '/home/sbp29/R/Projects/adaptivefitness/output/translatome/tr_oe_val_data.csv')

data.val %>%
  ggplot(aes(x = Time, y = OD)) +
  geom_line(aes(col = condition, group = Sample)) +
  facet_wrap(.~orf_name * val_id) +
  coord_cartesian(xlim = c(0,2000))

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 8
txt <- 8
lbls <- 9

##### GROWTH CHARACTERISTICS
res.val1 <- SummarizeGrowthByPlate(data.val1[1:(dim(data.val1)[1]-1),1:(dim(data.val1)[2]-1)])
res.val1$val_id <- 'PR1'
res.val1 <- merge(res.val1, platemap, by.x = c('sample','val_id'), by.y = c('pos_id','val_id'))
res.val2 <- SummarizeGrowthByPlate(data.val2[1:(dim(data.val2)[2]-1)])
res.val2$val_id <- 'PR2'
res.val2 <- merge(res.val2, platemap, by.x = c('sample','val_id'), by.y = c('pos_id','val_id'))
res.val3 <- SummarizeGrowthByPlate(data.val3[1:(dim(data.val3)[2]-1)])
res.val3$val_id <- 'PR3'
res.val3 <- merge(res.val3, platemap, by.x = c('sample','val_id'), by.y = c('pos_id','val_id'))

res.val <- rbind(res.val1, res.val2, res.val3)

res.val %>%
  ggplot(aes(x = orf_name, y = auc_e)) +
  geom_point() +
  facet_wrap(.~condition*val_id) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 90),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt,
                                              margin = margin(0.1,0,0.1,0, "mm"))) #+
  #coord_cartesian(ylim = c(0,0.04))
write.csv(res.val %>% filter(orf_name %notin% c('BLANK','REAL BLANK')), 
                             file = '/home/sbp29/R/Projects/adaptivefitness/output/translatome/tr_oe_val_results.csv')

merge(res.biotek, platemap, by = c('pos_id','val_id')) %>%
  ggplot(aes(x = orf_name, y = gr)) +
  geom_point() +
  facet_wrap(.~condition*val_id) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 90),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt,
                                              margin = margin(0.1,0,0.1,0, "mm")))


##### DELETION HITS
platemap.del <- read_xlsx('/home/sbp29/R/Projects/adaptivefitness/rawdata/translatome/validation/TR_DEL_HIT_PLATEMAP.xlsx') %>% data.frame()
data.del.val1 <- read_xlsx('/home/sbp29/R/Projects/adaptivefitness/rawdata/translatome/validation/TR_DEL_HIT_PR1_VALUES.xlsx') %>% data.frame()
data.del.val2 <- read_xlsx('/home/sbp29/R/Projects/adaptivefitness/rawdata/translatome/validation/TR_DEL_HIT_PR2_VALUES.xlsx') %>% data.frame()

data.del.val <- merge(data.frame(rbind(cbind(data.del.val1 %>% melt(id.vars = 'Time', variable.name = 'Sample', value.name = 'OD'), val_id = 'PR1'),
                                   cbind(data.del.val2 %>% melt(id.vars = 'Time', variable.name = 'Sample', value.name = 'OD'), val_id = 'PR2'))),
                  platemap.del,
                  by.x = c('val_id','Sample'), by.y = c('val_id','pos_id')) %>%
  filter(orf_name %notin% c('BLANK','REAL BLANK'),
         Time <= 2500)

write.csv(data.del.val, file = '/home/sbp29/R/Projects/adaptivefitness/output/translatome/tr_del_val_data.csv')

data.del.val %>%
  ggplot(aes(x = Time, y = OD)) +
  geom_line(aes(col = orf_name, group = Sample)) +
  # stat_summary(aes(col = orf_name),
  #              geom = 'line', fun = 'mean') +
  # stat_summary(aes(fill = orf_name, group = Sample),
  #              geom = 'ribbon', fun.data = 'mean_sdl') +
  facet_wrap(.~condition)

##### GROWTH CHARACTERISTICS
res.del.val1 <- SummarizeGrowthByPlate(data.del.val1)
res.del.val1$val_id <- 'PR1'
res.del.val1 <- merge(res.del.val1, platemap.del, by.x = c('sample','val_id'), by.y = c('pos_id','val_id'))
res.del.val2 <- SummarizeGrowthByPlate(data.del.val2)
res.del.val2$val_id <- 'PR2'
res.del.val2 <- merge(res.del.val2, platemap.del, by.x = c('sample','val_id'), by.y = c('pos_id','val_id'))

res.del.val <- rbind(res.del.val1, res.del.val2)

res.del.val %>%
  filter(condition != 'YPDA', orf_name != 'BLANK') %>%
  ggplot(aes(x = orf_name, y = r)) +
  geom_point() +
  stat_summary(geom = 'point', fun = 'mean', col = 'red', size = 1) +
  # stat_summary(geom = 'errorbar', fun.data = 'mean_sdl', col = 'red', size = 0.5) +
  facet_wrap(.~condition*val_id) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 90),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt,
                                              margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0,0.015))

write.csv(res.val %>% filter(orf_name %notin% c('BLANK','REAL BLANK')), 
          file = '/home/sbp29/R/Projects/adaptivefitness/output/translatome/tr_oe_val_results.csv')

merge(res.biotek, platemap, by = c('pos_id','val_id')) %>%
  ggplot(aes(x = orf_name, y = gr)) +
  geom_point() +
  facet_wrap(.~condition*val_id) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 90),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt,
                                              margin = margin(0.1,0,0.1,0, "mm")))