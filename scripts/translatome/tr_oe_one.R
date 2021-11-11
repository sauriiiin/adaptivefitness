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

expt.name <- "tr_oe_one"

source("~/R/Projects/adaptivefitness/R/functions/isoutlier.R")
source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

# load('output/translatome/AllData.RData')

controls <- read.csv(file = 'rawdata/translatome/condition_controls.csv', stringsAsFactors = F)

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 8
txt <- 8
lbls <- 9

##### GATHER DATA
tr.conds <- data.frame(arms = c('ONE','ONE','ONE','ONE','TWO','TWO','TWO','TWO'),
                       conds = c('GA','SA','HU','HO','GA','DM','TN','FL'))

# pgs <- dbGetQuery(conn, 'select orf_name from PROTOGENES where pg_2012 = 1')
pgs <- read.csv(file = '~/R/Projects/adaptivefitness/rawdata/translatome/oe_transient_translated.csv', stringsAsFactors = F)
tr_status <- read.csv(file = '~/R/Projects/adaptivefitness/rawdata/translatome/overexpression_collection_translation_status.csv',
                      stringsAsFactors = F)
colnames(tr_status) <- c('orf_name','is_tr','on_list')

data.fit <- NULL
# data.stats <- NULL
for (a in unique(tr.conds$arms)) {
  for (c in tr.conds$conds[tr.conds$arms == a]) {
    temp.fit <- dbGetQuery(conn, sprintf('select a.*, b.orf_name, c.strain_id, d.density, d.plate, d.col, d.row
                                         from TR_OE_%s_FS_%s_6144_NORM a, TR_OE_pos2orf_name b, TR_OE_pos2strainid c,
                                         TR_OE_pos2coor d
                                         where a.pos = b.pos and b.pos = c.pos and c.pos = d.pos
                                         order by a.hours, d.plate, d.col, d.row', a, c))
    temp.fit$arm <- a
    temp.fit$condition <- c
    temp.fit$saturation <- max(temp.fit$hours)
    data.fit <- rbind(data.fit, temp.fit)
    
  }
}
data.fit$orf_type[data.fit$orf_name %in% pgs$orf_name[pgs$translated == 0]] <- 'Not Translated'
data.fit$orf_type[data.fit$orf_name %in% pgs$orf_name[pgs$is_transient == 1]] <- 'Transient'
data.fit$orf_type[data.fit$orf_name == 'BF_control'] <- 'Reference'
data.fit$orf_type[is.na(data.fit$orf_type)] <- 'Not Transient'
data.fit <- merge(data.fit, tr_status, by = 'orf_name', all.x = T)
data.fit$rep <- as.numeric(str_trunc(as.character(data.fit$pos), 5, side = 'left', ellipsis = ''))

data.fit$saturation[data.fit$condition == 'FL'] <- 46
# data.stats$saturation[data.stats$condition == 'FL'] <- 46
time.sat <- data.fit %>%
  group_by(arm, condition) %>%
  summarize(saturation = max(saturation, na.rm = T), .groups = 'keep') %>%
  data.frame()

##### FITNESS PLATEMAPS
data.fit %>%
  filter(condition == 'TN', hours == saturation) %>%
  ggplot(aes(x = col, y = row, fill = fitness)) +
  geom_tile() +
  scale_y_reverse() +
  facet_wrap(.~plate)

## CLEANING DATA AFTER LOOKING AT FITNESS PLATEMAPS
data.fit$fitness[data.fit$condition == 'HO' &
                   data.fit$plate %in% c(2,3,10,18,19,23,26,32)] <- NA

##### REMOVE OUTLIERS
data.mad <- data.fit %>%
  group_by(arm, condition, hours, rep, orf_name, strain_id, orf_type, is_tr, on_list) %>%
  summarize(fitness.median = median(fitness, na.rm = T), cs.median = median(average, na.rm = T),
            fitness.mad = mad(fitness, na.rm = T), cs.mad = mad(average, na.rm = T),
            .groups = 'keep') %>%
  data.frame()

data.fit <- merge(data.fit, data.mad, by = c('arm','condition','hours','rep','orf_name','strain_id','orf_type','is_tr','on_list'))
data.fit$fitness[data.fit$fitness < (data.fit$fitness.median - 2*data.fit$fitness.mad) |
                   data.fit$fitness > (data.fit$fitness.median + 2*data.fit$fitness.mad)] <- NA
data.fit$average[data.fit$average < (data.fit$cs.median - 2*data.fit$cs.mad) |
                   data.fit$average > (data.fit$cs.median + 2*data.fit$cs.mad)] <- NA

data.mad <- merge(data.mad, time.sat, by = c('arm','condition'))

##### REFERENCE LIMITS
data.lim <- data.fit %>%
  filter(orf_name == 'BF_control') %>%
  group_by(arm, condition, hours, orf_name, rep) %>%
  summarize(average = median(average, na.rm = T),
            fitness = median(fitness, na.rm = T),
            .groups = 'keep') %>%
  group_by(arm, condition, hours, orf_name) %>%
  summarize(cs_ll = quantile(average, 0.025, na.rm = T),
            cs_m = median(average, na.rm = T),
            cs_ul = quantile(average, 0.975, na.rm = T),
            fitness_ll = quantile(fitness, 0.025, na.rm = T),
            fitness_m = median(fitness, na.rm = T),
            fitness_ul = quantile(fitness, 0.975, na.rm = T),
            .groups = 'keep') %>%
  data.frame()

# data.stats <- merge(data.stats, data.lim[,-4], by = c('arm','condition','hours'))
# data.stats$es <- abs(data.stats$cs_median - data.stats$fitness_m)/data.stats$fitness_m

## YBR196C-A
# write.csv(merge(data.mad %>%
#         filter(orf_name == 'YBR196C-A', hours == saturation), data.lim[,-4], by = c('arm','condition','hours')),
#       file = 'output/translatome/ybr_fitness_data.csv')

##### DIFFERENTIAL FITNESS ANALYSIS
data.ref.dist <- data.fit %>%
  filter(orf_name == 'BF_control', hours == saturation) %>%
  group_by(arm, condition, rep) %>%
  summarize(average = median(average, na.rm = T),
            fitness = median(fitness, na.rm = T),
            .groups = 'keep') %>%
  data.frame()
data.ref.dist$id <- paste(data.ref.dist$arm, data.ref.dist$condition, sep = '_')

data.mad <- merge(data.mad, data.lim[,-4], by = c('arm','condition','hours'))
data.mad$id <- paste(data.mad$arm, data.mad$condition, sep = '_')
data.diff <- NULL
diff.dist <- NULL
for (id1 in unique(data.mad$id)) {
  for (id2 in unique(data.mad$id[data.mad$id != id1])) {
    temp.data.diff <- merge(data.mad[data.mad$id == id1 & data.mad$hours == data.mad$saturation,-3],
                       data.mad[data.mad$id == id2 & data.mad$hours == data.mad$saturation,c(-1,-2,-3)],
                       by = c('rep','orf_name','strain_id','orf_type','is_tr','on_list'), suffixes = c('_cond','_ref'), all = T)
    temp.data.diff$fitness_diff <- temp.data.diff$fitness.median_cond - temp.data.diff$fitness.median_ref
    temp.data.diff$cs_diff <- temp.data.diff$cs.median_cond - temp.data.diff$cs.median_ref

    temp.diff.dist <- NULL
    temp1 <- data.ref.dist$fitness[data.ref.dist$id == id1]
    temp2 <- data.ref.dist$fitness[data.ref.dist$id == id2]
    temp3 <- NULL
    temp4 <- data.ref.dist$average[data.ref.dist$id == id1]
    temp5 <- data.ref.dist$average[data.ref.dist$id == id2]
    temp6 <- NULL
    for (i in 1:50000) {
      temp3 <- c(temp3, sample(temp1[!is.na(temp1)],1) - sample(temp2[!is.na(temp2)],1))
      temp6 <- c(temp6, sample(temp4[!is.na(temp4)],1) - sample(temp5[!is.na(temp5)],1))
    }
    temp.diff.dist <- rbind(temp.diff.dist, cbind(id_cond = id1, id_ref = id2,
                                        ll = quantile(temp3, 0.025)[[1]],
                                        m = quantile(temp3, 0.5)[[1]],
                                        ul = quantile(temp3, 0.975)[[1]],
                                        cs_ll = quantile(temp6, 0.025)[[1]],
                                        cs_m = quantile(temp6, 0.5)[[1]],
                                        cs_ul = quantile(temp6, 0.975)[[1]]))
    temp.diff.dist <- data.frame(temp.diff.dist)
    temp.diff.dist$ll <- as.numeric(as.character(temp.diff.dist$ll))
    temp.diff.dist$m <- as.numeric(as.character(temp.diff.dist$m))
    temp.diff.dist$ul <- as.numeric(as.character(temp.diff.dist$ul))
    temp.diff.dist$cs_ll <- as.numeric(as.character(temp.diff.dist$cs_ll))
    temp.diff.dist$cs_m <- as.numeric(as.character(temp.diff.dist$cs_m))
    temp.diff.dist$cs_ul <- as.numeric(as.character(temp.diff.dist$cs_ul))

    data.diff <- rbind(data.diff, temp.data.diff)
    diff.dist <- rbind(diff.dist, temp.diff.dist)
  }
}
data.diff <- merge(data.diff, diff.dist, by = c('id_cond','id_ref'), all = T)

data.diff$phenotype[data.diff$fitness_diff >= data.diff$ul] <- 'Beneficial'
data.diff$phenotype[data.diff$fitness_diff <= data.diff$ll] <- 'Deleterious'
data.diff$phenotype[is.na(data.diff$phenotype)] <- 'Neutral'

data.diff$phenotype_cs[data.diff$cs_diff >= data.diff$cs_ul] <- 'Beneficial'
data.diff$phenotype_cs[data.diff$cs_diff <= data.diff$cs_ll] <- 'Deleterious'
data.diff$phenotype_cs[is.na(data.diff$phenotype_cs)] <- 'Neutral'

temp1 <- str_split(data.diff$id_cond, '_', simplify = T)
colnames(temp1) <- c('arm1','cond1')
temp2 <- str_split(data.diff$id_ref, '_', simplify = T)
colnames(temp2) <- c('arm2','cond2')
data.diff <- cbind(data.diff, temp1, temp2)

diff.dist2 <- NULL
for (id1 in unique(diff.dist$id_cond)) {
  for (id2 in unique(diff.dist$id_ref[diff.dist$id_cond == id1])) {
    diff.dist2 <- rbind(diff.dist2, data.frame(id_cond = id1, id_ref = id2,
                                               x = seq(0,10,0.1),
                                               y.ll = seq(0,10,0.1) - diff.dist$ll[diff.dist$id_cond == id1 & diff.dist$id_ref == id2],
                                               y.m = seq(0,10,0.1) - diff.dist$m[diff.dist$id_cond == id1 & diff.dist$id_ref == id2],
                                               y.ul = seq(0,10,0.1) - diff.dist$ul[diff.dist$id_cond == id1 & diff.dist$id_ref == id2],
                                               x.cs = seq(-5000,5000,100),
                                               y.cs_ll = seq(-5000,5000,100) - diff.dist$cs_ll[diff.dist$id_cond == id1 & diff.dist$id_ref == id2],
                                               y.cs_m = seq(-5000,5000,100) - diff.dist$cs_m[diff.dist$id_cond == id1 & diff.dist$id_ref == id2],
                                               y.cs_ul = seq(-5000,5000,100) - diff.dist$cs_ul[diff.dist$id_cond == id1 & diff.dist$id_ref == id2]))
  }
}
diff.dist2 <- data.frame(diff.dist2)
temp1 <- str_split(diff.dist2$id_cond, '_', simplify = T)
colnames(temp1) <- c('arm1','cond1')
temp2 <- str_split(diff.dist2$id_ref, '_', simplify = T)
colnames(temp2) <- c('arm2','cond2')
diff.dist2 <- cbind(diff.dist2, temp1, temp2)

fig.diff <- data.diff[!(data.diff$strain_id %in% controls$strain_id),] %>%
  filter(arm1 == arm2, id_ref %in% c('ONE_GA','TWO_DM'), id_cond != 'TWO_GA') %>%
  filter(!is.na(id_cond), !is.na(id_ref), orf_type != 'Reference') %>%
  ggplot(aes(x = fitness.median_cond, y = fitness.median_ref)) +
  geom_point(aes(col = phenotype), size = 1) +
  geom_vline(aes(xintercept = fitness_ll_cond), linetype = 'dashed', size = 0.25, col = 'red') +
  geom_vline(aes(xintercept = fitness_ul_cond), linetype = 'dashed', size = 0.25, col = 'red') +
  geom_hline(aes(yintercept = fitness_ll_ref), linetype = 'dashed', size = 0.25, col = 'red') +
  geom_hline(aes(yintercept = fitness_ul_ref), linetype = 'dashed', size = 0.25, col = 'red') +
  geom_line(data = diff.dist2 %>% filter(arm1 == arm2, id_ref %in% c('ONE_GA','TWO_DM'), id_cond != 'TWO_GA'),
            aes(x = x , y = y.ul),
            linetype = 'dashed', size = 0.5) +
  geom_line(data = diff.dist2 %>% filter(arm1 == arm2, id_ref %in% c('ONE_GA','TWO_DM'), id_cond != 'TWO_GA'),
            aes(x = x , y = y.ll),
            linetype = 'dashed', size = 0.5) +
  geom_text_repel(data = data.diff %>%
                    filter(arm1 == arm2, id_ref %in% c('ONE_GA','TWO_DM'), id_cond != 'TWO_GA') %>%
                    filter(!is.na(id_cond), !is.na(id_ref), orf_name == 'YBR196C-A'),
                  aes(x = fitness.median_cond, y = fitness.median_ref, label = orf_name), size = 2,
                  force = 2, max.overlaps = 30) +
  geom_point(data = data.diff %>%
               filter(arm1 == arm2, id_ref %in% c('ONE_GA','TWO_DM'), id_cond != 'TWO_GA') %>%
               filter(!is.na(id_cond), !is.na(id_ref), orf_name == 'YBR196C-A'),
             aes(x = fitness.median_cond, y = fitness.median_ref, fill = phenotype),
             col = 'black', shape = 21, size = 1) +
  scale_color_manual(name = 'Differential Phenotype',
                     values = c('Beneficial' = '#FFC107',
                                'Deleterious' = '#3F51B5',
                                'Neutral' = '#9E9E9E'),
                     limits = c('Deleterious','Neutral','Beneficial')) +
  scale_fill_manual(name = 'Differential Phenotype',
                     values = c('Beneficial' = '#FFC107',
                                'Deleterious' = '#3F51B5',
                                'Neutral' = '#9E9E9E'),
                     limits = c('Deleterious','Neutral','Beneficial'),
                    guide = F) +
  labs(x = 'Relative Fitness in Stress Condition',
       y = 'Relative Fitness in Reference Condition') +
  coord_cartesian(xlim = c(0, 2),
                  ylim = c(0, 2)) +
  facet_wrap(.~id_ref * id_cond, nrow = 2,
             labeller = labeller(id_ref = c('ONE_GA' = 'Ref.: Galactose',
                                            'TWO_GA' = 'Ref.: Galactose',
                                            'TWO_DM' = 'Ref.: DMSO'),
                                 id_cond = c('ONE_HO' = 'Str.: Hydrogen Peroxide',
                                             'ONE_HU' = 'Str.: Hydroxyurea',
                                             'ONE_SA' = 'Str.: Salt',
                                             'TWO_GA' = 'Str.: Galactose',
                                             'TWO_DM' = 'Str.: DMSO',
                                             'TWO_FL' = 'Str.: Fluconazole',
                                             'TWO_TN' = 'Str.: Tunicamycin'))) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  guides(color = guide_legend(nrow=1, byrow=TRUE, order = 1,
                              override.aes = list(size = 2)))
ggsave(sprintf("%s/AllDifferentialFitness.jpg",fig_path), fig.diff,
       height = one.5c, width = two.c, units = 'mm',
       dpi = 600)

fig.cs <- data.diff[!(data.diff$strain_id %in% controls$strain_id),] %>%
  filter(arm1 == arm2, id_ref %in% c('ONE_GA','TWO_DM'), id_cond != 'TWO_GA') %>%
  filter(!is.na(id_cond), !is.na(id_ref), orf_type != 'Reference') %>%
  ggplot(aes(x = cs.median_cond, y = cs.median_ref)) +
  geom_point(aes(col = phenotype), size = 1) +
  geom_vline(aes(xintercept = cs_ll_cond), linetype = 'dashed', size = 0.25, col = 'red') +
  geom_vline(aes(xintercept = cs_ul_cond), linetype = 'dashed', size = 0.25, col = 'red') +
  geom_hline(aes(yintercept = cs_ll_ref), linetype = 'dashed', size = 0.25, col = 'red') +
  geom_hline(aes(yintercept = cs_ul_ref), linetype = 'dashed', size = 0.25, col = 'red') +
  geom_line(data = diff.dist2 %>% filter(arm1 == arm2, id_ref %in% c('ONE_GA','TWO_DM'), id_cond != 'TWO_GA'),
            aes(x = x.cs , y = y.cs_ul),
            linetype = 'dashed', size = 0.5) +
  geom_line(data = diff.dist2 %>% filter(arm1 == arm2, id_ref %in% c('ONE_GA','TWO_DM'), id_cond != 'TWO_GA'),
            aes(x = x.cs , y = y.cs_ll),
            linetype = 'dashed', size = 0.5) +
  geom_text_repel(data = data.diff %>%
                    filter(arm1 == arm2, id_ref %in% c('ONE_GA','TWO_DM'), id_cond != 'TWO_GA') %>%
                    filter(!is.na(id_cond), !is.na(id_ref), orf_name == 'YBR196C-A'),
                  aes(x = cs.median_cond, y = cs.median_ref, label = orf_name), size = 2,
                  force = 2, max.overlaps = 30) +
  geom_point(data = data.diff %>%
               filter(arm1 == arm2, id_ref %in% c('ONE_GA','TWO_DM'), id_cond != 'TWO_GA') %>%
               filter(!is.na(id_cond), !is.na(id_ref), orf_name == 'YBR196C-A'),
             aes(x = cs.median_cond, y = cs.median_ref, fill = phenotype),
             col = 'black', shape = 21, size = 1) +
  scale_color_manual(name = 'Differential Phenotype',
                     values = c('Beneficial' = '#FFC107',
                                'Deleterious' = '#3F51B5',
                                'Neutral' = '#9E9E9E'),
                     limits = c('Deleterious','Neutral','Beneficial')) +
  scale_fill_manual(name = 'Differential Phenotype',
                    values = c('Beneficial' = '#FFC107',
                               'Deleterious' = '#3F51B5',
                               'Neutral' = '#9E9E9E'),
                    limits = c('Deleterious','Neutral','Beneficial'),
                    guide = F) +
  labs(x = 'Colony Size in Stress Condition (pixels)',
       y = 'Colony Size in Reference Condition (pixels)') +
  coord_cartesian(xlim = c(0, 2500),
                  ylim = c(0, 2500)) +
  facet_wrap(.~id_ref * id_cond,
             labeller = labeller(id_ref = c('ONE_GA' = 'Ref.: Galactose',
                                            'TWO_GA' = 'Ref.: Galactose',
                                            'TWO_DM' = 'Ref.: DMSO'),
                                 id_cond = c('ONE_HO' = 'Str.: Hydrogen Peroxide',
                                             'ONE_HU' = 'Str.: Hydroxyurea',
                                             'ONE_SA' = 'Str.: Salt',
                                             'TWO_GA' = 'Str.: Galactose',
                                             'TWO_DM' = 'Str.: DMSO',
                                             'TWO_FL' = 'Str.: Fluconazole',
                                             'TWO_TN' = 'Str.: Tunicamycin'))) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  guides(color = guide_legend(nrow=1, byrow=TRUE, order = 1,
                              override.aes = list(size = 2)))
ggsave(sprintf("%s/AllColonySize.jpg",fig_path), fig.cs,
       height = one.5c, width = two.c, units = 'mm',
       dpi = 600)

data.mad$id <- factor(data.mad$id, levels = c('ONE_GA','ONE_HO','ONE_HU','ONE_SA',
                                              'TWO_GA','TWO_DM','TWO_FL','TWO_TN'))
fig.csvsfit <- data.mad[!(data.mad$strain_id %in% controls$strain_id),] %>%
  filter(orf_type != 'Reference', hours == saturation) %>%
  ggplot(aes(x = fitness.median, y = cs.median)) +
  geom_point(size = 1, col = '#9E9E9E') +
  geom_vline(aes(xintercept = fitness_ll), linetype = 'dashed', size = 0.25, col = 'red') +
  geom_vline(aes(xintercept = fitness_ul), linetype = 'dashed', size = 0.25, col = 'red') +
  geom_hline(aes(yintercept = cs_ll), linetype = 'dashed', size = 0.25, col = 'red') +
  geom_hline(aes(yintercept = cs_ul), linetype = 'dashed', size = 0.25, col = 'red') +
  geom_text_repel(data = data.mad %>% filter(orf_name == 'YBR196C-A', hours == saturation),
                  aes(x = fitness.median, y = cs.median, label = orf_name), size = 2,
                  force = 2, max.overlaps = 30) +
  geom_point(data = data.mad %>% filter(orf_name == 'YBR196C-A', hours == saturation),
             aes(x = fitness.median, y = cs.median),
             col = 'black', shape = 21, size = 1) +
  stat_cor(method = 'pearson', label.x = 0, label.y = 2200,
           # aes(label = paste('Pearson',..r.label.., sep = "~`,`~")),
           size = 2) +
  labs(x = 'Relative Fitness in Stress Condition',
       y = 'Colony Size in Stress Condition (pixels)') +
  coord_cartesian(xlim = c(0, 2),
                  ylim = c(0, 2500)) +
  facet_wrap(.~id, nrow = 2,
             labeller = labeller(id = c('ONE_GA' = 'Galactose (Arm 1)',
                                        'ONE_HO' = 'Hydrogen Peroxide',
                                        'ONE_HU' = 'Hydroxyurea',
                                        'ONE_SA' = 'Salt',
                                        'TWO_GA' = 'Galactose (Arm 2)',
                                        'TWO_DM' = 'DMSO',
                                        'TWO_FL' = 'Fluconazole',
                                        'TWO_TN' = 'Tunicamycin'))) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  guides(color = guide_legend(nrow=1, byrow=TRUE, order = 1,
                              override.aes = list(size = 2)))
ggsave(sprintf("%s/AllColonySizeVSFitness.jpg",fig_path), fig.csvsfit,
       height = one.5c, width = two.c, units = 'mm',
       dpi = 600)


fig.gal.corr <- data.diff %>%
  filter(id_ref == 'TWO_GA', id_cond == 'ONE_GA') %>%
  filter(!is.na(id_cond), !is.na(id_ref), orf_type != 'Reference') %>%
  ggplot(aes(x = fitness.median_cond, y = fitness.median_ref)) +
  geom_point(col = '#9E9E9E', size = 1) +
  geom_smooth(method = 'lm', size = 0.5, linetype = 'dashed', col = 'black') +
  stat_cor(method = 'pearson', size = 3) +
  labs(x = 'Relative Fitness in Arm 1',
       y = 'Relative Fitness in Arm 2') +
  coord_cartesian(xlim = c(0, 1.5),
                  ylim = c(0, 1.5)) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%s/GALCorr.jpg",fig_path), fig.gal.corr,
       height = one.c, width = one.c, units = 'mm',
       dpi = 600)

##### CONDITION CONTROLS
# merge(data.diff[data.diff$strain_id %in% controls$strain_id &
#                   data.diff$arm1 == data.diff$arm2 &
#                   !is.na(data.diff$id_cond) & !is.na(data.diff$id_ref) &
#                   data.diff$id_ref %in% c('ONE_GA','TWO_GA') &
#                   data.diff$id_cond != 'TWO_GA',],
#       melt(controls, id.vars = c('strain_id','standard_name','orf_name'),
#            variable.name = 'cond1', value.name = 'control_type'),
#       by = c('cond1','strain_id','orf_name'))

fig.conts <- merge(data.mad, melt(controls, id.vars = c('strain_id','standard_name','orf_name'),
                     variable.name = 'condition', value.name = 'control_type'), by = c('condition','strain_id','orf_name')) %>%
  filter(hours %in% c(141, 36, 20, 125, 26), control_type != '') %>%
  ggplot(aes(x = fitness.median, y = orf_name)) +
  geom_boxplot(aes(fill = control_type, col = control_type), outlier.shape = NA, size = 0.6) +
  labs(x = 'Fitness',
       y = 'Control Mutant') +
  scale_color_manual(name = 'Control Type',
                       values = c('Resistant' = '#FFA000',
                                  'Sensitive' = '#303F9F')) +
  scale_fill_manual(name = 'Control Type',
                     values = c('Resistant' = '#FFA000',
                                'Sensitive' = '#303F9F')) +
  # geom_jitter(size = 0.2) +
  facet_wrap(.~condition, nrow = 1) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%s/CONDITIONCONTROLS.jpg",fig_path), fig.conts,
       height = one.c, width = two.c, units = 'mm',
       dpi = 600)


##### GO/KEGG ENRICHMENT FOR THE DELETERIOUS AND BENEFICIAL MUTANTS
data.diff2 <- data.diff[!(data.diff$strain_id %in% controls$strain_id) &
                          data.diff$arm1 == data.diff$arm2 &
                          !is.na(data.diff$id_cond) & !is.na(data.diff$id_ref) &
                          data.diff$id_ref %in% c('ONE_GA','TWO_GA','TWO_DM'),]

goe <- data.frame()
kegg <- data.frame()
for (id_ref in c('ONE_GA','TWO_GA','TWO_DM')) {
  allgenes <- unique(data.diff2$orf_name[data.diff2$id_ref == id_ref])
  allgenes <- bitr(allgenes, fromType = "ORF",
                   toType = c("ENTREZID","GENENAME","ENSEMBL"),
                   OrgDb = org.Sc.sgd.db)
  allgenes <- allgenes[!is.na(allgenes$ENSEMBL),]
  for (id_cond in unique(data.diff2$id_cond[data.diff2$id_ref == id_ref])) {
    for (p in unique(data.diff2$phenotype[data.diff2$id_cond == id_cond & data.diff2$id_ref == id_ref])) {
      temp.deg <- data.diff2$orf_name[data.diff2$id_ref == id_ref & data.diff2$id_cond == id_cond &
                                       data.diff2$phenotype == p]
      temp.deg <- bitr(temp.deg, fromType = "ORF",
                       toType = c("ENTREZID","GENENAME","ENSEMBL"),
                       OrgDb = org.Sc.sgd.db)
      temp.deg <- temp.deg[!is.na(temp.deg$ENSEMBL),]

      temp.goe <- enrichGO(gene          = temp.deg$ENSEMBL,
                           universe      = allgenes$ENSEMBL,
                           OrgDb         = org.Sc.sgd.db,
                           keyType       = "ENSEMBL",
                           ont           = "ALL",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.01,
                           qvalueCutoff  = 0.05)
      if (length(temp.goe) == 1) {
        if (dim(temp.goe)[1] == 0) {
          cat(sprintf('There are no GO term enrichment for %s ORFs in REF: %s | COND: %s comparison.\n',
                      p,id_ref,id_cond))
        } else{
          cat(sprintf('GO term enrichment for %s ORFs in REF: %s | COND: %s comparison are:\n%s\n',
                      p,id_ref,id_cond,
                      paste(temp.goe$Description,collapse = ', ')))
          goe <- rbind(goe, data.frame(temp.goe, id_ref = id_ref, id_cond = id_cond, phenotype = p))
        }
      }

      temp.kegg <- enrichKEGG(gene         = temp.deg$ENSEMBL,
                              universe     = allgenes$ENSEMBL,
                              organism     = 'sce',
                              pvalueCutoff = 0.05)
      if (length(temp.kegg) == 1) {
        if (dim(temp.kegg)[1] == 0) {
          cat(sprintf('There are no KEGG pathway enriched for %s ORFs in REF: %s | COND: %s comparison.\n',
                      p,id_ref,id_cond))
        } else{
          cat(sprintf('KEGG pathway enrichment for %s ORFs in REF: %s | COND: %s comparison are:\n%s\n',
                      p,id_ref,id_cond,
                      paste(temp.kegg$Description,collapse = ', ')))
          kegg <- rbind(kegg, data.frame(temp.kegg, id_ref = id_ref, id_cond = id_cond, phenotype = p))
        }
      }
    }
  }
}
goe$GeneRatio <- as.numeric(str_split(goe$GeneRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(goe$GeneRatio,'/',simplify = T)[,2])
goe$BgRatio <- as.numeric(str_split(goe$BgRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(goe$BgRatio,'/',simplify = T)[,2])
goe$GO <- paste0(goe$ONTOLOGY, '_', goe$Description)

kegg$GeneRatio <- as.numeric(str_split(kegg$GeneRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(kegg$GeneRatio,'/',simplify = T)[,2])
kegg$BgRatio <- as.numeric(str_split(kegg$BgRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(kegg$BgRatio,'/',simplify = T)[,2])

goe <- goe[order(goe$id_ref,goe$id_cond,goe$phenotype,-goe$GeneRatio,-goe$Count,goe$qvalue),]
kegg <- kegg[order(kegg$id_ref,kegg$id_cond,kegg$phenotype,-kegg$GeneRatio,-kegg$Count,kegg$qvalue),]

write.csv(goe, file = 'output/translatome/go_enrichments.csv')
write.csv(kegg, file = 'output/translatome/kegg_enrichments.csv')


##### ORF_TYPE ENRICHMENT
## NOT TRANSLATED VS TRANSIENT VS NOT TRANSIENT
# head(data.diff2)
# 
# data.cnt.all <- merge(data.diff2 %>%
#   group_by(id_ref, id_cond, cond2, cond1, orf_type, phenotype) %>%
#   count() %>%
#   data.frame(), data.diff2 %>%
#     group_by(id_ref, id_cond, cond2, cond1, orf_type) %>%
#     count() %>%
#     data.frame(), by = c('id_ref','id_cond','cond1','cond2','orf_type'),
#   suffixes = c('','_total'))
# data.cnt.all$percentage <- data.cnt.all$n/data.cnt.all$n_total * 100 
# 
# fig.pheno.prop <- data.cnt.all %>%
#   filter(orf_type != 'Reference', id_ref != 'TWO_GA', id_cond != 'TWO_GA') %>%
#   ggplot(aes(x = orf_type, y = percentage, fill = phenotype)) +
#   geom_bar(stat="identity") +
#   geom_text(aes(label = sprintf('%0.2f%%',percentage)),
#             col = 'black', size = 2,
#             position = position_stack(vjust = 0.5)) +
#   facet_wrap(.~cond2*cond1, nrow = 1,
#              labeller = labeller(cond2 = c('DM' = 'Ref: DMSO',
#                                            'GA' = 'Ref: Galactose'),
#                                  cond1 = c('FL' = 'Str: Fluconazole',
#                                            'TN' = 'Str: Tunicamycin',
#                                            'HO' = 'Str: Hydrogen Peroxide',
#                                            'HU' = 'Str: Hydroxyurea',
#                                            'SA' = 'Str: Salt'))) +
#   labs(x = 'ORF Type',
#        y = 'Phenotype Proportion') +
#   scale_fill_manual(name = 'Phenotype',
#                     values = c('Beneficial' = '#FFA000',
#                                'Neutral' = '#757575',
#                                'Deleterious' = '#303F9F')) +
#   theme_linedraw() +
#   theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
#         axis.title = element_text(size = titles),
#         axis.text = element_text(size = txt),
#         legend.title = element_text(size = titles),
#         legend.text = element_text(size = txt),
#         legend.position = 'bottom',
#         legend.key.size = unit(3, "mm"),
#         legend.box.spacing = unit(0.5,"mm"),
#         strip.text = element_text(size = txt,
#                                   face = 'bold',
#                                   margin = margin(0.1,0,0.1,0, "mm")))
# ggsave(sprintf("%s/PhenotypeProportions.jpg",fig_path), fig.pheno.prop,
#        height = one.c, width = two.c, units = 'mm',
#        dpi = 600)

## ANNOTATED VS NOT ANNOTATED
head(data.diff2)
data.diff2$annotation[str_detect(data.diff2$orf_name, 'sm')] <- 'Unannotated'
data.diff2$annotation[is.na(data.diff2$annotation)] <- 'Annotated'

length(unique(data.diff2$orf_name[data.diff2$orf_type == 'Not Translated' & data.diff2$annotation == 'Annotated']))
length(unique(data.diff2$orf_name[data.diff2$orf_type == 'Not Translated' & data.diff2$annotation == 'Unannotated']))
length(unique(data.diff2$orf_name[data.diff2$orf_type == 'Transient' & data.diff2$annotation == 'Annotated']))
length(unique(data.diff2$orf_name[data.diff2$orf_type == 'Transient' & data.diff2$annotation == 'Unannotated']))
length(unique(data.diff2$orf_name[data.diff2$orf_type == 'Not Transient' & data.diff2$annotation == 'Annotated']))
length(unique(data.diff2$orf_name[data.diff2$orf_type == 'Not Transient' & data.diff2$annotation == 'Unannotated']))

data.cnt.ann <- merge(data.diff2 %>%
                        group_by(id_ref, id_cond, cond2, cond1, orf_type, annotation, phenotype) %>%
                        count() %>%
                        data.frame(), data.diff2 %>%
                        group_by(id_ref, id_cond, cond2, cond1, orf_type, annotation) %>%
                        count() %>%
                        data.frame(), by = c('id_ref','id_cond','cond1','cond2','orf_type','annotation'),
                      suffixes = c('','_total'))
data.cnt.ann$percentage <- data.cnt.ann$n/data.cnt.ann$n_total * 100 

fig.pheno.prop.ann <- data.cnt.ann %>%
  filter(orf_type != 'Reference', id_ref != 'TWO_GA', id_cond != 'TWO_GA') %>%
  ggplot(aes(x = orf_type, y = percentage, fill = phenotype)) +
  geom_bar(stat="identity") +
  geom_label_repel(aes(label = sprintf('%0.2f%%\n%d',percentage,n), fill = phenotype),
            col = 'white', size = 2, direction = 'y', label.size = 0.15,
            force = 2, seed = 111,
            position = position_stack(vjust = 0.5)) +
  scale_x_discrete(limits = c('Not Translated','Not Transient','Transient')) +
  scale_y_continuous(breaks = seq(0,200,10)) +
  facet_wrap(.~annotation*cond2*cond1, nrow = 2,
             labeller = labeller(cond2 = c('DM' = 'Ref: DMSO',
                                           'GA' = 'Ref: Galactose'),
                                 cond1 = c('FL' = 'Str: Fluconazole',
                                           'TN' = 'Str: Tunicamycin',
                                           'HO' = 'Str: Hydrogen Peroxide',
                                           'HU' = 'Str: Hydroxyurea',
                                           'SA' = 'Str: Salt'))) +
  labs(x = 'ORF Type',
       y = 'Phenotype Proportion') +
  scale_fill_manual(name = 'Phenotype',
                    values = c('Beneficial' = '#FFA000',
                               'Neutral' = '#757575',
                               'Deleterious' = '#303F9F')) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt, angle = 30, hjust = 1, vjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0,110))
ggsave(sprintf("%s/PhenotypeProportionsAnnotation.jpg",fig_path), fig.pheno.prop.ann,
       height = two.c, width = two.c, units = 'mm',
       dpi = 600)


##### DIFF DENSITY PLOT
fig.den <- data.diff2 %>%
  filter(orf_type != 'Reference',
         id_ref != 'TWO_GA', id_cond != 'TWO_GA') %>%
  ggplot(aes(x = fitness_diff, y = orf_type)) +
  geom_density_ridges() +
  geom_vline(aes(xintercept = ll), linetype = 'dashed', col = 'red') +
  geom_vline(aes(xintercept = ul), linetype = 'dashed', col = 'red') +
  labs(x = 'Differential Fitness',
       y = '') +
  facet_wrap(.~annotation*cond2*cond1, nrow = 2,
             labeller = labeller(cond2 = c('DM' = 'Ref: DMSO',
                                           'GA' = 'Ref: Galactose'),
                                 cond1 = c('FL' = 'Str: Fluconazole',
                                           'TN' = 'Str: Tunicamycin',
                                           'HO' = 'Str: Hydrogen Peroxide',
                                           'HU' = 'Str: Hydroxyurea',
                                           'SA' = 'Str: Salt'))) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.title.y = element_blank(),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(xlim = c(-0.5,0.5))
ggsave(sprintf("%s/DifferentialFitnessDensity.jpg",fig_path), fig.den,
       height = one.5c, width = two.c, units = 'mm',
       dpi = 600)

##### PHENOTYPE COUNTS DENSITY
data.diff2 %>%
  filter(orf_type != 'Reference') %>%
  group_by(id_ref, orf_type, annotation, orf_name, strain_id, phenotype) %>%
  count() %>% data.frame() %>%
  filter(phenotype == 'Beneficial') %>%
  ggplot(aes(x = n, y = orf_type)) +
  geom_density_ridges() +
  facet_wrap(.~annotation)

# write.csv(data.diff2 %>%
#             filter(on_list == 1) %>%
#             group_by(id_cond, condition, rep, orf_name, orf_type, fitness.median_cond, cs.median_cond, phenotype) %>%
#             data.frame(), file = 'output/translatome/orfs_on_list_phenotype.csv')

# ##### PHENOTYPE ODDS RATIO
# data.or <- NULL
# for (id_ref in unique(data.cnt.all$cond2)) {
#   for (id_cond in unique(data.cnt.all$cond1[data.cnt.all$cond2 == id_ref])) {
#     for (p in data.cnt.all$phenotype[data.cnt.all$cond2 == id_ref & data.cnt.all$cond1 == id_cond]) {
#       n_pg_effect <- data.cnt.all$n[data.cnt.all$cond2 == id_ref & data.cnt.all$cond1 == id_cond & data.cnt.all$phenotype == p & data.cnt.all$orf_type == 'Transient']
#       if (length(n_pg_effect) == 0) {
#         n_pg_effect <- 0
#         n_pg_noeffect <- data.cnt.all$n_total[data.cnt.all$cond2 == id_ref & data.cnt.all$cond1 == id_cond & data.cnt.all$phenotype == 'Neutral' & data.cnt.all$orf_type == 'Transient']
#       } else {
#         n_pg_noeffect <- data.cnt.all$n_total[data.cnt.all$cond2 == id_ref & data.cnt.all$cond1 == id_cond & data.cnt.all$phenotype == p & data.cnt.all$orf_type == 'Transient'] -
#           data.cnt.all$n[data.cnt.all$cond2 == id_ref & data.cnt.all$cond1 == id_cond & data.cnt.all$phenotype == p & data.cnt.all$orf_type == 'Transient']
#       }
#      
#       n_g_effect <- data.cnt.all$n[data.cnt.all$cond2 == id_ref & data.cnt.all$cond1 == id_cond & data.cnt.all$phenotype == p & data.cnt.all$orf_type == 'Not Transient']
#       if (length(n_g_effect) == 0) {
#         n_g_effect <- 0
#         n_g_noeffect <- data.cnt.all$n_total[data.cnt.all$cond2 == id_ref & data.cnt.all$cond1 == id_cond & data.cnt.all$phenotype == 'Neutral' & data.cnt.all$orf_type == 'Not Transient']
#       } else {
#         n_g_noeffect <- data.cnt.all$n_total[data.cnt.all$cond2 == id_ref & data.cnt.all$cond1 == id_cond & data.cnt.all$phenotype == p & data.cnt.all$orf_type == 'Not Transient'] -
#           data.cnt.all$n[data.cnt.all$cond2 == id_ref & data.cnt.all$cond1 == id_cond & data.cnt.all$phenotype == p & data.cnt.all$orf_type == 'Not Transient']
#       }
#       
#       ftest <- fisher.test(matrix(c(n_pg_effect, n_pg_noeffect, n_g_effect, n_g_noeffect),2,2))
#       data.or <- rbind(data.or, c(id_ref,id_cond,p,n_pg_effect, n_pg_noeffect, n_g_effect, n_g_noeffect,
#                                   ftest$conf.int[1],ftest$estimate[[1]],ftest$conf.int[2],ftest$p.value))
#     }
#   }
# }
# data.or <- data.frame(data.or, stringsAsFactors = F)
# colnames(data.or) <- c('id_ref','id_cond','phenotype','pg_w_effect','pg_wo_effect','g_w_effect','g_wo_effect','top','or','bottom','p')
# head(data.or)
# data.or$significant[data.or$p <= 0.05] <- 'Yes'
# data.or$significant[is.na(data.or$significant)] <- 'No'
# 
# data.or$pg_w_effect <- as.numeric(data.or$pg_w_effect)
# data.or$pg_wo_effect <- as.numeric(data.or$pg_wo_effect)
# data.or$g_w_effect <- as.numeric(data.or$g_w_effect)
# data.or$g_wo_effect <- as.numeric(data.or$g_wo_effect)
# data.or$top <- as.numeric(data.or$top)
# data.or$or <- as.numeric(data.or$or)
# data.or$bottom <- as.numeric(data.or$bottom)
# data.or$p <- as.numeric(data.or$p)
# 
# data.or$label[data.or$p > 0.05] <- 'ns'
# data.or$label[data.or$p <= 0.05] <- '*'
# data.or$label[data.or$p <= 0.01] <- '**'
# data.or$label[data.or$p <= 0.001] <- '***'
# data.or$label[data.or$p <= 0.0001] <- '****'
# 
# plot.or <- ggplot(data.or,aes(x=phenotype,y=or,ymin=bottom,ymax=top))+
#   geom_point(stat="identity",  shape=21, size=2, stroke=1, fill = "white")+
#   geom_errorbar(width=0.25)+
#   geom_hline(yintercept=1,linetype="dashed",color="red")+
#   geom_text(aes(x = phenotype, y = 28, label = label), size = 2, col = 'red') +
#   labs(x='Phenotype',y='Odds Ratio') +
#   # scale_y_continuous(trans="log10", breaks=c(0.2,0.5,1,2,5,10)) + 
#   # annotation_logticks(sides="l") +
#   facet_wrap(.~id_ref * id_cond) +
#   theme_linedraw() +
#   theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
#         axis.title = element_text(size = titles),
#         axis.text = element_text(size = txt),
#         legend.title = element_text(size = titles),
#         legend.text = element_text(size = txt),
#         legend.position = 'bottom',
#         legend.key.size = unit(3, "mm"),
#         legend.box.spacing = unit(0.5,"mm"),
#         strip.text = element_text(size = txt,
#                                   face = 'bold',
#                                   margin = margin(0.1,0,0.1,0, "mm")))
# ggsave(sprintf("%s/PHENOTYPE_ODDS.jpg",fig_path), plot.or,
#        height = one.5c, width = one.5c, units = 'mm',
#        dpi = 600) 


##### SAVE ALL
# save(list = ls(.GlobalEnv), file = "output/translatome/AllData.RData")


##### UPLOAD TO MySQL
temp <- unique(data.diff$orf_name[!(data.diff$strain_id %in% controls$strain_id) & data.diff$orf_name != 'BF_control'])
temp <- bitr(temp, fromType = "ORF",
     toType = c("DESCRIPTION"),
     OrgDb = org.Sc.sgd.db)

dbWriteTable(conn, 'TR_OE_DIFF_FITNESS',
             merge(data.diff[!(data.diff$strain_id %in% controls$strain_id),] %>% 
               filter(orf_name != 'BF_control', arm1 == arm2, 
                      id_ref %in% c('ONE_GA','TWO_DM'), id_cond != 'TWO_GA'), temp, by.x = 'orf_name', by.y = 'ORF', all = T), 
             overwrite = T)




