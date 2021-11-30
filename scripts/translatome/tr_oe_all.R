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
controls <- read.csv(file = 'rawdata/translatome/condition_controls.csv', stringsAsFactors = F)
pgs <- dbGetQuery(conn, 'select orf_name from PROTOGENES where pg_2012 = 1')

orf_types <- read.csv(file = 'rawdata/translatome/oe_transient')
orf_types$category[orf_types$is_transient + orf_types$translated + orf_types$is_candidate == 3] <- 'Transient'
orf_types$category[orf_types$is_conserved + orf_types$translated + orf_types$is_candidate == 3] <- 'Conserved'
orf_types$category[is.na(orf_types$category)] <- 'Others'

tr.conds <- data.frame(arms = c('ONE','ONE','ONE','ONE','TWO','TWO','TWO'),
                       conds = c('GA','SA','HO','HU','DM','FL','TN'))

borders <- dbGetQuery(conn, 'select * from TR_OE_borderpos')
smudge <- dbGetQuery(conn, 'select * from TR_OE_ONE_FS_GA_smudgebox')


##### GATHER DATA
data.fit.all <- NULL
for (a in unique(tr.conds$arms)) {
  for (c in tr.conds$conds[tr.conds$arms == a]) {
    temp.fit <- dbGetQuery(conn, sprintf('select a.*, b.density, b.plate, b.col, b.row
                                         from TR_OE_ALL_FS_%s_%s_6144_FITNESS a, TR_OE_pos2coor b
                                         where a.pos = b.pos
                                         order by a.hours, b.plate, b.col, b.row', a, c))
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

data.fit.all <- merge(data.fit.all, orf_types, by = 'orf_name', all.x = T)
data.fit.all$category[data.fit.all$orf_name == 'BF_control'] <- 'Reference'
# data.fit.all[is.na(data.fit.all$category),] %>%
#   group_by(orf_name) %>%
#   count() %>% data.frame()
data.fit.all$rep <- as.numeric(str_trunc(as.character(data.fit.all$pos), 5, side = 'left', ellipsis = ''))

data.fit.all <- data.fit.all[!(data.fit.all$condition == 'FL' & data.fit.all$hours == 192) &
                             !(data.fit.all$condition == 'GA' & data.fit.all$hours == 67) &
                             !(data.fit.all$condition == 'SA' & data.fit.all$hours %in% c(51,63,72)),]

data.fit.sum <- data.fit.all %>%
  group_by(arm, condition, hours, strain_id, orf_name, category) %>%
  summarize(avg.median = median(average, na.rm = T), avg.mad = mad(average, na.rm = T),
            fitness.median = median(fitness, na.rm = T), fitness.mad = mad(fitness, na.rm = T),
            .groups = 'keep') %>%
  data.frame()

##### REMOVE OUTLIERS
data.fit.all <- merge(data.fit.all, data.fit.sum,
                      by = c('arm','condition','hours','strain_id','orf_name','category'), all = T)

data.fit.all$average[data.fit.all$average < (data.fit.all$avg.median - 2*data.fit.all$avg.mad) |
                       data.fit.all$average > (data.fit.all$avg.median + 2*data.fit.all$avg.mad)] <- NA
data.fit.all$fitness[data.fit.all$fitness < (data.fit.all$fitness.median - 2*data.fit.all$fitness.mad) |
                       data.fit.all$fitness > (data.fit.all$fitness.median + 2*data.fit.all$fitness.mad)] <- NA

data.fit.sum <- data.fit.all %>%
  group_by(arm, condition, hours, strain_id, orf_name, category) %>%
  summarize(avg.median = median(average, na.rm = T), avg.mad = mad(average, na.rm = T),
            fitness.median = median(fitness, na.rm = T), fitness.mad = mad(fitness, na.rm = T),
            .groups = 'keep') %>%
  data.frame()
head(data.fit.sum)

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
data.fit.sum$norm_cs <- data.fit.sum$fitness.median * data.fit.sum$avg_m

##### PLOT GCS
data.fit.sum$condition <- factor(data.fit.sum$condition, levels = c('GA','HO','SA','FL','TN'))

fig.gc <- data.fit.sum %>% #[!(data.cs.sum$condition == 'SA' & data.cs.sum$orf_name == 'YGR170W'),] %>%
  filter(category != 'Reference') %>%
  ggplot(aes(x = hours, y = norm_cs)) +
  # geom_smooth(aes(group = strain_id, col = orf_type), method = 'loess') +
  # geom_line(aes(group = strain_id, col = category)) +
  # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
  #              aes(group = orf_type, fill = orf_type), geom="ribbon", alpha = 0.4) +
  # stat_summary(aes(group = orf_type, col = orf_type), fun=mean, geom="line", lwd =0.7) +
  geom_line(data = merge(data.fit.sum, melt(controls, id.vars = c('strain_id','standard_name','orf_name'),
                                           variable.name = 'condition', value.name = 'control_type'),
                         by = c('condition','strain_id','orf_name')) %>%
              filter(control_type != '', standard_name != 'ALD3') %>%
              group_by(arm, condition, hours, strain_id, orf_name, standard_name, control_type) %>%
              summarize(norm_cs = median(norm_cs, na.rm =T), .groups = 'keep') %>% data.frame(),
            aes(group = strain_id, col = control_type), size = 1, linetype = 'solid') +
  geom_line(aes(y = avg_ll), size = 1, col = 'red', linetype = 'dashed') +
  geom_line(aes(y = avg_ul), size = 1, col = 'red', linetype = 'dashed') +
  geom_label_repel(data = merge(merge(data.fit.sum, melt(controls, id.vars = c('strain_id','standard_name','orf_name'),
                                                        variable.name = 'condition', value.name = 'control_type'),
                                      by = c('condition','strain_id','orf_name')) %>%
                                  filter(control_type != '', standard_name != 'ALD3') %>%
                                  group_by(arm, condition, hours, strain_id, orf_name, standard_name, control_type) %>%
                                  summarize(norm_cs = median(norm_cs, na.rm =T), .groups = 'keep') %>% data.frame(),
                                data.fit.sum %>%
                                  group_by(arm, condition) %>%
                                  summarize(max_hrs = max(hours, na.rm = T), .groups = 'keep') %>%
                                  data.frame(), by = c('arm','condition')) %>%
                     filter(hours == max_hrs),
                   aes(x = hours, y = norm_cs, label = standard_name, fill = control_type), alpha = 0.6,
                   size = 1.2, min.segment.length = 0, max.overlaps = 20) +
  # scale_color_manual(name = '',
  #                    breaks = c('Resistant',
  #                               'Sensitive',
  #                               'Not Transient',
  #                               'Transient',
  #                               'Not Translated'),
  #                    values = c('Resistant' = '#FFA000',
  #                               'Sensitive' = '#303F9F',
  #                               'Not Translated' = '#009688',
  #                               'Not Transient' = '#FF5722',
  #                               'Transient' = '#E040FB')) +
  # scale_fill_manual(name = '',
  #                    breaks = c('Resistant',
  #                               'Sensitive',
  #                               'Not Transient',
  #                               'Transient',
  #                               'Not Translated'),
  #                    values = c('Resistant' = '#FFA000',
  #                               'Sensitive' = '#303F9F',
  #                               'Not Translated' = '#009688',
  #                               'Not Transient' = '#FF5722',
  #                               'Transient' = '#E040FB'),
  #                   guide = F) +
  labs(x = 'Colony Size (pixels)',
       y = 'Time (hours)') +
  facet_wrap(.~condition, scales = 'free_x', nrow = 3,
             labeller = labeller(condition = c('GA' = 'Galactose (Arm 1)',
                                               'HO' = 'Hydrogen Peroxide',
                                               'HU' = 'Hydroxyurea',
                                               'SA' = 'Salt',
                                               # 'TWO_GA' = 'Galactose (Arm 2)',
                                               # 'TWO_DM' = 'DMSO',
                                               'FL' = 'Fluconazole',
                                               'TN' = 'Tunicamycin'))) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
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
                              override.aes = list(size = 2))) +
  coord_cartesian(ylim = c(0,4500))
ggsave(sprintf("%s/GROWTHCURVES_CONTROLS.jpg",fig_path), fig.gc,
       height = two.c, width = two.c, units = 'mm',
       dpi = 600)


##### TOP 100 ENRICHMENT ANALYSIS
head(data.fit.sum)

goe.100 <- data.frame()
for (a in unique(data.fit.sum$arm)) {
  for (c in unique(data.fit.sum$condition[data.fit.sum$arm == a])) {
    for (h in unique(data.fit.sum$hours[data.fit.sum$arm == a & data.fit.sum$condition == c])) {
      temp <- data.fit.sum[data.fit.sum$arm == a & data.fit.sum$condition == c & data.fit.sum$hours == h &
                           data.fit.sum$orf_name != 'BF_control' & data.fit.sum$strain_id < 1000000,]
      temp.top <- head(temp[order(-temp$fitness.median),],100)
      temp.bot <- head(temp[order(temp$fitness.median),],100)
      
      allgenes <- unique(temp$orf_name)
      allgenes <- bitr(allgenes, fromType = "ORF",
                       toType = c("ENTREZID","GENENAME","ENSEMBL"),
                       OrgDb = org.Sc.sgd.db)
      allgenes <- allgenes[!is.na(allgenes$ENSEMBL),]
      
      temp.top <- bitr(temp.top$orf_name, fromType = "ORF",
                       toType = c("ENTREZID","GENENAME","ENSEMBL"),
                       OrgDb = org.Sc.sgd.db)
      temp.top <- temp.top[!is.na(temp.top$ENSEMBL),]
      temp.top.goe <- enrichGO(gene          = temp.top$ENSEMBL,
                               universe      = allgenes$ENSEMBL,
                               OrgDb         = org.Sc.sgd.db,
                               keyType       = "ENSEMBL",
                               ont           = "ALL",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.01,
                               qvalueCutoff  = 0.05)
      if (length(temp.top.goe) == 1) {
        if (dim(temp.top.goe)[1] == 0) {
        } else{
          goe.100 <- rbind(goe.100, data.frame(arm = a, condition = c, hours = h, category = 'top100', temp.top.goe))
        }
      }
      
      temp.bot <- bitr(temp.bot$orf_name, fromType = "ORF",
                       toType = c("ENTREZID","GENENAME","ENSEMBL"),
                       OrgDb = org.Sc.sgd.db)
      temp.bot <- temp.bot[!is.na(temp.bot$ENSEMBL),]
      temp.bot.goe <- enrichGO(gene          = temp.bot$ENSEMBL,
                               universe      = allgenes$ENSEMBL,
                               OrgDb         = org.Sc.sgd.db,
                               keyType       = "ENSEMBL",
                               ont           = "ALL",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.01,
                               qvalueCutoff  = 0.05)
      if (length(temp.bot.goe) == 1) {
        if (dim(temp.bot.goe)[1] == 0) {
        } else{
          goe.100 <- rbind(goe.100, data.frame(arm = a, condition = c, hours = h, category = 'bottom100', temp.bot.goe))
        }
      }
    }
  }
}

##### SAVE DATA
write.csv(data.fit.sum, file = 'output/translatome/tr_oe_fitandcs_highres.csv')
write.csv(data.fit.all, file = 'output/translatome/tr_oe_fitandcs_highres_all.csv')




