##### SAN DIEGO ANALYSIS
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 09/10/2019

##### INITIALIZE
library(RMariaDB)
library(ggplot2)
library(ggExtra)
library(gridExtra)
library(ggpubr)
library(reshape2)
library(cluster)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)
source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

out_path = 'figs/sandiego/';

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 7
txt <- 5
lbls <- 9

# ##### 5 CONDITION DATA
# cond5 <- dbGetQuery(conn, 'select distinct a.orf_name, hours, n, colony_size normalized_cs, q_cs, effect_cs,exp_id,
#                     AGE, chromosome, translation_media, category
#                     from brian_031918.DATASET_6 a, ANNE.SUMMARY_2012 b where a.orf_name = b.orf_name and exp_id in (28,31,91,33, 93) and a.orf_name !="BF_control"')
# save(cond5,
#      file = sprintf('%sCONDITION5.RData',out_path))
# 
# ##### FITNESS RESULTS
# load(sprintf('%sCONDITION5.RData',out_path))
# cond5$effect_cs <- factor(cond5$effect_cs, levels = c('beneficial','neutral','deleterious'))
# 
# ggplot(cond5,
#        aes(x = category, y = normalized_cs)) +
#   # geom_boxplot() +
#   geom_jitter(aes(col = effect_cs), alpha = 0.5, size = 0.2) +
#   geom_violin(col = '#757575', fill = 'transparent') +
#   # stat_compare_means() +
#   facet_wrap(.~exp_id, nrow = 1) +
#   scale_color_manual(breaks = c('beneficial','neutral','deleterious'),
#                      values = c('deleterious'='#3F51B5',
#                                 'neutral'='#212121',
#                                 'beneficial'='#FFC107')) +
#   theme_linedraw() +
#   theme(axis.title = element_text(size = txt),
#         axis.text = element_text(size = txt),
#         legend.title = element_text(size = titles),
#         legend.text = element_text(size = txt),
#         legend.position = 'right',
#         strip.text = element_text(size = txt,
#                                   margin = margin(0.1,0,0.1,0, "mm")),
#         legend.key.size = unit(3, "mm"),
#         legend.box.spacing = unit(0.5,"mm"))
# ggsave(sprintf("%sfitness.jpg",out_path),
#        height = 60, width = two.c, units = 'mm',
#        dpi = 300)


##### FITNESS HEATMAP
# ggplot(cond5) +
#   geom_tile(aes(x = as.factor(exp_id), y = orf_name, col = normalized_cs)) +
#   scale_color_distiller(palette = 'Set1')

# cond5.hm <- with(cond5, tapply(normalized_cs, list(orf_name, exp_id), c))
# cond5.hm[is.na(cond5.hm)] <- 0
# 
# heatmap(cond5.hm, scale = "none")
# 
# cond5.hm <- data.frame(cond5.hm)
# orf_type <- NULL
# age <- NULL
# for (orf in rownames(cond5.hm)) {
#   orf_type <- c(orf_type, cond5$category[cond5$orf_name == orf][1])
#   age <- c(age, cond5$AGE[cond5$orf_name == orf][1])
# }
# cond5.hm$orf_type <- orf_type
# cond5.hm$age <- age
# 
# save(cond5.hm,
#      file = sprintf('%sCONDITION5_HEATMAP.RData',out_path))
# load(sprintf('%sCONDITION5_HEATMAP.RData',out_path))
# 
# annot_df <- data.frame(orf_type = cond5.hm$orf_type,
#                        age = cond5.hm$age)
# col = list(orf_type = c("gene" = "green", "proto-gene" = "black"),
#            age = circlize::colorRamp2(c(0, 10), 
#                                       c("lightblue", "purple")))
# ha <- HeatmapAnnotation(df = annot_df, col = col,
#                         which = 'row')
# 
# Heatmap(cond5.hm[1:5], name = "fitness",
#         left_annotation  = ha,
#         show_row_names = F)

# ##### LID ANALYZED RESULTS
# fit.sd <- dbGetQuery(conn, 'select * from SDS_LI_6144_FITNESS
#                       where orf_name in
#                       (select distinct orf_name from BARFLEX_SPACE_AGAR
#                       where 384plate in (11,22) or orf_name = "BF_control")
#                       order by orf_name')
# 
# fit.oesp <- dbGetQuery(conn, 'select * from OESP1_FS_6144_FITNESS
#                       where hours = 22 and orf_name in
#                       (select distinct orf_name from BARFLEX_SPACE_AGAR
#                       where 384plate in (11,22) or orf_name = "BF_control")
#                       order by orf_name')
# 
# fit.all <- NULL
# i <- 1
# for (orf in unique(fit.sd$orf_name)) {
#   fit.all$orf_name[i] <- orf
#   fit.all$sd[i] <- median(fit.sd$fitness[fit.sd$orf_name == orf], na.rm = T)
#   fit.all$pitt[i] <- median(fit.oesp$fitness[fit.sd$orf_name == orf], na.rm = T)
#   i <- i + 1
# }
# fit.all <- data.frame(fit.all)
# 
# ggplot(fit.all) +
#   geom_abline() +
#   geom_point(aes(x = sd, y = pitt)) +
#   coord_cartesian(xlim = c(0,1.2),
#                   ylim = c(0,1.2))
# 
# ggplot() +
#   geom_line(data = fit.sd[fit.sd$orf_name == 'BF_control',],
#             aes(x = fitness, col = 'sd'), stat = 'density', trim = T) +
#   geom_line(data = fit.oesp[fit.oesp$orf_name == 'BF_control',],
#             aes(x = fitness, col = 'pitt'), stat = 'density', trim = T) +
#   coord_cartesian(xlim = c(0.8,1.2))

##### PLATEMAPS
platemap.sd <- dbGetQuery(conn, "select a.*, c.exp_id, c.orf_name, c.hours, c.average, c.fitness, d.effect_cs
                          from brian_data.POS2COOR a, brian_031918.POS2ORF_NAME_new b, brian_031918.FITNESS_ALL6 c,
                          brian_031918.DATASET_6 d
                          where a.pos = b.pid and a.pos = c.pid and c.exp_id in (28,31,91,33,93)
                          and c.orf_name = d.orf_name and c.exp_id = d.exp_id
                          order by c.exp_id, plate, col, row")
platemap.sd$colony[platemap.sd$orf_name == "BF_control"] <- 'reference'
platemap.sd$colony[is.na(platemap.sd$colony)] <- 'query'

platemap.sd$source[platemap.sd$row%%2==1 & platemap.sd$col%%2==1] = 'Top Left'
platemap.sd$source[platemap.sd$row%%2==0 & platemap.sd$col%%2==1] = 'Bottom Left'
platemap.sd$source[platemap.sd$row%%2==1 & platemap.sd$col%%2==0] = 'Top Right'
platemap.sd$source[platemap.sd$row%%2==0 & platemap.sd$col%%2==0] = 'Bottom Right'

platemap.sd$source <- factor(platemap.sd$source, levels = c('Top Left','Top Right','Bottom Left','Bottom Right'))

ggplot(platemap.sd) +
  geom_point(aes(x= average, y= fitness, col = effect_cs, shape = colony), size = 0.5) +
  scale_color_manual(breaks = c('beneficial','neutral','deleterious'),
                     values = c('deleterious'='#3F51B5',
                                'neutral'='#212121',
                                'beneficial'='#FFC107')) +
  facet_wrap(.~exp_id*plate*source, ncol = 4) +
  theme_linedraw() +
  theme(axis.title = element_text(size = txt),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'right',
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")),
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"))
ggsave(sprintf("%saverageVSfitness.jpg",out_path),
       height = two.c*5, width = two.c, units = 'mm',
       dpi = 300)

platemap.sd <- dbGetQuery(conn, "select a.*, c.exp_id, c.orf_name, c.hours, c.average, c.fitness
                          from brian_data.POS2COOR a, brian_031918.POS2ORF_NAME_new b, brian_031918.FITNESS_ALL6 c
                          where a.pos = b.pid and a.pos = c.pid and c.exp_id in (28,31,91,33,93)
                          order by c.exp_id, plate, col, row")
platemap.sd$colony[platemap.sd$orf_name == "BF_control"] <- 'reference'
platemap.sd$colony[is.na(platemap.sd$colony)] <- 'query'

platemap.sd$source[platemap.sd$row%%2==1 & platemap.sd$col%%2==1] = 'Top Left'
platemap.sd$source[platemap.sd$row%%2==0 & platemap.sd$col%%2==1] = 'Bottom Left'
platemap.sd$source[platemap.sd$row%%2==1 & platemap.sd$col%%2==0] = 'Top Right'
platemap.sd$source[platemap.sd$row%%2==0 & platemap.sd$col%%2==0] = 'Bottom Right'

platemap.sd$source <- factor(platemap.sd$source, levels = c('Top Left','Top Right','Bottom Left','Bottom Right'))

ggplot(platemap.sd) +
  geom_point(aes(x= col, y=row, col = fitness, shape = colony), size = 0.5) +
  scale_color_distiller(palette = 'PRGn', na.value = 'white',
                        limits = c(0,1.3), oob = squish) +
  facet_wrap(.~exp_id*plate*source, ncol = 4) +
  theme_linedraw() +
  theme(axis.title = element_text(size = txt),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'right',
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")),
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"))
ggsave(sprintf("%sFitnessLandscape.jpg",out_path),
       height = 150*5, width = two.c, units = 'mm',
       dpi = 300)

ggplot(platemap.sd) +
  # geom_boxplot(aes(x = source, y = fitness)) +
  geom_violin(aes(x = source, y = fitness), fill = 'transparent', trim = T) +
  facet_wrap(.~exp_id*plate, ncol = 5) +
  theme_linedraw() +
  theme(axis.title = element_text(size = txt),
        axis.text = element_text(angle = 45, hjust = 1, size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'right',
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")),
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"))
ggsave(sprintf("%sFitnessDensity.jpg",out_path),
       height = two.c, width = two.c, units = 'mm',
       dpi = 300)

ggplot(platemap.sd) +
  geom_point(aes(x= col, y=row, col = colony), size = 0.3) +
  facet_wrap(.~exp_id*plate*source, ncol = 4) +
  theme_linedraw() +
  theme(axis.title = element_text(size = txt),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'right',
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")),
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"))
ggsave(sprintf("%splatemap.jpg",out_path),
       height = 150*5, width = two.c, units = 'mm',
       dpi = 300)

##### 
plate.colony <- plyr::count(platemap.sd, vars = c("plate"," source", "colony"))
plate.colony[plate.colony$colony == 'reference',]

for (pl in unique(platemap.sd$plate)) {
  for (sr in unique(platemap.sd$source[platemap.sd$plate == pl])) {
    if (sum(plate.colony$colony[plate.colony$plate == pl & plate.colony$source == sr] == 'reference') == 0) {
      platemap.sd$orf_name[platemap.sd$plate == pl & platemap.sd$source == sr & platemap.sd$colony == 'query'] <- NA
    }
  }
}

platemap.sd$colony <- NA
platemap.sd$colony[platemap.sd$orf_name == "BF_control"] <- 'reference'
platemap.sd$colony[is.na(platemap.sd$colony)] <- 'query'

dbWriteTable(conn, "SDS_LI_pos2orf_name2", platemap.sd[c("pos","source","orf_name")], overwrite = T)

##### REFERENCE REPLICATES
ref.pos <- dbGetQuery(conn, 'select *
                      from SDS_LI_pos2orf_name2 a, SDS_LI_pos2coor b
                      where a.pos = b.pos and a.orf_name = "BF_control"
                      order by plate, col, row')

ref.reps <- NULL
i = 1
for (pl in unique(ref.pos$plate)) {
  if (i == 1) {
    ref.reps <- as.numeric(ref.pos$pos[ref.pos$plate == pl])
    i = i + 1
  } else {
    ref.reps <- rbind(ref.reps,as.numeric(ref.pos$pos[ref.pos$plate == pl]))
  }
}
ref.reps <- data.frame(t(ref.reps))
colnames(ref.reps) <- c('one','two','three','four')

dbWriteTable(conn, "SDS_LI_REFREPS", ref.reps, overwrite = T)

##### SDS vs OESP1
sdsoesp1 <- dbGetQuery(conn, "select a.orf_name, a.cs_mean oesp_cs, b.cs_mean sds_cs,
                       a.cs_median oesp_cs_median, b.cs_median sds_cs_median,
                       c.p oesp_p, c.stat oesp_stat, d.p sds_p, d.stat sds_stat
                       from OESP1_FS_6144_FITNESS_STATS a, SDS_LI_6144_FITNESS_STATS b,
                       OESP1_FS_6144_PVALUE c, SDS_LI_6144_PVALUE d
                       where a.orf_name = b.orf_name and a.hours = 22
                       and a.orf_name = c.orf_name and a.hours = c.hours
                       and a.orf_name = d.orf_name")

sdsoesp1$effect[sdsoesp1$oesp_p <= 0.05 & sdsoesp1$oesp_stat < 0 &
                  sdsoesp1$sds_p <= 0.05 & sdsoesp1$sds_stat < 0] <- 'D,D'
sdsoesp1$effect[sdsoesp1$oesp_p > 0.05 & sdsoesp1$sds_p > 0.05] <- 'N,N'
sdsoesp1$effect[sdsoesp1$oesp_p <= 0.05 & sdsoesp1$oesp_stat < 0 &
                  sdsoesp1$sds_p > 0.05] <- 'D,N'

ggplot(sdsoesp1) +
  geom_abline() +
  geom_point(aes(x = oesp_cs, y = sds_cs, col = effect)) +
  scale_color_discrete(name = 'EFFECT\n(PITT,SD)') +
  labs(x = 'PITT', y = 'SD') +
  coord_cartesian(xlim = c(0.75,1),
                  ylim = c(0.75,1)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = txt),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'right',
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")),
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"))
ggsave(sprintf("%ssdVSpitt.jpg",out_path),
       height = 80, width = one.c, units = 'mm',
       dpi = 300)

