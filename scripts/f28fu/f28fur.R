##### F28FU - CONT EXP
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 02/19/2021

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
library(wesanderson)
library(ggrepel)
library(reshape2)
library(effsize)
library(qvalue)
library(RMariaDB)

source("R/functions/isoutlier.R")
source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")
out_path = 'figs/f28fu/repeat2/'

##### FIGURE SIZE
one.c <- 90 #single column 
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 7
txt <- 7
lbls <- 9

##### 5FOA EXPERIMENT
##### GATHER DATA

tblname_dat1 <- 'F28FUR2_FS2_6144_FITNESS'
tblname_dat2 <- 'F28FUR2_FS2_R1_1536_FITNESS'
tblname_dat3 <- 'F28FUR2_FS2_R2_1536_FITNESS'
tblname_p2o <- 'F28FUR_pos2orf_name'
tblname_p2c <- 'F28FUR_pos2coor'

temp_6144 <- dbGetQuery(conn, sprintf('select c.*, b.orf_name, a.hours, a.average, a.fitness
                                     from %s a, %s b, %s c
                                     where a.pos = b.pos and a.pos = c.pos
                                     order by density, hours, plate, col, row',
                                     tblname_dat1, tblname_p2o, tblname_p2c))
temp_1536_r1 <- dbGetQuery(conn, sprintf('select c.*, b.orf_name, a.hours, a.average, a.fitness
                                     from %s a, %s b, %s c
                                      where a.pos = b.pos and a.pos = c.pos
                                      order by density, hours, plate, col, row',
                                      tblname_dat2, tblname_p2o, tblname_p2c))
temp_1536_r2 <- dbGetQuery(conn, sprintf('select c.*, b.orf_name, a.hours, a.average, a.fitness
                                     from %s a, %s b, %s c
                                      where a.pos = b.pos and a.pos = c.pos
                                      order by density, hours, plate, col, row',
                                      tblname_dat3, tblname_p2o, tblname_p2c))
temp_6144$replicate <- 'R1'
temp_1536_r1$replicate <- 'R1'
temp_1536_r2$replicate <- 'R2'

foa.data <- rbind(temp_6144, temp_1536_r1, temp_1536_r2)

foa.data$condition[foa.data$plate == 1] <- 'GLU'
foa.data$condition[foa.data$plate == 2] <- 'GAL'
foa.data$condition[foa.data$plate == 3] <- 'CAS'
foa.data$condition[foa.data$plate == 4] <- 'SDA'

foa.data$source[foa.data$row%%2==1 & foa.data$col%%2==1] = 'TL'
foa.data$source[foa.data$row%%2==0 & foa.data$col%%2==1] = 'BL'
foa.data$source[foa.data$row%%2==1 & foa.data$col%%2==0] = 'TR'
foa.data$source[foa.data$row%%2==0 & foa.data$col%%2==0] = 'BR'
foa.data$source <- factor(unique(foa.data$source), levels = c('TL','TR','BL','BR'))

foa.data$rep <- as.numeric(str_trunc(as.character(foa.data$pos), 4, side = 'left', ellipsis = ''))

##### CLEAN DATA
head(foa.data)

for (d in unique(foa.data$density)) {
  for (rr in unique(foa.data$replicate[foa.data$density == d])) {
    for (h in unique(foa.data$hours[foa.data$density == d & foa.data$replicate == rr])) {
      for (p in unique(foa.data$plate[foa.data$hours == h & foa.data$density == d & foa.data$replicate == rr])) {
        # temp <- foa.data[foa.data$hours == h & foa.data$density == d & foa.data$plate == p & foa.data$replicate == rr,]
        # temp$average[isoutlier(temp$average, 2)] <- NA
        # foa.data$average_clean[foa.data$hours == h & foa.data$density == d & foa.data$plate == p & foa.data$replicate == rr] <- temp$average
        # 
        # temp$fitness[isoutlier(temp$fitness, 2)] <- NA
        # foa.data$fitness_clean[foa.data$hours == h & foa.data$density == d & foa.data$plate == p & foa.data$replicate == rr] <- temp$fitness
        # 
        for (o in unique(foa.data$orf_name[foa.data$hours == h & foa.data$density == d & foa.data$replicate == rr & foa.data$plate == p & foa.data$orf_name != 'BOR' & !is.na(foa.data$orf_name)])) {
          temp <- foa.data[foa.data$hours == h & foa.data$density == d & foa.data$replicate == rr & foa.data$plate == p & foa.data$orf_name == o & !is.na(foa.data$orf_name),]
          temp$average[isoutlier(temp$average, 2)] <- NA
          foa.data$average_clean[foa.data$hours == h & foa.data$density == d & foa.data$replicate == rr & foa.data$plate == p & foa.data$orf_name == o & !is.na(foa.data$orf_name)] <- temp$average
          
          temp$fitness[isoutlier(temp$fitness, 2)] <- NA
          foa.data$fitness_clean[foa.data$hours == h & foa.data$density == d & foa.data$replicate == rr & foa.data$plate == p & foa.data$orf_name == o & !is.na(foa.data$orf_name)] <- temp$fitness
          
          # for (r in unique(foa.data$rep[foa.data$hours == h & foa.data$density == d & foa.data$replicate == rr & foa.data$plate == p & foa.data$orf_name == o & !is.na(foa.data$orf_name)])) {
          #   temp <- foa.data[foa.data$hours == h & foa.data$density == d & foa.data$replicate == rr & foa.data$plate == p & foa.data$orf_name == o & !is.na(foa.data$orf_name) & foa.data$rep == r,]
          #   temp$average_clean[isoutlier(temp$average_clean, 2)] <- NA
          #   foa.data$average_clean[foa.data$hours == h & foa.data$density == d & foa.data$replicate == rr & foa.data$plate == p & foa.data$orf_name == o & !is.na(foa.data$orf_name) & foa.data$rep == r] <- temp$average_clean
          #   
          #   temp$fitness_clean[isoutlier(temp$fitness_clean, 2)] <- NA
          #   foa.data$fitness_clean[foa.data$hours == h & foa.data$density == d & foa.data$replicate == rr & foa.data$plate == p & foa.data$orf_name == o & !is.na(foa.data$orf_name) & foa.data$rep == r] <- temp$fitness_clean
          # }
        }
      }
    }
  }
}
save(foa.data, file = sprintf('%sf28fur2_fs2_data.RData',out_path))

# ##### PLOT COLONY SIZE
# plot.den.cs <- ggplot(foa.data[!(foa.data$orf_name %in% c('BOR','REF')) & !is.na(foa.data$orf_name) & foa.data$hours > 0,],
#                       aes(x = average_clean, y = orf_name, fill = orf_name)) +
#   geom_density_ridges(quantile_lines = TRUE,
#                       quantiles = c(0.25, 0.5, 0.75),
#                       scale = 3, alpha = 0.8, size = 0.3,
#                       vline_size = 0.2, vline_color = "black",
#                       na.rm = T) +
#   # geom_vline(xintercept = c(0.99,1.01), linetype = 'dashed', col = 'red', lwd = 0.5) +
#   # scale_x_continuous(breaks = seq(0,2,0.1),
#   #                    minor_breaks = seq(0,2,0.01)) +
#   labs(y = 'Strains', x = 'Colony Size (pix.)') +
#   scale_fill_discrete(guide = F) +
#   facet_wrap(.~density*replicate*hours*condition,
#              ncol = 4, scales = 'free_x') +
#   theme_linedraw() +
#   theme(plot.title = element_blank(),
#         axis.title = element_text(size = titles),
#         axis.text = element_text(size = txt),
#         axis.text.y = element_text(angle = -45, vjust = 1),
#         legend.title = element_text(size = titles),
#         legend.text = element_text(size = txt),
#         legend.key.size = unit(3, "mm"),
#         legend.position = "bottom",
#         legend.box.spacing = unit(0.5,"mm"),
#         strip.text = element_text(size = txt,
#                                   margin = margin(0.1,0,0.1,0, "mm")))
# 
# plot.den.fit <- ggplot(foa.data[!(foa.data$orf_name %in% c('BOR','REF')) & !is.na(foa.data$orf_name) & foa.data$hours > 0,],
#                        aes(x = fitness_clean, y = orf_name, fill = orf_name)) +
#   geom_density_ridges(quantile_lines = TRUE,
#                       quantiles = c(0.25, 0.5, 0.75),
#                       scale = 3, alpha = 0.8, size = 0.3,
#                       vline_size = 0.2, vline_color = "black",
#                       na.rm = T) +
#   # geom_vline(xintercept = c(0.99,1.01), linetype = 'dashed', col = 'red', lwd = 0.5) +
#   scale_x_continuous(breaks = seq(0,2,0.1),
#                      minor_breaks = seq(0,2,0.01)) +
#   labs(y = 'Strains', x = 'Fitness') +
#   scale_fill_discrete(guide = F) +
#   facet_wrap(.~density*replicate*hours*condition,
#              ncol = 4) +
#   theme_linedraw() +
#   theme(plot.title = element_blank(),
#         axis.title = element_text(size = titles),
#         axis.text = element_text(size = txt),
#         axis.text.y = element_text(angle = -45, vjust = 1),
#         legend.title = element_text(size = titles),
#         legend.text = element_text(size = txt),
#         legend.key.size = unit(3, "mm"),
#         legend.position = "bottom",
#         legend.box.spacing = unit(0.5,"mm"),
#         strip.text = element_text(size = txt,
#                                   margin = margin(0.1,0,0.1,0, "mm")))
# 
# plot.box.fit <- ggplot(foa.data[!(foa.data$orf_name %in% c('BOR','REF')) & !is.na(foa.data$orf_name) & foa.data$hours > 0,],
#                        aes(y = fitness_clean, x = orf_name, fill = orf_name)) +
#   geom_boxplot(outlier.shape = NA) +
#   stat_compare_means(method = 'wilcox.test', ref.group = "BF_control",
#                      label = "p.signif", label.y = 1.2,
#                      size = 1.5, angle = 90, vjust = 0.5) +
#   scale_y_continuous(breaks = seq(0,2,0.05),
#                      minor_breaks = seq(0,2,0.01)) +
#   coord_cartesian(ylim = c(0.8,1.2)) +
#   scale_fill_discrete(guide = F) +
#   labs(x = 'Strains', y = 'Fitness') +
#   facet_wrap(.~density*replicate*hours*condition,
#              ncol = 4) +
#   theme_linedraw() +
#   theme(plot.title = element_blank(),
#         axis.title = element_text(size = titles),
#         axis.text = element_text(size = txt),
#         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#         legend.title = element_text(size = titles),
#         legend.text = element_text(size = txt),
#         legend.key.size = unit(3, "mm"),
#         legend.position = "bottom",
#         legend.box.spacing = unit(0.5,"mm"),
#         strip.text = element_text(size = txt,
#                                   margin = margin(0.1,0,0.1,0, "mm")))
# 
# ggsave(sprintf("%sDENSITY_COLONYSIZES.jpg",out_path), plot.den.cs,
#        height = one.c*14, width = two.c, units = 'mm',
#        dpi = 300)
# ggsave(sprintf("%sDENSITY_FITNESS.jpg",out_path), plot.den.fit,
#        height = one.c*14, width = two.c, units = 'mm',
#        dpi = 300)
# ggsave(sprintf("%sBOXPLOT_FITNESS.jpg",out_path), plot.box.fit,
#        height = one.c*14, width = two.c*1.5, units = 'mm',
#        dpi = 300)
# 
# ##### CALCULATING EFFECT SIZE (CLIFF'S DELTA)
# foa.es <- NULL
# i <- 1
# for (d in unique(foa.data$density)) {
#   for (rr in unique(foa.data$replicate[foa.data$density == d])) {
#     for (h in unique(foa.data$hours[foa.data$density == d & foa.data$replicate == rr])){
#       for (c in unique(foa.data$condition[foa.data$density == d & foa.data$replicate == rr & foa.data$hours == h])) {
#         cont_fit <- foa.data$fitness_clean[foa.data$density == d & foa.data$replicate == rr & foa.data$hours == h &
#                                              foa.data$condition == c & foa.data$orf_name == 'BF_control']
#         cont_mean <- median(foa.data$fitness_clean[foa.data$density == d & foa.data$replicate == rr & foa.data$hours == h &
#                                                      foa.data$condition == c & foa.data$orf_name == 'BF_control'], na.rm = T)
#         for (o in unique(foa.data$orf_name[foa.data$density == d & foa.data$replicate == rr & foa.data$hours == h & foa.data$condition == c &
#                                            foa.data$orf_name != 'BF_control' & foa.data$orf_name != 'BOR' & !is.na(foa.data$orf_name)])){
#           orf_fit <- foa.data$fitness_clean[foa.data$density == d & foa.data$replicate == rr & foa.data$hours == h &
#                                               foa.data$condition == c & foa.data$orf_name == o]
#           orf_mean <- median(foa.data$fitness_clean[foa.data$density == d & foa.data$replicate == rr & foa.data$hours == h &
#                                                       foa.data$condition == c & foa.data$orf_name == o], na.rm = T)
#           if (!is.na(orf_mean)) {
#             temp_cd <- cliff.delta(orf_fit, cont_fit, conf.level=.95,
#                                    use.unbiased=TRUE, use.normal=FALSE, return.dm=FALSE)
#             temp_p <- wilcox.test(orf_fit, cont_fit, alternative = 'two.sided', conf.int = 0.95)
#             foa.es$density[i] <- d
#             foa.es$replicate[i] <- rr
#             foa.es$hours[i] <- h
#             foa.es$condition[i] <- c
#             foa.es$group1[i] <- 'BF_control'
#             foa.es$group2[i] <- o
#             foa.es$cont_fit[i] <- cont_mean
#             foa.es$orf_fit[i] <- orf_mean
#             foa.es$cliff.delta[i] <- temp_cd$estimate
#             foa.es$magnitude[i] <- as.character(temp_cd$magnitude)
#             foa.es$effect_size[i] <- (orf_mean - cont_mean)/cont_mean * 100
#             foa.es$wilcox.p[i] <- temp_p$p.value
#             i <- i + 1
#           }
#         }
#       }
#       foa.es$p.adj[foa.es$density == d & foa.es$replicate == rr & foa.es$hours == h & foa.es$condition == c] <-
#         p.adjust(foa.es$wilcox.p[foa.es$density == d & foa.es$replicate == rr & foa.es$hours == h & foa.es$condition == c], method = 'BH')
#     } 
#   }
# }
# foa.es <- data.frame(foa.es)
# head(foa.es)
# 
# plot.es.mag <- ggplot(foa.es[foa.es$condition != 'GLY' & foa.es$hours > 0,]) +
#   geom_tile(aes(x = group1, y = group2, fill = magnitude), col = 'black') +
#   facet_wrap(.~density*replicate*hours*condition, ncol = 4) +
#   theme_linedraw() +
#   scale_fill_discrete(name = 'Magnitude') +
#   theme(plot.title = element_blank(),
#         axis.title = element_text(size = titles),
#         axis.text = element_text(size = txt),
#         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#         legend.title = element_text(size = titles),
#         legend.text = element_text(size = txt),
#         legend.key.size = unit(3, "mm"),
#         legend.position = "bottom",
#         legend.box.spacing = unit(0.5,"mm"),
#         strip.text = element_text(size = txt,
#                                   margin = margin(0.1,0,0.1,0, "mm")))
# 
# plot.es.val <- ggplot(foa.es[foa.es$condition != 'GLY' & foa.es$hours > 0,]) +
#   geom_tile(aes(x = group1, y = group2, fill = cliff.delta), col = 'black') +
#   facet_wrap(.~density*replicate*hours*condition, ncol = 4) +
#   theme_linedraw() +
#   scale_fill_continuous(name = "Cliff's Delta") +
#   theme(plot.title = element_blank(),
#         axis.title = element_text(size = titles),
#         axis.text = element_text(size = txt),
#         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#         legend.title = element_text(size = titles),
#         legend.text = element_text(size = txt),
#         legend.key.size = unit(3, "mm"),
#         legend.position = "bottom",
#         legend.box.spacing = unit(0.5,"mm"),
#         strip.text = element_text(size = txt,
#                                   margin = margin(0.1,0,0.1,0, "mm")))
# 
# ggsave(sprintf("%sES_MAGNITUDE.jpg",out_path), plot.es.mag,
#        height = one.c*3, width = two.c, units = 'mm',
#        dpi = 300)
# ggsave(sprintf("%sES_VALUE.jpg",out_path), plot.es.val,
#        height = one.c*3, width = two.c, units = 'mm',
#        dpi = 300)
# 
# ##### SAN DIEGO FOA
# foa.sd.es <- dbGetQuery(conn, "select a.orf_name, a.exp_id, a.hours, a.p, a.stat,
#                      b.mean, b.median, c.perc5, c.perc95
#                      from QVALUES_v2_CND1_SC a, STATS_v2_CND1_SC b, PERC_v2_CND1_SC c
#                      where a.orf_name = b.orf_name and a.exp_id = b.exp_id 
#                      and b.exp_id = c.exp_id
#                      and a.orf_name in
#                      (select orf_name from F28FUR_strainid2orf_name)
#                      and a.orf_name not in ('YHR021W-A')")
# 
# foa.sd.data <- dbGetQuery(conn, "select *
#                           from FITNESS_v2_CND1_SC a
#                           where a.orf_name in
#                           (select orf_name from F28FUR_strainid2orf_name)
#                           and a.orf_name not in ('YHR021W-A')
#                           order by orf_name desc") 
# 
# 
# ##### SAVING DATA
# save(foa.data, foa.es, foa.sd.data, foa.sd.es, file = sprintf('%sf28fu_data.RData',out_path))
# 

