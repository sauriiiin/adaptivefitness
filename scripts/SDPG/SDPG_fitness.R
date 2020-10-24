##### SDPG FITNESS RESULTS
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 10/05/2020

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
library(ggrepel)
library(reshape2)
library(RMariaDB)

source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")
out_path = 'figs/SDPG/';

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 7
txt <- 7
lbls <- 9

# ##### GATHER DATA
# expt <- c('GLU','CAS','SDA')
# reps <- c('R1','R2','R3')
# dens <- c(1536, 6144)
# 
# lidres <- NULL
# mcatres <- NULL
# for (d in dens) {
#   for (e in expt) {
#     for (r in reps) {
#       temp1 <- dbGetQuery(conn,
#                               sprintf('select c.*, d.p, d.es from
#                                       (select e.*, f.N, f.cs_mean, f.cs_median, f.cs_std
#                                       from 
#                                       (select b.*, a.strain_id, a.orf_name, a.hours, a.average, a.fitness
#                                       from
#                                       SDPG_%s_FS_%s_%d_FITNESS a, SDPG_pos2coor b 
#                                       where a.pos = b.pos) e
#                                       left join
#                                       SDPG_%s_FS_%s_%d_FITNESS_STATS f
#                                       on
#                                       e.orf_name = f.orf_name and e.hours = f.hours) c
#                                       left join
#                                       SDPG_%s_FS_%s_%d_PVALUE d
#                                       on c.hours = d.hours and c.strain_id = d.strain_id
#                                       order by c.hours, c.density, c.plate, c.col, c.row',
#                                       e,r,d,e,r,d,e,r,d))
#       temp1$arm <- e
#       temp1$rep <- r
#       for (h in unique(sort(temp1$hours))) {
#         temp1$ccs[temp1$hours == h] <- temp1$fitness[temp1$hours == h] *
#           median(temp1$average[temp1$orf_name == 'BF_control' & temp1$hours == h],na.rm = T)
#       }
# 
#       temp2 <- dbGetQuery(conn,
#                           sprintf('select c.*, d.p, d.es from
#                                   (select e.*, f.N, f.cs_mean, f.cs_median, f.cs_std
#                                   from 
#                                   (select b.*, a.strain_id, a.orf_name, a.hours, a.average, a.fitness
#                                   from
#                                   SDPG_%s_FS_%s_%d_MCAT_FITNESS a, SDPG_pos2coor b 
#                                   where a.pos = b.pos) e
#                                   left join
#                                   SDPG_%s_FS_%s_%d_MCAT_FITNESS_STATS f
#                                   on
#                                   e.orf_name = f.orf_name and e.hours = f.hours) c
#                                   left join
#                                   SDPG_%s_FS_%s_%d_MCAT_PVALUE d
#                                   on c.hours = d.hours and c.strain_id = d.strain_id
#                                   order by c.hours, c.density, c.plate, c.col, c.row',
#                                   e,r,d,e,r,d,e,r,d))
#       temp2$arm <- e
#       temp2$rep <- r
#       for (h in unique(sort(temp2$hours))) {
#         temp2$ccs[temp2$hours == h] <- temp2$fitness[temp2$hours == h] *
#           median(temp2$average[temp2$orf_name == 'BF_control' & temp2$hours == h],na.rm = T)
#       }
#       
#       lidres <- rbind(lidres, temp1)
#       mcatres <- rbind(mcatres, temp2)
#     }
#   }
# }
# save(lidres, mcatres, file = sprintf('%s/results.RDATA',out_path))

##### PLATE HEATMAP
load(sprintf('%s/results.RDATA',out_path))

# ggplot(lidres[lidres$hours == 24 & lidres$density == 1536 & lidres$arm == 'GLU',]) +
#   geom_violin(aes(x = orf_name, y = fitness),
#               trim = T) +
#   coord_cartesian(ylim = c(0,2)) +
#   facet_wrap(.~density*arm*rep*hours)

ggplot(lidres[lidres$hours == 68 & lidres$density == 1536 & lidres$arm == 'SDA',]) +
  geom_tile(aes(x = col, y = row, fill = fitness)) +
  scale_y_reverse() +
  # scale_fill_continuous(limits = c(0,2)) +
  facet_wrap(.~density*arm*rep*hours)

ggplot(mcatres[mcatres$hours == 68 & mcatres$density == 1536 & mcatres$arm == 'SDA',]) +
  geom_tile(aes(x = col, y = row, fill = fitness)) +
  scale_y_reverse() +
  # scale_fill_continuous(limits = c(0,2)) +
  facet_wrap(.~density*arm*rep*hours)


##### GROWTH CURVES
lidsum <- lidres %>%
  group_by(density, arm, rep, hours, orf_name) %>%
  summarise(cs = mean(ccs),fitness = mean(cs_mean), p = mean(p), es = mean(es))

ggplot(lidsum,
       aes(x = hours, y = cs, col = orf_name)) +
  # geom_point() +
  stat_summary(fun = mean, geom = 'line') +
  facet_wrap(.~density*arm*rep)

##### CORRELATION HEATMAPS
library(Hmisc)
library(ComplexHeatmap)
library(circlize)

head(lidsum)
lidsum <- lidsum[!(lidsum$density == 6144 & lidsum$arm == 'GLU' & lidsum$rep != 'R3' & lidsum$hours >= 20),]
lidsum <- lidsum[!(lidsum$density == 6144 & lidsum$arm == 'CAS' & lidsum$rep == 'R2' & lidsum$hours >= 20),]

cor_dat <- NULL
cols <- data.frame()
for (d in unique(lidsum$density)) {
  for (a in unique(lidsum$arm[lidsum$density == d])) {
    for (r in unique(lidsum$rep[lidsum$density == d & 
                                lidsum$arm == a])) {
      for (h in unique(lidsum$hours[lidsum$density == d & 
                                    lidsum$arm == a & 
                                    lidsum$rep == r])) {
        if (h > 10) {
          cor_dat <- cbind(cor_dat, t(t(lidsum$fitness[lidsum$density == d & 
                                                         lidsum$arm == a & 
                                                         lidsum$rep == r & 
                                                         lidsum$hours == h])))
          cols <- rbind(cols,
                        data.frame(title = sprintf('%d_%s_%s_%0.0f',d,a,r,h),
                          density = d,
                          arm = a,
                          replicate = r,
                          hours = h))
        }
      }
    }
  }
}
colnames(cor_dat) <- cols$title
cols$time[cols$density == 1536 & cols$hours == 48 & cols$arm != 'SDA'] <- "Saturated"
cols$time[cols$density == 1536 & cols$hours == 68 & cols$arm == 'SDA'] <- "Saturated"
cols$time[cols$density == 6144 & cols$hours == 24 & cols$arm != 'SDA'] <- "Saturated"
cols$time[cols$density == 6144 & cols$hours == 36 & cols$arm == 'SDA'] <- "Saturated"
cols$time[is.na(cols$time)] <- "Growing"

ha <- HeatmapAnnotation(
  Density = cols$density,
  Condition = cols$arm,
  Replicate = cols$replicate,
  Time = cols$time,
  col = list(Density = c("1536" = "#607D8B", "6144" = "#FFA000"),
             Condition = c("GLU" = "#B3E5FC", "CAS" = "#448AFF", "SDA" = "#303F9F"),
             Replicate = c("R1" = "#E1BEE7", "R2" = "#E040FB", "R3" = "#7B1FA2"),
             Time = c("Saturated" = "#FF5252", "Growing" = "#512DA8")),
  gp = gpar(col = "black", lwd = 0.2)
)
col_fun = colorRamp2(c(0, 5, 10), c("blue", "white", "red"))

mat <- cor(cor_dat, use = 'complete.obs', method = 'pearson')

ht <- Heatmap(mat,
              name = "Correlation",
              column_title = "The Famous 28 Heatmap",
              column_title_gp = gpar(fontsize = 12, fontface = "bold"),
              show_row_dend = F,
              rect_gp = gpar(col = "black", lwd = 0.2),
              show_row_names = F,
              show_column_names = F,
              # row_km = 3,
              # column_km = 3,
              # show_parent_dend_line = T,
              border = T,
              top_annotation = ha)
draw(ht)

##### FITNESS VIOLIN PLOTS
ggplot(lidres[lidres$hours > 10,]) +
  geom_violin(aes(x = orf_name, y = fitness)) +
  facet_wrap(.~density*arm*rep*hours)
