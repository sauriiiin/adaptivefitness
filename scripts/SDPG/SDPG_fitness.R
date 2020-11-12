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

# ggplot(lidres[lidres$hours == 68 & lidres$density == 1536 & lidres$arm == 'SDA',]) +
#   geom_tile(aes(x = col, y = row, fill = fitness)) +
#   scale_y_reverse() +
#   # scale_fill_continuous(limits = c(0,2)) +
#   facet_wrap(.~density*arm*rep*hours)
# 
# ggplot(mcatres[mcatres$hours == 68 & mcatres$density == 1536 & mcatres$arm == 'SDA',]) +
#   geom_tile(aes(x = col, y = row, fill = fitness)) +
#   scale_y_reverse() +
#   # scale_fill_continuous(limits = c(0,2)) +
#   facet_wrap(.~density*arm*rep*hours)

##### CLEAN DATA
lidres <- lidres[lidres$orf_name != 'YHR021W-A' &
                   lidres$orf_name != 'BOR' &
                   lidres$orf_name != 'REF' &
                   !is.na(lidres$orf_name),]
# Removing plates with bad colony grids
lidres <- lidres[!(lidres$density == 6144 & lidres$arm == 'GLU' & lidres$rep != 'R3' & lidres$hours >= 20),]
lidres <- lidres[!(lidres$density == 6144 & lidres$arm == 'CAS' & lidres$rep == 'R2' & lidres$hours >= 20),]

# Removing outliers
for (d in unique(lidres$density)) {
  for (a in unique(lidres$arm[lidres$density == d])) {
    for (r in unique(lidres$rep[lidres$density == d & 
                                lidres$arm == a])) {
      for (h in unique(lidres$hours[lidres$density == d & 
                                    lidres$arm == a & 
                                    lidres$rep == r])) {
        for (o in unique(lidres$orf_name[lidres$density == d & 
                                     lidres$arm == a & 
                                     lidres$rep == r &
                                     lidres$hours == h])) {
          temp <- lidres$fitness[lidres$density == d & 
                           lidres$arm == a & 
                           lidres$rep == r &
                           lidres$hours == h &
                           lidres$orf_name == o]
          m <- median(temp, na.rm = T)
          madev <- mad(temp)
          ul <- m + 3*madev
          ll <- m - 3*madev
          
          lidres$fitness[lidres$density == d & 
                           lidres$arm == a & 
                           lidres$rep == r &
                           lidres$hours == h &
                           lidres$orf_name == o &
                           lidres$fitness > ul] <- NA
          lidres$fitness[lidres$density == d & 
                           lidres$arm == a & 
                           lidres$rep == r &
                           lidres$hours == h &
                           lidres$orf_name == o &
                           lidres$fitness < ll] <- NA
          
           
        }
        cont_median <- median(lidres$average[lidres$density == d & 
                                      lidres$arm == a & 
                                      lidres$rep == r &
                                      lidres$hours == h &
                                      lidres$orf_name == 'BF_control'], na.rm = T)
        lidres$ccs[lidres$density == d & 
                         lidres$arm == a & 
                         lidres$rep == r &
                         lidres$hours == h] <- lidres$cs_median[lidres$density == d & 
                                                            lidres$arm == a & 
                                                            lidres$rep == r &
                                                            lidres$hours == h] * cont_median
      }
    }
  }
}

# Saturation information
lidres$time <- NA
lidres$time[lidres$density == 1536 & lidres$hours == 48 & lidres$arm != 'SDA'] <- "Saturated"
lidres$time[lidres$density == 1536 & lidres$hours %in% c(72, 78) & lidres$arm == 'SDA'] <- "Saturated"
lidres$time[lidres$density == 6144 & lidres$hours == 24 & lidres$arm != 'SDA'] <- "Saturated"
lidres$time[lidres$density == 6144 & lidres$hours == 36 & lidres$arm == 'SDA'] <- "Saturated"
lidres$time[is.na(lidres$time)] <- "Other"

# Factorize the Experimental Arms
lidres$arm <- factor(lidres$arm, levels = c('GLU','CAS','SDA'))

##### GROWTH CURVES
lidsum <- lidres %>%
  group_by(density, arm, rep, hours, orf_name) %>%
  summarise(cs = mean(ccs, na.rm = T),fitness = mean(cs_mean, na.rm = T), p = mean(p, na.rm = T), es = mean(es, na.rm = T))
lidsum$arm <- factor(lidsum$arm, levels = c('GLU','CAS','SDA'))

fit_gc <- ggplot(lidsum[lidsum$hours > 0,],
       aes(x = hours, y = cs, col = orf_name)) +
  # geom_point() +
  stat_summary(fun = mean, geom = 'line') +
  labs(x = "Time (hours)",
       y = 'Normalized Colony Size') +
  scale_color_discrete(name = 'ORFs') +
  # scale_y_log10() +
  facet_wrap(.~density*arm*rep,
             scales = 'free',
             nrow = 3) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%sFITNESS_GC.jpg",out_path), fit_gc,
       height = 250, width = 500, units = 'mm',
       dpi = 300)

##### CORRELATION HEATMAPS
library(Hmisc)
library(ComplexHeatmap)
library(circlize)

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
cols$time[cols$density == 1536 & cols$hours >= 40 & cols$arm != 'SDA'] <- "Saturated"
cols$time[cols$density == 1536 & cols$hours >= 60 & cols$arm == 'SDA'] <- "Saturated"
cols$time[cols$density == 6144 & cols$hours >= 20 & cols$arm != 'SDA'] <- "Saturated"
cols$time[cols$density == 6144 & cols$hours >= 30 & cols$arm == 'SDA'] <- "Saturated"
cols$time[is.na(cols$time)] <- "Other"

cols$saturation <- cols$hours
cols$saturation[cols$density == 1536 & cols$arm != 'SDA'] <- cols$saturation[cols$density == 1536 & cols$arm != 'SDA']/48
cols$saturation[cols$density == 1536 & cols$arm == 'SDA'] <- cols$saturation[cols$density == 1536 & cols$arm == 'SDA']/72
cols$saturation[cols$density == 6144 & cols$arm != 'SDA'] <- cols$saturation[cols$density == 6144 & cols$arm != 'SDA']/24
cols$saturation[cols$density == 6144 & cols$arm == 'SDA'] <- cols$saturation[cols$density == 6144 & cols$arm == 'SDA']/36

col_fun = colorRamp2(c(0.3, 0.5, 1), c("white", "#FF5252", "#512DA8"))
ha <- HeatmapAnnotation(
  Density = cols$density,
  Condition = cols$arm,
  Replicate = cols$replicate,
  Saturation = cols$saturation,
  col = list(Density = c("1536" = "#607D8B", "6144" = "#FFA000"),
             Condition = c("GLU" = "#B3E5FC", "CAS" = "#448AFF", "SDA" = "#303F9F"),
             Replicate = c("R1" = "#E1BEE7", "R2" = "#E040FB", "R3" = "#7B1FA2"),
             Saturation = col_fun),
  gp = gpar(col = "black", lwd = 0.2)
)

ra <- rowAnnotation(
  Density = cols$density,
  Condition = cols$arm,
  Replicate = cols$replicate,
  Saturation = cols$saturation,
  col = list(Density = c("1536" = "#607D8B", "6144" = "#FFA000"),
             Condition = c("GLU" = "#B3E5FC", "CAS" = "#448AFF", "SDA" = "#303F9F"),
             Replicate = c("R1" = "#E1BEE7", "R2" = "#E040FB", "R3" = "#7B1FA2"),
             Saturation = col_fun),
  gp = gpar(col = "black", lwd = 0.2),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

mat <- cor(cor_dat, use = 'complete.obs', method = 'pearson')

ht <- Heatmap(mat,
              name = "Correlation",
              column_title = "The Famous 28 Heatmap",
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%0.2f", mat[i, j]), x, y, gp = gpar(fontsize = 2))
                },
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
draw(ra + ht)

##### FITNESS BOX PLOTS
fit_box <- ggplot(lidres[lidres$time == 'Saturated' &
                              lidres$orf_name != 'YHR021W-A' &
                              lidres$orf_name != 'BOR' &
                              !is.na(lidres$orf_name),],
                     aes(x = orf_name, y = fitness, fill = orf_name)) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(method = 'wilcox.test', ref.group = "BF_control",
                     label = "p.signif", label.y = 1.2,
                     size = 1.5, angle = 90, vjust = 0.5) +
  labs(x = "ORFs",
       y = 'Fitness') +
  scale_fill_discrete(name = 'ORFs') +
  facet_wrap(.~arm*density*rep*hours,
             nrow = 3) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 90),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0.8, 1.2))
ggsave(sprintf("%sFITNESS_BOX.jpg",out_path), fit_box,
       height = 250, width = 500, units = 'mm',
       dpi = 300)

fit_box2 <- ggplot(lidres[lidres$time == "Saturated",],
                   aes(x = orf_name, y = fitness)) +
  geom_boxplot(aes(fill = rep),
               outlier.shape = NA) +
  # stat_compare_means(method = 'wilcox',
  #                    ref.group = "BF_control",
  #                    label = "p.signif", label.y = 1.2,
  #                    size = 1.5, angle = 90, vjust = 0.5) +
  stat_compare_means(aes(group = rep),
                     # comparisons = my_comparisons,
                     # method = 'wilcox',
                     label = "p.signif",
                     label.y = 1.18,
                     size = 1.5, hide.ns = T) +
  stat_compare_means(method = 'wilcox',
                     ref.group = 'BF_control',
                     label = "p.signif",
                     label.y = 1.2,
                     size = 1.5, hide.ns = T) +
  coord_cartesian(ylim = c(0.8, 1.2)) +
  facet_wrap(.~arm*density,
             nrow = 3) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 90),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%sFITNESS_BOX2.jpg",out_path), fit_box2,
       height = 250, width = 300, units = 'mm',
       dpi = 300)

##### FITNESS DENSITY PLOTS
fit_den <- ggplot(lidres[lidres$time == 'Saturated' &
                lidres$orf_name != 'YHR021W-A' &
                lidres$orf_name != 'BOR' &
                !is.na(lidres$orf_name),],
       aes(x = fitness, y = orf_name, group = orf_name, fill = orf_name)) +
  geom_density_ridges(quantile_lines = TRUE,
                      scale = 3, alpha = 0.8, size = 0.3,
                      vline_size = 0.2, vline_color = "black") +
  labs(x = "Fitness",
       y = 'ORFs') +
  scale_fill_discrete(name = 'ORFs') +
  facet_wrap(.~arm*density*rep*hours,
             nrow = 3) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(xlim = c(0.8, 1.2))
ggsave(sprintf("%sFITNESS_DEN.jpg",out_path), fit_den,
       height = 300, width = 500, units = 'mm',
       dpi = 300)

##### STATS
### ES AND PVALUES
library(effsize)
fit_stat2 <- data.frame(compare_means(fitness ~ orf_name, lidres[lidres$hours >= 8,],
                                      method = "wilcox.test", paired = FALSE,
                                      p.adjust.method = 'bonferroni',
                                      group.by = c("density", "arm", "rep", "hours"), ref.group = "BF_control"))
head(fit_stat2)

for (d in unique(fit_stat2$density)) {
  for (a in unique(fit_stat2$arm[fit_stat2$density == d])) {
    for (r in unique(fit_stat2$rep[fit_stat2$density == d & fit_stat2$arm == a])) {
      for (h in unique(fit_stat2$hours[fit_stat2$density == d & fit_stat2$arm == a & fit_stat2$rep == r])) {
        cont_fit <- lidres$fitness[lidres$density == d & lidres$arm == a & lidres$rep == r &
                                     lidres$hours == h & lidres$orf_name == 'BF_control']
        cont_mean <- median(lidres$fitness[lidres$density == d & lidres$arm == a & lidres$rep == r &
                                             lidres$hours == h & lidres$orf_name == 'BF_control'], na.rm = T)
        for (o in unique(fit_stat2$group2[fit_stat2$density == d & fit_stat2$arm == a & fit_stat2$rep == r &
                                          fit_stat2$hours == h & fit_stat2$group2 != 'BF_control'])) {
          orf_fit <- lidres$fitness[lidres$density == d & lidres$arm == a & lidres$rep == r &
                                      lidres$hours == h & lidres$orf_name == o]
          orf_mean <- median(lidres$fitness[lidres$density == d & lidres$arm == a & lidres$rep == r &
                                              lidres$hours == h & lidres$orf_name == o], na.rm = T)
          temp_cd <- cliff.delta(orf_fit, cont_fit, conf.level=.95,
                                 use.unbiased=TRUE, use.normal=FALSE, return.dm=FALSE)
          fit_stat2$empirical_p[fit_stat2$density == d & fit_stat2$arm == a & fit_stat2$rep == r &
                                  fit_stat2$hours == h & fit_stat2$group2 == o] <- lidres$p[lidres$density == d & lidres$arm == a & lidres$rep == r &
                                                                                              lidres$hours == h & lidres$orf_name == o][1]
          fit_stat2$cliff.delta[fit_stat2$density == d & fit_stat2$arm == a & fit_stat2$rep == r &
                                  fit_stat2$hours == h & fit_stat2$group2 == o] <- temp_cd$estimate
          fit_stat2$magnitude[fit_stat2$density == d & fit_stat2$arm == a & fit_stat2$rep == r &
                                fit_stat2$hours == h & fit_stat2$group2 == o] <- as.character(temp_cd$magnitude)
          fit_stat2$effect_size[fit_stat2$density == d & fit_stat2$arm == a & fit_stat2$rep == r &
                                  fit_stat2$hours == h & fit_stat2$group2 == o] <- (orf_mean - cont_mean)/cont_mean * 100
          
        }
      }
    }
  }
}

fit_stat2$phenotype[fit_stat2$p.signif == 'ns'] <- 'Neutral'
fit_stat2$phenotype[fit_stat2$p.signif != 'ns' & fit_stat2$cliff.delta < 0] <- 'Deleterious'
fit_stat2$phenotype[fit_stat2$p.signif != 'ns' & fit_stat2$cliff.delta > 0] <- 'Beneficial'

fit_stat2$empirical_phenotype[fit_stat2$empirical_p > 0.05] <- 'Neutral'
fit_stat2$empirical_phenotype[fit_stat2$empirical_p <= 0.05 & fit_stat2$cliff.delta < 0] <- 'Deleterious'
fit_stat2$empirical_phenotype[fit_stat2$empirical_p <= 0.05 & fit_stat2$cliff.delta > 0] <- 'Beneficial'

fit_stat2$saturation[fit_stat2$density == 1536 & fit_stat2$hours >= 40 & fit_stat2$arm != 'SDA'] <- "Saturated"
fit_stat2$saturation[fit_stat2$density == 1536 & fit_stat2$hours >= 60 & fit_stat2$arm == 'SDA'] <- "Saturated"
fit_stat2$saturation[fit_stat2$density == 6144 & fit_stat2$hours >= 20 & fit_stat2$arm != 'SDA'] <- "Saturated"
fit_stat2$saturation[fit_stat2$density == 6144 & fit_stat2$hours >= 30 & fit_stat2$arm == 'SDA'] <- "Saturated"
fit_stat2$saturation[is.na(fit_stat2$saturation)] <- "Other"

fit_stat2$name <- paste(fit_stat2$density, fit_stat2$arm, fit_stat2$rep, fit_stat2$hours, sep = '_')
fit_stat <- fit_stat2[fit_stat2$saturation == 'Saturated',]

# plyr::count(fit_stat, vars = c("arm", "density", "hours", "group2", "phenotype"))
head(fit_stat)

fit_pheno <- ggplot(fit_stat,
       aes(x = group2, fill = phenotype)) +
  geom_bar() +
  # geom_text(stat = 'count', aes(label = cliff.delta)) +
  labs(x = "ORFs",
       y = "Count") +
  facet_wrap(.~arm*density, ncol = 2) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Deleterious','Neutral','Beneficial'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107')) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 90),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%sFITNESS_PHENO.jpg",out_path), fit_pheno,
       height = 250, width = 300, units = 'mm',
       dpi = 300)


fit_es <- ggplot(fit_stat,
       aes(x = effect_size, y = group2, col = phenotype)) +
  geom_vline(xintercept = c(-10,-5,-1,1,5,10),
             linetype = 'dashed', lwd = 0.5) +
  geom_point(size = 3) +
  labs(x = 'Effect Size (%)\n(mean ORF - mean BF_control Fitness)/mean BF_control Fitness',
       y = 'ORFs') +
  scale_x_continuous(breaks = seq(-20,20,2.5)) +
  geom_text(aes(label = rep), col = 'white', size = 1) +
  facet_wrap(.~arm*density, ncol = 2) +
  scale_color_manual(name = 'Phenotype',
                    breaks = c('Deleterious','Neutral','Beneficial'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107')) +
  coord_cartesian(xlim = c(-15, 15)) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%sFITNESS_ES.jpg",out_path), fit_es,
       height = 250, width = 300, units = 'mm',
       dpi = 300)

# Romano 2006 == |d|<0.147 "negligible", |d|<0.33 "small", |d|<0.474 "medium", otherwise "large"
fit_cd <- ggplot(fit_stat,
                 aes(x = cliff.delta, y = group2, col = phenotype)) +
  geom_vline(xintercept = c(-0.474,-0.33,-0.147,0.147,0.33,0.474),
             linetype = 'dashed', lwd = 0.5) +
  geom_point(size = 3) +
  labs(x = "Cliff's Delta",
       y = "ORFs") +
  scale_x_continuous(breaks = seq(-1,1,0.1)) +
  geom_text(aes(label = rep), col = 'white', size = 1) +
  facet_wrap(.~arm*density, ncol = 2) +
  scale_color_manual(name = 'Phenotype',
                     breaks = c('Deleterious','Neutral','Beneficial'),
                     values = c('Deleterious'='#3F51B5',
                                'Neutral'='#212121',
                                'Beneficial'='#FFC107')) +
  coord_cartesian(xlim = c(-1, 1)) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%sFITNESS_CLIFF.jpg",out_path), fit_cd,
       height = 250, width = 300, units = 'mm',
       dpi = 300)

# EMPIRICAL PHENOTYPE
fit_ep <- ggplot(fit_stat,
                 aes(x = effect_size, y = group2, col = empirical_phenotype)) +
  geom_vline(xintercept = c(-10,-5,-1,1,5,10),
             linetype = 'dashed', lwd = 0.5) +
  geom_point(size = 3) +
  labs(x = 'Effect Size (%)\n(mean ORF - mean BF_control Fitness)/mean BF_control Fitness',
       y = 'ORFs') +
  scale_x_continuous(breaks = seq(-20,20,2.5)) +
  geom_text(aes(label = rep), col = 'white', size = 1) +
  facet_wrap(.~arm*density, ncol = 2) +
  scale_color_manual(name = 'Empirical Phenotype',
                     breaks = c('Deleterious','Neutral','Beneficial'),
                     values = c('Deleterious'='#3F51B5',
                                'Neutral'='#212121',
                                'Beneficial'='#FFC107')) +
  coord_cartesian(xlim = c(-15, 15)) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%sFITNESS_ES_EP.jpg",out_path), fit_ep,
       height = 250, width = 300, units = 'mm',
       dpi = 300)

##### ANOVA
# library(AICcmodavg)
# 
# res.aov <- aov(fitness ~ density*arm*rep*orf_name, data = lidres[lidres$time == 'Saturated',])
# summary(res.aov)
# 
# a1 <- aov(fitness ~ density + arm + rep + orf_name, data = lidres[lidres$time == 'Saturated',])
# a2 <- aov(fitness ~ density + arm + rep * orf_name, data = lidres[lidres$time == 'Saturated',])
# a3 <- aov(fitness ~ density * arm + rep * orf_name, data = lidres[lidres$time == 'Saturated',])
# a4 <- aov(fitness ~ density * arm * rep * orf_name, data = lidres[lidres$time == 'Saturated',])
# 
# 
# model.set <- list(a1, a2, a3, a4)
# model.names <- c("no int", "RO", "DA+RO", "DARO")
# 
# aictab(model.set, modnames = model.names)
# 
# summary(a4)

##### CLIFFS DELTA ALL THE WAY
### SDPG_pairwise.R

### CLEAN RESULTS FROM ABOVE SCRIPT
# library(qvalue)
# load(sprintf('%spair_comp.RDATA',out_path))
# pair_comp <- data.frame(pair_comp)
# 
# pair_comp$g1_density <- as.numeric(as.character(pair_comp$g1_density))
# pair_comp$g1_density <- factor(pair_comp$g1_density, levels = c(1536, 6144))
# pair_comp$g2_density <- as.numeric(as.character(pair_comp$g2_density))
# pair_comp$g2_density <- factor(pair_comp$g2_density, levels = c(1536, 6144))
# pair_comp$pvalue <- as.numeric(as.character(pair_comp$pvalue))
# pair_comp$cliffdelta <- as.numeric(as.character(pair_comp$cliffdelta))
# 
# pair_comp$g1 <- sprintf('%d_%s_%s_%s',pair_comp$g1_density,pair_comp$g1_arm,pair_comp$g1_rep,pair_comp$g1_orf_name)
# pair_comp$g2 <- sprintf('%d_%s_%s_%s',pair_comp$g2_density,pair_comp$g2_arm,pair_comp$g2_rep,pair_comp$g2_orf_name)
# 
# temp <- qvalue(pair_comp$pvalue)
# pair_comp$qvalue <- temp$qvalues
# 
# pair_comp$phenotype[pair_comp$g1_orf_name == 'BF_control' & pair_comp$pvalue > 0.05] <- 'Neutral'
# pair_comp$phenotype[pair_comp$g1_orf_name == 'BF_control' & pair_comp$pvalue <= 0.05 & pair_comp$cliffdelta < 0] <- 'Deleterious'
# pair_comp$phenotype[pair_comp$g1_orf_name == 'BF_control' & pair_comp$pvalue <= 0.05 & pair_comp$cliffdelta > 0] <- 'Beneficial'
# pair_comp$phenotype[is.na(pair_comp$phenotype)] <- 'NA'
# pair_comp$phenotype <- factor(pair_comp$phenotype, levels = c('Deleterious','Neutral','Beneficial','NA'))
# save(pair_comp, file = sprintf('%spair_comp.RDATA',out_path))

##### USE CLEANED PAIR_COMP ONLY
# load(sprintf('%spair_comp.RDATA',out_path))
# source('/home/sbp29/R/Projects/adaptivefitness/scripts/SDPG/make_comp_plot.R')
# head(pair_comp)
# 
# g1_orf_name <- c('BF_control')
# g1_density <- c(1536)
# g1_arm <- c('GLU')
# g1_rep <- c('R1','R2','R3')
# 
# g2_orf_name <- c('BF_control')
# g2_density <- c(6144)
# g2_arm <- c('GLU')
# g2_rep <- c('R1','R2','R3')
# 
# make_comp_plot(pair_comp,
#                g1_orf_name, g1_density, g1_arm, g1_rep,
#                g2_orf_name, g2_density, g2_arm, g2_rep)



#####
# my_comparisons <- list(c(1:1))

#####
# library(rstatix)
# 
# res.aov <- anova_test(
#   data = lidres[lidres$time == 'Saturated' & lidres$density == 1536,],
#   dv = fitness, wid = orf_name, within = rep
# )
# get_anova_table(res.aov)


