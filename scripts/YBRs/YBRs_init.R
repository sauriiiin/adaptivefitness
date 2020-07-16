library(ggplot2)
library(gridExtra)
library(grid)
library(tidyverse)
library(ggpubr)
library(stringr)
library(egg)
library(zoo)
library(ggrepel)
library(plotly)
library(scales)
library(reshape2)
library(RMariaDB)

source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

out_path <- 'figs/ybr/';

remove_outliers <- function(x) {
  qnt <- quantile(x, probs=c(.05, .95), na.rm = TRUE)
  y <- x
  y[x < qnt[1]] <- NA
  y[x > qnt[2]] <- NA
  y
}

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 7
txt <- 5
lbls <- 9

##### PLATEMAPS
ybrs_pm <- dbGetQuery(conn, 'select * from YBRs_pos2coor a, YBRs_pos2orf_name b
                  where a.pos = b.pos')

ggplot(ybrs_pm[ybrs_pm$density == 1536,],
       aes(x = col, y = row, fill = orf_name)) +
  geom_tile(col = 'black') +
  scale_y_reverse() +
  facet_wrap(.~density*plate, 
             scale = 'free',
             ncol = 2) +
  theme_linedraw()

plyr::count(ybrs_pm, vars = c('density','orf_name'))

##### FITNESS RESULTS
fitness <- dbGetQuery(conn, 'select *
                   from YBRs_YPD_1536_FITNESS
                   union
                   select *
                   from YBRs_YPD_6144_FITNESS')
fitness$density[fitness$pos <= 999999] <- 1536
fitness$density[fitness$pos >= 999999] <- 6144
plyr::count(fitness, vars = 'density')

fitness$orf_name <- factor(fitness$orf_name, levels = c('REF',
                                                        'FY4',
                                                        'BY4741',
                                                        'BY4741_KAN',
                                                        'BY4742',
                                                        'BY4742_KAN',
                                                        'CRISPY',
                                                        'ATG',
                                                        'PAM',
                                                        'TATA',
                                                        'C_STOP',
                                                        'NC_STOP',
                                                        'NC_STOP2',
                                                        'BOR'))


# fit.stats <- dbGetQuery(conn, 'select a.orf_name, a.cs_mean 1536_mean, a.cs_median 1536_median, a.cs_std 1536_std,
#                         b.cs_mean 6144_mean, b.cs_median 6144_median, b.cs_std 6144_std
#                         from YBRs_YPD_1536_FITNESS_STATS a, YBRs_YPD_6144_FITNESS_STATS b
#                         where a.orf_name = b.orf_name')
# 
# fit.stats$orf_name <- factor(fit.stats$orf_name, levels = c('REF',
#                                                         'FY4',
#                                                         'BY4741',
#                                                         'BY4741_KAN',
#                                                         'BY4742',
#                                                         'BY4742_KAN',
#                                                         'CRISPY',
#                                                         'ATG',
#                                                         'PAM',
#                                                         'TATA',
#                                                         'C_STOP',
#                                                         'NC_STOP',
#                                                         'NC_STOP2',
#                                                         'BOR'))

## CLEAN fitness
for (d in unique(fitness$density)) {
  for (h in unique(fitness$hours[fitness$density == d])) {
    for (o in unique(fitness$orf_name[fitness$density == d & fitness$hours == h])) {
      fitness$fitness[fitness$density == d & fitness$hours == h & fitness$orf_name == o] <-
        remove_outliers(fitness$fitness[fitness$density == d & fitness$hours == h & fitness$orf_name == o])
    }
  }
}

## PLOT fitness
ggplot(fitness[fitness$orf_name != 'BOR' & fitness$hours > 5,],
       aes(x = orf_name, y = fitness, fill = orf_name)) +
  geom_hline(yintercept = 1, col = 'red') +
  geom_jitter(size = 0.05, alpha = 0.7) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5) +
  labs(title = 'YPD') +
  scale_fill_discrete(guide = F) +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.5)) +
  facet_wrap(.~density*hours) +
  theme(plot.title = element_text(size = titles,
                                  face = 'bold',
                                  hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        strip.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"))

ggsave(sprintf("%sFITNESS_YPD.jpg",out_path),
       height = 150, width = two.c, units = 'mm',
       dpi = 300)

# ggplot(fit.stats[fit.stats$orf_name != 'BOR',],
#        aes(x = `1536_median`, y = `6144_median`, col = orf_name)) +
#   geom_point() +
#   geom_text(aes(label = orf_name))

##### PHENOTYPE


