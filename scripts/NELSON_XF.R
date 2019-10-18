##### NELSON_XF
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 10/17/2019

##### INITIALIZE
library(RMariaDB)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggExtra)
library(grid)
library(tidyverse)
library(ggpubr)
library(stringr)
source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

out_path = 'figs/nelson/';
expt_name <- 'XF'


alldat = dbGetQuery(conn, 'select *
              from XF_1536_FITNESS_STATS')

ref.pos <- dbGetQuery(conn, 'select pos from XF_pos2orf_name
                      where orf_name = "BY4741" and pos < 1000')
ref.all <- dbGetQuery(conn, 'select * from XF_1536_FITNESS
                      where orf_name = "BY4741"')

temp <- c(14000,24000,34000,44000)
ref.fit <- NULL
i = 1
for (pos in ref.pos$pos) {
  ref.fit$cs_mean[i] <- mean(ref.all$fitness[ref.all$pos %in% c(temp + pos)], na.rm = T)
  ref.fit$cs_median[i] <- median(ref.all$fitness[ref.all$pos %in% c(temp + pos)], na.rm = T)
  i = i + 1
}  
ref.fit <- data.frame(ref.fit)

ggplot(alldat) +
  geom_line(aes(x = cs_median, col = orf_name2), stat = 'density') +
  theme_linedraw()

#              draw_quantiles = c(0.25, 0.5, 0.75)
ggplot() +
  geom_boxplot(data = ref.fit,
              aes(x = 'BY4741', y = cs_median, fill = 'BY4741'),
              width = 0.4) +
  geom_violin(data = ref.fit,
              aes(x = 'BY4741', y = cs_median, fill = 'BY4741'), alpha = 0.6) +
  geom_boxplot(data = alldat[alldat$orf_name != 'BY4741',],
              aes(x = orf_name2, y = cs_median, fill = orf_name2),
              width = 0.4) +
  geom_violin(data = alldat[alldat$orf_name != 'BY4741',],
              aes(x = orf_name2, y = cs_median, fill = orf_name2), alpha = 0.6) +
  labs(x = '',
       y = 'Relative Fitness') +
  scale_fill_discrete(name = 'Strains',
                      breaks = c('BY4741','dYBR','dbck','A25','A26','A27','A28'),
                      labels = c('BY4741','dYBR','dBCK','A25','A26','A27','A28')) +
  scale_x_discrete(limits = c('BY4741','dYBR','dbck','A25','A26','A27','A28'),
                   labels = c('BY4741','dYBR','dBCK','A25','A26','A27','A28')) +
  theme_linedraw() +
  theme(axis.title.x = element_blank())
ggsave(sprintf("%sfitness.jpg",out_path),
       width = 10, height = 10,
       dpi = 300)


