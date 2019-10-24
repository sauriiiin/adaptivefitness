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
              from XF_1536_FITNESS_STATS where hours = 26')

ref.pos <- dbGetQuery(conn, 'select pos from XF_pos2orf_name
                      where orf_name = "BY4741" and pos < 1000')
ref.all <- dbGetQuery(conn, 'select * from XF_1536_FITNESS
                      where orf_name = "BY4741" and hours = 26')

temp <- c(14000,24000,34000,44000)
ref.fit <- NULL
i = 1
for (pos in ref.pos$pos) {
  ref.fit$cs_mean[i] <- mean(ref.all$fitness[ref.all$pos %in% c(temp + pos)], na.rm = T)
  ref.fit$cs_median[i] <- median(ref.all$fitness[ref.all$pos %in% c(temp + pos)], na.rm = T)
  i = i + 1
}  
ref.fit <- data.frame(ref.fit)
ref.fit$orf_name <- 'BY4741'

for (orf in unique(alldat$orf_name2)) {
  temp <- NULL
  temp$cs_mean <- alldat$cs_mean[alldat$orf_name2 == orf]
  temp$cs_median <- alldat$cs_median[alldat$orf_name2 == orf]
  temp$orf_name <- orf
  temp <- data.frame(temp)
  
  ref.fit <- rbind(ref.fit, temp)
}

ref.med <- median(ref.fit$cs_median[ref.fit$orf_name == 'BY4741'], na.rm = T)
#              draw_quantiles = c(0.25, 0.5, 0.75)
ggplot() +
  geom_hline(yintercept = ref.med, linetype = 'dashed', col = 'red') +
  geom_boxplot(data = ref.fit,
              aes(x = orf_name, y = cs_median, fill = orf_name),
              width = 0.4) +
  geom_violin(data = ref.fit,
              aes(x = orf_name, y = cs_median, fill = orf_name), alpha = 0.6) +
  # geom_boxplot(data = alldat[alldat$orf_name != 'BY4741',],
  #             aes(x = orf_name2, y = cs_median, fill = orf_name2),
  #             width = 0.4) +
  # geom_violin(data = alldat[alldat$orf_name != 'BY4741',],
  #             aes(x = orf_name2, y = cs_median, fill = orf_name2), alpha = 0.6) +
  labs(x = '',
       y = 'Relative Fitness') +
  scale_fill_discrete(name = '',
                      breaks = c('A25','A26','A27','A28','dbck','dYBR','BY4741'),
                      labels = c('A25','A26','A27','A28','dBCK1','dYBR','BY4741')) +
  scale_x_discrete(limits = c('A25','A26','A27','A28','dbck','dYBR','BY4741'),
                   labels = c('A25','A26','A27','A28','dBCK1','dYBR','BY4741')) +
  theme_linedraw() +
  theme(axis.title = element_blank()) +
  coord_cartesian(ylim = c(0.75,1.25))
ggsave(sprintf("%sfitness.jpg",out_path),
       width = 7, height = 7,
       dpi = 300)


