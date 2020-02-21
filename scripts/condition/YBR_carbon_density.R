##### YBR SCREENS
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 01/21/2020
##### Screens were done in November 2018

##### INITIALIZE
library(RMariaDB)
library(ggplot2)
library(ggExtra)
library(gridExtra)
library(ggpubr)
library(reshape2)
source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

out_path <- 'figs/ybr/'

expt <- 'YBR'
c.src <- c('CAS','GAL','RAF')
den <- c(96,384,1536,6144)

data <- NULL
for (s in c.src) {
  for (d in den) {
    temp <- dbGetQuery(conn, sprintf('select * from %s_%s_%d_FITNESS',
                                     expt,s,d))
    temp$carbon <- s
    temp$density <- d
    data <- rbind(data, temp)
  }
}

data$colony[data$orf_name == 'BF_control'] <- 'BF_control'
data$colony[data$orf_name != 'BF_control'] <- 'YBR196C-A'

ggplot(data) +
  geom_boxplot(aes(x = colony, y = fitness)) +
  geom_violin(aes(x = colony, y = fitness), fill = 'transparent') +
  stat_compare_means(aes(x = colony, y = fitness)) +
  labs(y = 'Fitness') +
  theme_linedraw() +
  theme(axis.title.x = element_blank()) +
  facet_wrap(.~density * carbon,
             nrow = 4, ncol = 3,
             scales = 'free')
ggsave(sprintf("%s%s_FITNESS.jpg",out_path,expt),
       height = 10, width = 10,
       dpi = 1000)


ggplot(data) +
  geom_boxplot(aes(x = colony, y = average)) +
  geom_violin(aes(x = colony, y = average), fill = 'transparent') +
  stat_compare_means(aes(x = colony, y = fitness)) +
  labs(y = 'Colony Size (pix. count)') +
  theme_linedraw() +
  theme(axis.title.x = element_blank()) +
  facet_wrap(.~density * carbon,
             nrow = 4, ncol = 3,
             scales = 'free')
ggsave(sprintf("%s%s_CS.jpg",out_path,expt),
       height = 10, width = 10,
       dpi = 1000)
