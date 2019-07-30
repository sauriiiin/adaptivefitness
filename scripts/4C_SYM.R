##### SYMMETRY ISSUE
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 07/30/2019

##### INITIALIZE
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyverse)
# library(egg)
library(ggpubr)
library(stringr)

##### LOAD DATA
cfit.dat <- read.csv("/home/sbp29/R/Projects/proto_plots/rawdata/4C3_GA3_CONTFITALL_8.csv",
                     na.strings = "NaN")
colnames(cfit.dat) <- c("cont_hrs","hours","fitness")

q = 8
r = 18

# for (q in unique(stats.all$hours)) {
#   for (r in unique(stats.all$cont_hrs)) {
ggplot(cfit.dat[cfit.dat$hours == q & cfit.dat$cont_hrs == r,]) +
  geom_line(aes(x = fitness, col = "Reference"), stat = "density", lwd = 2) +
  geom_line(data = stats.all[stats.all$hours == q & stats.all$cont_hrs == r,],
            aes(x = cs_mean, col = "Query"), stat = "density", lwd = 2) +
  scale_color_discrete(name = "Strain") +
  labs(title = sprintf('BEAN | t(R) = %d | t(Q) = %d', r, q),
       subtitle = sprintf('D = %s | N = %s | B = %s',
                          sprintf("%s",dat.cnt[dat.cnt$hours == q & dat.cnt$cont_hrs == r,][4:6])[1],
                          sprintf("%s",dat.cnt[dat.cnt$hours == q & dat.cnt$cont_hrs == r,][4:6])[2],
                          sprintf("%s",dat.cnt[dat.cnt$hours == q & dat.cnt$cont_hrs == r,][4:6])[3]),
       x = "Fitness", y = "Density") +
  theme_linedraw() +
  coord_cartesian(xlim = c(0.8,1.2))
ggsave(sprintf("figs/%s_FDIS_%d_%d.jpg",expt_name,r,q),
       height = 20, width = 20, units = "cm",
       dpi = 300)
#   }
# }
