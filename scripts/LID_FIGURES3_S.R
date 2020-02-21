##### LID PAPER SUPPLEMENTARY FIGURES #3
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 02/16/2020

##### INITIALIZE
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyverse)
library(ggpubr)
library(stringr)
library(scales)
library(egg)
library(zoo)
library(ggrepel)
library(RMariaDB)

source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")
out_path = 'figs/paper/';

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 7
txt <- 5
lbls <- 9

##### FIGURE S1. Prescreens Plate Effect
load(sprintf('%s4C4PSDATA.RData',out_path))

ps.plateeffect <- ggplot(data.ps,
       aes(x = source, y = average-median_cs)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, col = 'red', linetype = 'dashed') +
  labs(x = 'Source Plate',
       y = 'Median Substracted Pixel Counts') +
  scale_x_discrete(breaks=c("1TL","2TR","3BL","4BR"),
                   labels=c("Top Left","Top Right","Bottom Left","Bottom Right")) +
  stat_compare_means(label.x = 2.5, label.y = 700, hjust = 0.5, size = 1.5) +
  facet_wrap(.~stage) +
  theme_linedraw()+
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = txt,
                                  margin = margin(2,0,2,0, "mm")))
ggsave(sprintf("%sFIGURES1.jpg",out_path), ps.plateeffect,
       height = one.c, width = one.c, units = 'mm',
       dpi = 300)

##### FIGURE S2 PLATEMAPS
p2c.dat <- dbGetQuery(conn, 'select * from 4C4_pos2coor a, 4C4_pos2orf_name b
                      where a.pos = b.pos')
p2c.dat$colony[is.na(p2c.dat$orf_name)] <- 'Gap'
p2c.dat$colony[p2c.dat$orf_name == 'BF_control'] <- 'Reference'
p2c.dat$colony[p2c.dat$orf_name != 'BF_control'] <- 'Query'

stocks <- ggplot(p2c.dat[p2c.dat$density == 384,],
                 aes(x = col, y = row, col = colony)) +
  geom_point(size = 1, shape = 15) +
  scale_color_manual(name = 'Colony Type',
                     breaks = c("Reference","Query","Gap"),
                     values = c("Reference" = "#3F51B5",
                                "Query" = "#FFC107",
                                "Gap" = "#D32F2F"),
                     guide = F) +
  labs(title = 'Glycerol Stocks',
       x = 'Column', y = "Row") +
  facet_wrap(.~density*plate,
             scales = 'free',
             nrow = 1) +
  theme_linedraw() +
  theme(title = element_text(size = titles),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.title.x = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        strip.text = element_text(size = txt,
                                  margin = margin(0.5,0,0.5,0, "mm")))

wctp1 <- ggplot(p2c.dat[p2c.dat$density == 384,],
                 aes(x = col, y = row, col = colony)) +
  geom_point(size = 1, shape = 16) +
  scale_color_manual(name = 'Colony Type',
                     breaks = c("Reference","Query","Gap"),
                     values = c("Reference" = "#3F51B5",
                                "Query" = "#FFC107",
                                "Gap" = "#D32F2F"),
                     guide = F) +
  labs(title = 'Working Copies & Transition Plates (#1)',
       x = 'Column', y = "Row") +
  facet_wrap(.~density*plate,
             scales = 'free',
             nrow = 1) +
  theme_linedraw() +
  theme(title = element_text(size = titles),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.title.x = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        strip.text = element_text(size = txt,
                                  margin = margin(0.5,0,0.5,0, "mm")))

up1tp2 <- ggplot(p2c.dat[p2c.dat$density == 1536,],
                aes(x = col, y = row, col = colony)) +
  geom_point(size = 0.5, shape = 16) +
  scale_color_manual(name = 'Colony Type',
                     breaks = c("Reference","Query","Gap"),
                     values = c("Reference" = "#3F51B5",
                                "Query" = "#FFC107",
                                "Gap" = "#D32F2F"),
                     guide = F) +
  labs(title = 'Upscale Plates (#1) & Transition Plates (#2)',
       x = 'Column', y = "Row") +
  facet_wrap(.~density*plate,
             scales = 'free',
             nrow = 1) +
  theme_linedraw() +
  theme(title = element_text(size = titles),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.title.x = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        strip.text = element_text(size = txt,
                                  margin = margin(0.5,0,0.5,0, "mm")))

up2 <- ggplot(p2c.dat[p2c.dat$density == 6144,],
                 aes(x = col, y = row, col = colony)) +
  geom_point(size = 0.2, shape = 16) +
  scale_color_manual(name = 'Colony Type',
                     breaks = c("Reference","Query","Gap"),
                     values = c("Reference" = "#3F51B5",
                                "Query" = "#FFC107",
                                "Gap" = "#D32F2F")) +
  labs(title = 'Upscale Plates (#2)',
       x = 'Column', y = "Row") +
  facet_wrap(.~density*plate,
             scales = 'free',
             nrow = 1) +
  theme_linedraw() +
  theme(title = element_text(size = titles),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        strip.text = element_text(size = txt,
                                  margin = margin(0.5,0,0.5,0, "mm"))) +
  guides(color = guide_legend(override.aes = list(size = 3)))

plt.maps <- ggarrange(stocks, wctp1, up1tp2, up2,
                      ncol = 1, nrow = 4)

ggsave(sprintf("%sFIGURES2.jpg",out_path), plt.maps,
       height = two.c, width = two.c, units = 'mm',
       dpi = 300)
  

##### FIGURE S2: VIRTUAL PLATE 1 EXAMPLE
load(sprintf('%sVP1EGDATA.RData',out_path))

vir1plt <- ggplot(vir1.eg) +
  geom_point(aes(x = col, y = row, size = average, col = colony)) +
  scale_color_manual(name = 'Colony Type',
                     breaks = c("Reference","Query"),
                     values = c("Reference" = "#3F51B5",
                                "Query" = "#FFC107")) +
  scale_size_continuous(range = c(0,2.5), guide = F) +
  scale_x_continuous(breaks = seq(0,96,4)) +
  scale_y_continuous(breaks = seq(0,64,4),trans = 'reverse') +
  labs(x = 'Column', y = 'Row') +
  coord_cartesian(xlim = c(10,33), ylim = c(10,25)) +
  facet_wrap(.~kind, nrow = 1) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        strip.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"))

ggsave(sprintf("%sFIGURES3.jpg",out_path), vir1plt,
       height = 70, width = two.c, units = 'mm',
       dpi = 300)

#####  FIGURE S3: VIRTUAL PLATES CS DISTRIBUTION
load(sprintf("/home/sbp29/R/Projects/adaptivefitness/output/4C4_FS_CC_VIR_PLATE1.Rdata"))

vir.cs <- ggplot(fit.all, aes(x = average, fill = colony)) +
  # geom_line(stat = 'density', trim = T, lwd = 1) +
  geom_density(stat = 'density', alpha = 0.8) +
  labs(x = 'Colony Size (pixel)',
       y = 'Density') +
  scale_fill_manual(name = 'Colony Type',
                    breaks = c("Reference","Query"),
                    values = c("Reference" = "#3F51B5",
                               "Query" = "#FFC107")) +
  facet_wrap(.~vir_plate) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")),
        legend.key.size = unit(3, "mm"))

ggsave(sprintf("%sFIGURES4.jpg",out_path), vir.cs,
       height = two.c, width = two.c, units = 'mm',
       dpi = 300)

##### FIGURE S5: sensitivity when p <= 0.05
## LID_FIGURES3.R
  
##### FIGURE S6: REF REP WITH P CUT OFF
load(sprintf('%sREFREPDATA.RData',out_path))
refrep$power[refrep$hours < refrep$cont_hrs] <- refrep$Deleterious_p[refrep$hours < refrep$cont_hrs]/910 * 100
refrep$power[refrep$hours > refrep$cont_hrs] <- refrep$Beneficial_p[refrep$hours > refrep$cont_hrs]/910 * 100
refrep$abs_cen <- abs(1-refrep$cen) * 100

refrep$rep <- as.factor(refrep$rep)
refrep$ref_prop <- as.factor(refrep$ref_prop)

sen.rep <- ggplot(refrep[round(refrep$abs_cen) <= 5,],
                  aes(x = rep, y = power, fill = ref_prop)) +
  geom_boxplot() +
  # geom_smooth(aes(x = rep, y = power), method = 'loess', se = F, lwd = 1.2) +
  labs(x = 'No. of Replicates',
       y = 'Sensitivity') +
  scale_x_discrete(breaks = seq(0,16,2)) +
  scale_y_continuous(breaks = seq(0,100,10),
                     minor_breaks = seq(0,105,5),
                     labels = paste(seq(0,100,10),'%',sep='')) +
  scale_fill_manual(name = 'Reference\nProportion',
                    breaks = c(0.0625,0.125,0.1875,0.25),
                    values = c('#C5CAE9','#448AFF','#3F51B5','#212121'),
                    labels = c(sprintf('%0.2f%%',c(0.0625,0.125,0.1875,0.25)*100))) +
  coord_cartesian(ylim = c(0,100)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt))
ggsave(sprintf("%sFIGURES6.jpg",out_path), sen.rep,
       height = 70, width = one.c, units = 'mm',
       dpi = 300)

