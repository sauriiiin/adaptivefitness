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
library(egg)
library(zoo)
library(ggrepel)
library(plotly)
library(scales)
library(reshape2)
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
  

##### FIGURE S3: VIRTUAL PLATE 1 EXAMPLE
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

#####  FIGURE S4: VIRTUAL PLATES CS DISTRIBUTION
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


##### PLATE SURFACE AND SPATIAL BIAS
load(sprintf('%s4C4SURFACE.RData',out_path))
dat.raw$average[dat.raw$average == 0] <- NA
load(sprintf('%sSOURCENORMALIZATIONDATA.RData',out_path))
dat.sn$bg[is.na(dat.sn$average)] <- NA

dat.sn$source[dat.sn$source == '1TL'] <- 'Top Left'
dat.sn$source[dat.sn$source == '2TR'] <- 'Top Right'
dat.sn$source[dat.sn$source == '3BL'] <- 'Bottom Left'
dat.sn$source[dat.sn$source == '4BR'] <- 'Bottom Right'

dat.sn$source <- factor(dat.sn$source, levels = c('Top Left','Top Right','Bottom Left','Bottom Right'))
dat.raw$source <- factor(dat.raw$source, levels = c('Top Left','Top Right','Bottom Left','Bottom Right'))

ul <- quantile(dat.sn$average, 0.9999, na.rm = T)
ll <- quantile(dat.sn$average, 0.0001, na.rm = T)
m <- quantile(dat.sn$average, 0.5, na.rm = T)

plt.avg <- ggplot(dat.raw[dat.raw$plate == 3,],
                  aes(x = col, y = row)) +
  geom_tile(aes(fill = average)) +
  scale_size_continuous(range = c(0,2.5), guide = F) +
  scale_x_continuous(breaks = seq(0,96,4)) +
  scale_y_continuous(breaks = seq(0,64,4),trans = 'reverse') +
  # scale_fill_gradient2(low = "#F57C00", high = "black",
  #                      mid = "grey80", midpoint = m,
  #                      guide = F) +
  scale_fill_gradient(low = "#FF9800", high = "black", na.value = "white",
                      limits = c(ll,ul), oob = squish,
                      guide = F) +
  # scale_fill_distiller(name = 'Colony Size (pix.)',
  #                      palette = "RdYlGn",
  #                      guide = F) +
  labs(x = 'Column', y = 'Row') +
  coord_cartesian(xlim = c(1,96), ylim = c(64,1)) +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")),
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"))

plt.avg.cont <- ggplot(dat.sn[dat.sn$method == 3 & dat.sn$plate == 3,],
                  aes(x = col, y = row)) +
  geom_tile(aes(fill = average)) +
  geom_point(data = dat.sn[dat.sn$method == 3 & dat.sn$plate == 3 & dat.sn$orf_name == 'BF_control',],
             aes(x = col, y = row, col = 'black'), col = 'grey60', shape = 4, size = 0.4) +
  scale_size_continuous(range = c(0,2.5), guide = F) +
  scale_x_continuous(breaks = seq(0,96,4)) +
  scale_y_continuous(breaks = seq(0,64,4),trans = 'reverse') +
  scale_fill_gradient(low = "#FFA000", high = "black",
                      limits = c(ll,ul),
                      guide = F) +
  # scale_fill_distiller(name = 'Colony Size (pix.)',
  #                      palette = "Set1",
  #                      guide = F) +
  scale_color_manual(name = '',
                      label = '= Ref.',
                     values = 'black',
                     guide = F) +
  labs(x = 'Column', y = 'Row') +
  coord_cartesian(xlim = c(1,96), ylim = c(64,1)) +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")),
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"))

plt.avg.sr <- ggplot(dat.raw[dat.raw$plate == 3,],
                     aes(x = col, y = row)) +
  geom_point(aes(col = average), shape = 15, size = 0.81) +
  # scale_size_continuous(range = c(10,20), guide = F) +
  scale_x_continuous(breaks = seq(0,96,4)) +
  scale_y_continuous(breaks = seq(0,64,4),trans = 'reverse') +
  scale_color_gradient(name = 'Colony\nSize (pix.)',
                       low = "#FFA000", high = "black",na.value = "white",
                       limits = c(ll,ul), oob = squish) +
  # scale_color_distiller(name = 'Colony\nSize (pix.)',
  #                      palette = "Set1") +
  labs(x = 'Column', y = 'Row') +
  coord_cartesian(xlim = c(1,96), ylim = c(64,1)) +
  facet_wrap(.~source, ncol = 2) +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'right',
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")),
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"))

plt.avg.sr.cont <- ggplot(dat.sn[dat.sn$method == 3 & dat.sn$plate == 3,],
                     aes(x = col, y = row)) +
  geom_point(aes(col = average), shape = 15, size = 0.81) +
  geom_point(data = dat.sn[dat.sn$method == 3 & dat.sn$plate == 3 & dat.sn$orf_name == 'BF_control',],
             aes(x = col, y = row, fill = 'black'), col = 'grey60', shape = 4, size = 0.4) +
  # scale_size_continuous(range = c(10,20), guide = F) +
  scale_x_continuous(breaks = seq(0,96,4)) +
  scale_y_continuous(breaks = seq(0,64,4),trans = 'reverse') +
  scale_color_gradient(name = 'Colony\nSize (pix.)',
                      low = "#FFA000", high = "black",
                      limits = c(ll,ul)) +
  # scale_color_distiller(name = 'Colony\nSize (pix.)',
  #                       palette = "Set1") +
  scale_fill_discrete(name = '',
                      label = '= Ref.') +
  labs(x = 'Column', y = 'Row') +
  coord_cartesian(xlim = c(1,96), ylim = c(64,1)) +
  facet_wrap(.~source, ncol = 2) +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'right',
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")),
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"))


plt.srf.nocont <- ggpubr::ggarrange(plt.avg, plt.avg.sr,
                     ncol = 2, nrow = 1,
                     common.legend = T,
                     legend = 'right')
# # ggsave(sprintf("%sFIGURES7_1.jpg",out_path), plt.srf.nocont,
# #        height = 60, width = two.c, units = 'mm',
# #        dpi = 300)

plt.srf.cont <- ggpubr::ggarrange(plt.avg.cont, plt.avg.sr.cont,
                             ncol = 2, nrow = 1,
                             common.legend = T,
                             legend = 'right')
# # ggsave(sprintf("%sFIGURES7_2.jpg",out_path), plt.srf.cont,
# #        height = 60, width = two.c, units = 'mm',
# #        dpi = 300)

# plt.srf <- ggarrange(plt.srf.nocont, plt.srf.cont,
#                      labels = c('a','b'),
#                      label.args = list(gp=gpar(font = 2, fontsize = lbls),
#                                        hjust=0),
#                      nrow = 2, ncol = 1)
# # ggsave(sprintf("%sFIGURES7.jpg",out_path), plt.srf,
# #        height = 120, width = two.c, units = 'mm',
# #        dpi = 300)


ul <- quantile(dat.sn$fitness, 0.9999, na.rm = T)
ll <- quantile(dat.sn$fitness, 0.0001, na.rm = T)
m <- quantile(dat.sn$fitness, 0.5, na.rm = T)

plt.fit <- ggplot(dat.sn[dat.sn$method == 3 & dat.sn$plate == 3,],
                  aes(x = col, y = row)) +
  geom_tile(aes(fill = fitness)) +
  scale_size_continuous(range = c(0,2.5), guide = F) +
  scale_x_continuous(breaks = seq(0,96,4)) +
  scale_y_continuous(breaks = seq(0,64,4),trans = 'reverse') +
  scale_fill_gradient2(name = 'Fitness    ',
                       low = "#D32F2F", high = "#FFA000",
                       mid = "#303F9F", midpoint = m,
                       limits = c(0.7,1.3),
                       breaks = seq(0.8,1.2,0.2)) +
  # scale_fill_distiller(name = 'Fitness',
  #                      palette = "Set1",
  #                      guide = F) +
  labs(x = 'Column', y = 'Row') +
  coord_cartesian(xlim = c(1,96), ylim = c(64,1)) +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")),
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"))


plt.fit.sr <- ggplot(dat.sn[dat.sn$method == 3 & dat.sn$plate == 3,],
                     aes(x = col, y = row)) +
  geom_point(aes(col = fitness), shape = 15, size = 0.81) +
  # scale_size_continuous(range = c(10,20), guide = F) +
  scale_x_continuous(breaks = seq(0,96,4)) +
  scale_y_continuous(breaks = seq(0,64,4),trans = 'reverse') +
  scale_color_gradient2(name = 'Fitness    ',
                       low = "#D32F2F", high = "#FFA000",
                       mid = "#303F9F", midpoint = m,
                       limits = c(0.7,1.3),
                       breaks = seq(0.8,1.2,0.2),
                       guide = F) +
  # scale_color_distiller(name = 'Fitness',
  #                       palette = "Set1") +
  labs(x = 'Column', y = 'Row') +
  coord_cartesian(xlim = c(1,96), ylim = c(64,1)) +
  facet_wrap(.~source, ncol = 2) +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'right',
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")),
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"))

plt.srf.fit <- ggpubr::ggarrange(plt.fit.sr, plt.fit,
                                  ncol = 2, nrow = 1,
                                  common.legend = T,
                                  legend = 'right')

plt.srf.all <- ggarrange(plt.srf.nocont, plt.srf.cont,plt.srf.fit,
                         labels = c('a','b','c'),
                         label.args = list(gp=gpar(font = 2, fontsize = lbls),
                                           hjust=0),
                         nrow = 3, ncol = 1)

ggsave(sprintf("%sFIGURES7.jpg",out_path), plt.srf.all,
       height = two.c, width = two.c, units = 'mm',
       dpi = 300)

##### REFERENCE FITNESS DISTRIBUTION
load(sprintf('%sCONTFIT.RData',out_path))

ggplot(contfit) +
  geom_line(aes(x = fitness), stat = 'density', trim = T) +
  geom_vline(xintercept = quantile(contfit$fitness, c(0.025,0.975)),
             col = 'red', linetype = 'dashed') +
  labs(x = 'Reference Fitness Distribution',
       y = 'Density') +
  coord_cartesian(xlim = c(0.90,1.10)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt))
ggsave(sprintf("%sFIGURES8.jpg",out_path),
       height = one.c, width = one.c, units = 'mm',
       dpi = 300)

##### BEAN SPECIFICITY AND SENSITIVITY
load(sprintf('%sBEANSPEDATA.RData',out_path))
load(sprintf('%sBEANSENDATA.RData',out_path))

t <- dat.cnt2$Beneficial[1]
sen.fdr <- ggplot(dat.cnt2) +
  geom_area(aes(x = (cen-1)*100, y = Neutral, fill = 'Neutral'), alpha = 1) +
  geom_area(aes(x = (cen-1)*100, y = Beneficial, fill = 'Beneficial'), alpha = 1) +
  geom_area(aes(x = (cen-1)*100, y = Deleterious, fill = 'Deleterious'), alpha = 1) +
  geom_vline(xintercept = seq(-2,2,0.025)*100, col = '#757575', lwd = 0.5, alpha =0.5) +
  geom_hline(yintercept = c(t * seq(0,1,0.05)), col = '#757575', lwd = 0.5, alpha =0.5) +
  labs(x = 'Fitness Effect',
       y = 'Sensitivity') +
  scale_y_continuous(breaks = c(t * seq(0,1,0.1)),
                     minor_breaks = c(t * seq(0,1,0.05)),
                     labels = c('0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')) +
  scale_x_continuous(breaks = seq(-2,2,0.05)*100,
                     minor_breaks = seq(-2,2,0.025)*100) +
  scale_color_manual(name = 'Effects',
                     breaks = c('Beneficial','Neutral','Deleterious'),
                     values = c('Deleterious'='#3F51B5',
                                'Neutral'='#212121',
                                'Beneficial'='#FFC107'),
                     guide = F) +
  scale_fill_manual(name = 'Effects',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107')) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(5, "mm")) +
  coord_cartesian(xlim = c(-10,10),
                  ylim = c(0,t))


sen.p <- ggplot(dat.cnt2) +
  geom_area(aes(x = (cen-1)*100, y = Neutral_p, fill = 'Neutral'), alpha = 1) +
  geom_area(aes(x = (cen-1)*100, y = Beneficial_p, fill = 'Beneficial'), alpha = 1) +
  geom_area(aes(x = (cen-1)*100, y = Deleterious_p, fill = 'Deleterious'), alpha = 1) +
  geom_vline(xintercept = seq(-2,2,0.025)*100, col = '#757575', lwd = 0.5, alpha =0.5) +
  geom_hline(yintercept = c(t * seq(0,1,0.05)), col = '#757575', lwd = 0.5, alpha =0.5) +
  labs(x = 'Fitness Effect',
       y = 'Sensitivity') +
  scale_y_continuous(breaks = c(t * seq(0,1,0.1)),
                     minor_breaks = c(t * seq(0,1,0.05)),
                     labels = paste(sprintf('%0.0f',seq(0,1,0.1)*100),'%', sep = '')) +
  scale_x_continuous(breaks = seq(-2,2,0.05)*100,
                     minor_breaks = seq(-2,2,0.025)*100,
                     labels = paste(sprintf('%0.0f',seq(-2,2,0.05)*100),'%', sep = '')) +
  scale_color_manual(name = 'Effects',
                     breaks = c('Beneficial','Neutral','Deleterious'),
                     values = c('Deleterious'='#3F51B5',
                                'Neutral'='#212121',
                                'Beneficial'='#FFC107')) +
  scale_fill_manual(name = 'Effects',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107')) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(5, "mm")) +
  coord_cartesian(xlim = c(-10,10),
                  ylim = c(0,t))


spe.1 <- ggplot(stats.tmp[stats.tmp$hours == stats.tmp$cont_hrs,],
                aes(x = p, y = (1-fpr)*100, col = as.factor(hours))) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), aes(group = 1),
               geom="ribbon", color="blue", alpha = 0.4) +
  stat_summary(fun=mean, geom="line", color="blue", lwd =0.5) +
  # stat_summary(fun=mean, geom="point", color="#FFC107", size = 2) +
  labs(x = "p",
       y = "Specificity") +
  scale_x_continuous(breaks = seq(-1,1,0.025),
                     minor_breaks = seq(-1,1,0.0125)) +
  scale_y_continuous(breaks = seq(0,200,1),
                     minor_breaks = seq(0,200,0.5),
                     labels = paste(sprintf('%d',seq(0,200,1)),'%', sep = '')) +
  coord_cartesian(xlim = c(0, 0.1), ylim = c(90, 100)) +
  scale_color_discrete(name = "Hours",
                       breaks=as.character(c(2.9,4.02,4.89,6.14,6.88,7.85,8.96,9.97,11.04))) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(5, "mm")) #+
# coord_cartesian(xlim = c(0, 0.1), ylim = c(95, 100))

bean.perf <- ggarrange(spe.1, sen.fdr,
                  nrow = 1, ncol = 2,
                  labels = c('a','b'),
                  label.args = list(gp=gpar(font = 2, fontsize = lbls),
                                    hjust=-1))
ggsave(sprintf("%sFIGURES9.jpg",out_path), bean.perf,
       height = 80, width = two.c, units = 'mm',
       dpi = 300)

##### FIGURE S10
## CHANGE IN SPECIFICITY WITH REF. PROPORTIONS
load(sprintf('%sSPECIFICITY.RData',out_path))

spe.data$abs_cen <- abs(1-spe.data$cen) * 100
spe.data$rep <- as.factor(spe.data$rep)
spe.data$ref <- as.factor(spe.data$ref)

spe.ref <- spe.data[spe.data$cont_hrs == spe.data$hours &
                    round(spe.data$abs_cen) <= 5 &
                    spe.data$p <= 0.05,]

spe.r <- ggplot(spe.ref,
       aes(x = rep, y = round((1-fpr),4)*100,
           fill = ref)) +
  geom_boxplot() +
  labs(x = 'No. of Replicates',
       y = 'Specificity') +
  scale_x_discrete(breaks = seq(0,16,2)) +
  scale_y_continuous(breaks = seq(0,100,2),
                     minor_breaks = seq(0,105,1),
                     labels = paste(seq(0,100,2),'%',sep='')) +
  scale_fill_manual(name = 'Reference\nProportion',
                    breaks = c(0.0625,0.125,0.1875,0.25),
                    values = c('#C5CAE9','#448AFF','#3F51B5','#212121'),
                    labels = c(sprintf('%0.2f%%',c(0.0625,0.125,0.1875,0.25)*100)),
                    guide = F) +
  coord_cartesian(ylim = c(90,100)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt))

# ggsave(sprintf("%sFIGURES10.jpg",out_path),
#        height = 70, width = one.c, units = 'mm',
#        dpi = 300)


## RMSE WITH CHANGE IN REF PROP
# No. of replicates do not affect rmse
load(sprintf('%sRMSE_REFREP.RData',out_path))
rmse.rr$rep <- as.factor(rmse.rr$rep)
rmse.rr$ref <- as.factor(rmse.rr$ref)

rmse.r <- ggplot(rmse.rr[rmse.rr$rep == 16 & rmse.rr$hours > 2,],
       aes(x = hours, y = per, col = ref)) + 
  geom_point(size = 2, alpha = 0.9) +
  labs(x = 'Time (hour)', y = 'RMSE %') +
  scale_x_continuous(breaks = seq(0,12,2)) +
  scale_y_continuous(breaks = seq(0,100,2),
                     minor_breaks = seq(0,105,1),
                     labels = paste(seq(0,100,2),'%',sep='')) +
  scale_color_manual(name = 'Reference\nProportion',
                     breaks = c(0.0625,0.125,0.1875,0.25),
                     values = c('#C5CAE9','#448AFF','#3F51B5','#212121'),
                     labels = c(sprintf('%0.2f%%',c(0.0625,0.125,0.1875,0.25)*100))) +
  # coord_cartesian(ylim = c(90,100)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt))

figs10 <- ggarrange(spe.r, rmse.r,
                       nrow = 1, ncol = 2,
                       labels = c('a','b'),
                       label.args = list(gp=gpar(font = 2, fontsize = lbls),
                                         hjust=-1))
ggsave(sprintf("%sFIGURES10.jpg",out_path), figs10,
       height = 80, width = two.c, units = 'mm',
       dpi = 300)


##### VP2 Hourly and PIX ES wise results
load(sprintf('%sVP2PIXES.RData',out_path))

rnd.es$lid_result <- factor(rnd.es$lid_result, levels = c('Ben/Ben',
                                                          'Del/Del',
                                                          'Neu/Ben',
                                                          'Neu/Del',
                                                          'Ben/Del',
                                                          'Del/Ben'))

mcatbar <- ggplot(rnd.es) +
  geom_bar(aes(x = as.factor(hours), fill = mcat_result)) +
  labs(title = 'MCAT',
       x = 'Reference Population Time Point (hours)',
       y = 'Count') +
  scale_fill_manual(name = 'Effects',
                    breaks = c('Neu/Ben',
                               'Ben/Ben',
                               'Del/Ben',
                               'Ben/Del',
                               'Del/Del',
                               'Neu/Del'),
                    values = c('Neu/Del'='#536DFE',
                               'Del/Del'='#303F9F',
                               'Ben/Del'='#C5CAE9',
                               'Del/Ben'='#FFECB3',
                               'Ben/Ben'='#FFA000',
                               'Neu/Ben'='#FFC107'),
                    labels = c('Neu/Del'='False-Neutral Deleterious',
                               'Del/Del'='True Deleterious',
                               'Ben/Del'='False-Beneficial Deleterious',
                               'Del/Ben'='False-Deleterious Beneficial',
                               'Ben/Ben'='True Beneficial',
                               'Neu/Ben'='False-Neutral Beneficial'),
                    guide = F) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt))

lidbar <- ggplot(rnd.es) +
  geom_bar(aes(x = as.factor(hours), fill = lid_result)) +
  labs(title = 'LID',
       x = 'Reference Population Time Point (hours)',
       y = 'Count') +
  scale_fill_manual(name = 'Effects',
                    breaks = c('Neu/Ben',
                               'Ben/Ben',
                               'Del/Ben',
                               'Ben/Del',
                               'Del/Del',
                               'Neu/Del'),
                    values = c('Neu/Del'='#536DFE',
                               'Del/Del'='#303F9F',
                               'Ben/Del'='#C5CAE9',
                               'Del/Ben'='#FFECB3',
                               'Ben/Ben'='#FFA000',
                               'Neu/Ben'='#FFC107'),
                    labels = c('Neu/Del'='False-Neutral Deleterious',
                               'Del/Del'='True Deleterious',
                               'Ben/Del'='False-Beneficial Deleterious',
                               'Del/Ben'='False-Deleterious Beneficial',
                               'Ben/Ben'='True Beneficial',
                               'Neu/Ben'='False-Neutral Beneficial'),
                    guide = F) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt))

truthbar <- ggplot(rnd.es,
       aes(x = as.factor(hours))) +
  geom_bar(aes(fill = truth)) +
  labs(title = 'TRUTH',
       x = 'Reference Population Time Point (hours)',
       y = 'Count') +
  scale_fill_manual(name = 'Effects',
                    breaks = c('Beneficial',
                               'Deleterious'),
                    values = c('Deleterious'='#303F9F',
                               'Beneficial'='#FFA000'),
                    labels = c('Deleterious'='True Deleterious',
                               'Beneficial'='True Beneficial'),
                    guide = F) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt))


lidpixes <- ggplot(rnd.es,
       aes(x = lid_pix_es * 100, fill = lid_result)) +
  geom_histogram(binwidth = 10) +
  labs(title = 'LID',
       x = 'Fitness Effect',
       y = 'Count') +
  scale_x_continuous(breaks = seq(-1000,1000,100),
                     minor_breaks = seq(-1000,1000,10),
                     labels = paste(seq(-1000,1000,100),'%', sep = '')) +
  scale_fill_manual(name = 'Effects',
                    breaks = c('Neu/Ben',
                               'Ben/Ben',
                               'Del/Ben',
                               'Ben/Del',
                               'Del/Del',
                               'Neu/Del'),
                    values = c('Neu/Del'='#536DFE',
                               'Del/Del'='#303F9F',
                               'Ben/Del'='#C5CAE9',
                               'Del/Ben'='#FFECB3',
                               'Ben/Ben'='#FFA000',
                               'Neu/Ben'='#FFC107'),
                    labels = c('Neu/Del'='False-Neutral Deleterious',
                               'Del/Del'='True Deleterious',
                               'Ben/Del'='False-Beneficial Deleterious',
                               'Del/Ben'='False-Deleterious Beneficial',
                               'Ben/Ben'='True Beneficial',
                               'Neu/Ben'='False-Neutral Beneficial'),
                    drop = F) +
  theme_linedraw() + 
  # facet_wrap(.~hours) +
  theme(plot.title = element_text(size = titles),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt))

mcatpixes <- ggplot(rnd.es,
       aes(x = lid_pix_es*100, fill = mcat_result)) +
  geom_histogram(binwidth = 10) +
  labs(title = 'MCAT',
       x = 'Fitness Effect',
       y = 'Count') +
  scale_x_continuous(breaks = seq(-1000,1000,100),
                     minor_breaks = seq(-1000,1000,10),
                     labels = paste(seq(-1000,1000,100),'%', sep = '')) +
  scale_fill_manual(name = 'Effects',
                    breaks = c('Neu/Ben',
                               'Ben/Ben',
                               'Del/Ben',
                               'Ben/Del',
                               'Del/Del',
                               'Neu/Del'),
                    values = c('Neu/Del'='#536DFE',
                               'Del/Del'='#303F9F',
                               'Ben/Del'='#C5CAE9',
                               'Del/Ben'='#FFECB3',
                               'Ben/Ben'='#FFA000',
                               'Neu/Ben'='#FFC107'),
                    labels = c('Neu/Del'='False-Neutral Deleterious',
                               'Del/Del'='True Deleterious',
                               'Ben/Del'='False-Beneficial Deleterious',
                               'Del/Ben'='False-Deleterious Beneficial',
                               'Ben/Ben'='True Beneficial',
                               'Neu/Ben'='False-Neutral Beneficial'),
                    drop = F) +
  theme_linedraw() +
  # facet_wrap(.~hours) +
  theme(plot.title = element_text(size = titles),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt))

hrlybars <- ggpubr::ggarrange(lidbar, truthbar, mcatbar,
                              nrow = 1)
pixes <- ggpubr::ggarrange(lidpixes, mcatpixes,
                           common.legend = T,
                           legend = 'bottom')
pixesbars <- ggarrange(hrlybars, pixes,
                       nrow = 2,
                       heights = c(1.4,1),
                       labels = c('a','b'),
                       label.args = list(gp=gpar(font = 2, fontsize = lbls),
                                         hjust=-1))
ggsave(sprintf("%sFIGURES11.jpg",out_path), pixesbars,
       height = 120, width = two.c, units = 'mm',
       dpi = 300)
