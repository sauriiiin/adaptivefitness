##### LID MAIN MANUSCRIPT FIGURES
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 07/06/2020

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

##### FUNCTIONS
lm_eqn <- function(df){
  m <- lm(average ~ bg, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}


##### FIGURE 2
load("/home/sbp29/R/Projects/adaptivefitness/output/4C4_FS_CC_VIR_PLATE1_ES.Rdata")
hhours <- c(0,1.02,1.38,2.9,4.02,4.89,6.14,6.88,7.85,8.96,9.97,11.04)

load(sprintf('%sVP2RNDDATA.RData',out_path))
# rnd_data$hours <- factor(rnd_data$hours, levels = c(sort(unique(rnd_data$hours))))
rnd_data$virtual <- sprintf('Ref. CS Time = %0.2f hr', rnd_data$hours)
rnd_data$virtual <- factor(rnd_data$virtual, levels = sprintf('Ref. CS Time = %0.2f hr', c(unique(sort(rnd_data$hours)))))


virP1$text <- '12x12 Fitness Effect Matrix'
es.heatmap <- ggplot(virP1) +
  geom_tile(aes(x = sprintf('%0.1f',que_hrs), y = sprintf('%0.1f',ref_hrs), fill = round((es-1)*100,0)),
            col = 'black') +
  geom_text(aes(x = sprintf('%0.1f',que_hrs), y = sprintf('%0.1f',ref_hrs), label = round((es-1)*100,0)),
            col = 'black', size = 2) +
  scale_x_discrete(limits = sprintf('%0.1f',hhours)) +
  scale_y_discrete(limits = sprintf('%0.1f',sort(hhours, decreasing = T))) +
  scale_fill_gradient2(name = 'Fitness Effect %',
                       low = "#3F51B5", high = "#FFC107", mid = "white",
                       trans = 'pseudo_log',
                       midpoint = 0,
                       breaks = c(-50,-5,0,5,50,500)) +
  labs(title = 'Bimodal Colony-Size Distribution',
       x = 'Mutant Colony-Size Time (hour)',
       y = 'Reference Colony-Size Time (hour)') +
  facet_wrap(.~text) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        strip.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"))

rnd_den <- ggplot(rnd_data[!is.na(rnd_data$average) & !is.na(rnd_data$colony)
                           & rnd_data$hours > 0,],
                  aes(x = average, fill = colony)) +
  # geom_line(stat = 'density', trim = T, lwd = 1.2) +
  geom_density(stat = 'density', alpha = 0.8) +
  labs(title = 'Random Colony-Size Distribution',
       x = 'Colony Size (pixel)',
       y = 'Density') +
  scale_fill_manual(name = 'Colony Type',
                    breaks = c("Reference","Query"),
                    labels = c("Reference","Mutant"),
                    values = c("Reference" = "#212121",
                               "Query" = "#9E9E9E")) +
  facet_wrap(.~virtual) +
  coord_cartesian(ylim = c(0,0.015)) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        strip.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"))

fig2 <- ggarrange(es.heatmap, rnd_den,
                  nrow = 1, widths = c(0.5,1))

ggsave(sprintf("%sFIGURE2_lower.jpg",out_path), fig2,
       height = 80, width = two.c, units = 'mm',
       dpi = 300)


##### FIGURE 3
#### a-e. BACKGROUND PREDICTION
load(sprintf('%sBACKGROUNDDATA.RData',out_path))

## LID BACKGROUND
pred.lid <- ggplot(dat.bg[!is.na(dat.bg$fitness) &
                            dat.bg$method == 3 &
                            dat.bg$hours > 0,],
                   aes(x = average, y = bg)) +
  geom_point(col = '#1976D2', alpha = 0.5) +
  geom_abline(col = 'red', linetype = 'dashed') +
  # geom_density_2d(col = 'red', lwd = 0.1) +
  annotate("text", y = 0, x = 580,
           label = lm_eqn(dat.bg[!is.na(dat.bg$fitness) & dat.bg$method == 3,]),
           size = 2, parse = T,
           hjust = 0) +
  labs(title = 'LID',
       x = 'Observed Colony Size (pixel)',
       y = 'Background Colony Size (pixel)') +
  coord_flip(xlim = c(0,650),
             ylim = c(0,650)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        plot.title = element_text(size = titles,
                                  hjust = 0.5))

## NOSN BACKGROUND
pred.nosn <- ggplot(dat.bg[!is.na(dat.bg$fitness) &
                             dat.bg$method == 2 &
                             dat.bg$hours > 0,],
                    aes(x = average, y = bg)) +
  geom_point(col = '#03A9F4', alpha = 0.5) +
  geom_abline(col = 'red', linetype = 'dashed') +
  # geom_density_2d(col = 'red', lwd = 0.1) +
  annotate("text", y = 0, x = 580,
           label = lm_eqn(dat.bg[!is.na(dat.bg$fitness) & dat.bg$method == 2,]),
           size = 2, parse = T,
           hjust = 0) +
  labs(title = 'LID - SN',
       x = 'Observed Colony Size (pixel)',
       y = 'Background Colony Size (pixel)') +
  coord_flip(xlim = c(0,650),
             ylim = c(0,650)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        plot.title = element_text(size = titles,
                                  hjust = 0.5))

## NOCC BACKGROUND
pred.nocc <- ggplot(dat.bg[!is.na(dat.bg$fitness) &
                             dat.bg$method == 5 &
                             dat.bg$hours > 0,],
                    aes(x = average, y = bg)) +
  geom_point(col = '#BBDEFB', alpha = 0.5) +
  geom_abline(col = 'red', linetype = 'dashed') +
  # geom_density_2d(col = 'red', lwd = 0.1) +
  annotate("text", y = 0, x = 580,
           label = lm_eqn(dat.bg[!is.na(dat.bg$fitness) & dat.bg$method == 5,]),
           size = 2, parse = T,
           hjust = 0) +
  labs(title = 'LID - AC',
       x = 'Observed Colony Size (pixel)',
       y = 'Background Colony Size (pixel)') +
  coord_flip(xlim = c(0,650),
             ylim = c(0,650)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        plot.title = element_text(size = titles,
                                  hjust = 0.5))

## BEAN BACKGROUND
pred.bean <- ggplot(dat.bg[!is.na(dat.bg$fitness) &
                             dat.bg$method == 4 &
                             dat.bg$hours > 0,],
                    aes(x = average, y = bg)) +
  geom_point(col = '#009688', alpha = 0.5) +
  geom_abline(col = 'red', linetype = 'dashed') +
  # geom_density_2d(col = 'red', lwd = 0.1) +
  annotate("text", y = 0, x = 580,
           label = lm_eqn(dat.bg[!is.na(dat.bg$fitness) & dat.bg$method == 4,]),
           size = 2, parse = T,
           hjust = 0) +
  labs(title = 'MCAT',
       x = 'Observed Colony Size (pixel)',
       y = 'Background Colony Size (pixel)') +
  coord_flip(xlim = c(0,650),
             ylim = c(0,650)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        plot.title = element_text(size = titles,
                                  hjust = 0.5))


## RND BACKGROUND
pred.rnd <- ggplot(dat.bg[!is.na(dat.bg$fitness) &
                            dat.bg$method == 6 &
                            dat.bg$hours > 0,],
                   aes(x = average, y = bg)) +
  geom_point(col = '#795548', alpha = 0.5) +
  geom_abline(col = 'red', linetype = 'dashed') +
  # geom_density_2d(col = 'red', lwd = 0.1) +
  annotate("text", y = 0, x = 580,
           label = lm_eqn(dat.bg[!is.na(dat.bg$fitness) & dat.bg$method == 6,]),
           size = 2, parse = T,
           hjust = 0) +
  labs(title = 'RND',
       x = 'Observed Colony Size (pixel)',
       y = 'Randomly Picked Colony Size (pixel)') +
  coord_flip(xlim = c(0,650),
             ylim = c(0,650)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        plot.title = element_text(size = titles,
                                  hjust = 0.5))

#### f. RMSE
load(sprintf('%sRMSEDATA.RData',out_path))
rmse <- rmse[rmse$hour > 0,]

pred.rmse <- ggplot(rmse[order(rmse$method, decreasing = T),]) +
  geom_point(aes(x = hour, y = per, fill = as.factor(method)),
             size = 1.2, alpha = 0.9, shape = 21, col = 'black') +
  labs(title = 'RMSE',
       x = 'Time (hour)', y = 'RMSE %') +
  scale_x_continuous(breaks = seq(0,12,2)) +
  scale_fill_manual(name = 'Method',
                     breaks = c('3','2','5','4','6'),
                     limits = c('3','2','5','4','6'),
                     labels = c('LID','LID-SN','LID-AC','MCAT','RND'),
                     values = c('3'='#1976D2','2'='#03A9F4','5'='#BBDEFB',
                                '4'='#009688','6'='#795548')) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles,
                                  hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = c(0.6, 0.75),
        legend.background = element_rect(fill = 'transparent'),
        legend.key = element_rect(size = 0.5),
        legend.key.size = unit(.8, 'lines')) +
  guides(fill = guide_legend(ncol=2))


#### g-k. SOURCE NORMALIZATION
load(sprintf('%sSOURCENORMALIZATIONDATA.RData',out_path))
dat.sn$source <- factor(dat.sn$source, levels = c("4BR","3BL","2TR","1TL"))

dat.sn$name <- substr(dat.sn$name, 4, 100)

sn.raw <- ggplot(dat.sn[!is.na(dat.sn$fitness) & dat.sn$method == 1,],
                 aes(x = source, y = average)) +
  # geom_boxplot(width = 0.5, outlier.size = 0.5, fill = '#C5CAE9', lwd = 0.15) +
  geom_violin(fill = 'grey90', draw_quantiles = c(0.25, 0.5, 0.75), lwd = 0.4) +
  stat_compare_means(label.x = 4.5, label.y = mean(dat.sn$average[!is.na(dat.sn$fitness) & dat.sn$method == 1], na.rm = T),
                     hjust = 0.5, size = 1.5) +
  scale_x_discrete(name="Source",
                   limits=c("4BR","3BL","2TR","1TL"),
                   labels=c("Bottom\nRight","Bottom\nLeft","Top\nRight","Top\nLeft")) +
  # scale_y_continuous(limits = c(0.7,1.3)) +
  facet_wrap(.~name, nrow = 1,
             scales = 'free') +
  labs(x = 'Source',
       y = 'Colony Size (pixel)') +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.y = element_text(angle = 90, hjust = .5),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        strip.text = element_text(size = txt, margin = margin(0.5,0,0.5,0, "mm")),
        legend.key.size = unit(5, "mm")) +
  coord_flip(xlim = c(1,4.1))

sn.sn <- ggplot(dat.sn[!is.na(dat.sn$fitness) & dat.sn$method == 4,],
                aes(x = source, y = fitness)) +
  # geom_boxplot(width = 0.5, outlier.size = 0.5, fill = '#303F9F', lwd = 0.15) +
  geom_violin(fill = '#1976D2', draw_quantiles = c(0.25, 0.5, 0.75), lwd = 0.4) +
  stat_compare_means(label.x = 4.5, label.y = 1, hjust = 0.5, size = 1.5) +
  scale_x_discrete(name="Source",
                   limits=c("4BR","3BL","2TR","1TL"),
                   labels=c("Bottom Right","Bottom Left","Top Right","Top Left")) +
  scale_y_continuous(limits = c(0.4,1.6)) +
  facet_wrap(.~name, nrow = 1,
             scales = 'free') +
  labs(x = 'Source',
       y = 'Fitness') +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        strip.text = element_text(size = txt, margin = margin(0.5,0,0.5,0, "mm")),
        legend.key.size = unit(5, "mm")) +
  coord_flip(xlim = c(1,4.1))

sn.nosn <- ggplot(dat.sn[!is.na(dat.sn$fitness) & dat.sn$method == 2,],
                  aes(x = source, y = fitness)) +
  # geom_boxplot(width = 0.5, outlier.size = 0.5, fill = '#536DFE', lwd = 0.15) +
  geom_violin(fill = '#03A9F4', draw_quantiles = c(0.25, 0.5, 0.75), lwd = 0.4) +
  stat_compare_means(label.x = 4.5, label.y = 1, hjust = 0.5, size = 1.5) +
  scale_x_discrete(name="Source",
                   limits=c("4BR","3BL","2TR","1TL"),
                   labels=c("Bottom Right","Bottom Left","Top Right","Top Left")) +
  scale_y_continuous(limits = c(0.4,1.6)) +
  facet_wrap(.~name, nrow = 1,
             scales = 'free') +
  labs(x = 'Source',
       y = 'Fitness') +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        strip.text = element_text(size = txt, margin = margin(0.5,0,0.5,0, "mm")),
        legend.key.size = unit(5, "mm")) +
  coord_flip(xlim = c(1,4.1))

sn.nocc <- ggplot(dat.sn[!is.na(dat.sn$fitness) & dat.sn$method == 3,],
                  aes(x = source, y = fitness)) +
  # geom_boxplot(width = 0.5, outlier.size = 0.5, fill = '#FFA000', lwd = 0.15) +
  geom_violin(fill = '#BBDEFB', draw_quantiles = c(0.25, 0.5, 0.75), lwd = 0.4) +
  stat_compare_means(label.x = 4.5, label.y = 1, hjust = 0.5, size = 1.5) +
  scale_x_discrete(name="Source",
                   limits=c("4BR","3BL","2TR","1TL"),
                   labels=c("Bottom Right","Bottom Left","Top Right","Top Left")) +
  scale_y_continuous(limits = c(0.4,1.6)) +
  facet_wrap(.~name, nrow = 1,
             scales = 'free') +
  labs(x = 'Source',
       y = 'Fitness') +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        strip.text = element_text(size = txt, margin = margin(0.5,0,0.5,0, "mm")),
        legend.key.size = unit(5, "mm")) +
  coord_flip(xlim = c(1,4.1))

sn.bean <- ggplot(dat.sn[!is.na(dat.sn$fitness) & dat.sn$method == 5,],
                  aes(x = source, y = fitness)) +
  # geom_boxplot(width = 0.5, outlier.size = 0.5, fill = '#FFEB3B', lwd = 0.15) +
  geom_violin(fill = '#009688', draw_quantiles = c(0.25, 0.5, 0.75), lwd = 0.4) +
  stat_compare_means(label.x = 4.5, label.y = 1, hjust = 0.5, size = 1.5) +
  scale_x_discrete(name="Source",
                   limits=c("4BR","3BL","2TR","1TL"),
                   labels=c("Bottom Right","Bottom Left","Top Right","Top Left")) +
  scale_y_continuous(limits = c(0.4,1.6)) +
  facet_wrap(.~name, nrow = 1,
             scales = 'free') +
  labs(x = 'Source',
       y = 'Fitness') +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        strip.text = element_text(size = txt, margin = margin(0.5,0,0.5,0, "mm")),
        legend.key.size = unit(5, "mm")) +
  coord_flip(xlim = c(1,4.1))


#### COMPILING THE FIGURE
fig3a_f <- ggpubr::ggarrange(pred.lid, pred.nosn, pred.nocc,
                             pred.bean, pred.rnd, pred.rmse,
                             ncol = 3, nrow = 2,
                             labels = c('a','b','c',
                                        'd','e','f'),
                             font.label = list(face = 'bold', size = lbls),
                             hjust=-1)

fig3g_k <- ggarrange(sn.raw,sn.sn,sn.nosn,sn.nocc,sn.bean,
                     nrow = 1, ncol = 5,
                     labels = c('g','h','i','j','k'),
                     label.args = list(gp=gpar(font = 2, fontsize = lbls),
                                       hjust = -1))

fig3 <- ggpubr::ggarrange(fig3a_f, fig3g_k, ncol = 1, nrow = 2,
                          heights = c(2,1))
ggsave(sprintf("%sFIGURE3.jpg",out_path), fig3,
       height = 180, width = two.c, units = 'mm',
       dpi = 300)


###### FIGURE 4: OVERALL PERFORMACE
load(sprintf('%sLID_SENSITIVITY_DATA.RData',out_path))
lid.sen.data <- dat.cnt2
load(sprintf('%sLID_SPECIFICITY_DATA.RData',out_path))
lid.spe.data <- stats.tmp

load(sprintf('%sREFREPDATA.RData',out_path))
refrep$power[refrep$hours < refrep$cont_hrs] <- refrep$Deleterious[refrep$hours < refrep$cont_hrs]/910 * 100
refrep$power[refrep$hours > refrep$cont_hrs] <- refrep$Beneficial[refrep$hours > refrep$cont_hrs]/910 * 100
refrep$abs_cen <- abs(1-refrep$cen) * 100
refrep$rep <- as.factor(refrep$rep)
refrep$ref_prop <- as.factor(refrep$ref_prop)

load(sprintf('%sMCAT_SENSITIVITY_DATA.RData',out_path))
mcat.sen.data <- dat.cnt2
load(sprintf('%sMCAT_SPECIFICITY_DATA.RData',out_path))
mcat.spe.data <- stats.tmp

t <- lid.sen.data$Beneficial[1]
lid.sen.fdr <- ggplot(lid.sen.data) +
  geom_area(aes(x = (cen-1)*100, y = Neutral, fill = 'Neutral'), alpha = 1) +
  geom_area(aes(x = (cen-1)*100, y = Beneficial, fill = 'Beneficial'), alpha = 1) +
  geom_area(aes(x = (cen-1)*100, y = Deleterious, fill = 'Deleterious'), alpha = 1) +
  geom_vline(xintercept = seq(-2,2,0.025)*100, col = '#757575', lwd = 0.5, alpha =0.5) +
  geom_hline(yintercept = c(t * seq(0,1,0.05)), col = '#757575', lwd = 0.5, alpha =0.5) +
  geom_vline(xintercept = c(-5,5), col = 'red', lwd = 1, linetype = 'dashed') +
  geom_text(x=0, y=t*0.9, label="LID", col = 'white', size = 3) +
  labs(x = 'Fitness Effect',
       y = 'Sensitivity') +
  scale_y_continuous(breaks = c(t * seq(0,1,0.1)),
                     minor_breaks = c(t * seq(0,1,0.05)),
                     labels = c('0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')) +
  scale_x_continuous(breaks = seq(-2,2,0.05)*100,
                     minor_breaks = seq(-2,2,0.025)*100,
                     labels = paste0(round(seq(-2,2,0.05)*100), '%', sep = '')) +
  scale_color_manual(name = 'Phenotype',
                     breaks = c('Beneficial','Neutral','Deleterious'),
                     values = c('Deleterious'='#3F51B5',
                                'Neutral'='#212121',
                                'Beneficial'='#FFC107'),
                     guide = F) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107')) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(4, "mm"),
        legend.position = c(0.87,0.2),
        legend.margin = margin(0.5,0.5,0.5,0.5, "mm")) +
  coord_cartesian(xlim = c(-10,10),
                  ylim = c(0,t))

sen.ref.rep <- ggplot(refrep[round(refrep$abs_cen) <= 5,],
                  aes(x = rep, y = power, fill = ref_prop)) +
  geom_boxplot(outlier.shape = NA) +
  # geom_smooth(aes(x = rep, y = power), method = 'loess', se = F, lwd = 1.2) +
  labs(x = 'No. of Replicates',
       y = 'Sensitivity') +
  scale_x_discrete(breaks = seq(0,16,2)) +
  scale_y_continuous(breaks = seq(0,100,10),
                     minor_breaks = seq(0,105,5),
                     labels = paste(seq(0,100,10),'%',sep='')) +
  scale_fill_manual(name = 'Reference Proportion',
                    breaks = c(0.0625,0.125,0.1875,0.25),
                    values = c('#607D8B','#757575','#BDBDBD','#FFFFFF'),
                    labels = c(sprintf('%0.2f%%',c(0.0625,0.125,0.1875,0.25)*100))) +
  coord_cartesian(ylim = c(0,100)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key = element_rect(size = 0.5),
        legend.key.size = unit(.8, 'lines'),
        legend.position = c(0.5,0.1),
        legend.margin = margin(0.5,0.5,0.5,0.5, "mm")) +
  guides(fill = guide_legend(nrow=1))


t <- mcat.sen.data$Beneficial[1]
mcat.sen.fdr <- ggplot(mcat.sen.data) +
  geom_area(aes(x = (cen-1)*100, y = Neutral, fill = 'Neutral'), alpha = 1) +
  geom_area(aes(x = (cen-1)*100, y = Beneficial, fill = 'Beneficial'), alpha = 1) +
  geom_area(aes(x = (cen-1)*100, y = Deleterious, fill = 'Deleterious'), alpha = 1) +
  geom_vline(xintercept = seq(-2,2,0.025)*100, col = '#757575', lwd = 0.5, alpha =0.5) +
  geom_hline(yintercept = c(t * seq(0,1,0.05)), col = '#757575', lwd = 0.5, alpha =0.5) +
  geom_vline(xintercept = c(-5,5), col = 'red', lwd = 1, linetype = 'dashed') +
  geom_text(x=0, y=t*0.9, label="MCAT", col = 'white', size = 3) +
  labs(x = 'Fitness Effect',
       y = 'Sensitivity') +
  scale_y_continuous(breaks = c(t * seq(0,1,0.1)),
                     minor_breaks = c(t * seq(0,1,0.05)),
                     labels = c('0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')) +
  scale_x_continuous(breaks = seq(-2,2,0.05)*100,
                     minor_breaks = seq(-2,2,0.025)*100,
                     labels = paste0(round(seq(-2,2,0.05)*100), '%', sep = '')) +
  scale_color_manual(name = 'Phenotype',
                     breaks = c('Beneficial','Neutral','Deleterious'),
                     values = c('Deleterious'='#3F51B5',
                                'Neutral'='#212121',
                                'Beneficial'='#FFC107'),
                     guide = F) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107')) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(4, "mm"),
        legend.position = c(0.87,0.2),
        legend.margin = margin(0.5,0.5,0.5,0.5, "mm")) +
  coord_cartesian(xlim = c(-10,10),
                  ylim = c(0,t))

spe.res <- ggplot(mcat.spe.data[mcat.spe.data$hours == mcat.spe.data$cont_hrs,],
                aes(x = p, y = (1-fpr)*100)) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), aes(group = 1),
               geom="ribbon", alpha = 0.4) +
  stat_summary(aes(color="MCAT"), fun=mean, geom="line", lwd =0.7) +
  stat_summary(data = lid.spe.data[lid.spe.data$hours == lid.spe.data$cont_hrs,],
               fun.data=mean_sdl, fun.args = list(mult=1),
               aes(x = p, y = (1-fpr)*100, group = 1),
               geom="ribbon", alpha = 0.4) +
  stat_summary(data = lid.spe.data[lid.spe.data$hours == lid.spe.data$cont_hrs,],
               aes(color="LID"), fun=mean, geom="line", lwd =0.7) +
  labs(x = "p-value",
       y = "Specificity") +
  scale_x_continuous(breaks = seq(-1,1,0.025),
                     minor_breaks = seq(-1,1,0.0125)) +
  scale_y_continuous(breaks = seq(0,200,1),
                     minor_breaks = seq(0,200,0.5),
                     labels = paste(sprintf('%d',seq(0,200,1)),'%', sep = '')) +
  coord_cartesian(xlim = c(0, 0.1), ylim = c(90, 100)) +
  scale_color_manual(name = "",
                     breaks = c('LID','MCAT'),
                     values = c('LID' = '#1976D2', 'MCAT' = '#009688')) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(4, "mm"),
        legend.position = c(0.2,0.2),
        legend.margin = margin(0.5,0.5,0.5,0.5, "mm"))

fig4 <- ggarrange(spe.res, lid.sen.fdr, mcat.sen.fdr, sen.ref.rep,
                  nrow = 2, ncol = 2,
                  labels = c('a','b','c','d'),
                  label.args = list(gp=gpar(font = 2, fontsize = lbls),
                                    hjust=-1))
ggsave(sprintf("%sFIGURE4.jpg",out_path), fig4,
       height = 180, width = two.c, units = 'mm',
       dpi = 300)


###### FIGURE 5: RANDOM CS DISTRIBUTION
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
       y = '') +
  scale_fill_manual(name = 'Phenotype',
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
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm"))

lidbar <- ggplot(rnd.es) +
  geom_bar(aes(x = as.factor(hours), fill = lid_result)) +
  labs(title = 'LID',
       x = 'Reference Population Time Point (hours)',
       y = 'Number of Mutant Strains') +
  scale_fill_manual(name = 'Phenotype',
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
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm"))

truthbar <- ggplot(rnd.es,
                   aes(x = as.factor(hours))) +
  geom_bar(aes(fill = truth)) +
  labs(title = 'TRUTH',
       x = 'Reference Population Time Point (hours)',
       y = '') +
  scale_fill_manual(name = 'Phenotype',
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
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm"))


lidpixes <- ggplot(rnd.es[!(rnd.es$lid_pix_es > 0 & rnd.es$truth == 'Deleterious') &
                            !(rnd.es$lid_pix_es < 0 & rnd.es$truth == 'Beneficial'),],
                   aes(x = lid_pix_es * 100, fill = lid_result)) +
  geom_histogram(binwidth = 10) +
  labs(title = 'LID',
       x = 'Fitness Effect',
       y = 'Number of Mutant Strains') +
  scale_x_continuous(breaks = seq(-1000,1000,100),
                     minor_breaks = seq(-1000,1000,10),
                     labels = paste(seq(-1000,1000,100),'%', sep = '')) +
  scale_fill_manual(name = 'Phenotype',
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
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm"))

mcatpixes <- ggplot(rnd.es[!(rnd.es$lid_pix_es > 0 & rnd.es$truth == 'Deleterious') &
                             !(rnd.es$lid_pix_es < 0 & rnd.es$truth == 'Beneficial'),],
                    aes(x = lid_pix_es*100, fill = mcat_result)) +
  geom_histogram(binwidth = 10) +
  labs(title = 'MCAT',
       x = 'Fitness Effect',
       y = '') +
  scale_x_continuous(breaks = seq(-1000,1000,100),
                     minor_breaks = seq(-1000,1000,10),
                     labels = paste(seq(-1000,1000,100),'%', sep = '')) +
  scale_fill_manual(name = 'Phenotype',
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
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm"))

hrlybars <- ggpubr::ggarrange(lidbar, truthbar, mcatbar,
                              nrow = 1)
pixes <- ggpubr::ggarrange(lidpixes, mcatpixes,
                           common.legend = T,
                           legend = 'bottom')
pixesbars <- ggarrange(hrlybars, pixes,
                       nrow = 2,
                       heights = c(1,1.4),
                       labels = c('a','b'),
                       label.args = list(gp=gpar(font = 2, fontsize = lbls),
                                         hjust=-1))

load(sprintf('%sVP2PIEDATA.RData',out_path))

rnd_pie.lid <- ggplot(data = pie.per[pie.per$method == 'LID',], aes(x = "", y = per, fill = value)) +
  geom_bar(stat = 'identity') +
  coord_polar("y", start = 0) +
  geom_label_repel(aes(label = sprintf("%0.2f%%",per)),
                   position = position_stack(vjust = 0.5),
                   label.size = 0.15,
                   show.legend = F) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('6BenNeutral',
                               '5TrueBeneficial',
                               '4FalseDeleterious',
                               '3FalseBeneficial',
                               '2TrueDeleterious',
                               '1DelNeutral'),
                    values = c('1DelNeutral'='#536DFE',
                               '2TrueDeleterious'='#303F9F',
                               '3FalseBeneficial'='#C5CAE9',
                               '4FalseDeleterious'='#FFECB3',
                               '5TrueBeneficial'='#FFA000',
                               '6BenNeutral'='#FFC107'),
                    labels = c('1DelNeutral'='False-Neutral Deleterious',
                               '2TrueDeleterious'='True Deleterious',
                               '3FalseBeneficial'='False-Beneficial Deleterious',
                               '4FalseDeleterious'='False-Deleterious Beneficial',
                               '5TrueBeneficial'='True Beneficial',
                               '6BenNeutral'='False-Neutral Beneficial'),
                    guide = F) +
  labs(title = 'LID') +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(size = titles),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = 'bottom',
        strip.text = element_text(size = titles),
        legend.box.spacing = unit(0.5,"mm"))

rnd_pie.truth <- ggplot(data = pie.per[pie.per$method == 'TRUTH',], aes(x = "", y = per, fill = value)) +
  geom_bar(stat = 'identity') +
  coord_polar("y", start = 0) +
  geom_label_repel(aes(label = sprintf("%0.2f%%",per)),
                   position = position_stack(vjust = 0.5),
                   label.size = 0.15,
                   show.legend = F) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('6BenNeutral',
                               '5TrueBeneficial',
                               '4FalseDeleterious',
                               '3FalseBeneficial',
                               '2TrueDeleterious',
                               '1DelNeutral'),
                    values = c('1DelNeutral'='#536DFE',
                               '2TrueDeleterious'='#303F9F',
                               '3FalseBeneficial'='#C5CAE9',
                               '4FalseDeleterious'='#FFECB3',
                               '5TrueBeneficial'='#FFA000',
                               '6BenNeutral'='#FFC107'),
                    labels = c('1DelNeutral'='False-Neutral Deleterious',
                               '2TrueDeleterious'='True Deleterious',
                               '3FalseBeneficial'='False-Beneficial Deleterious',
                               '4FalseDeleterious'='False-Deleterious Beneficial',
                               '5TrueBeneficial'='True Beneficial',
                               '6BenNeutral'='False-Neutral Beneficial'),
                    guide = F) +
  labs(title = 'TRUTH') +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(size = titles),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = 'bottom',
        strip.text = element_text(size = titles),
        legend.box.spacing = unit(0.5,"mm"))

rnd_pie.mcat <- ggplot(data = pie.per[pie.per$method == 'MCAT',], aes(x = "", y = per, fill = value)) +
  geom_bar(stat = 'identity') +
  coord_polar("y", start = 0) +
  geom_label_repel(aes(label = sprintf("%0.2f%%",per)),
                   position = position_stack(vjust = 0.5),
                   label.size = 0.15,
                   show.legend = F) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('6BenNeutral',
                               '5TrueBeneficial',
                               '4FalseDeleterious',
                               '3FalseBeneficial',
                               '2TrueDeleterious',
                               '1DelNeutral'),
                    values = c('1DelNeutral'='#536DFE',
                               '2TrueDeleterious'='#303F9F',
                               '3FalseBeneficial'='#C5CAE9',
                               '4FalseDeleterious'='#FFECB3',
                               '5TrueBeneficial'='#FFA000',
                               '6BenNeutral'='#FFC107'),
                    labels = c('1DelNeutral'='False-Neutral Deleterious',
                               '2TrueDeleterious'='True Deleterious',
                               '3FalseBeneficial'='False-Beneficial Deleterious',
                               '4FalseDeleterious'='False-Deleterious Beneficial',
                               '5TrueBeneficial'='True Beneficial',
                               '6BenNeutral'='False-Neutral Beneficial'),
                    guide = F) +
  labs(title = 'MCAT') +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(size = titles),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = 'bottom',
        strip.text = element_text(size = titles),
        legend.box.spacing = unit(0.5,"mm"))


fig5.a <- ggpubr::ggarrange(rnd_pie.lid, rnd_pie.truth, rnd_pie.mcat,
                            nrow = 1)

fig5 <- ggarrange(fig5.a, hrlybars, pixes,
                  nrow = 3,
                  ncol = 1,
                  labels = c('a','b','c'),
                  label.args = list(gp=gpar(font = 2, fontsize = lbls),
                                    hjust=-1))
ggsave(sprintf("%sFIGURE5.jpg",out_path), fig5,
       height = two.c, width = two.c, units = 'mm',
       dpi = 300)

