##### LID MAIN MANUSCRIPT FIGURES
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 09/25/2020

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
library(png)
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
txt <- 7
lbls <- 9

##### FIGURE 2. EMPIRICAL STRATEGY FOR LI DETECTOR EVALUATION
load(file = sprintf('%sCS_DENSITY_DATA.RData', out_path))
density.data$colony[density.data$orf_name == 'BF_control'] <- 'Reference'
density.data$colony[density.data$orf_name != 'BF_control' & !is.na(density.data$orf_name)] <- 'Mutant'
density.data$title <- sprintf('"t"["R"]=="%0.1fhr"', density.data$hours)
density.data$title <- factor(density.data$title, levels = sprintf('"t"["R"]=="%0.1fhr"',
                                                                  c(unique(sort(density.data$hours)))))

load("/home/sbp29/R/Projects/adaptivefitness/output/4C4_FS_CC_VIR_PLATE1_ES.Rdata")
hhours <- c(0,1.02,1.38,2.9,4.02,4.89,6.14,6.88,7.85,8.96,9.97,11.04)
virP1$text <- '11x11 Fitness Effect Matrix'

##### 2A
illustration <- readPNG('figs/paper/FIGURE2A.png')

fig.2a <- ggplot() + 
  background_image(illustration) +
  theme(plot.margin = margin(t=6, l=2, r=2, b=6, unit = "mm"),
        plot.background = element_blank())

##### 2B
fig.2b <- ggplot(density.data[!is.na(density.data$colony) &
                                density.data$hours > 0,],
                 aes(x = average, y = colony, group = colony, fill = colony)) +
  geom_density_ridges(quantile_lines = TRUE,
                      scale = 3, alpha = 0.8, size = 0.3,
                      vline_size = 0.2, vline_color = "black") +
  labs(title = '',
       x = 'Colony Size (pixel counts)',
       y = 'Colony Type') +
  scale_fill_manual(name = 'Colony Type',
                    breaks = c("Reference","Mutant"),
                    labels = c("Reference","Mutant"),
                    values = c("Reference" = "#9E9E9E",
                               "Mutant" = "#7B1FA2"),
                    guide = F) +
  coord_cartesian(xlim = c(0,650)) +
  facet_wrap(.~title, labeller = label_parsed) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.y = element_text(angle = -45, vjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))

##### 2C
fig.2c <- ggplot(virP1[virP1$ref_hrs > 0 & virP1$que_hrs > 0,]) +
  geom_tile(aes(x = sprintf('%0.1f',que_hrs), y = sprintf('%0.1f',ref_hrs), fill = round((es-1)*100,0)),
            col = 'black') +
  geom_text(aes(x = sprintf('%0.1f',que_hrs), y = sprintf('%0.1f',ref_hrs), label = round((es-1)*100,0)),
            col = 'black', size = 2) +
  scale_x_discrete(limits = sprintf('%0.1f',hhours[2:12])) +
  scale_y_discrete(limits = sprintf('%0.1f',sort(hhours[2:12], decreasing = T))) +
  scale_fill_gradient2(name = 'Fitness\nEffect %',
                       low = "#3F51B5", high = "#FFC107", mid = "white",
                       trans = 'pseudo_log',
                       midpoint = 0,
                       breaks = c(-50,-5,0,5,50,500)) +
  labs(title = 'Bimodal Colony-Size Distribution',
       x = expression(t["M"]~~(hours)),
       y = expression(t["R"]~~(hours))) +
  facet_wrap(.~text) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'right',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))

##### 2D
fig.2d <- ggplot(density.data[!is.na(density.data$colony) &
                                     density.data$hours > 0,],
                      aes(x = rnd_avg, y = colony, group = colony, fill = colony)) +
  geom_density_ridges(quantile_lines = TRUE,
                      scale = 3, alpha = 0.8, size = 0.3,
                      vline_size = 0.2, vline_color = "black") +
  labs(title = '',
       x = 'Colony Size (pixel counts)',
       y = 'Colony Type') +
  scale_fill_manual(name = 'Colony Type',
                    breaks = c("Reference","Mutant"),
                    labels = c("Reference","Mutant"),
                    values = c("Reference" = "#9E9E9E",
                               "Mutant" = "#7B1FA2"),
                    guide = F) +
  coord_cartesian(xlim = c(0,650)) +
  facet_wrap(.~title, labeller = label_parsed) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.y = element_text(angle = -45, vjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))

##### FINAL FIGURE
fig2 <- annotate_figure(ggpubr::ggarrange(fig.2a, fig.2b,
                                          fig.2c, fig.2d, 
                                          nrow = 2, ncol = 2,
                                          labels = c('A','B','C','D'),
                                          heights = c(1,1),
                                          font.label = list(face = 'bold', size = lbls, family = "sans"),
                                          hjust=-1),
                        bottom = text_grob(expression("t"["R"]~"= Reference colony size time ,"~"t"["M"]~"= Mutant colony size time"),
                                           face = "bold", family = "sans", size = titles))
ggsave(sprintf("%sFIGURE2.jpg",out_path), fig2,
       height = two.c*0.9, width = two.c, units = 'mm',
       dpi = 300)


##### FIGURE 3. SPECIFITICTY AND SENSITIVITY OF LID AND MCAT
load(sprintf('%sLID_SENSITIVITY_DATA.RData',out_path))
lid.sen.data <- dat.cnt2
load(sprintf('%sLID_SPECIFICITY_DATA.RData',out_path))
lid.spe.data <- stats.tmp

load(sprintf('%sMCAT_SENSITIVITY_DATA.RData',out_path))
mcat.sen.data <- dat.cnt2
load(sprintf('%sMCAT_SPECIFICITY_DATA.RData',out_path))
mcat.spe.data <- stats.tmp

##### 3A
fig.3a <- ggplot(mcat.spe.data[mcat.spe.data$hours == mcat.spe.data$cont_hrs,],
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

##### 3B
t <- lid.sen.data$Beneficial[1]
fig.3b <- ggplot(lid.sen.data) +
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

##### 3C
t <- mcat.sen.data$Beneficial[1]
fig.3c <- ggplot(mcat.sen.data) +
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

##### FINAL FIGURE
fig3 <- ggpubr::ggarrange(fig.3a,
                          ggpubr::ggarrange(fig.3b, fig.3c,
                                            ncol = 2,
                                            labels = c('B','C'),
                                            font.label = list(face = 'bold', size = lbls, family = "sans"),
                                            hjust=-1),
                          nrow = 2, heights = c(1.2,2),
                          labels = c('A',''),
                          font.label = list(face = 'bold', size = lbls, family = "sans"),
                          hjust=-1)

ggsave(sprintf("%sFIGURE3.jpg",out_path), fig3,
       height = 140, width = two.c, units = 'mm',
       dpi = 300)


##### FIGURE 4. SENSITIVITY WHEN UNDERLYING DFE IS RANDOM
load(file = sprintf('%sRND2_RESULTS.RData',out_path))

mcat.rnd2.res$effect <- factor(mcat.rnd2.res$effect, levels = c('Beneficial','Neutral','Deleterious'))
lid.rnd2.res$effect <- factor(lid.rnd2.res$effect, levels = c('Beneficial','Neutral','Deleterious'))
mcat.rnd2.res$truth <- factor(mcat.rnd2.res$truth, levels = c('Beneficial','Neutral','Deleterious'))
lid.rnd2.res$truth <- factor(lid.rnd2.res$truth, levels = c('Beneficial','Neutral','Deleterious'))

mcat.pie <- plyr::count(mcat.rnd2.res, vars = 'effect')
mcat.pie$freq <- mcat.pie$freq/sum(mcat.pie$freq) * 100
truth.pie <- plyr::count(mcat.rnd2.res, vars = 'truth')
truth.pie$freq <- truth.pie$freq/sum(truth.pie$freq) * 100
lid.pie <- plyr::count(lid.rnd2.res[lid.rnd2.res$ref == 0.25 & lid.rnd2.res$rep == 16 & lid.rnd2.res$attempt == 1,], vars = 'effect')
lid.pie$freq <- lid.pie$freq/sum(lid.pie$freq) * 100

##### 4A
fig.4a <- ggplot(truth.pie, aes(x = "", y = freq, fill = truth)) +
  geom_bar(stat = 'identity') +
  coord_polar("y", start = 0) +
  geom_label_repel(aes(label = sprintf("%0.2f%%",freq)),
                   position = position_stack(vjust = 0.5),
                   colour ='white',
                   label.size = 0.15,
                   show.legend = F) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107'),
                    guide = F) +
  labs(title = 'ACTUAL CLASSIFICATION') +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(size = titles, hjust = 0.5),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = 'bottom',
        strip.text = element_text(size = titles),
        legend.box.spacing = unit(0.5,"mm"))

##### 4B
fig.4b <- ggplot(lid.pie, aes(x = "", y = freq, fill = effect)) +
  geom_bar(stat = 'identity') +
  coord_polar("y", start = 0) +
  geom_label_repel(aes(label = sprintf("%0.2f%%",freq)),
                   position = position_stack(vjust = 0.5),
                   colour ='white',
                   label.size = 0.15,
                   show.legend = F) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107'),
                    guide = F) +
  labs(title = 'LID CLASSIFICATION') +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(size = titles, hjust = 0.5),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = 'bottom',
        strip.text = element_text(size = titles),
        legend.box.spacing = unit(0.5,"mm"))

##### 4C
fig.4c <- ggplot(mcat.pie, aes(x = "", y = freq, fill = effect)) +
  geom_bar(stat = 'identity') +
  coord_polar("y", start = 0) +
  geom_label_repel(aes(label = sprintf("%0.2f%%",freq)),
                   position = position_stack(vjust = 0.5),
                   colour ='white',
                   label.size = 0.15,
                   show.legend = F) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107'),
                    guide = F) +
  labs(title = 'MCAT CLASSIFICATION') +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(size = titles, hjust = 0.5),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = 'bottom',
        strip.text = element_text(size = titles),
        legend.box.spacing = unit(0.5,"mm"))

##### 4D
fig.4d <- ggplot(lid.rnd2.res[lid.rnd2.res$ref == 0.25 & lid.rnd2.res$rep == 16 & lid.rnd2.res$attempt == 1,],
                 aes(x = as.factor(hours), fill = truth)) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..),
            col = '#F5F5F5', size = 2,
            position = position_stack(vjust = 0.5)) +
  labs(title = 'TRUTH',
       x = 'Reference Population Time Point (hours)',
       y = 'Number of Mutants') +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107'),
                    guide = F) +
  scale_x_discrete(labels = sprintf('%0.1f',sort(unique(lid.rnd2.res$hours)))) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm"))

##### 4E
fig.4e <- ggplot(lid.rnd2.res[lid.rnd2.res$ref == 0.25 & lid.rnd2.res$rep == 16 & lid.rnd2.res$attempt == 1,],
                 aes(x = as.factor(hours), fill = effect)) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..),
            col = '#F5F5F5', size = 2,
            position = position_stack(vjust = 0.5)) +
  labs(title = 'LID',
       x = 'Reference Population Time Point (hours)',
       y = '') +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107'),
                    guide = F) +
  scale_x_discrete(labels = sprintf('%0.1f',sort(unique(lid.rnd2.res$hours)))) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm"))

##### 4F
fig.4f <- ggplot(mcat.rnd2.res,
                  aes(x = as.factor(hours), fill = effect)) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..),
            col = '#F5F5F5', size = 2,
            position = position_stack(vjust = 0.5)) +
  labs(title = 'MCAT',
       x = 'Reference Population Time Point (hours)',
       y = '') +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107'),
                    guide = F) +
  scale_x_discrete(labels = sprintf('%0.1f',sort(unique(mcat.rnd2.res$hours)))) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm"))

##### 4G
mcat.rnd2.res$effect <- factor(mcat.rnd2.res$effect, levels = c('Beneficial','Deleterious','Neutral'))
lid.rnd2.res$effect <- factor(lid.rnd2.res$effect, levels = c('Beneficial','Deleterious','Neutral'))
lid.rnd2.res$truth <- factor(lid.rnd2.res$truth, levels = c('Beneficial','Deleterious','Neutral'))

fig.4g <- ggplot(lid.rnd2.res[lid.rnd2.res$ref == 0.25 & lid.rnd2.res$rep == 16 & lid.rnd2.res$attempt == 1,]) +
  geom_histogram(aes(x = (es - 1) * 100, y=..count../5, fill = truth), binwidth = 10) +
  labs(title = 'TRUTH',
       x = 'Fitness Effect',
       y = 'Number of Mutants') +
  scale_x_continuous(breaks = seq(-1000,1000,100),
                     minor_breaks = seq(-1000,1000,10),
                     labels = paste(seq(-1000,1000,100),'%', sep = '')) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107')) +
  # facet_wrap(.~ref*rep, ncol = 8) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm"),
        legend.position = 'bottom',
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))

##### 4H
fig.4h <- ggplot(lid.rnd2.res[lid.rnd2.res$ref == 0.25 & lid.rnd2.res$rep == 16 & lid.rnd2.res$attempt == 1,]) +
  geom_histogram(aes(x = (es - 1) * 100, y=..count../5, fill = effect), binwidth = 10) +
  labs(title = 'LID',
       x = 'Fitness Effect',
       y = '') +
  scale_x_continuous(breaks = seq(-1000,1000,100),
                     minor_breaks = seq(-1000,1000,10),
                     labels = paste(seq(-1000,1000,100),'%', sep = '')) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107')) +
  # facet_wrap(.~ref*rep, ncol = 8) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm"),
        legend.position = 'bottom',
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))

##### 4I
fig.4i <- ggplot(mcat.rnd2.res) +
  geom_histogram(aes(x = (es - 1) * 100, y=..count../5, fill = effect), binwidth = 10) +
  labs(title = 'MCAT',
       x = 'Fitness Effect',
       y = '') +
  scale_x_continuous(breaks = seq(-1000,1000,100),
                     minor_breaks = seq(-1000,1000,10),
                     labels = paste(seq(-1000,1000,100),'%', sep = '')) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107')) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm"),
        legend.position = 'bottom',
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))

##### FINAL FIGURE
fig4 <- ggpubr::ggarrange(fig.4a, fig.4b, fig.4c,
                          fig.4d, fig.4e, fig.4f,
                          fig.4g, fig.4h, fig.4i,
                          nrow = 3,
                          ncol = 3,
                          common.legend = T,
                          legend = 'bottom',
                          labels = c('A','B','C','D','E','F','G','H','I'),
                          font.label = list(face = 'bold', size = lbls, family = "sans"),
                          hjust=-1)
ggsave(sprintf("%sFIGURE4.jpg",out_path), fig4,
       height = two.c, width = two.c, units = 'mm',
       dpi = 300)


##### FIGURE 5. SENSITIVITY'S RELATIONSHIP WITH REPLICATES AND REFERENCES
load(file = sprintf('%sLID_RNDVDATA.RData',out_path))
load(sprintf('%sREFREPDATA.RData',out_path))
refrep$power[refrep$hours < refrep$cont_hrs] <- refrep$Deleterious[refrep$hours < refrep$cont_hrs]/910 * 100
refrep$power[refrep$hours > refrep$cont_hrs] <- refrep$Beneficial[refrep$hours > refrep$cont_hrs]/910 * 100
refrep$abs_cen <- abs(1-refrep$cen) * 100
refrep$rep <- as.factor(refrep$rep)
# refrep$ref_prop <- as.factor(refrep$ref_prop)
refrep$ref <- refrep$ref_prop
refrep$ref_title <- sprintf('References/Plate=%0.2f%%', refrep$ref*100)
refrep$ref_title <- factor(refrep$ref_title, levels = unique(refrep$ref_title)[order(unique(refrep$ref))])

es = 5
fig5 <- ggplot(rnd.v.data[round(rnd.v.data$abs_cen) == es,],
               aes(x = rep,
                   y = power)) +
  geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1),
                size = 0.3, width = 0.5) +
  geom_line(stat = 'summary', fun = 'mean', group = 1,
            size = 0.3) +
  geom_point(stat = 'summary', fun = 'mean', group = 1,
             aes(col = 'Random CS Distribution'),
             size = 1) +
  geom_errorbar(data = refrep[round(refrep$abs_cen) == es,],
                aes(x = rep,
                    y = power),
                stat="summary", fun.data="mean_se", fun.args = list(mult = 1),
                size = 0.3, width = 0.5) +
  geom_line(data = refrep[round(refrep$abs_cen) == es,],
            aes(x = rep,
                y = power),
            stat = 'summary', fun = 'mean', group = 1,
            size = 0.3) +
  geom_point(data = refrep[round(refrep$abs_cen) == es,],
             stat = 'summary', fun = 'mean', group = 1,
             aes(x = rep,
                 y = power,
                 col = 'Bimodal CS Distribution'),
             size = 1) +
  scale_color_manual(name = 'Virtual Plates',
                     breaks = c('Bimodal CS Distribution', 'Random CS Distribution'),
                     values = c('#7C4DFF', '#E64A19')) +
  scale_y_continuous(breaks = seq(0,100,10),
                     labels = paste0(seq(0,100,10),'%',sep = '')) +
  labs(x = 'Replicates per strain',
       y = 'Sensitivity') +
  facet_wrap(.~ref_title, nrow = 1) +
  theme_linedraw() +
  theme( plot.title = element_text(size = titles),
         axis.title = element_text(size = titles),
         axis.text = element_text(size = txt),
         legend.title = element_text(size = titles),
         legend.text = element_text(size = txt),
         legend.key = element_rect(size = 0.5),
         legend.key.size = unit(.8, 'lines'),
         legend.position = 'bottom',
         legend.margin = margin(0.5,0.5,0.5,0.5, "mm"),
         strip.text = element_text(size = txt,
                                   margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%sFIGURE5.jpg",out_path), fig5,
       height = 70, width = two.c, units = 'mm',
       dpi = 300)
