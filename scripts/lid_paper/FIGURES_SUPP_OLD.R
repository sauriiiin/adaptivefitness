##### LID SUPPLEMENTARY FIGURES
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 07/14/2020

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

##### FIGURE S91. Prescreens Plate Effect
load(sprintf('%s4C4PSDATA.RData',out_path))

ps.plateeffect <- ggplot(data.ps,
                         aes(x = source, y = average-median_cs, fill = source)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, col = 'red', linetype = 'dashed') +
  labs(x = 'Source Plate',
       y = 'Median Substracted Pixel Counts') +
  scale_x_discrete(breaks=c("1TL","2TR","3BL","4BR"),
                   labels=c("Top Left","Top Right","Bottom Left","Bottom Right")) +
  scale_fill_manual(breaks = c('1TL','2TR','3BL','4BR'),
                    values = c('1TL' = '#303F9F',
                               '2TR' = '#D32F2F',
                               '3BL' = '#FFEB3B',
                               '4BR' = '#388E3C')) +
  stat_compare_means(label.x = 2.5, label.y = 700, hjust = 0.5, size = 1.5) +
  facet_wrap(.~stage) +
  theme_linedraw()+
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = txt,
                                  margin = margin(2,0,2,0, "mm")))
ggsave(sprintf("%sFIGURE_S91.jpg",out_path), ps.plateeffect,
       height = one.c, width = one.c, units = 'mm',
       dpi = 300)


##### FIGURE S2. PLATEMAPS
p2c.dat <- dbGetQuery(conn, 'select * from 4C4_pos2coor a, 4C4_pos2orf_name b
                      where a.pos = b.pos')
p2c.dat$colony[is.na(p2c.dat$orf_name)] <- 'Gap'
p2c.dat$colony[p2c.dat$orf_name == 'BF_control'] <- 'Reference'
p2c.dat$colony[p2c.dat$orf_name != 'BF_control'] <- 'Mutant'

stocks <- ggplot(p2c.dat[p2c.dat$density == 384,],
                 # aes(x = col, y = row, col = orf_name)) +
                 aes(x = col, y = row, col = colony)) +
  geom_point(size = 1, shape = 15) +
  # scale_color_discrete(guide = F) +
  scale_y_reverse() +
  scale_color_manual(name = 'Colony Type',
                     breaks = c("Reference","Mutant"),
                     labels = c("Reference","Mutant"),
                     values = c("Reference" = "#9E9E9E",
                                "Mutant" = "#7B1FA2",
                                "Gap" = 'transparent'),
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
                # aes(x = col, y = row, col = orf_name)) +
                aes(x = col, y = row, col = colony)) +
  geom_point(size = 1, shape = 16) +
  # scale_color_discrete(guide = F) +
  scale_y_reverse() +
  scale_color_manual(name = 'Colony Type',
                     breaks = c("Reference","Mutant"),
                     labels = c("Reference","Mutant"),
                     values = c("Reference" = "#9E9E9E",
                                "Mutant" = "#7B1FA2",
                                "Gap" = 'transparent'),
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
                 # aes(x = col, y = row, col = orf_name)) +
                 aes(x = col, y = row, col = colony)) +
  geom_point(size = 0.5, shape = 16) +
  # scale_color_discrete(guide = F) +
  scale_y_reverse() +
  scale_color_manual(name = 'Colony Type',
                     breaks = c("Reference","Mutant"),
                     labels = c("Reference","Mutant"),
                     values = c("Reference" = "#9E9E9E",
                                "Mutant" = "#7B1FA2",
                                "Gap" = 'transparent'),
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
              # aes(x = col, y = row, col = orf_name)) +
              aes(x = col, y = row, col = colony)) +
  geom_point(size = 0.2, shape = 16) +
  # scale_color_discrete(guide = F) +
  scale_y_reverse() +
  scale_color_manual(name = 'Colony Type',
                     breaks = c("Reference","Mutant"),
                     labels = c("Reference","Mutant"),
                     values = c("Reference" = "#9E9E9E",
                                "Mutant" = "#7B1FA2",
                                "Gap" = 'transparent')) +
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

ggsave(sprintf("%sFIGURE_S2.jpg",out_path), plt.maps,
       height = two.c, width = two.c, units = 'mm',
       dpi = 300)


#####  FIGURE_S4: VIRTUAL PLATES CS DISTRIBUTION
load(sprintf("/home/sbp29/R/Projects/adaptivefitness/output/4C4_FS_CC_VIR_PLATE1.Rdata"))
fit.all$colony[fit.all$colony == 'Query'] <- 'Mutant'

vir.cs <- ggplot(fit.all, aes(x = average, fill = colony)) +
  geom_density(stat = 'density', alpha = 0.8) +
  labs(x = 'Colony Size (pixel)',
       y = 'Density') +
  scale_fill_manual(name = 'Colony Type',
                    breaks = c("Reference","Mutant"),
                    labels = c("Reference","Mutant"),
                    values = c("Reference" = "#9E9E9E",
                               "Mutant" = "#7B1FA2")) +
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

ggsave(sprintf("%sFIGURE_S4.jpg",out_path), vir.cs,
       height = two.c, width = two.c, units = 'mm',
       dpi = 300)

##### FIGURE S5: VIRTUAL PLATE 1 EXAMPLE
load(sprintf('%sVP1EGDATA.RData',out_path))
vir1.eg$colony[vir1.eg$colony == 'Query'] <- 'Mutant'
vir1.eg$colony[vir1.eg$orf_name == ''] <- 'Gap'

vir1plt.a <- ggplot(vir1.eg[vir1.eg$kind == "1. Time = 2.9 hr",]) +
  geom_point(aes(x = col, y = row, size = average, col = colony)) +
  geom_point(data = vir1.eg[vir1.eg$colony == 'Gap' & vir1.eg$kind == "1. Time = 2.9 hr",],
             aes(x = col, y =row, col = colony), shape = 4) +
  scale_color_manual(name = 'Colony Type',
                     breaks = c("Reference","Mutant","Gap"),
                     values = c("Reference" = "#9E9E9E",
                                "Mutant" = "#7B1FA2",
                                "Gap" = "#D32F2F"),
                     guide = F) +
  scale_size_continuous(range = c(0,2.5), guide = F) +
  scale_x_continuous(breaks = seq(0,96,4)) +
  scale_y_continuous(breaks = seq(0,64,4),trans = 'reverse') +
  labs(x = 'Column', y = 'Row') +
  coord_cartesian(xlim = c(10,33), ylim = c(10,25)) +
  facet_wrap(.~kind, nrow = 1) +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        strip.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"))

vir1plt.b <- ggplot(vir1.eg[vir1.eg$kind == "2. Time = 11.04 hr",]) +
  geom_point(aes(x = col, y = row, size = average, col = colony)) +
  geom_point(data = vir1.eg[vir1.eg$colony == 'Gap' & vir1.eg$kind == "2. Time = 11.04 hr",],
             aes(x = col, y =row, col = colony), shape = 4) +
  scale_color_manual(name = 'Colony Type',
                     breaks = c("Reference","Mutant","Gap"),
                     values = c("Reference" = "#9E9E9E",
                                "Mutant" = "#7B1FA2",
                                "Gap" = "#D32F2F")) +
  scale_size_continuous(name = 'Colony Size (pixels)',
                        range = c(0,2.5)) +
  scale_x_continuous(breaks = seq(0,96,4)) +
  scale_y_continuous(breaks = seq(0,64,4),trans = 'reverse') +
  labs(x = 'Column', y = 'Row') +
  coord_cartesian(xlim = c(10,33), ylim = c(10,25)) +
  facet_wrap(.~kind, nrow = 1) +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        strip.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"))

vir1plt.c <- ggplot(vir1.eg[vir1.eg$kind == "3. Virtual Plate",]) +
  geom_point(aes(x = col, y = row, size = average, col = colony)) +
  geom_point(data = vir1.eg[vir1.eg$colony == 'Gap' & vir1.eg$kind == "3. Virtual Plate",],
             aes(x = col, y =row, col = colony), shape = 4) +
  scale_color_manual(name = 'Colony Type',
                     breaks = c("Reference","Mutant","Gap"),
                     values = c("Reference" = "#9E9E9E",
                                "Mutant" = "#7B1FA2",
                                "Gap" = "#D32F2F"),
                     guide = F) +
  scale_size_continuous(range = c(0,2.5), guide = F) +
  scale_x_continuous(breaks = seq(0,96,4)) +
  scale_y_continuous(breaks = seq(0,64,4),trans = 'reverse') +
  labs(x = 'Column', y = 'Row') +
  coord_cartesian(xlim = c(10,33), ylim = c(10,25)) +
  facet_wrap(.~kind, nrow = 1) +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        strip.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"))

vir1plt <- ggarrange(vir1plt.a, vir1plt.b, vir1plt.c,
                     labels = c('a','b','c'),
                     label.args = list(gp=gpar(font = 2, fontsize = lbls),
                                       hjust=0),
                     nrow = 1, ncol = 3)

ggsave(sprintf("%sFIGURE_S5.jpg",out_path), vir1plt,
       height = 60, width = two.c, units = 'mm',
       dpi = 300)


##### FIGURE S10 - RMSE WITH CHANGE IN REF PROP
# No. of replicates do not affect rmse
load(sprintf('%sRMSE_REFREP.RData',out_path))
rmse.rr$rep <- as.factor(rmse.rr$rep)
rmse.rr$ref <- as.factor(rmse.rr$ref)

rmse.r <- ggplot(rmse.rr[rmse.rr$rep == 16 & rmse.rr$hours > 2,],
                 aes(x = hours, y = per, fill = ref), col = 'black') + 
  geom_point(size = 2, alpha = 0.9, shape = 21) +
  labs(x = 'Time (hour)', y = 'RMSE %') +
  scale_x_continuous(breaks = seq(0,12,2)) +
  scale_y_continuous(breaks = seq(0,100,2),
                     minor_breaks = seq(0,105,1),
                     labels = paste(seq(0,100,2),'%',sep='')) +
  scale_fill_manual(name = 'Reference\nProportion',
                     breaks = c(0.0625,0.125,0.1875,0.25),
                     values = c('#607D8B','#757575','#BDBDBD','#FFFFFF'),
                     labels = c(sprintf('%0.2f%%',c(0.0625,0.125,0.1875,0.25)*100))) +
  # coord_cartesian(ylim = c(90,100)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt))

ggsave(sprintf("%sFIGURE_S10.jpg",out_path), rmse.r,
       height = 70, width = one.c, units = 'mm',
       dpi = 300)


#### FIGURE S11. SPECIFICITY
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
  geom_boxplot(outlier.shape = NA) +
  labs(x = 'Number of Replicates',
       y = 'Specificity') +
  scale_x_discrete(breaks = seq(0,16,2)) +
  scale_y_continuous(breaks = seq(0,100,2),
                     minor_breaks = seq(0,105,1),
                     labels = paste(seq(0,100,2),'%',sep='')) +
  scale_fill_manual(name = 'Reference\nProportion',
                    breaks = c(0.0625,0.125,0.1875,0.25),
                    values = c('#607D8B','#757575','#BDBDBD','#FFFFFF'),
                    labels = c(sprintf('%0.2f%%',c(0.0625,0.125,0.1875,0.25)*100))) +
  coord_cartesian(ylim = c(90,100)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt))

## SPECIFICITY OF ALL TYPES OF NORMALIZATION
load(sprintf('%sSPE_METHODS.RData',out_path))

spe.methods$abs_cen <- abs(1-spe.methods$cen) * 100

spe.met <- spe.methods[spe.methods$cont_hrs == spe.methods$hours &
                         # round(spe.methods$abs_cen) <= 5 &
                         spe.methods$p <= 0.05,]

spe.m <- ggplot(spe.met,
                aes(x = method, y = round((1-fpr),4)*100)) +
  geom_boxplot(fill = 'white', outlier.shape = NA) +
  labs(x = 'Bias Correction Method',
       y = 'Specificity') +
  scale_x_discrete(limits = c('No Normalization',
                              'LID-SN',
                              'LID-AC',
                              'LID')) +
  scale_y_continuous(breaks = seq(0,100,2),
                     minor_breaks = seq(0,105,1),
                     labels = paste(seq(0,100,2),'%',sep='')) +
  coord_cartesian(ylim = c(90,100)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt))

fig_s11 <- ggarrange(spe.m,spe.r,
                    nrow = 1, ncol = 2,
                    labels = c('a','b'),
                    label.args = list(gp=gpar(font = 2, fontsize = lbls),
                                      hjust=-1))
ggsave(sprintf("%sFIGURE_S11.jpg",out_path), fig_s11,
       height = 80, width = two.c, units = 'mm',
       dpi = 300)
