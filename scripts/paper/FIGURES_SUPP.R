##### LID SUPPLEMENTARY FIGURES
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

#####  FIGURE_S4: COLONY SIZE DISTRIBUTIONS OF VIRTUAL PLATES WITH BIMODAL DISTRIBUTION
load(sprintf("/home/sbp29/R/Projects/adaptivefitness/output/4C4_FS_CC_VIR_PLATE1.Rdata"))
fit.all$colony[fit.all$colony == 'Query'] <- 'Mutant'

# fit.all$vir_plate <- NULL
fit.all <- fit.all[fit.all$hours != 0,]
fit.all <- fit.all[fit.all$cont_hrs != 0,]

# fit.all$vir_plate <- sprintf('"t"["R"]=="%0.1f"~"\n"~"t"["M"]=="%0.1f"', fit.all$cont_hrs, fit.all$hours)
# fit.all$vir_plate <- factor(fit.all$vir_plate, levels=c(sprintf('"t"["R"]=="%0.1f"~"\n"~"t"["M"]=="%0.1f"', sort(unique(fit.all$cont_hrs)), sort(unique(fit.all$hours)))))

fit.all$title1 <- sprintf('"t"["R"]=="%0.1fhrs"', fit.all$cont_hrs)
fit.all$title1 <- factor(fit.all$title1, levels = sprintf('"t"["R"]=="%0.1fhrs"', sort(unique(fit.all$cont_hrs))))
fit.all$title2 <- sprintf('"t"["M"]=="%0.1fhrs"', fit.all$hours)
fit.all$title2 <- factor(fit.all$title2, levels = sprintf('"t"["M"]=="%0.1fhrs"', sort(unique(fit.all$hours))))

fig.s4 <- ggplot(fit.all, aes(x = average, fill = colony)) +
  geom_density(stat = 'density', alpha = 0.8) +
  labs(x = 'Colony Size (pixel)',
       y = 'Density') +
  scale_fill_manual(name = 'Colony Type',
                    breaks = c("Reference","Mutant"),
                    labels = c("Reference","Mutant"),
                    values = c("Reference" = "#9E9E9E",
                               "Mutant" = "#7B1FA2")) +
  facet_wrap(~title1 + title2, labeller = labeller(.cols = label_parsed, .multi_line = T),
             nrow = 11) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt-2),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        strip.text = element_text(size = txt-1,
                                  margin = margin(0.1,0,0.1,0, "mm")),
        legend.key.size = unit(3, "mm"))

ggsave(sprintf("%sFIGURE_S4.jpg",out_path), fig.s4,
       height = 220, width = two.c, units = 'mm',
       dpi = 300)


##### FIGURE S5. ACCURACY OF BACKGROUND COLONY SIZES AND SOURCE NORMALIZATION
load(sprintf('%sSPATIAL.RData',out_path))
load(sprintf('%sRMSEDATA.RData',out_path))
load(sprintf('%sSOURCENORMALIZATIONDATA.RData',out_path))

##### S5A
spatial$plate <- as.factor(spatial$plate)
spatial$hours <- as.factor(spatial$hours)
my_comparisons <- list(c("NO-NORM", "RND"), c("NO-NORM", "MCAT"), c("NO-NORM", "LID-SN"),
                       c("NO-NORM", "LID-AC"), c("NO-NORM", "LID") )

fig.s5a <- ggplot(spatial[spatial$hours == 11.04 &
                            spatial$plate == 999,],
                  aes(x = name, y = cv, fill = plate)) +
  geom_hline(yintercept = median(spatial$cv[spatial$name == "LID" & spatial$plate == 999], na.rm = T),
             col = 'red', linetype = 'dashed') +
  geom_boxplot(outlier.colour = "black", outlier.shape = NA)  +
  labs(x = 'Fitness Estimation Strategy',
       y = 'CV%\n(SD/Mean * 100)') +
  scale_fill_discrete(name = '',
                      breaks = c(1,2,3,4,999),
                      labels = c('Plate 1', 'Plate 2', 'Plate 3',
                                 'Plate 4', 'All Plates'),
                      guide = F) +
  stat_compare_means(comparisons = my_comparisons,
                     size = 1.5, label.y = c(32,35,38,41,44)) +
  # facet_wrap(.~plate) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt))

##### S5B
rmse <- rmse[rmse$hour > 0,]
fig.s5b <- ggplot(rmse[order(rmse$method, decreasing = T),]) +
  geom_point(aes(x = hour, y = per, fill = as.factor(method)),
             size = 2, alpha = 0.9, shape = 21, col = 'black') +
  labs(title = 'RMSE',
       x = 'Time (hour)', y = 'RMSE %') +
  scale_x_continuous(breaks = seq(0,12,2)) +
  scale_fill_manual(name = 'Strategy',
                    breaks = c('3','2','5','4','6'),
                    limits = c('3','2','5','4','6'),
                    labels = c('LID','LID-SN','LID-AC','MCAT','RND'),
                    values = c('3'='#1976D2','2'='#03A9F4','5'='#BBDEFB',
                               '4'='#009688','6'='#795548')) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        # plot.title = element_text(size = titles,
        #                           hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = c(0.8, 0.75),
        legend.background = element_rect(fill = 'transparent'),
        legend.key = element_rect(size = 0.5),
        legend.key.size = unit(.8, 'lines')) +
  guides(fill = guide_legend(ncol=2))

##### S5C
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

##### FINAL FIGURE
fig.s5ab <- ggarrange(fig.s5a, fig.s5b, 
                     nrow = 1, ncol = 2,
                     labels = c('A','B'),
                     label.args = list(gp=gpar(font = 2, fontsize = lbls, family = "sans"),
                                       hjust = -1))
fig.s5c <- ggarrange(sn.raw,sn.sn,sn.nosn,sn.nocc,sn.bean,
                   nrow = 1, ncol = 5,
                   labels = c('C','','','',''),
                   label.args = list(gp=gpar(font = 2, fontsize = lbls, family = "sans"),
                                     hjust = -1))
fig.s5 <- ggpubr::ggarrange(fig.s5ab, fig.s5c,
                          ncol = 1, nrow = 2,
                          heights = c(1.5,2))
ggsave(sprintf("%sFIGURE_S5.jpg",out_path), fig.s5,
       height = 120, width = two.c, units = 'mm',
       dpi = 300)


##### FIGURE S6. SENSITIVITY AND DIFFERENT SCHEMES OF NORMALIZATION
load(sprintf('%sSEN_METHODS.RData',out_path))
sen.methods$method[sen.methods$method == 'No Normalization'] <- 'NO-NORM'

t <- sen.methods$Beneficial[1]
fig.s6 <- ggplot(sen.methods) +
  geom_area(aes(x = (cen-1)*100, y = Neutral, fill = 'Neutral'), alpha = 1) +
  geom_area(aes(x = (cen-1)*100, y = Beneficial, fill = 'Beneficial'), alpha = 1) +
  geom_area(aes(x = (cen-1)*100, y = Deleterious, fill = 'Deleterious'), alpha = 1) +
  geom_vline(xintercept = seq(-2,2,0.025)*100, col = '#757575', lwd = 0.5, alpha =0.5) +
  geom_hline(yintercept = c(t * seq(0,1,0.05)), col = '#757575', lwd = 0.5, alpha =0.5) +
  geom_vline(xintercept = c(-5,5), col = 'red', linetype = 'dashed', alpha = 0.8) +
  labs(x = 'Fitness Effect',
       y = 'Sensitivity') +
  scale_y_continuous(breaks = c(t * seq(0,1,0.1)),
                     minor_breaks = c(t * seq(0,1,0.05)),
                     labels = c('0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')) +
  scale_x_continuous(breaks = seq(-2,2,0.05)*100,
                     minor_breaks = seq(-2,2,0.025)*100) +
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
  facet_wrap(.~method) +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(5, "mm")) +
  coord_cartesian(xlim = c(-10,10),
                  ylim = c(0,t))

ggsave(sprintf("%sFIGURE_S6.jpg",out_path), fig.s6,
       height = 160, width = two.c, units = 'mm',
       dpi = 300)


##### FIGURE S7. SENSITIVITY'S RELATIONSHIP WITH REFERENCES AND REPLICATES
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

##### S7A
es = 7
fig.s7a <- ggplot(rnd.v.data[round(rnd.v.data$abs_cen) == es,],
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
         legend.text = element_text(size = txt+1),
         legend.key = element_rect(size = 0.5),
         legend.key.size = unit(.8, 'lines'),
         legend.position = 'bottom',
         legend.margin = margin(0.5,0.5,0.5,0.5, "mm"),
         strip.text = element_text(size = txt+2,
                                   margin = margin(0.1,0,0.1,0, "mm")))

##### S7B
t <- 100
fig.s7b <- ggplot(rnd.v.data) +
  geom_area(aes(x = (es-1)*100, y = Neutral, fill = 'Neutral'), alpha = 1) +
  geom_area(aes(x = (es-1)*100, y = Beneficial, fill = 'Beneficial'), alpha = 1) +
  geom_area(aes(x = (es-1)*100, y = Deleterious, fill = 'Deleterious'), alpha = 1) +
  geom_vline(xintercept = seq(-2,2,0.025)*100, col = '#757575', lwd = 0.25, alpha =0.5) +
  geom_hline(yintercept = c(t * seq(0,1,0.05)), col = '#757575', lwd = 0.25, alpha =0.5) +
  geom_vline(xintercept = c(-5,5), col = 'red', lwd = 0.25, linetype = 'dashed') +
  # geom_text(x=0, y=t*0.9, label="LID", col = 'white', size = 3) +
  labs(x = 'Fitness Effect',
       y = 'Sensitivity') +
  scale_y_continuous(breaks = c(t * seq(0,1,0.2)),
                     minor_breaks = c(t * seq(0,1,0.05)),
                     labels = paste(t * seq(0,1,0.2), '%', sep = "")) +
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
  facet_wrap(.~ref_title*rep_title, ncol = 8) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt-2),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(4, "mm"),
        legend.position = 'bottom',
        legend.margin = margin(0.5,0.5,0.5,0.5, "mm"),
        strip.text = element_text(size = txt-2,
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(xlim = c(-10,10),
                  ylim = c(0,t))

##### FINAL FIGURE
fig.s7 <- ggarrange(fig.s7a, fig.s7b,
                    labels = c('A','B'),
                    label.args = list(gp=gpar(font = 2, fontsize = lbls, family = "sans"),
                                      hjust=-1),
                    heights = c(1,2))
ggsave(sprintf("%sFIGURE_S7.jpg",out_path), fig.s7,
       height = two.c, width = two.c, units = 'mm',
       dpi = 300)


##### FIGURE S8. RMSE'S RELATIONSHIP WITH REFERENCES AND REPLICATES
load(sprintf('%sRMSE_REFREP.RData',out_path))
rmse.rr$rep <- as.factor(rmse.rr$rep)
rmse.rr$ref <- as.factor(rmse.rr$ref)

fig.s8 <- ggplot(rmse.rr[rmse.rr$rep == 16 & rmse.rr$hours > 2,],
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

ggsave(sprintf("%sFIGURE_S8.jpg",out_path), fig.s8,
       height = 70, width = one.c, units = 'mm',
       dpi = 300)


##### FIGURE S9. SPECIFICITY OF DIFFERENT WAYS OF USING LID
load(sprintf('%sSPE_METHODS.RData',out_path))
load(sprintf('%sSPECIFICITY.RData',out_path))

##### S9A
spe.methods$abs_cen <- abs(1-spe.methods$cen) * 100
spe.met <- spe.methods[spe.methods$cont_hrs == spe.methods$hours & spe.methods$p <= 0.05,]

fig.s9a <- ggplot(spe.met,
                  aes(x = method, y = round((1-fpr),4)*100)) +
  geom_boxplot(fill = 'white', outlier.shape = NA) +
  labs(x = 'Fitness Estimation Strategy',
       y = 'Specificity') +
  scale_x_discrete(limits = c('LID',
                              'LID-AC',
                              'LID-SN',
                              'No Normalization'),
                   labels = c('LID',
                              'LID-AC',
                              'LID-SN',
                              'NO-NORM')) +
  scale_y_continuous(breaks = seq(0,100,2),
                     minor_breaks = seq(0,105,1),
                     labels = paste(seq(0,100,2),'%',sep='')) +
  coord_cartesian(ylim = c(90,100)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt))

##### S9B
spe.data$abs_cen <- abs(1-spe.data$cen) * 100
spe.data$rep <- as.factor(spe.data$rep)
spe.data$ref <- as.factor(spe.data$ref)

spe.ref <- spe.data[spe.data$cont_hrs == spe.data$hours &
                      round(spe.data$abs_cen) <= 5 &
                      spe.data$p <= 0.05,]

fig.s9b <- ggplot(spe.ref,
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

##### FINAL FIGURE
fig.s9 <- ggarrange(fig.s9a, fig.s9b,
                     nrow = 1, ncol = 2,
                     labels = c('A','B'),
                     label.args = list(gp=gpar(font = 2, fontsize = lbls, family = "sans"),
                                       hjust=-1))
ggsave(sprintf("%sFIGURE_S9.jpg",out_path), fig.s9,
       height = 80, width = two.c, units = 'mm',
       dpi = 300)


##### FIGURE S10. PLATEMAPS
p2c.dat <- dbGetQuery(conn, 'select * from 4C4_pos2coor a, 4C4_pos2orf_name b
                      where a.pos = b.pos')
p2c.dat$colony[is.na(p2c.dat$orf_name)] <- 'Gap'
p2c.dat$colony[p2c.dat$orf_name == 'BF_control'] <- 'Reference'
p2c.dat$colony[p2c.dat$orf_name != 'BF_control'] <- 'Mutant'
save(p2c.dat, file = sprintf('%sPLATEMAPS.RData',out_path))

load(sprintf('%sPLATEMAPS.RData',out_path))

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

##### FINAL FIGURE
fig.s10 <- ggarrange(stocks, wctp1, up1tp2, up2,
                      ncol = 1, nrow = 4)
ggsave(sprintf("%sFIGURE_S10.jpg",out_path), fig.s10,
       height = two.c, width = two.c, units = 'mm',
       dpi = 300)

##### FIGURE S11. VIRTUAL PLATE WITH BIMODAL DISTRIBUTION EXAMPLE
load(sprintf('%sVP1EGDATA.RData',out_path))
vir1.eg$colony[vir1.eg$colony == 'Query'] <- 'Mutant'
vir1.eg$colony[vir1.eg$orf_name == ''] <- 'Gap'

unique(vir1.eg$kind)

##### S11A
fig.s11a <- ggplot(vir1.eg[vir1.eg$kind == "1. Time = 2.9 hr",]) +
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

##### S11B
vir1.eg$kind[vir1.eg$kind == "2. Time = 11.04 hr"] <- "2. Time = 11.0 hr"
fig.s11b <- ggplot(vir1.eg[vir1.eg$kind == "2. Time = 11.0 hr",]) +
  geom_point(aes(x = col, y = row, size = average, col = colony)) +
  geom_point(data = vir1.eg[vir1.eg$colony == 'Gap' & vir1.eg$kind == "2. Time = 11.0 hr",],
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

##### S11C
fig.s11c <- ggplot(vir1.eg[vir1.eg$kind == "3. Virtual Plate",]) +
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

##### FINAL FIGURE
fig.s11 <- ggarrange(fig.s11a, fig.s11b, fig.s11c,
                     labels = c('A','B','C'),
                     label.args = list(gp=gpar(font = 2, fontsize = lbls, family = "sans"),
                                       hjust=0),
                     nrow = 1, ncol = 3)

ggsave(sprintf("%sFIGURE_S11.jpg",out_path), fig.s11,
       height = 60, width = two.c, units = 'mm',
       dpi = 300)
