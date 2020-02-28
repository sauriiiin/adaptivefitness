##### LID PAPER FIGURES NEW
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 09/28/2019

##### INITIALIZE
library(RMariaDB)
library(readxl)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(reshape2)
library(grid)
library(tidyverse)
# library(egg)
source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

load("figs/lid_paper/fitness.RData")

##### INITIALIZE
expt_name = '4C4_FS_CC'
expt = 'VAL2'
out_path = 'figs/lid_paper/4C4/';

density = 6144;
# tbl_pval = sprintf('%s_%d_PVALUE',expt_name,density)
# hours = dbGetQuery(conn, sprintf('select distinct hours from %s order by hours asc', tbl_pval))
# pvals = seq(0,1,0.01)

# tbl_fit_nonrm = sprintf('%s_NONORM_%d_FITNESS',substr(expt_name,1,7),density);
# tbl_fit_nil = sprintf('%s_NIL_%d_FITNESS',substr(expt_name,1,7),density);
# tbl_fit_ncc = sprintf('%s_%d_FITNESS',substr(expt_name,1,7),density);
# tbl_fit = sprintf('%s_%d_FITNESS',expt_name,density);
# tbl_fit_bean = sprintf('%s_BEAN_%d_FITNESS',substr(expt_name,1,7),density);
# tbl_p2o = '4C3_pos2orf_name3';
# tbl_bpos = '4C3_borderpos';
# 
# p2c_info = NULL
# p2c_info[1] = '4C3_pos2coor6144'
# p2c_info[2] = '6144plate'
# p2c_info[3] = '6144col'
# p2c_info[4] = '6144row'
# 
# p2c = dbGetQuery(conn, sprintf('select * from %s a order by a.%s, a.%s, a.%s',
#                                p2c_info[1],
#                                p2c_info[2],
#                                p2c_info[3],
#                                p2c_info[4]))
# 
# n_plates = dbGetQuery(conn, sprintf('select distinct %s from %s a order by %s asc',
#                                     p2c_info[2],
#                                     p2c_info[1],
#                                     p2c_info[2]))
# 
# dat.nonrm <- dbGetQuery(conn, sprintf('select b.*, a.orf_name, a.hours, a.bg, a.average, a.fitness
#                                     from %s a, %s b
#                                     where a.pos = b.pos and average is not NULL and bg is not NULL and orf_name is not NULL
#                                     order by a.hours, b.%s, b.%s, b.%s',
#                                     tbl_fit_nonrm,
#                                     p2c_info[1],
#                                     p2c_info[2],p2c_info[3],p2c_info[4]))
# dat.nonrm$method <- 'nonorm'
# 
# dat.nil <- dbGetQuery(conn, sprintf('select b.*, a.orf_name, a.hours, a.bg, a.average, a.fitness
#                                     from %s a, %s b
#                                     where a.pos = b.pos and average is not NULL and bg is not NULL and orf_name is not NULL
#                                     order by a.hours, b.%s, b.%s, b.%s',
#                                       tbl_fit_nil,
#                                       p2c_info[1],
#                                       p2c_info[2],p2c_info[3],p2c_info[4]))
# dat.nil$method <- 'nil'
# 
# dat.ncc <- dbGetQuery(conn, sprintf('select b.*, a.orf_name, a.hours, a.bg, a.average, a.fitness
#                                     from %s a, %s b
#                                     where a.pos = b.pos and average is not NULL and bg is not NULL and orf_name is not NULL
#                                     order by a.hours, b.%s, b.%s, b.%s',
#                                     tbl_fit_ncc,
#                                     p2c_info[1],
#                                     p2c_info[2],p2c_info[3],p2c_info[4]))
# dat.ncc$method <- 'ncc'
# 
# dat.lid <- dbGetQuery(conn, sprintf('select b.*, a.orf_name, a.hours, a.bg, a.average, a.fitness
#                                     from %s a, %s b
#                                     where a.pos = b.pos and average is not NULL and bg is not NULL and orf_name is not NULL
#                                     order by a.hours, b.%s, b.%s, b.%s',
#                                     tbl_fit,
#                                     p2c_info[1],
#                                     p2c_info[2],p2c_info[3],p2c_info[4]))
# dat.lid$method <- 'lid'
# 
# dat.bean <- dbGetQuery(conn, sprintf('select b.*, a.orf_name, a.hours, a.bg, a.average, a.fitness
#                                     from %s a, %s b
#                                     where a.pos = b.pos and average is not NULL and bg is not NULL and orf_name is not NULL
#                                     order by a.hours, b.%s, b.%s, b.%s',
#                                    tbl_fit_bean,
#                                    p2c_info[1],
#                                    p2c_info[2],p2c_info[3],p2c_info[4]))
# dat.bean$method <- 'bean'
# 
# dat.all <- rbind(dat.nonrm, dat.bean, dat.nil, dat.ncc, dat.lid)
# 
# dat.all$source[dat.all$row%%2==1 & dat.all$col%%2==1] = '1TL'
# dat.all$source[dat.all$row%%2==0 & dat.all$col%%2==1] = '3BL'
# dat.all$source[dat.all$row%%2==1 & dat.all$col%%2==0] = '2TR'
# dat.all$source[dat.all$row%%2==0 & dat.all$col%%2==0] = '4BR'
# 
# dat.all$colony[dat.all$orf_name == 'BF_control'] = 'Reference'
# dat.all$colony[dat.all$orf_name != 'BF_control'] = 'Query'
# dat.all$colony[is.na(dat.all$orf_name)] = 'Gap'
# save(dat.all, file = "figs/lid_paper/fitness.RData")

##### RMSE OF PREDICTION
# dat.all$average[dat.all$average == 0 & !is.na(dat.all$average)] <- NA
rmse <- NULL
i = 1
for (m in unique(dat.all$method)) {
  if (m == 'nonorm') {
    for (hr in unique(dat.all$hours[dat.all$method == m])) {
      dat.all$bg[dat.all$hours == hr & dat.all$method == m] <- sample(dat.all$average[dat.all$hours == hr & dat.all$method == m])
      rmse$method[i] <- m
      rmse$hour[i] <- hr
      rmse$avg[i] <- mean(dat.all$average[dat.all$hours == hr & dat.all$method == m], na.rm = T)
      rmse$rmse[i] <- sqrt(mean((dat.all$average[dat.all$hours == hr & dat.all$method == m] - dat.all$bg[dat.all$hours == hr & dat.all$method == m])^2,na.rm = T))
      i = i + 1
    }
  } else {
    for (hr in unique(dat.all$hours[dat.all$method == m])) {
      rmse$method[i] <- m
      rmse$hour[i] <- hr
      rmse$avg[i] <- mean(dat.all$average[dat.all$hours == hr & dat.all$method == m], na.rm = T)
      rmse$rmse[i] <- sqrt(mean((dat.all$average[dat.all$hours == hr & dat.all$method == m] - dat.all$bg[dat.all$hours == hr & dat.all$method == m])^2,na.rm = T))
      i = i + 1
    }
  }
}
rmse <- data.frame(rmse)
rmse$per <- rmse$rmse/rmse$avg * 100

ggplot(rmse[rmse$hour > 0,]) +
  geom_point(aes(x = avg, y = rmse, fill = method),
             shape = 21, size = 5, alpha = 0.8, col = 'black') +
  scale_x_continuous(breaks = seq(0,500,50), minor_breaks = seq(0,500,25)) +
  scale_y_continuous(breaks = seq(0,100,10), minor_breaks = seq(0,100,5)) +
  scale_fill_manual(name = 'Method',
                    breaks = c('nonorm','bean','nil','ncc','lid'),
                    values = c('nonorm' = '#212121',
                               'nil' = '#448AFF',
                               'ncc' = '#7C4DFF',
                               'lid' = '#388E3C',
                               'bean' = '#FFA000'),
                    labels = c('RND', 'MCAT', 'LID-SN', 'LID-CC', 'LID')) +
  labs(title = 'Compairing RMSE',
       x = 'Mean Pixel Count',
       y = 'RMSE') +
  theme_linedraw() +
  coord_cartesian(xlim = c(150, 400),
                  ylim = c(0, 80))
ggsave(sprintf("%srmse.jpg",out_path),
       width = 6, height = 5,
       dpi = 300)

ggplot(rmse[rmse$hour > 0,]) +
  geom_point(aes(x = avg, y = rmse/avg * 100, fill = method),
             shape = 21, size = 5, alpha = 0.8, col = 'black') +
  scale_x_continuous(breaks = seq(0,500,50), minor_breaks = seq(0,500,25)) +
  scale_y_continuous(breaks = seq(0,100,10), minor_breaks = seq(0,100,5)) +
  scale_fill_manual(name = 'Method',
                    breaks = c('nonorm','bean','nil','ncc','lid'),
                    values = c('nonorm' = '#212121',
                               'nil' = '#448AFF',
                               'ncc' = '#7C4DFF',
                               'lid' = '#388E3C',
                               'bean' = '#FFA000'),
                    labels = c('RND', 'MCAT', 'LID-SN', 'LID-CC', 'LID')) +
  labs(title = 'Compairing RMSE',
       x = 'Mean Pixel Count',
       y = 'RMSE %') +
  theme_linedraw() +
  coord_cartesian(xlim = c(150, 400),
                  ylim = c(0, 40))
ggsave(sprintf("%srmse2.jpg",out_path),
       width = 6, height = 5,
       dpi = 300)

##### SOURCE NORMALIZATION

sn.nonorm <- ggplot() +
  geom_line(data = dat.all[dat.all$method == 'nonorm' & dat.all$hours == 18,],
            aes(x = fitness, col = source), stat = 'density', lwd = 1.2) +
  scale_colour_manual(name="Source",
                      values=c("1TL"="#D32F2F","2TR"="#536DFE","3BL"="#388E3C","4BR"="#795548"),
                      breaks=c("1TL","2TR","3BL","4BR"),
                      labels=c("Top Left","Top Right","Bottom Left","Bottom Right")) +
  labs(title = 'Raw Colony Size',
       x = 'Colony Size (Pix)',
       y = 'Density') +
  theme_linedraw() +
  coord_cartesian(xlim = c(200, 600))

sn.nil <- ggplot() +
  geom_line(data = dat.all[dat.all$method == 'nil' & dat.all$hours == 18,],
            aes(x = fitness, col = source), stat = 'density', lwd = 1.2) +
  scale_colour_manual(name="Source",
                      values=c("1TL"="#D32F2F","2TR"="#536DFE","3BL"="#388E3C","4BR"="#795548"),
                      breaks=c("1TL","2TR","3BL","4BR"),
                      labels=c("Top Left","Top Right","Bottom Left","Bottom Right")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  labs(title = 'LID w/o Source Normalization',
       x = 'Fitness',
       y = '') +
  theme_linedraw() +
  coord_cartesian(xlim = c(0.7, 1.3))

sn.lid <- ggplot() +
  geom_line(data = dat.all[dat.all$method == 'lid' & dat.all$hours == 18,],
            aes(x = fitness, col = source), stat = 'density', lwd = 1.2) +
  scale_colour_manual(name="Source",
                      values=c("1TL"="#D32F2F","2TR"="#536DFE","3BL"="#388E3C","4BR"="#795548"),
                      breaks=c("1TL","2TR","3BL","4BR"),
                      labels=c("Top Left","Top Right","Bottom Left","Bottom Right")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  labs(title = 'LID',
       x = 'Fitness',
       y = '') +
  theme_linedraw() +
  coord_cartesian(xlim = c(0.7, 1.3))

annotate_figure(ggarrange(sn.nonorm, sn.nil, sn.lid,
                          nrow = 1,
                          common.legend = T,
                          legend = 'bottom'))#,
          # top = text_grob('Source Normalization', size = 15))

# ggsave(sprintf("%ssn.jpg",out_path),
#        width = 16, height = 5,
#        dpi = 300)

##### CHANGE IN PWR WITH METHOD
pwr.change <- NULL
i = 1
pwr.change$method[i] <- 'nonorm'
pwr.change$es[i] <- .975
pwr.change$sen[i] <- 6
pwr.change$method[i+1] <- 'nonorm'
pwr.change$es[i+1] <- 1.025
pwr.change$sen[i+1] <- 6
i = i + 2
pwr.change$method[i] <- 'nonorm'
pwr.change$es[i] <- .95
pwr.change$sen[i] <- 10
pwr.change$method[i+1] <- 'nonorm'
pwr.change$es[i+1] <- 1.05
pwr.change$sen[i+1] <- 7
i = i + 2
pwr.change$method[i] <- 'nil'
pwr.change$es[i] <- .975
pwr.change$sen[i] <- 23
pwr.change$method[i + 1] <- 'nil'
pwr.change$es[i + 1] <- 1.025
pwr.change$sen[i + 1] <- 10
i = i + 2
pwr.change$method[i] <- 'nil'
pwr.change$es[i] <- .95
pwr.change$sen[i] <- 55
pwr.change$method[i + 1] <- 'nil'
pwr.change$es[i + 1] <- 1.05
pwr.change$sen[i + 1] <- 26
i = i + 2
pwr.change$method[i] <- 'ncc'
pwr.change$es[i] <- .975
pwr.change$sen[i] <- 34
pwr.change$method[i+1] <- 'ncc'
pwr.change$es[i+1] <- 1.025
pwr.change$sen[i+1] <- 26
i = i + 2
pwr.change$method[i] <- 'ncc'
pwr.change$es[i] <- .95
pwr.change$sen[i] <- 70
pwr.change$method[i+1] <- 'ncc'
pwr.change$es[i+1] <- 1.05
pwr.change$sen[i+1] <- 66
i = i + 2
pwr.change$method[i] <- 'lid'
pwr.change$es[i] <- .975
pwr.change$sen[i] <- 35
pwr.change$method[i+1] <- 'lid'
pwr.change$es[i+1] <- 1.025
pwr.change$sen[i+1] <- 33
i = i + 2
pwr.change$method[i] <- 'lid'
pwr.change$es[i] <- .95
pwr.change$sen[i] <- 70
pwr.change$method[i+1] <- 'lid'
pwr.change$es[i+1] <- 1.05
pwr.change$sen[i+1] <- 70

pwr.change$effect[pwr.change$es > 1] <- 'Beneficial'
pwr.change$effect[pwr.change$es < 1] <- 'Deleterious'

pwr.change <- data.frame(pwr.change)
pwr.temp <- aggregate(pwr.change$sen[pwr.change$es %in% c(0.95,1.05)], list(pwr.change$method[pwr.change$es %in% c(0.95,1.05)]), mean)
colnames(pwr.temp) <- c('method', 'sen')

ggplot(pwr.temp) +
  geom_line(aes(y = sen, x = 4:1), lwd = 2, col = '#303F9F', alpha = 0.9) +
  geom_point(aes(y = sen, x = 4:1), size = 5, col = '#FFA000') +
  geom_point(aes(y = sen, x = 4:1), size = 1) +
  scale_y_continuous(breaks = seq(0,100,10), minor_breaks = seq(0,100,5)) +
  scale_x_continuous(breaks = 1:4,
                     labels = c('No\nNormalization','Ref. based\nLinear Interpolant',
                                '+ Source\nDeconstruction','+ Local\nArtefact Correction')) +
  labs(title = 'Change in Sensitivity for 5% Fitness Effect',
       subtitle = 'At 5% False Discovery Rate',
       x = '',
       y = 'Sensitivity') +
  theme_linedraw() +
  # theme(axis.title.x = element_blank()) +
  coord_cartesian(xlim = c(0.5,4.5),
                  ylim = c(0,80))
ggsave(sprintf("%spwr_chng.jpg",out_path),
       width = 5, height = 5,
       dpi = 300)

##### FPR FOLLOWS EXPECTATION

ggplot(stats.tmp[stats.tmp$hours == stats.tmp$cont_hrs,]) +
  geom_abline(col = 'red', linetype = 'dashed', lwd = 1) +
  geom_line(aes(x = p, y = fpr, col = hours), lwd = 1.2) +
  # geom_segment(aes(x=0,xend=0.2,y=0,yend=0), col = 'black', lwd = 0.8) +
  # geom_segment(aes(x=0.2,xend=0.2,y=0,yend=0.2), col = 'black', lwd = 0.8) +
  # geom_segment(aes(x=0.2,xend=0,y=0.2,yend=0.2), col = 'black', lwd = 0.8) +
  # geom_segment(aes(x=0,xend=0,y=0.2,yend=0), col = 'black', lwd = 0.8) +
  labs(title = "FPR Follows Random Expectation",
       # subtitle = "for p between 0 and 1",
       x = "P-Value",
       y = "False Positive Rate") +
  scale_x_continuous(breaks = seq(0,1,0.2),
                     minor_breaks = seq(0,1,0.05),
                     limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0,1,0.2),
                     minor_breaks = seq(0,1,0.05),
                     limits = c(0,1)) +
  scale_color_manual(name = "Avg. CS\n(Pix)",
                     breaks=c("13","14","16","17","18"),
                     values=c("13"="#D32F2F","14"="#536DFE","16"="#388E3C","17"="#795548","18"="#00BCD4",
                              "0"="transparent","8"="transparent","9"="transparent","10"="transparent","11"="transparent"),
                     labels=sprintf('%0.1f',rmse$avg[rmse$method == 'lid' & rmse$hour %in% c(13,14,16,17,18)])) +
  theme_linedraw() +
  # theme(axis.text.x = element_text(size=8),
  #       axis.title.x = element_text(size=10),
  #       axis.text.y = element_text(size=8),
  #       axis.title.y = element_text(size=10),
  #       # axis.title = element_blank(),
  #       legend.text = element_text(size=8),
  #       legend.title = element_text(size=10),
  #       legend.position = "bottom",
  #       plot.title = element_text(size=12),
  #       plot.subtitle = element_text(size=10)) +
  guides(color = guide_legend(override.aes = list(size=4)),
         shape = guide_legend(override.aes = list(size=4))) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

# ggsave(sprintf("%sfpr.jpg",out_path),
#        width = 6, height = 5,
#        dpi = 300)

##### SENSITIVITY WITH REPLICATE
pwr.rep <- NULL
pwr.rep$rep <- c(1:8)
pwr.rep$sen <- c((15+15)/2, (35+37)/2, (45+50)/2, (58+58)/2, (58+58)/2, (65+65)/2, (65+65)/2, (70+70)/2)
pwr.rep <- data.frame(pwr.rep)

ggplot(pwr.rep) +
  geom_line(aes(y = sen, x = rep), lwd = 2, col = '#303F9F', alpha = 0.9) +
  geom_point(aes(y = sen, x = rep), size = 5, col = '#FFA000') +
  scale_y_continuous(breaks = seq(0,100,10), minor_breaks = seq(0,100,5)) +
  scale_x_continuous(breaks = seq(0,10,1), minor_breaks = seq(0,10,1)) +
  labs(title = 'Change in Sensitivity for 5% Fitness Effects',
       subtitle = 'at p <= 0.05 with number of replicates per mutant',
       x = 'No. of Replicates',
       y = 'Sensitivity') +
  theme_linedraw() +
  coord_cartesian(ylim = c(0,80))
ggsave(sprintf("%spwr_rep.jpg",out_path),
       width = 5, height = 5,
       dpi = 300)

##### DETECTING SMALL FITNESS EFFECTS
ggplot(dat.cnt2) +
  geom_area(aes(x = cen, y = Neutral_p, fill = 'Neutral'), alpha = 1) +
  geom_area(aes(x = cen, y = Beneficial_p, fill = 'Beneficial'), alpha = 1) +
  geom_area(aes(x = cen, y = Deleterious_p, fill = 'Deleterious'), alpha = 1) +
  geom_vline(xintercept = seq(0,2,0.01), col = '#757575', lwd = 0.5, alpha =0.5) +
  geom_hline(yintercept = c(t * seq(0,1,0.05)), col = '#757575', lwd = 0.5, alpha =0.5) +
  labs(title = "Small Fitness Effects Can Be Detected",
       # subtitle = sprintf('with %d technical replicate (p = 0.05)',rep),
       x = 'Fitness Effect',
       y = 'Sensitivity') +
  scale_y_continuous(breaks = c(t * seq(0,1,0.1)),
                     minor_breaks = c(t * seq(0,1,0.05)),
                     labels = c('0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')) +
  scale_x_continuous(breaks = seq(0,2,0.05),
                     minor_breaks = seq(0,2,0.01),
                     labels = sprintf('%.1f%%',(seq(0,2,0.05)-1)*100)) +
  scale_color_manual(name = 'Effects',
                     breaks = c('Beneficial','Neutral','Deleterious'),
                     values = c('Deleterious'='#F44336',
                                'Neutral'='#303F9F',
                                'Beneficial'='#4CAF50')) +
  scale_fill_manual(name = 'Effects',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#F44336',
                               'Neutral'='#303F9F',
                               'Beneficial'='#4CAF50')) +
  theme_linedraw() +
  # theme(axis.text.x = element_text(size=8),
  #       axis.title.x = element_text(size=10),
  #       axis.text.y = element_text(size=8),
  #       axis.title.y = element_blank(),
  #       # axis.title = element_blank(),
  #       legend.text = element_text(size=8),
  #       legend.title = element_text(size=10),
  #       legend.position = "bottom",
  #       plot.title = element_text(size=12),
  #       plot.subtitle = element_text(size=10)) +
  guides(color = guide_legend(override.aes = list(size=2)),
         shape = guide_legend(override.aes = list(size=2))) +
  coord_cartesian(xlim = c(0.9,1.1),
                  ylim = c(0,t))
ggsave(sprintf("%ssml_fit.jpg",out_path),
       width = 6, height = 5,
       dpi = 300)

##### CHANGE IN SENSITIVITY WITH NO OF REFERENCES & REPLICATES
load(file = "/home/sbp29/R/Projects/adaptivefitness/figs/lid_paper/ref_rep_sen.RData")

ggplot(sen[sen$x %in% c(0.95,1.05),]) +
  geom_line(aes(x = rep, y = y, col = as.character(ref),
                group = ref), stat = 'summary', lwd = 1.2,
            linetype = 'dotted') +
  geom_point(aes(x = rep, y = y, col = as.character(ref),
                 group = ref), stat = 'summary', size = 4) +
  labs(title = 'Change in Sensitivity for 5% Fitness Effect',
       subtitle = 'at p < 0.05 with change in number of replicates and references',
       x = 'No. of Replicates',
       y = 'Sensitivity') +
  scale_y_continuous(breaks = seq(0,100,10), minor_breaks = seq(0,100,5)) +
  scale_x_continuous(breaks = seq(0,10,1), minor_breaks = seq(0,10,1)) +
  scale_color_manual(name = 'Reference\nProportion',
                     breaks = c('6.25','12.5','18.75','25'),
                     values = c('6.25'  ='#D32F2F',
                                '12.5'  ='#FFA000',
                                '18.75' ='#388E3C',
                                '25'    ='#303F9F'),
                     labels = c('6.25'  ='06.25 %',
                                '12.5'  ='12.50 %',
                                '18.75' ='18.75 %',
                                '25'    ='25.00 %')) +
  theme_linedraw() +
  guides(color = guide_legend(override.aes = list(size=4)),
         shape = guide_legend(override.aes = list(size=4))) +
  coord_cartesian(ylim = c(0,80))
ggsave(sprintf("%sref_prop.jpg",out_path),
       width = 6, height = 5,
       dpi = 300)

pie <- data.frame(
  group = c("Reference", "Query"),
  value = cbind(c(25, 75),
                c(18.75, 100 - 18.75),
                c(12.5, 100 - 12.5),
                c(6.25, 100 - 6.25))
)

ggplot(pie) +
  geom_bar(aes(x = "", y = value.4, fill = group), stat = 'identity') +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = c('Reference' = '#FFC107',
                                 'Query' = '#3F51B5'),
                      guide = F) +
  theme_minimal() +
  theme(axis.ticks = element_blank(),
        axis.title = element_blank(),
        rect = element_rect(fill = "transparent"))
ggsave(sprintf("%spie625.png",out_path),  bg = "transparent")

##### EFFECT SIZE
es.dat <- data.frame('hours'=c(8,9,10,11,13,14,16,17,18),
                     'P1'=double(length=9),
                     'P2'=double(length=9),
                     'P3'=double(length=9),
                     'P4'=double(length=9),
                     'P5'=double(length=9),
                     'P6'=double(length=9),
                     'P7'=double(length=9),
                     'P8'=double(length=9),
                     'P9'=double(length=9),
                     row.names = c('P1','P2','P3','P4','P5','P6','P7','P8','P9'))
es.dat <- rbind(hours=c(0,8,9,10,11,13,14,16,17,18),es.dat)

for (r in unique(dat.all$hours[dat.all$method == 'lid'])) {
  if (r > 0) {
    # ref.mean <- mean(dat.all$average[dat.all$hours == r & dat.all$orf_name == 'BF_control'], na.rm = T)
    ref.mean <- mean(dat.all$average[dat.all$hours == r], na.rm = T)
    for (q in unique(dat.all$hours[dat.all$method == 'lid'])) {
      if ( q > 0) {
        # que.mean <- mean(dat.all$average[dat.all$hours == q & !is.na(dat.all$orf_name) & dat.all$orf_name != 'BF_control'], na.rm = T)
        que.mean <- mean(dat.all$average[dat.all$hours == q], na.rm = T)
        es.dat[es.dat[,1] == r,es.dat[1,] == q] <- que.mean/ref.mean
      }
    }
  }
}

write.csv(es.dat, file = sprintf('%seffect_sizes.csv',out_path))

##### PLATE LAYOUT
lay.384$density <- 384
lay.1536$density <- 1536
lay.6144$density <- 6144
colnames(lay.384) <- c('pos','strain_id','orf_name','pos..4','plate','row','col','colony','density')
colnames(lay.1536) <- c('pos','strain_id','orf_name','pos..4','plate','row','col','colony','density')
colnames(lay.6144) <- c('pos','strain_id','orf_name','pos..4','plate','row','col','colony','density')

lay <- rbind(lay.384, lay.1536, lay.6144) 

layout <- ggplot(lay) +
  geom_point(aes(x = col, y = row, col = colony)) +
  scale_color_manual(name = 'Colony Type',
                     breaks = c('Stock','Reference', 'Query', 'Gap'),
                     values = c("Stock" = "#689F38",
                                "Reference" = "#FFC107",
                                "Query" = "#303F9F",
                                "Gap" = "#FF5252"),
                     guide = F) +
  facet_wrap(.~density*plate, nrow = 3, scales = "free") +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())
ggsave(sprintf("%slayout.jpg",out_path),
       layout,
       width = 12, height = 6,
       dpi = 1000)  

for (d in unique(lay$density)) {
  if (d == 384) {
    ggplot(lay[lay$density == d,]) +
      geom_point(aes(x = col, y = row, col = colony)) +
      scale_color_manual(name = 'Colony Type',
                         breaks = c('Stock','Reference', 'Query', 'Gap'),
                         values = c("Stock" = "#689F38",
                                    "Reference" = "#FFC107",
                                    "Query" = "#303F9F",
                                    "Gap" = "#FF5252"),
                         guide = F) +
      scale_x_continuous(limits = c(1,) )
      facet_wrap(.~plate, nrow = 1, scales = "free") +
      theme_linedraw() +
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            strip.background = element_blank(),
            strip.text.x = element_blank())
  }
}

##### POS2COOR

p2c.data <- dbGetQuery(conn, 'select * from 4C4_pos2coor a, 4C4_pos2orf_name b
                       where a.pos = b.pos')
p2c.data$colony[p2c.data$orf_name == 'BF_control'] <- 'Reference'
p2c.data$colony[p2c.data$orf_name != 'BF_control'] <- 'Query'
p2c.data$colony[is.na(p2c.data$orf_name)] <- 'Gap'

plate.map <- ggplot(p2c.data) +
  geom_point(aes(x = col, y = row, col = colony, size = as.character(density))) +
  scale_color_manual(name = 'Colony\nType',
                     values = c("Gap" = "#FFFFFF",
                                "Reference" = "#3F51B5",
                                "Query" = "#FFC107")) +
  scale_size_manual(guide = 'none',
                      values = c("384" = 4,
                                 "1536" = 1.5,
                                 "6144" = 0.4)) +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank()) +
  facet_wrap(.~density * plate,
             nrow = 3,
             scales = 'free')
ggsave(sprintf("%s%s_PLATE_MAP.jpg",out_path,expt_name),plate.map,
       height = 9, width = 16,
       dpi = 1000)

##### SOURCE NORMALIZATION
dat.avg <- dbGetQuery(conn, 'select b.*, orf_name, hours, average fitness
                      from
                      4C4_FS_NOSN_6144_FITNESS a, 4C4_pos2coor b
                      where hours = 11.04 and a.pos = b.pos
                      order by plate, col, row')
dat.avg$method <- '1. Raw Colony Size'

dat.nosn <- dbGetQuery(conn, 'select b.*, orf_name, hours, fitness
                        from
                        4C4_FS_NOSN_6144_FITNESS a, 4C4_pos2coor b
                        where hours = 11.04 and a.pos = b.pos
                        and fitness between 0 and 2
                        order by plate, col, row')
dat.nosn$method <- '2. LID w/o Source Normalization'

dat.lid <- dbGetQuery(conn, 'select b.*, orf_name, hours, fitness
                      from
                      4C4_FS_CC_6144_FITNESS a, 4C4_pos2coor b
                      where hours = 11.04 and a.pos = b.pos
                      and fitness between 0 and 2
                      order by plate, col, row')
dat.lid$method <- '3. LID'

dat.bean <- dbGetQuery(conn, 'select b.*, orf_name, hours, fitness
                       from
                       4C4_FS_BEAN_6144_FITNESS a, 4C4_pos2coor b
                       where hours = 11.04 and a.pos = b.pos
                       and fitness between 0 and 2
                       order by plate, col, row')
dat.bean$method <- '4. MCAT'

dat.sn <- rbind(dat.avg, dat.nosn, dat.lid, dat.bean)
dat.sn$source[dat.sn$row%%2==1 & dat.sn$col%%2==1] = '1TL'
dat.sn$source[dat.sn$row%%2==0 & dat.sn$col%%2==1] = '3BL'
dat.sn$source[dat.sn$row%%2==1 & dat.sn$col%%2==0] = '2TR'
dat.sn$source[dat.sn$row%%2==0 & dat.sn$col%%2==0] = '4BR'

ggplot(dat.sn[!is.na(dat.sn$fitness),]) +
  geom_line(aes(x = fitness, col = source),
            stat = 'density', trim = T,
            lwd = 1.2) +
  labs(title = 'Effect of Source Normalization') +
  scale_colour_manual(name="Source",
                      breaks=c("1TL","2TR","3BL","4BR"),
                      values=c("1TL"="#303F9F","2TR"="#536DFE","3BL"="#FFA000","4BR"="#FFECB3"),
                      labels=c("Top Left","Top Right","Bottom Left","Bottom Right")) +
  facet_wrap(.~method, nrow = 1,
             scales = 'free') +
  theme_linedraw() +
  theme(axis.title = element_blank())
ggsave(sprintf("%s%s_SN.jpg",out_path,expt_name),
       height = 3, width = 10,
       dpi = 1000)



