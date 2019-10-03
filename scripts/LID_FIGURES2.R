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
expt_name = '4C3_GA3_CC2'
expt = 'VAL'
out_path = 'figs/lid_paper/';
density = 6144;

tbl_pval = sprintf('%s_%d_PVALUE',expt_name,density)
hours = dbGetQuery(conn, sprintf('select distinct hours from %s order by hours asc', tbl_pval))
pvals = seq(0,1,0.01)

# tbl_fit_nonrm = sprintf('%s_NONORM_%d_FITNESS',substr(expt_name,1,7),density);
# tbl_fit_nil = sprintf('%s_NIL_%d_FITNESS',substr(expt_name,1,7),density);
# tbl_fit_ncc = sprintf('%s_%d_FITNESS',substr(expt_name,1,7),density);
# tbl_fit = sprintf('%s_%d_FITNESS',expt_name,density);
# tbl_fit_bean = sprintf('%s_BEAN_%d_FITNESS',substr(expt_name,1,7),density);
# tbl_p2o = '4C3_pos2orf_name3';
# tbl_bpos = '4C3_borderpos';

p2c_info = NULL
p2c_info[1] = '4C3_pos2coor6144'
p2c_info[2] = '6144plate'
p2c_info[3] = '6144col'
p2c_info[4] = '6144row'

p2c = dbGetQuery(conn, sprintf('select * from %s a order by a.%s, a.%s, a.%s',
                               p2c_info[1],
                               p2c_info[2],
                               p2c_info[3],
                               p2c_info[4]))

n_plates = dbGetQuery(conn, sprintf('select distinct %s from %s a order by %s asc',
                                    p2c_info[2],
                                    p2c_info[1],
                                    p2c_info[2]))


# dat.nonrm <- dbGetQuery(conn, sprintf('select b.*, a.orf_name, a.hours, a.bg, a.average, a.fitness
#                                     from %s a, %s b
#                                     where a.pos = b.pos
#                                     order by a.hours, b.%s, b.%s, b.%s',
#                                     tbl_fit_nonrm,
#                                     p2c_info[1],
#                                     p2c_info[2],p2c_info[3],p2c_info[4]))
# dat.nonrm$method <- 'nonorm'
# 
# dat.nil <- dbGetQuery(conn, sprintf('select b.*, a.orf_name, a.hours, a.bg, a.average, a.fitness
#                                     from %s a, %s b
#                                     where a.pos = b.pos
#                                     order by a.hours, b.%s, b.%s, b.%s',
#                                       tbl_fit_nil,
#                                       p2c_info[1],
#                                       p2c_info[2],p2c_info[3],p2c_info[4]))
# dat.nil$method <- 'nil'
# 
# dat.ncc <- dbGetQuery(conn, sprintf('select b.*, a.orf_name, a.hours, a.bg, a.average, a.fitness
#                                     from %s a, %s b
#                                     where a.pos = b.pos
#                                     order by a.hours, b.%s, b.%s, b.%s',
#                                     tbl_fit_ncc,
#                                     p2c_info[1],
#                                     p2c_info[2],p2c_info[3],p2c_info[4]))
# dat.ncc$method <- 'ncc'
# 
# dat.lid <- dbGetQuery(conn, sprintf('select b.*, a.orf_name, a.hours, a.bg, a.average, a.fitness
#                                     from %s a, %s b
#                                     where a.pos = b.pos
#                                     order by a.hours, b.%s, b.%s, b.%s',
#                                     tbl_fit,
#                                     p2c_info[1],
#                                     p2c_info[2],p2c_info[3],p2c_info[4]))
# dat.lid$method <- 'lid'
# 
# dat.bean <- dbGetQuery(conn, sprintf('select b.*, a.orf_name, a.hours, a.bg, a.average, a.fitness
#                                     from %s a, %s b
#                                     where a.pos = b.pos
#                                     order by a.hours, b.%s, b.%s, b.%s',
#                                    tbl_fit_bean,
#                                    p2c_info[1],
#                                    p2c_info[2],p2c_info[3],p2c_info[4]))
# dat.bean$method <- 'bean'
# 
# dat.all <- rbind(dat.nonrm, dat.nil, dat.ncc, dat.lid, dat.bean)

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
rmse$avg[rmse$method == 'bean'] <- rmse$avg[rmse$method == 'nonorm' & rmse$hour >0]
rmse$per <- rmse$rmse/rmse$avg * 100

?t.test
t.test(rmse$rmse[rmse$method == 'lid'], rmse$rmse[rmse$method == 'ncc'], alternative = 'less')

ggplot(rmse) +
  geom_point(aes(x = avg, y = rmse, fill = method),
             shape = 21, size = 3, alpha = 0.8, col = 'black') +
  scale_x_continuous(breaks = seq(0,500,50), minor_breaks = seq(0,500,25)) +
  scale_y_continuous(breaks = seq(0,100,10), minor_breaks = seq(0,100,5)) +
  scale_fill_manual(name = 'Method',
                    breaks = c('nonorm','nil','ncc','lid','bean'),
                    values = c('nonorm' = '#212121',
                               'nil' = '#448AFF',
                               'ncc' = '#7C4DFF',
                               'lid' = '#388E3C',
                               'bean' = '#FFA000'),
                    labels = c('RND', 'LID-SN', 'LID-CC', 'LID', 'MCAT')) +
  labs(title = 'Compairing RMSE',
       x = 'Mean Pixel Count',
       y = 'RMSE') +
  theme_linedraw() +
  coord_cartesian(xlim = c(0, 400),
                  ylim = c(0, 80))
ggsave(sprintf("%srmse.jpg",out_path),
       width = 6, height = 5,
       dpi = 300)

ggplot(rmse) +
  geom_point(aes(x = avg, y = rmse/avg * 100, fill = method),
             shape = 21, size = 3, alpha = 0.8, col = 'black') +
  scale_x_continuous(breaks = seq(0,500,50), minor_breaks = seq(0,500,25)) +
  scale_y_continuous(breaks = seq(0,100,10), minor_breaks = seq(0,100,5)) +
  scale_fill_manual(name = 'Method',
                    breaks = c('nonorm','nil','ncc','lid','bean'),
                    values = c('nonorm' = '#212121',
                               'nil' = '#448AFF',
                               'ncc' = '#7C4DFF',
                               'lid' = '#388E3C',
                               'bean' = '#FFA000'),
                    labels = c('RND', 'LID-SN', 'LID-CC', 'LID', 'MCAT')) +
  labs(title = 'Compairing RMSE',
       x = 'Mean Pixel Count',
       y = 'RMSE %') +
  theme_linedraw() +
  coord_cartesian(xlim = c(0, 400),
                  ylim = c(0, 60))
ggsave(sprintf("%srmse2.jpg",out_path),
       width = 6, height = 5,
       dpi = 300)

##### SOURCE NORMALIZATION
# dat.all$source[dat.all$`6144row`%%2==1 & dat.all$`6144col`%%2==1] = '1TL'
# dat.all$source[dat.all$`6144row`%%2==0 & dat.all$`6144col`%%2==1] = '3BL'
# dat.all$source[dat.all$`6144row`%%2==1 & dat.all$`6144col`%%2==0] = '2TR'
# dat.all$source[dat.all$`6144row`%%2==0 & dat.all$`6144col`%%2==0] = '4BR'
# 
# dat.all$colony[dat.all$orf_name == 'BF_control'] = 'Reference'
# dat.all$colony[dat.all$orf_name != 'BF_control'] = 'Query'
# dat.all$colony[is.na(dat.all$orf_name)] = 'Gap'

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

# save(dat.all, file = "rawdata/fitness.RData")

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
  geom_line(aes(x = sen, y = 4:1), lwd = 1.6, col = '#303F9F', alpha = 0.9) +
  geom_point(aes(x = sen, y = 4:1), size = 3, col = '#FFA000') +
  scale_x_continuous(breaks = seq(0,100,10), minor_breaks = seq(0,100,5)) +
  scale_y_continuous(labels = c('No. Norm','LID - SN','LID - CC','LID')) +
  labs(title = 'Change in Sensitivity for 5% Fitness Effect',
       subtitle = 'At 5% False Discovery Rate',
       x = 'Sensitivity',
       y = '') +
  theme_linedraw() +
  theme(axis.title.y = element_blank()) +
  coord_cartesian(xlim = c(0,80))

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

ggsave(sprintf("%sfpr.jpg",out_path),
       width = 6, height = 5,
       dpi = 300)

##### SENSITIVITY WITH REPLICATE
pwr.rep <- NULL
pwr.rep$rep <- c(1:8)
pwr.rep$sen <- c((15+15)/2, (35+37)/2, (45+50)/2, (58+58)/2, (58+58)/2, (65+65)/2, (65+65)/2, (70+70)/2)
pwr.rep <- data.frame(pwr.rep)

ggplot(pwr.rep) +
  geom_line(aes(x = sen, y = rep), lwd = 1.6, col = '#303F9F', alpha = 0.9) +
  geom_point(aes(x = sen, y = rep), size = 3, col = '#FFA000') +
  scale_x_continuous(breaks = seq(0,100,10), minor_breaks = seq(0,100,5)) +
  scale_y_continuous(breaks = seq(0,10,2), minor_breaks = seq(0,10,1)) +
  labs(title = 'Change in Sensitivity for 5% Fitness Effects',
       subtitle = 'at p <= 0.05 with number of replicates per mutant',
       x = 'Sensitivity',
       y = 'No. of Replicates') +
  theme_linedraw() +
  coord_cartesian(xlim = c(0,80))
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
ref.change <- NULL
i = 1
ref.change$ref[i] <- '4'
ref.change$es[i] <- 5
ref.change$rep[i] <- 8
ref.change$sen[i] <- mean(c(70,70))
i = i + 1
ref.change$ref[i] <- '4'
ref.change$es[i] <- 5
ref.change$rep[i] <- 6
ref.change$sen[i] <- mean(c(65,65))
i = i + 1
ref.change$ref[i] <- '4'
ref.change$es[i] <- 5
ref.change$rep[i] <- 4
ref.change$sen[i] <- mean(c(58,58))
i = i + 1
ref.change$ref[i] <- '4'
ref.change$es[i] <- 5
ref.change$rep[i] <- 2
ref.change$sen[i] <- mean(c(35,37))
i = i + 1
ref.change$ref[i] <- '3'
ref.change$es[i] <- 5
ref.change$rep[i] <- 8
ref.change$sen[i] <- mean(c(63,57))
i = i + 1
ref.change$ref[i] <- '3'
ref.change$es[i] <- 5
ref.change$rep[i] <- 6
ref.change$sen[i] <- mean(c(55,50))
i = i + 1
ref.change$ref[i] <- '3'
ref.change$es[i] <- 5
ref.change$rep[i] <- 4
ref.change$sen[i] <- mean(c(54,50))
i = i + 1
ref.change$ref[i] <- '3'
ref.change$es[i] <- 5
ref.change$rep[i] <- 2
ref.change$sen[i] <- mean(c(35,40))
i = i + 1
ref.change$ref[i] <- '2'
ref.change$es[i] <- 5
ref.change$rep[i] <- 8
ref.change$sen[i] <- mean(c(30,38))
i = i + 1
ref.change$ref[i] <- '2'
ref.change$es[i] <- 5
ref.change$rep[i] <- 6
ref.change$sen[i] <- mean(c(25,34))
i = i + 1
ref.change$ref[i] <- '2'
ref.change$es[i] <- 5
ref.change$rep[i] <- 4
ref.change$sen[i] <- mean(c(17,17))
i = i + 1
ref.change$ref[i] <- '2'
ref.change$es[i] <- 5
ref.change$rep[i] <- 2
ref.change$sen[i] <- mean(c(7,15))
i = i + 1
ref.change$ref[i] <- '1'
ref.change$es[i] <- 5
ref.change$rep[i] <- 8
ref.change$sen[i] <- mean(c(32,35))
i = i + 1
ref.change$ref[i] <- '1'
ref.change$es[i] <- 5
ref.change$rep[i] <- 6
ref.change$sen[i] <- mean(c(30,17))
i = i + 1
ref.change$ref[i] <- '1'
ref.change$es[i] <- 5
ref.change$rep[i] <- 4
ref.change$sen[i] <- mean(c(19,26))
i = i + 1
ref.change$ref[i] <- '1'
ref.change$es[i] <- 5
ref.change$rep[i] <- 2
ref.change$sen[i] <- mean(c(2,2))
i = i + 1

ref.change <- data.frame(ref.change)

ggplot(ref.change) +
  geom_line(aes(x = sen, y = rep, col = ref),
            lwd = 1.6, alpha = 0.9) +
  geom_point(aes(x = sen, y = rep),
             col = 'black', size = 3) +
  labs(title = 'Change in Sensitivity for 5% Fitness Effect',
       subtitle = 'at p < 0.05 with change in number of replicates and references',
       x = 'Sensitivity',
       y = 'No. of Replicates') +
  scale_x_continuous(breaks = seq(0,100,10), minor_breaks = seq(0,100,5)) +
  scale_color_manual(name = 'Reference\nProportion',
                     breaks = c('1','2','3','4'),
                     values = c('1'='#D32F2F',
                                '2'='#303F9F',
                                '3'='#388E3C',
                                '4'='#FFA000'),
                     labels = c('1'='06.25%','2'='12.50%','3'='18.75%','4'='25.00%')) +
  theme_linedraw() +
  guides(color = guide_legend(override.aes = list(size=4)),
         shape = guide_legend(override.aes = list(size=4))) +
  coord_cartesian(xlim = c(0,80))
ggsave(sprintf("%sref_prop.jpg",out_path),
       width = 6, height = 5,
       dpi = 300)

