##### 4C_MESSUP
##### MESSING UP THE FINAL PLATES
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 06/26/2019

##### INITIALIZE
library(RMariaDB)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggExtra)
library(grid)
library(tidyverse)
# library(egg)
library(ggpubr)
library(stringr)
source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")
out_path = 'figs/lid_paper/4C4/';

##### MAKING THE MESS
# alldat = dbGetQuery(conn, 'select c.*, b.orf_name, a.hours, a.average
#                   from 4C4_FS_6144_RAW a, 4C4_pos2orf_name b, 4C4_pos2coor c
#                   where a.pos = b.pos and b.pos = c.pos')
# alldat$rnd_hrs <- alldat$hours
# alldat$rnd_avg <- alldat$average
# 
# ggplot(alldat) +
#   geom_point(aes(x = rnd_hrs, y = rnd_avg))
# 
# for (hr in unique(alldat$hours)) {
#   for (orf in unique(alldat$orf_name[alldat$hours == hr & !is.na(alldat$orf_name) & alldat$orf_name != "BF_control"])) {
#     hr.tmp <- sample(unique(alldat$hours[alldat$hours != hr]),1)
#     alldat$rnd_avg[!is.na(alldat$orf_name) & alldat$orf_name == orf & alldat$hours == hr] <-
#       alldat$average[!is.na(alldat$orf_name) & alldat$orf_name == orf & alldat$hours == hr.tmp]
#     alldat$rnd_hrs[!is.na(alldat$orf_name) & alldat$orf_name == orf & alldat$hours == hr] <-
#       alldat$hours[!is.na(alldat$orf_name) & alldat$orf_name == orf & alldat$hours == hr.tmp]
#   }
# }
# 
# ggplot(alldat) +
#   geom_point(aes(x = rnd_hrs, y = rnd_avg))
# 
# dbWriteTable(conn, "4C4_FS_RND_6144_DATA", alldat, overwrite = T)
# 
# dbExecute(conn, 'drop table 4C4_FS_RND_BEAN_6144_JPEG')
# dbExecute(conn, 'create table 4C4_FS_RND_BEAN_6144_JPEG
#     select pos, hours, rnd_hrs, rnd_avg average
#     from 4C4_FS_RND_6144_DATA')
# 
# dbExecute(conn, 'drop table 4C4_FS_RND_6144_JPEG')
# dbExecute(conn, 'create table 4C4_FS_RND_6144_JPEG
#     select pos, hours, rnd_hrs, rnd_avg average
#     from 4C4_FS_RND_6144_DATA')
# 
# dbExecute(conn, 'update 4C4_FS_RND_6144_JPEG
#     set average = NULL
#     where pos in
#     (select pos from 4C4_borderpos)')

##### ANALYZING THE MESS
dat.dir <- "/home/sbp29/R/Projects/adaptivefitness/rawdata/4C4_FS_RND_BEAN/"
expt_name <- '4C4_FS_RND_BEAN'
pvals = seq(0,1,0.005)

stats.files <- list.files(path = dat.dir,
                          pattern = "P.csv", recursive = TRUE)
fit.files <- list.files(path = dat.dir,
                        pattern = "S.csv", recursive = TRUE)

hours <- NULL
reps <- NULL
refs <- NULL
atmpt <- NULL
for (s in strsplit(stats.files,'_')) {
  # refs <- c(refs, as.numeric(s[4]))
  reps <- c(reps, as.numeric(s[5]))
  # reps <- 8
  hours <- c(hours, as.numeric(s[6]))
  atmpt <- c(atmpt, as.numeric(s[7]))
}
# refs <- unique(refs)
reps <- sort(unique(reps))
hours <- sort(unique(hours))
atmpt <- sort(unique(atmpt))

hhours <- c(0,1.02,1.38,2.9,4.02,4.89,6.14,6.88,7.85,8.96,9.97,11.04)

stats.all <- NULL
fit.all <- NULL
fpr.all <- NULL
##### PUTTING IT TOGETHER
for (iii in 1:length(atmpt)) {
  for (ii in 1:length(reps)) {
    rep <- reps[ii]
    for (i in 1:length(hours)) {
      hr <- hours[i]
      
      dat.stats <- read.csv(paste0(dat.dir,
                                   sprintf('%s_%d_%d_%d_STATS_P.csv',expt_name,rep,hr,iii)),
                            na.strings = "NaN")

      dat.stats$cont_hrs <- hhours[i]
      dat.stats$rep <- rep
      dat.fit <- read.csv(paste0(dat.dir,
                                 sprintf('%s_%d_%d_%d_FITNESS.csv',expt_name,rep,hr,iii)),
                          na.strings = "NaN")

      dat.fit$cont_hrs <- hhours[i]
      dat.fit$rep <- rep
      dat.fit$se <- dat.fit$average - dat.fit$bg
      
      dat.stats$attempt <- iii
      dat.stats$replicates <- rep
      dat.stats$pthresh <- quantile(sort(dat.stats$p[dat.stats$hours == hhours[i]]),.05)
      for (h in unique(dat.stats$hours)) {
        cont.mean <- mean(dat.fit$fitness[dat.fit$hours == h & dat.fit$orf_name == 'BF_control' & !is.na(dat.fit$fitness)])
        # dat.stats$es[dat.stats$hours == h] <- round(dat.stats$cs_mean[dat.stats$hours == h]/cont.mean,4)
        
        ref.pix.m <- mean(dat.fit$average[dat.fit$hours == h & dat.fit$orf_name == 'BF_control' & !is.na(dat.fit$fitness)])
        que.pix.m <- mean(dat.fit$average[dat.fit$hours == h & dat.fit$orf_name != 'BF_control' & !is.na(dat.fit$fitness)])
        dat.stats$es[dat.stats$hours == h] <- que.pix.m/ref.pix.m
        
        dat.stats$cen[dat.stats$hours == h] <- mean(dat.stats$es[dat.stats$hours == h])
        
        dat.stats$effect[dat.stats$hours == h & dat.stats$p <= dat.stats$pthresh & dat.stats$cs_mean > cont.mean & !is.na(dat.stats$cs_mean)] <- 'Beneficial'
        dat.stats$effect[dat.stats$hours == h & dat.stats$p <= dat.stats$pthresh & dat.stats$cs_mean < cont.mean & !is.na(dat.stats$cs_mean)] <- 'Deleterious'
        dat.stats$effect[dat.stats$hours == h & is.na(dat.stats$effect) & !is.na(dat.stats$cs_mean)] <- 'Neutral'
        
        dat.stats$effect_p[dat.stats$hours == h & dat.stats$p <= 0.05 & dat.stats$cs_mean > cont.mean & !is.na(dat.stats$cs_mean)] <- 'Beneficial'
        dat.stats$effect_p[dat.stats$hours == h & dat.stats$p <= 0.05 & dat.stats$cs_mean < cont.mean & !is.na(dat.stats$cs_mean)] <- 'Deleterious'
        dat.stats$effect_p[dat.stats$hours == h & is.na(dat.stats$effect_p) & !is.na(dat.stats$cs_mean)] <- 'Neutral'
      }
      
      stats.all <- rbind(stats.all,dat.stats)
      fit.all <- rbind(fit.all,dat.fit)
    }
  }
}

fit.all$strain[fit.all$orf_name == "BF_control"] = "Reference"
fit.all$strain[fit.all$orf_name != "BF_control"] = "Query"

# for (i in 1:length(fit.all$pos)) {
#   fit.all$average[i]
# }

for (hr in sort(unique(stats.all$hours))) {
  for (orf in unique(stats.all$orf_name)) {
    stats.all$pix_mean[stats.all$hours == hr & stats.all$orf_name == orf] <-
      mean(fit.all$average[fit.all$hours == hr & fit.all$orf_name == orf], na.rm = T)
  }
}

ref.all <- NULL
i = 0
for (hr in sort(unique(fit.all$hours))) {
  for (pos in 1153:1536) {
    i <- i + 1
    ref.all$pix_mean[i] <-
      mean(fit.all$average[fit.all$pos %in% c(c(110000,120000,130000,140000,210000,220000,230000,240000,
                                                310000,320000,330000,340000,410000,420000,430000,440000) + pos) & fit.all$hours == hr], na.rm = T)
    ref.all$cs_mean[i] <-
      mean(fit.all$fitness[fit.all$pos %in% c(c(110000,120000,130000,140000,210000,220000,230000,240000,
                                                310000,320000,330000,340000,410000,420000,430000,440000) + pos) & fit.all$hours == hr], na.rm = T)
    ref.all$hours[i] <- hr
  }
}
ref.all <- data.frame(ref.all)

rnd_hrs <- dbGetQuery(conn, 'select orf_name, hours, rnd_hrs
                      from 4C4_FS_RND_6144_DATA
                      where orf_name is not NULL
                      group by hours, orf_name, rnd_hrs
                      order by hours, orf_name')

for (hr in sort(unique(stats.all$hours))) {
  for (orf in unique(stats.all$orf_name[stats.all$hours == hr])) {
    stats.all$rnd_hrs[stats.all$hours == hr & stats.all$orf_name == orf] <-
      rnd_hrs$rnd_hrs[rnd_hrs$hours == hr & rnd_hrs$orf_name == orf]
  }
}

stats.all$from[stats.all$hours > stats.all$rnd_hrs] = 'less'
stats.all$from[stats.all$hours < stats.all$rnd_hrs] = 'more'

fit.all$method <- 'BEAN'
ref.all$method <- 'BEAN'
stats.all$method <- 'BEAN'

fit.data <- NULL
ref.data <- NULL
sta.data <- NULL

fit.data <- rbind(fit.data, fit.all)
ref.data <- rbind(ref.data, ref.all)
sta.data <- rbind(sta.data, stats.all)

save(fit.data, file = "output/4C4_rnd_fit_data.RData")
save(ref.data, file = "output/4C4_rnd_ref_data.RData")
save(sta.data, file = "output/4C4_rnd_sta_data.RData")

load("output/4C4_rnd_fit_data.RData")
load("output/4C4_rnd_ref_data.RData")
load("output/4C4_rnd_sta_data.RData")
rnd_hrs <- dbGetQuery(conn, 'select hours, orf_name, round(avg(rnd_hrs),2) rnd_hrs
                              from 4C4_FS_RND_6144_DATA
                              where orf_name != "BF_control"
                              group by hours, orf_name')
# rnd_hrs <- read.csv("figs/lid_paper/rnd_hrs.csv")

rnd_hrs$from[rnd_hrs$rnd_hrs > rnd_hrs$hours] <- 'more'
rnd_hrs$from[rnd_hrs$rnd_hrs < rnd_hrs$hours] <- 'less'

# hours <- c(8,9,10,11,13,14,16,17,18)
# hours <- c(0,1.02,1.38,2.9,4.02,4.89,6.14,6.88,7.85,8.96,9.97,11.04)
hours <- c(1.38,2.9,4.02,4.89,6.14,6.88,7.85,8.96,9.97,11.04)

# sta.data <- sta.data[sort(hours),]

for (h in sort(hours)) {
  for (o in unique(sta.data$orf_name[sta.data$hours == h])) {
    sta.data$from[sta.data$hours == h & sta.data$orf_name == o] <- rnd_hrs$from[rnd_hrs$hours == h & rnd_hrs$orf_name == o]
  }
}

ggplot(ref.data) +
  geom_jitter(data = sta.data,
              aes(x = hours, y = pix_mean, col = effect_p),
              size = 0.3, alpha = 0.7) +
  geom_violin(aes(x = hours, y = pix_mean, group = hours),
              draw_quantiles = c(0.25, 0.5, 0.75),
              fill = 'transparent') +
  labs(x = 'Hours',
       y = 'Colony Size (pix)') +
  # geom_jitter(aes(x = hours, y = pix_mean),
  #             size = 0.2, alpha = 0.7) +
  facet_wrap(~method, nrow = 2) +
  scale_color_manual(name = 'Effects',
                    breaks = c('Beneficial', 'Neutral', 'Deleterious'),
                    values = c('Beneficial' = '#388E3C',
                               'Deleterious' = '#FF5252',
                               'Neutral' = '#303F9F')) +
  scale_x_continuous(breaks = 1:20) +
  theme_linedraw() +
  guides(color = guide_legend(override.aes = list(size=3)),
         shape = guide_legend(override.aes = list(size=3)))
# ggsave(sprintf("%srnd_jitter.jpg",out_path),
#        width = 6, height = 5,
#        dpi = 300)

den.dat <- data.frame(rbind(ref.data[ref.data$method == 'BEAN',c(1,3)],
                            sta.data[sta.data$method == 'BEAN',c(15,2)]))
i = 1
for (hr in unique(den.dat$hours)) {
  den.dat$plate[den.dat$hours == hr] <- sprintf('Plate_%d',i)
  i = i + 1
}
ggplot(den.dat) +
  geom_line(aes(x = pix_mean),
              stat = 'density') +
  labs(title = 'Virtual Plate Colony Size Distribution',
       y = 'Density',
       x = 'Colony Size (pix)') +
  facet_wrap(.~plate) +
  theme_linedraw()
# ggsave(sprintf("%srnd_den.jpg",out_path),
#        width = 5, height = 5,
#        dpi = 300)

sta.data$res <- sta.data$effect_p
sta.data$res[sta.data$effect_p == 'Beneficial' & sta.data$from == 'more'] <- '5TrueBeneficial'
sta.data$res[sta.data$effect_p == 'Beneficial' & sta.data$from == 'less'] <- '3FalseBeneficial'
sta.data$res[sta.data$effect_p == 'Deleterious' & sta.data$from == 'more'] <- '4FalseDeleterious'
sta.data$res[sta.data$effect_p == 'Deleterious' & sta.data$from == 'less'] <- '2TrueDeleterious'
sta.data$res[sta.data$effect_p == 'Neutral' & sta.data$from == 'more'] <- '6BenNeutral'
sta.data$res[sta.data$effect_p == 'Neutral' & sta.data$from == 'less'] <- '1DelNeutral'
# sta.data$res[sta.data$effect_p == 'Neutral' & sta.data$from == 'less'] <- 'DelNeutral'
# sta.data$res[sta.data$effect_p == 'Deleterious' & sta.data$from == 'less'] <- 'TrueDeleterious'
# sta.data$res[sta.data$effect_p == 'Beneficial' & sta.data$from == 'less'] <- 'FalseBeneficial'
# sta.data$res[sta.data$effect_p == 'Deleterious' & sta.data$from == 'more'] <- 'FalseDeleterious'
# sta.data$res[sta.data$effect_p == 'Beneficial' & sta.data$from == 'more'] <- 'TrueBeneficial'
# sta.data$res[sta.data$effect_p == 'Neutral' & sta.data$from == 'more'] <- 'BenNeutral'

pie.dat <- data.frame(rbind(cbind(hours = sta.data$hours[sta.data$method == 'BEAN'],
                                  value = sta.data$from[sta.data$method == 'BEAN'],
                                  method = 'TRUTH'),
                            cbind(hours = sta.data$hours[sta.data$method == 'BEAN'],
                                  value = sta.data$res[sta.data$method == 'BEAN'],
                                  method = 'BEAN'),
                            cbind(hours = sta.data$hours[sta.data$method == 'BEAN'],
                                  value = sta.data$res[sta.data$method == 'LID'],
                                  method = 'LID')))

pie.dat$value[pie.dat$value == 'less'] <- '2TrueDeleterious'
pie.dat$value[pie.dat$value == 'more'] <- '5TrueBeneficial'

# pie.dat <- arrange(transform(pie.dat,
#                              value=factor(value,levels=c('DelNeutral','TrueDeleterious','FalseBeneficial',
#                                                          'FalseDeleterious','TrueBeneficial','BenNeutral'))),
#                    value)

pie.dat <- arrange(transform(pie.dat,
                           method=factor(method,levels=c('LID','TRUTH','BEAN'))),
                   method)
# pie.dat$value <- as.character(pie.dat$value)

pie.per <- plyr::count(pie.dat, vars = c('method','value'))
pie.per$per <- pie.per$freq/(sum(pie.per$freq)/3) * 100

ggplot(pie.dat) +
  geom_bar(aes(x = "", y = "", fill = value), stat = 'identity') +
  coord_polar("y", start = 0) +
  scale_fill_manual(name = 'Effects',
                    breaks = c('6BenNeutral',
                               '5TrueBeneficial',
                               '4FalseDeleterious',
                               '3FalseBeneficial',
                               '2TrueDeleterious',
                               '1DelNeutral'),
                    values = c('1DelNeutral'='#FF9800',
                               '2TrueDeleterious'='#D32F2F',
                               '3FalseBeneficial'='#FFEB3B',
                               '4FalseDeleterious'='#DCEDC8',
                               '5TrueBeneficial'='#388E3C',
                               '6BenNeutral'='#8BC34A'),
                    labels = c('1DelNeutral'='False Neutral Deleterious',
                               '2TrueDeleterious'='True Deleterious',
                               '3FalseBeneficial'='False Beneficial Deleterious',
                               '4FalseDeleterious'='False Deleterious Beneficial',
                               '5TrueBeneficial'='True Beneficial',
                               '6BenNeutral'='False Neutral Beneficial')) +
                      
  facet_wrap(.~method) +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'bottom')
ggsave(sprintf("%srnd_pie.jpg",out_path),
       width = 6, height = 3,
       dpi = 1000)  



tp <- sta.data[sta.data$effect_p == 'Beneficial' & sta.data$from == 'more' |
                             sta.data$effect_p == 'Deleterious' & sta.data$from == 'less',]
fp <- sta.data[sta.data$effect_p == 'Deleterious' & sta.data$from == 'more' |
                 sta.data$effect_p == 'Beneficial' & sta.data$from == 'less',]

lid.pow <- dim(tp[tp$method == 'LID',])[1]/dim(sta.data[sta.data$method == 'LID',])[1] * 100
mcat.pow <- dim(tp[tp$method == 'BEAN',])[1]/dim(sta.data[sta.data$method == 'BEAN',])[1] * 100

lid.fp <- dim(fp[fp$method == 'LID',])[1]/dim(sta.data[sta.data$method == 'LID',])[1] * 100
mcat.fp <- dim(fp[fp$method == 'BEAN',])[1]/dim(sta.data[sta.data$method == 'BEAN',])[1] * 100

lid.tp <- dim(tp[tp$method == 'LID',])[1]/dim(sta.data[sta.data$method == 'LID' &
                                                         (sta.data$effect_p == 'Beneficial' | sta.data$effect_p == 'Deleterious'),])[1] * 100
mcat.tp <- dim(tp[tp$method == 'BEAN',])[1]/dim(sta.data[sta.data$method == 'BEAN' &
                                                         (sta.data$effect_p == 'Beneficial' | sta.data$effect_p == 'Deleterious'),])[1] * 100

mean(ref.data$pix_mean[ref.data$hours == 18], na.rm = T)

mess.tmp <- NULL
for (h in sort(hours)) {
  if (h > 0) {
    
    tmp.ref.pix <- mean(ref.data$pix_mean[ref.data$hours == h], na.rm = T)
    
    if (dim(count(sta.data[sta.data$hours == h,], vars = "from"))[1] == 2) {
      tmp.ratio <- count(sta.data[sta.data$hours == h,], vars = "from")[[2,2]]/1818 * 100 #count(sta.data[sta.data$hours == h,], from)[[1,2]]
    } else {
      tmp.ratio <- 0
    }
    
    
    tmp.lid.pow <- dim(tp[tp$method == 'LID' & tp$hours == h,])[1]/
      dim(sta.data[sta.data$method == 'LID' & sta.data$hours == h,])[1] * 100
    tmp.mcat.pow <- dim(tp[tp$method == 'BEAN' & tp$hours == h,])[1]/
      dim(sta.data[sta.data$method == 'BEAN' & sta.data$hours == h,])[1] * 100
    
    tmp.lid.fp <- dim(fp[fp$method == 'LID' & fp$hours == h,])[1]/
      dim(sta.data[sta.data$method == 'LID' & sta.data$hours == h,])[1] * 100
    tmp.mcat.fp <- dim(fp[fp$method == 'BEAN' & fp$hours == h,])[1]/
      dim(sta.data[sta.data$method == 'BEAN' & sta.data$hours == h,])[1] * 100
    
    tmp.lid.tp <- dim(tp[tp$method == 'LID' & tp$hours == h,])[1]/
      dim(sta.data[sta.data$method == 'LID' & sta.data$hours == h &
                     (sta.data$effect_p == 'Beneficial' | sta.data$effect_p == 'Deleterious'),])[1] * 100
    tmp.mcat.tp <- dim(tp[tp$method == 'BEAN' & tp$hours == h,])[1]/
      dim(sta.data[sta.data$method == 'BEAN' & sta.data$hours == h &
                     (sta.data$effect_p == 'Beneficial' | sta.data$effect_p == 'Deleterious'),])[1] * 100
    
    mess.tmp <- rbind(mess.tmp,
                      cbind(h, tmp.ref.pix, tmp.ratio,
                            tmp.lid.tp, tmp.lid.fp, tmp.lid.pow,
                            tmp.mcat.tp, tmp.mcat.fp, tmp.mcat.pow))
  }
}
mess.tmp <- data.frame(mess.tmp)

mess.tmp.lid <- mess.tmp[,c(1:6)]
mess.tmp.lid$method <- 'LID'
mess.tmp.mcat <- mess.tmp[,c(1:3,7:9)]
mess.tmp.mcat$method <- 'BEAN'

colnames(mess.tmp.lid) <- c('hours','REF','BEN','TP','FP','POW','method')
colnames(mess.tmp.mcat) <- c('hours','REF','BEN','TP','FP','POW','method')

mess.res <- rbind(mess.tmp.lid,mess.tmp.mcat)
ggplot(mess.res) +
  geom_line(aes(x = POW, y = TP, col = method),lwd = 1.4) +
  labs(x = 'Sensitivity',
       y = 'Specificity')


mess.res <- melt(mess.res, id.vars = c('hours','method','REF','BEN'))
ggplot(mess.res[mess.res$variable != 'TP',]) +
  geom_line(aes(x = BEN, y = value, col = variable, linetype = method),
            lwd = 1.4) +
  scale_color_manual(name = 'Measure',
                     breaks = c('TP','FP','POW'),
                     values = c('TP'='#9C27B0','FP'='#009688','POW'='#1976D2'),
                     labels = c('True\nPos','False Pos.','Sensitivity')) +
  scale_linetype_discrete(name = 'Method') +
  scale_x_continuous(breaks = seq(0,100,10)) +
  scale_y_continuous(breaks = seq(0,100,10)) +
  labs(title = 'LID & MCAT Performance',
       subtitle = 'On Randomly Distributed Fitness Effects',
       x = 'True Beneficial Proportion (%)',
       y = 'Percentage') +
  theme_linedraw() +
  coord_cartesian(xlim = c(0,90)) +
  guides(color = guide_legend(override.aes = list(size=4)),
         linetype = guide_legend(override.aes = list(size=1)))
ggsave(sprintf("%sperf_rnd.jpg",out_path),
       width = 6, height = 5,
       dpi = 300)

ggplot(sta.data[sta.data$hours != 9,]) +
  geom_bar(aes(x = from, fill = effect_p)) +
  scale_fill_manual(name = 'Predicted\nEffects',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#F44336',
                               'Neutral'='#303F9F',
                               'Beneficial'='#4CAF50')) +
  scale_x_discrete(labels = c('True\nDel.','True\nBen')) +
  labs(title = 'LID & MCAT Performance',
       subtitle = 'On Randomly Distributed Fitness Effects',
       y = 'Count') +
  theme_linedraw() +
  theme(axis.title.x = element_blank()) +
  facet_wrap(.~hours*method, ncol = 4)
ggsave(sprintf("%sperf_rnd2.jpg",out_path),
       width = 10, height = 10,
       dpi = 300)

# hr = 18
# cs <- ggplot() +
#   geom_point(data = fit.all[fit.all$hours == hr,],
#              aes(x = x6144col_1, y = x6144row_1, col = average), shape = 15) +
#   scale_x_continuous(breaks = seq(1,96,1)) +
#   scale_y_continuous(breaks = seq(1,64,1),trans = 'reverse') +
#   scale_color_distiller(name = "PIX",
#                         limits = c(quantile(fit.all$average[fit.all$hours == hr],0.05,na.rm = T)[[1]],
#                                    quantile(fit.all$average[fit.all$hours == hr],0.95,na.rm = T)[[1]]),
#                         palette = "Spectral") +
#   labs(title = "",
#        subtitle = "Colony Sizes",
#        x = "Columns",
#        y = "Rows") +
#   theme_linedraw() +
#   theme(axis.ticks = element_blank(),
#         axis.text = element_blank())
# 
# f <- ggplot() +
#   geom_point(data = fit.all[fit.all$hours == hr,],
#              aes(x = x6144col_1, y = x6144row_1, col = fitness), shape = 15) +
#   scale_x_continuous(breaks = seq(1,96,1)) +
#   scale_y_continuous(breaks = seq(1,64,1),trans = 'reverse') +
#   scale_color_distiller(name = "FIT",
#                         limits = c(quantile(fit.all$fitness[fit.all$hours == hr],0.05,na.rm = T)[[1]],
#                                    quantile(fit.all$fitness[fit.all$hours == hr],0.95,na.rm = T)[[1]]),
#                         palette = "PiYG") +
#   labs(title = "",
#        subtitle = "Fitness",
#        x = "Columns",
#        y = "Rows") +
#   theme_linedraw() +
#   theme(axis.ticks = element_blank(),
#         axis.text = element_blank())
# 
# cs <- ggplot() +
#   geom_line(data = fit.all[fit.all$hours == hr,],
#                aes(x = average, col = strain), lwd = 1.2, stat = "density")
# 
# f <- ggplot() +
#   geom_line(data = fit.all[fit.all$hours == hr,],
#             aes(x = fitness, col = strain), lwd = 1.2, stat = "density")
# 
# csVf <- ggplot() +
#   geom_point(data = ref.all[ref.all$hours == hr,],
#              aes(x = pix_mean, y = cs_mean, col = "Reference"), alpha = 1) +
#   geom_smooth(data = ref.all[ref.all$hours == hr,],
#              aes(x = pix_mean, y = cs_mean),
#              col = "black",
#              method = lm, se = F) +
#   geom_point(data = stats.all[stats.all$hours == hr,],
#              aes(x = pix_mean, y = cs_mean, col = effect_p), alpha = 0.5) +
#   # geom_smooth(data = stats.all[stats.all$hours == hr,],
#   #            aes(x = pix_mean, y = cs_mean),
#   #            col = "black",
#   #            method = lm, se = F) +
#   scale_color_manual(name = "Strain Kind",
#                      breaks = c("Reference", "Beneficial", "Neutral", "Deleterious"),
#                      values = c("Reference" = "#607D8B",
#                                 "Beneficial" = "#388E3C",
#                                 "Neutral" = "#303F9F",
#                                 "Deleterious" = "#D32F2F"),
#                      labels = c("Ref","B","N","D")) +
#   scale_x_continuous(breaks = seq(0,1000,100), minor_breaks = seq(0,1000,25)) +
#   scale_y_continuous(breaks = seq(-10,10,0.2), minor_breaks = seq(-10,10,0.05)) +
#   labs(title = sprintf("LID at %d hours",hr),
#        subtitle = "CS and Fitness Correlation",
#        x = "Colony Size (pix)",
#        y = "Fitness") +
#   theme_linedraw() +
#   theme(legend.position = "right") +
#   coord_cartesian(xlim = c(100,500),
#                   ylim = c(0.2,2))
# 
# ggsave(sprintf("%s%s_CSVFIT_%d.png",
#                out_path,expt_name,hr),
#        csVf,
#        width = 7,height = 7)
# 
# ggarrange(csVf,cs,f,
#           nrow = 1, ncol = 3, widths = c(1.8,1.8,1.4))
# 
# ggsave(sprintf("%s%s_CSFIT_%d.png",
#                out_path,expt_name,hr),
#        width = 18,height = 4.5)

csVf.all <- ggplot() +
  geom_point(data = ref.all,
             aes(x = pix_mean, y = cs_mean, col = "Reference"), alpha = 1) +
  geom_smooth(data = ref.all,
              aes(x = pix_mean, y = cs_mean),
              col = "black",
              method = lm, se = F) +
  geom_point(data = stats.all,
             aes(x = pix_mean, y = cs_mean, col = effect_p), alpha = 0.5) +
  facet_wrap(.~hours, nrow = 3, ncol = 3) +
  scale_color_manual(name = "",
                     breaks = c("Reference", "Beneficial", "Neutral", "Deleterious"),
                     values = c("Reference" = "#607D8B",
                                "Beneficial" = "#388E3C",
                                "Neutral" = "#303F9F",
                                "Deleterious" = "#D32F2F")) +
  scale_x_continuous(breaks = seq(0,1000,100), minor_breaks = seq(0,1000,25)) +
  scale_y_continuous(breaks = seq(-10,10,0.2), minor_breaks = seq(-10,10,0.05)) +
  labs(title = "LID",
  # labs(title = "MCAT",
       subtitle = "CS and Fitness Correlation by Hour",
       x = "Colony Size (pix)",
       y = "Fitness") +
  theme_linedraw() +
  theme(legend.position = "right") +
  coord_cartesian(xlim = c(0,500),
                  ylim = c(0,2))

ggsave(sprintf("%s%s_CSVFIT_ALL.png",
               out_path,expt_name),
       csVf.all,
       width = 10, height = 10)

# csden.all <- ggplot() +
#   geom_line(data = ref.all,
#              aes(x = pix_mean, col = "Reference"), stat = "density", lwd = 1.2) +
#   geom_line(data = stats.all,
#              aes(x = pix_mean, col = effect_p), stat = "density", lwd = 1.2) +
#   facet_wrap(.~hours, nrow = 3, ncol = 3) +
#   scale_color_manual(name = "Strain Kind",
#                      breaks = c("Reference", "Beneficial", "Neutral", "Deleterious"),
#                      values = c("Reference" = "#607D8B",
#                                 "Beneficial" = "#388E3C",
#                                 "Neutral" = "#303F9F",
#                                 "Deleterious" = "#D32F2F"),
#                      labels = c("Ref","B","N","D")) +
#   scale_x_continuous(breaks = seq(0,1000,100), minor_breaks = seq(0,1000,25)) +
#   # scale_y_continuous(breaks = seq(-10,10,0.2), minor_breaks = seq(-10,10,0.05)) +
#   labs(title = "LID",
#        subtitle = "CS Density",
#        x = "Colony Size (pix)",
#        y = "Density") +
#   theme_linedraw() +
#   theme(legend.position = "right") +
#   coord_cartesian(xlim = c(100,500))
# 
# ggsave(sprintf("%s%s_CSDEN_ALL.png",
#                out_path,expt_name),
#        csden.all,
#        width = 15,height = 15)
# 
# fden.all <- ggplot() +
#   geom_line(data = ref.all,
#             aes(x = cs_mean, col = "Reference"), stat = "density", lwd = 1.2) +
#   geom_line(data = stats.all,
#             aes(x = cs_mean, col = effect_p), stat = "density", lwd = 1.2) +
#   facet_wrap(.~hours, nrow = 3, ncol = 3) +
#   scale_color_manual(name = "Strain Kind",
#                      breaks = c("Reference", "Beneficial", "Neutral", "Deleterious"),
#                      values = c("Reference" = "#607D8B",
#                                 "Beneficial" = "#388E3C",
#                                 "Neutral" = "#303F9F",
#                                 "Deleterious" = "#D32F2F"),
#                      labels = c("Ref","B","N","D")) +
#   scale_x_continuous(breaks = seq(-10,10,0.2), minor_breaks = seq(-10,10,0.05)) +
#   # scale_y_continuous(breaks = seq(-10,10,0.2), minor_breaks = seq(-10,10,0.05)) +
#   labs(title = "LID",
#        subtitle = "Fitness Density",
#        x = "Fitness",
#        y = "Density") +
#   theme_linedraw() +
#   theme(legend.position = "right") +
#   coord_cartesian(xlim = c(0.2,2))
# 
# ggsave(sprintf("%s%s_FDEN_ALL.png",
#                out_path,expt_name),
#        fden.all,
#        width = 15,height = 15)

edis.all <- ggplot() +
  geom_bar(data = stats.all,
           aes(x=effect_p, y=(..count..)/913*100, fill = effect_p)) +
  facet_wrap(.~hours, nrow = 3, ncol = 3) +
  scale_fill_manual(name = "Effect",
                     breaks = c("Beneficial", "Neutral", "Deleterious"),
                     values = c("Beneficial" = "#388E3C",
                                "Neutral" = "#303F9F",
                                "Deleterious" = "#D32F2F"),
                    guide = F) +
  scale_y_continuous(breaks = seq(-20,110,20), minor_breaks = seq(-20,110,5)) +
  scale_x_discrete(limits=c("Deleterious","Neutral","Beneficial")) +
  labs(title = "LID",
  # labs(title = "MCAT",
       subtitle = "Effect Distribution by Hour",
       x = "Effect",
       y = "Perc. (%)") +
  theme_linedraw() +
  theme(legend.position = "right") +
  coord_cartesian(ylim = c(0,100))

ggsave(sprintf("%s%s_EDIS_ALL.png",
               out_path,expt_name),
       edis.all,
       width = 10, height = 10)

##### LOOKING AT THE ORIGINAL DATA
oridata <- dbGetQuery(conn, "select * from 4C3_GA3_RND_BEAN_6144_DATA")
oridata <- oridata[oridata$orf_name != 'BF_control' & !is.na(oridata$orf_name) & oridata$hours != 0,]
oridata$from <- NULL
oridata$from[oridata$hours > oridata$rnd_hrs] <- 'less'
oridata$from[oridata$hours < oridata$rnd_hrs] <- 'more'

ggplot(oridata) +
  geom_bar(aes(x = from, y = (..count..)/9128*100, fill = from)) +
  facet_wrap(.~hours, nrow = 3, ncol = 3) +
  scale_fill_manual(name = "Colony Size",
                    breaks = c("less", "more"),
                    values = c("more" = "#388E3C",
                               "less" = "#D32F2F"),
                    labels = c("Less","More"),
                    guide = F) +
  scale_y_continuous(breaks = seq(-20,110,20), minor_breaks = seq(-20,110,5)) +
  scale_x_discrete(limits=c("less","more"),
                   labels = c("Less","More")) +
  labs(title = "Original Data",
       subtitle = "Is the query hour more or less than ref hour?",
       x = "Query Hour",
       y = "Perc. (%)") +
  theme_linedraw() +
  theme(legend.position = "right") +
  coord_cartesian(ylim = c(0,100))
# ggsave(sprintf("%s%s_OCS_ALL.png",
#                out_path,expt_name),
#        width = 10, height = 10)

##### T-TEST
sig.dat <- NULL
i = 0
for (hr in sort(unique(stats.all$hours))) {
  i = i + 1
  sig.dat$hours[i] <- hr
  sig.dat$con_ben[i] <- sum(oridata$hours == hr & oridata$from == 'more')/dim(oridata[oridata$hours == hr,])[1] * 100
  sig.dat$pre_ben[i] <- sum(stats.all$hours == hr & stats.all$effect_p == 'Beneficial')/dim(stats.all[stats.all$hours == hr,])[1] * 100
  sig.dat$con_del[i] <- sum(oridata$hours == hr & oridata$from == 'less')/dim(oridata[oridata$hours == hr,])[1] * 100
  sig.dat$pre_del[i] <- sum(stats.all$hours == hr & stats.all$effect_p == 'Deleterious')/dim(stats.all[stats.all$hours == hr,])[1] * 100
  sig.dat$p[i] <- t.test(c(sig.dat$con_ben[i],sig.dat$con_del[i]), c(sig.dat$pre_ben[i],sig.dat$pre_del[i]),paired=TRUE)$p.value
  sig.dat$p_unpaired[i] <- t.test(c(sig.dat$con_ben[i],sig.dat$con_del[i]), c(sig.dat$pre_ben[i],sig.dat$pre_del[i]))$p.value
}
sig.dat <- data.frame(sig.dat)
save(sig.dat, file = "figs/lid_paper/sig_bean.RData")

##### LOAD BACK LID AND BEAN DIS DATA
load("figs/lid_paper/sig_lid.RData")
sig.lid <- sig.dat
load("figs/lid_paper/sig_bean.RData")
sig.bean <- sig.dat
load("figs/lid_paper/rnd_pred.RData")

dis.dat <- rbind(cbind(sig.lid$hours, sig.lid$con_ben, sig.lid$con_del, sig.lid$pre_ben, sig.lid$pre_del, 100 -  sig.lid$pre_ben - sig.lid$pre_del),
                 cbind(sig.lid$hours, sig.lid$con_ben, sig.lid$con_del, sig.bean$pre_ben, sig.bean$pre_del, 100 - sig.bean$pre_ben - sig.bean$pre_del))
colnames(dis.dat) <- c('hours', 'con_ben', 'con_del',
                       'ben', 'del', 'neu')
dis.dat <- data.frame(dis.dat)
dis.dat$method[1:9] <- 'LID'
dis.dat$method[10:18] <- 'BEAN'

rmse.lid.ben <- sqrt(mean((dis.dat$con_ben[dis.dat$method == 'LID'] - dis.dat$ben[dis.dat$method == 'LID'])^2, na.rm = T))
rmse.lid.del <- sqrt(mean((dis.dat$con_del[dis.dat$method == 'LID'] - dis.dat$del[dis.dat$method == 'LID'])^2, na.rm = T))
rmse.bean.ben <- sqrt(mean((dis.dat$con_ben[dis.dat$method == 'BEAN'] - dis.dat$ben[dis.dat$method == 'BEAN'])^2, na.rm = T))
rmse.bean.del <- sqrt(mean((dis.dat$con_del[dis.dat$method == 'BEAN'] - dis.dat$del[dis.dat$method == 'BEAN'])^2, na.rm = T))

sig.lid$method <- 'lid'
sig.bean$method <- 'BEAN'

sig.dat <- rbind(sig.lid, sig.bean)
rnd.pred <- sig.dat[c(1:5,8)]
# save(rnd.pred, file = "figs/lid_paper/rnd_pred.RData")

# R2.lid.ben <- summary(lm(con_ben ~ ben, data = dis.dat))[8]
# R2.bean.ben <- summary(lm(con_ben ~ bean_ben, data = dis.dat))[8]
# R2.lid.del <- summary(lm(con_del ~ lid_del, data = dis.dat))[8]
# R2.bean.del <- summary(lm(con_del ~ bean_del, data = dis.dat))[8]

pre.ben <- ggplot(dis.dat) +
  geom_abline() +
  geom_point(aes(x = con_ben, y = ben, col = method), size = 5) +
  labs(title = 'LID does better than MCAT',
       subtitle = sprintf('Beneficial Fitness Effects\nLID RMSE = %.2f | MCAT RMSE = %.2f', rmse.lid.ben, rmse.bean.ben),
       x = 'Condition Positive (%)',
       y = 'Predicted Positive (%)') +
  scale_color_discrete(name = 'Method') +
  theme_linedraw() +
  coord_cartesian(xlim = c(0, 100),
                  ylim = c(0, 100))

pre.del <- ggplot(dis.dat) +
  geom_abline() +
  geom_point(aes(x = con_del, y = del, col = method), size = 5) +
  labs(title = '',
       subtitle = sprintf('Deleterious Fitness Effects\nLID RMSE = %.2f | MCAT RMSE = %.2f', rmse.lid.del, rmse.bean.del),
       x = 'Condition Positive (%)',
       y = '') +
  scale_color_discrete(name = 'Method') +
  theme_linedraw() +
  coord_cartesian(xlim = c(0, 100),
                  ylim = c(0, 100))

ggarrange(pre.ben, pre.del, nrow = 1,
          common.legend = T, legend = 'right')
ggsave(sprintf("%s%s_PREDICTIONS.png",
               out_path,expt_name),
       width = 10, height = 5)

ggplot(dis.dat) +
  geom_point(aes(x = con_ben, y = ben + del, col = method))

##### DETAILED LOOK AT THE PREDICTIONS
pred.ref = dbGetQuery(conn, sprintf('select hours, average
                                    from %s_6144_FITNESS
                                    where orf_name = "BF_control"',
                                    expt_name))

pred.bean = dbGetQuery(conn, sprintf('select b.*, a.p, a.stat
                from %s_BEAN_6144_PVALUE a, %s_6144_DATA b
                where a.orf_name = b.orf_name and a.hours = b.hours
                order by a.hours, b.6144plate, b.6144col, b.6144row',
                                     expt_name, expt_name))
pred.bean$effect <- NULL
pred.bean$effect[pred.bean$p < 0.05 & pred.bean$stat > 0] = 'Beneficial'
pred.bean$effect[pred.bean$p < 0.05 & pred.bean$stat < 0] = 'Deleterious'
pred.bean$effect[is.na(pred.bean$effect)] = 'Neutral'

pred.bean$effect_con <- NULL
pred.bean$effect_con[pred.bean$hours < pred.bean$rnd_hrs] = 'Beneficial'
pred.bean$effect_con[pred.bean$hours > pred.bean$rnd_hrs] = 'Deleterious'

pred.lid = dbGetQuery(conn, sprintf('select b.*, a.p, a.stat
                from %s_6144_PVALUE a, %s_6144_DATA b
                where a.orf_name = b.orf_name and a.hours = b.hours
                order by a.hours, b.6144plate, b.6144col, b.6144row',
                                    expt_name, expt_name))
pred.lid$effect <- NULL
pred.lid$effect[pred.lid$p < 0.05 & pred.lid$stat > 0] = 'Beneficial'
pred.lid$effect[pred.lid$p < 0.05 & pred.lid$stat < 0] = 'Deleterious'
pred.lid$effect[is.na(pred.lid$effect)] = 'Neutral'

for (hr in sort(unique(pred.ref$hours))) {
  pred.bean$ref[pred.bean$hours == hr] = median(pred.ref$average[pred.ref$hours == hr], na.rm = T)
  pred.lid$ref[pred.lid$hours == hr] = median(pred.ref$average[pred.ref$hours == hr], na.rm = T)
  for (rnd_hrs in sort(unique(pred.bean$rnd_hrs[pred.bean$hours == hr]))) {
    pred.bean$que[pred.bean$hours == hr & pred.bean$rnd_hrs == rnd_hrs] =
      median(pred.bean$average[pred.bean$hours == hr & pred.bean$rnd_hrs == rnd_hrs], na.rm = T)
    pred.lid$que[pred.lid$hours == hr & pred.lid$rnd_hrs == rnd_hrs] =
      median(pred.lid$average[pred.lid$hours == hr & pred.lid$rnd_hrs == rnd_hrs], na.rm = T) 
  }
}

# for (rnd_hrs in sort(unique(pred.bean$rnd_hrs))) {
#   pred.bean$que[pred.bean$rnd_hrs == rnd_hrs] =
#     median(pred.bean$average[pred.bean$rnd_hrs == rnd_hrs], na.rm = T)
#   pred.lid$que[pred.lid$rnd_hrs == rnd_hrs] =
#     median(pred.lid$average[pred.lid$rnd_hrs == rnd_hrs], na.rm = T) 
# }

pred.mcat.box <- ggplot(pred.bean[pred.bean$hours > 16,]) +
    geom_vline(xintercept = pred.bean$ref[pred.bean$hours == hr], lwd = 0.8,
               linetype = 'dashed', col = 'black') +
    geom_histogram(aes(x = que, fill = effect),
                   binwidth = 10,
                   position = 'stack',
                   alpha = 0.8) +
    scale_fill_manual(name = 'Effects',
                      breaks = c('Beneficial', 'Deleterious', 'Neutral'),
                      values = c('Beneficial' = '#388E3C',
                                 'Deleterious' = '#FF5252',
                                 'Neutral' = '#303F9F')) +
    labs(title = 'Effects Detected by MCAT',
         subtitle = 'per median cs (pix) of references',
         x = 'Median CS (pix)',
         y = 'Count') +
    facet_wrap(.~ref, nrow = 1, ncol = 2) +
    theme_linedraw() +
    coord_cartesian(xlim = c(0,400),
                    ylim = c(0, 1000))
  # ggsave(sprintf("%sBOX_PREDS_MCAT.png",
  #                out_path),
  #        width = 6, height = 5) 

ggplot(pred.lid[pred.lid$hours > 16,]) +
    geom_vline(xintercept = pred.lid$ref[pred.lid$hours == hr], lwd = 0.8,
               linetype = 'dashed', col = 'black') +
    geom_histogram(aes(x = que, fill = effect),
                   binwidth = 10,
                   position = 'stack',
                   alpha = 0.8) +
    scale_fill_manual(name = 'Effects',
                      breaks = c('Beneficial', 'Deleterious', 'Neutral'),
                      values = c('Beneficial' = '#388E3C',
                                 'Deleterious' = '#FF5252',
                                 'Neutral' = '#303F9F')) +
    labs(title = 'Effects Detected by LID',
         subtitle = 'per median cs (pix) of references',
         x = 'Median CS (pix)',
         y = 'Count') +
    facet_wrap(.~ref, nrow = 1, ncol = 2) +
    theme_linedraw() +
    coord_cartesian(xlim = c(0,400),
                    ylim = c(0, 1000))
  
  ggsave(sprintf("%sBOX_PREDS_LID_%d.png",
                 out_path, hr),
         width = 6, height = 5) 



