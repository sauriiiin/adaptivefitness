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
library(ggrepel)
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
expt_name <- '4C4_FS_CC'
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

# save(fit.data, file = "output/4C4_rnd_fit_data.RData")
# save(ref.data, file = "output/4C4_rnd_ref_data.RData")
# save(sta.data, file = "output/4C4_rnd_sta_data.RData")

load("output/4C4_rnd_fit_data.RData")
load("output/4C4_rnd_ref_data.RData")
load("output/4C4_rnd_sta_data.RData")

sta.data$res <- sta.data$effect_p
sta.data$res[sta.data$effect_p == 'Beneficial' & sta.data$from == 'more'] <- '5TrueBeneficial'
sta.data$res[sta.data$effect_p == 'Beneficial' & sta.data$from == 'less'] <- '3FalseBeneficial'
sta.data$res[sta.data$effect_p == 'Deleterious' & sta.data$from == 'more'] <- '4FalseDeleterious'
sta.data$res[sta.data$effect_p == 'Deleterious' & sta.data$from == 'less'] <- '2TrueDeleterious'
sta.data$res[sta.data$effect_p == 'Neutral' & sta.data$from == 'more'] <- '6BenNeutral'
sta.data$res[sta.data$effect_p == 'Neutral' & sta.data$from == 'less'] <- '1DelNeutral'

sta.data <- sta.data[!sta.data$hours == 0,]
sta.data$method[sta.data$method == "BEAN"] <- "MCAT"

pie.dat <- data.frame(rbind(cbind(hours = sta.data$hours[sta.data$method == 'MCAT'],
                                  value = sta.data$from[sta.data$method == 'MCAT'],
                                  method = 'TRUTH'),
                            cbind(hours = sta.data$hours[sta.data$method == 'MCAT'],
                                  value = sta.data$res[sta.data$method == 'MCAT'],
                                  method = 'MCAT'),
                            cbind(hours = sta.data$hours[sta.data$method == 'LID'],
                                  value = sta.data$res[sta.data$method == 'LID'],
                                  method = 'LID')))

pie.dat$value[pie.dat$value == 'less'] <- '2TrueDeleterious'
pie.dat$value[pie.dat$value == 'more'] <- '5TrueBeneficial'

pie.dat <- arrange(transform(pie.dat,
                           method=factor(method,levels=c('LID','TRUTH','MCAT'))),
                   method)

pie.per <- plyr::count(pie.dat, vars = c('method','value'))
pie.per$per <- pie.per$freq/(sum(pie.per$freq)/3) * 100

rnd_pie <- ggplot(data = pie.per, aes(x = "", y = per, fill = value)) +
  geom_bar(stat = 'identity') +
  coord_polar("y", start = 0) +
  # geom_label(aes(label = sprintf("%0.2f%%",per)),
  #           position = position_stack(vjust = 0.5),
  #           label.padding = unit(0.25, "lines")) +
  geom_label_repel(aes(label = sprintf("%0.2f%%",per)),
                  position = position_stack(vjust = 0.5),
                  show.legend = F) +
  scale_fill_manual(name = 'Effects',
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
                               '6BenNeutral'='False-Neutral Beneficial')) +
  
  facet_wrap(.~method) +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = 'bottom',
        text = element_text(size = 15))
ggsave(sprintf("%s%s_VIR_PLATE2_PIE.jpg",out_path,expt_name), rnd_pie,
       height = 5, width = 12,
       dpi = 300)


# rnd_pie <- ggplot(pie.dat) +
#   geom_bar(aes(x = "", y = "", fill = value), stat = 'identity') +
#   coord_polar("y", start = 0) +
#   scale_fill_manual(name = 'Effects',
#                     breaks = c('6BenNeutral',
#                                '5TrueBeneficial',
#                                '4FalseDeleterious',
#                                '3FalseBeneficial',
#                                '2TrueDeleterious',
#                                '1DelNeutral'),
#                     values = c('1DelNeutral'='#536DFE',
#                                '2TrueDeleterious'='#303F9F',
#                                '3FalseBeneficial'='#C5CAE9',
#                                '4FalseDeleterious'='#FFECB3',
#                                '5TrueBeneficial'='#FFA000',
#                                '6BenNeutral'='#FFC107'),
#                     labels = c('1DelNeutral'='False-Neutral Deleterious',
#                                '2TrueDeleterious'='True Deleterious',
#                                '3FalseBeneficial'='False-Beneficial Deleterious',
#                                '4FalseDeleterious'='False-Deleterious Beneficial',
#                                '5TrueBeneficial'='True Beneficial',
#                                '6BenNeutral'='False-Neutral Beneficial')) +
#                       
#   facet_wrap(.~method) +
#   theme_linedraw() +
#   theme(axis.title = element_blank(),
#         axis.ticks = element_blank(),
#         legend.position = 'bottom',
#         text = element_text(size = 10))
# ggsave(sprintf("%s%s_VIR_PLATE2_PIE.jpg",out_path,expt_name), rnd_pie,
#        height = 5, width = 12,
#        dpi = 300)

rnd_data <- dbGetQuery(conn, 'select * from 4C4_FS_RND_6144_FITNESS a, 4C4_pos2coor b
                                where a.pos = b.pos')
i <- 1
for (h in unique(rnd_data$hours)) {
  rnd_data$plate[rnd_data$hours == h] <- i
  i <- i + 1
}
rnd_data$colony[rnd_data$orf_name == 'BF_control'] <- 'Reference'
rnd_data$colony[rnd_data$orf_name != 'BF_control'] <- 'Query'
# rnd_data$colony[is.na(rnd_data$orf_name)] <- 'Gap'

vir_plates2 <- ggplot(rnd_data) +
  geom_line(aes(x = average, col = colony),
            stat = 'density', trim = T, lwd = 1.2) +
  labs(title = 'Virtual Plates',
       subtitle = 'With Random Colony Size Distribution',
       x = 'Colony Size (pix. count)',
       y = 'Density') +
  scale_color_manual(name = 'Colony\nType',
                     breaks = c("Reference","Query"),
                     values = c("Reference" = "#3F51B5",
                                "Query" = "#FFC107")) +
  facet_wrap(.~hours) +
  coord_cartesian(ylim = c(0,0.015)) +
  theme_linedraw()
ggsave(sprintf("%s%s_VIR_PLATE2.jpg",out_path,expt_name), vir_plates2,
       width = 7, height = 6,
       dpi = 1000)


