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
out_path = 'figs/mess/';
dat.dir <- "/home/sbp29/R/Projects/proto_plots/rawdata/4C3_GA3_RND_BEAN/"
expt_name <- '4C3_GA3_RND_BEAN'
pvals = seq(0,1,0.005)

##### MAKING THE MESS
# alldat = dbGetQuery(conn, 'select b.*, a.orf_name, a.hours, a.average
#               from 4C3_GA3_RAW_6144_FITNESS a, 4C3_pos2coor6144 b
#               where a.pos = b.pos
#               order by hours, 6144plate, 6144col, 6144row')
# alldat$rnd_hrs <- alldat$hours
# 
# for (orf in unique(alldat$orf_name[!is.na(alldat$orf_name) & alldat$orf_name != "BF_control"])) {
#   for (hr in unique(alldat$hours)) {
#     hr.tmp <- sample(unique(alldat$hours[alldat$hours != hr]),1)
#     alldat$average[!is.na(alldat$orf_name) & alldat$orf_name == orf & alldat$hours == hr] <-
#       alldat$average[!is.na(alldat$orf_name) & alldat$orf_name == orf & alldat$hours == hr.tmp]
#     alldat$rnd_hrs[!is.na(alldat$orf_name) & alldat$orf_name == orf & alldat$hours == hr] <-
#       alldat$hours[!is.na(alldat$orf_name) & alldat$orf_name == orf & alldat$hours == hr.tmp]
#   }
# }
# 
# # dbWriteTable(conn, "4C3_GA3_RND_6144_DATA", alldat, overwrite = T)
# dbWriteTable(conn, "4C3_GA3_RND_BEAN_6144_DATA", alldat, overwrite = T)
# 
# # dbExecute(conn, 'update 4C3_GA3_RND_6144_DATA
# #     set average = NULL
# #     where pos in
# #     (select pos from 4C3_borderpos)')
# # 
# dbExecute(conn, 'create table 4C3_GA3_RND_6144_JPEG
#     select pos, hours, rnd_hrs, average
#     from 4C3_GA3_RND_6144_DATA')
# 
# dbExecute(conn, 'create table 4C3_GA3_RND_BEAN_6144_JPEG
#     select pos, hours, rnd_hrs, average
#     from 4C3_GA3_RND_BEAN_6144_DATA')

##### ANALYZING THE MESS
stats.files <- list.files(path = dat.dir,
                          pattern = "P.csv", recursive = TRUE)
fit.files <- list.files(path = dat.dir,
                        pattern = "S.csv", recursive = TRUE)
hours <- NULL
reps <- NULL
for (s in strsplit(stats.files,'_')) {
  # reps <- c(reps, as.numeric(s[4]))
  reps <- 8
  hours <- c(hours, as.numeric(s[6]))
}
reps <- unique(reps)
hours <- unique(hours)

stats.all <- NULL
fit.all <- NULL
fpr.all <- NULL

# PUTTING IT TOGETHER
for (ii in 1:length(reps)) {
  rep <- reps[ii]
  for (i in 2:length(hours)) {
    hr <- hours[i]
    
    dat.stats <- read.csv(paste0(dat.dir,
                                 sprintf('%s_%d_%d_STATS_P.csv',expt_name,rep,hr)),
                          na.strings = "NaN")
    # dat.stats <- read.csv(paste0(dat.dir,
    #                              sprintf('%s_%d_STATS_P.csv',expt_name,hr)),
    #                       na.strings = "NaN")
    dat.stats <- dat.stats[dat.stats$hours != 0,]
    dat.stats$cont_hrs <- hr
    dat.stats$rep <- rep
    dat.fit <- read.csv(paste0(dat.dir,
                               sprintf('%s_%d_%d_FITNESS.csv',expt_name,rep,hr)),
                        na.strings = "NaN")
    # dat.fit <- read.csv(paste0(dat.dir,
    #                            sprintf('%s_%d_FITNESS.csv',expt_name,hr)),
    #                     na.strings = "NaN")
    dat.fit$cont_hrs <- hr
    dat.fit$rep <- rep
    dat.fit$se <- dat.fit$average - dat.fit$bg
    
    dat.stats$pthresh <- quantile(sort(dat.stats$p[dat.stats$hours == hr]),.05)
    for (h in unique(dat.stats$hours)) {
      cont.mean <- mean(dat.fit$fitness[dat.fit$hours == h & dat.fit$orf_name == 'BF_control' & !is.na(dat.fit$fitness)])
      dat.stats$es[dat.stats$hours == h] <- round(dat.stats$cs_mean[dat.stats$hours == h]/cont.mean,4)
      
      dat.stats$cen[dat.stats$hours == h] <- mean(dat.stats$es[dat.stats$hours == h])
      
      dat.stats$effect[dat.stats$hours == h & dat.stats$p <= dat.stats$pthresh & dat.stats$cs_mean > cont.mean] <- 'Beneficial'
      dat.stats$effect[dat.stats$hours == h & dat.stats$p <= dat.stats$pthresh & dat.stats$cs_mean < cont.mean] <- 'Deleterious'
      dat.stats$effect[dat.stats$hours == h & is.na(dat.stats$effect)] <- 'Neutral'
      
      dat.stats$effect_p[dat.stats$hours == h & dat.stats$p <= 0.05 & dat.stats$cs_mean > cont.mean] <- 'Beneficial'
      dat.stats$effect_p[dat.stats$hours == h & dat.stats$p <= 0.05 & dat.stats$cs_mean < cont.mean] <- 'Deleterious'
      dat.stats$effect_p[dat.stats$hours == h & is.na(dat.stats$effect_p)] <- 'Neutral'
    }
    
    stats.all <- rbind(stats.all,dat.stats)
    fit.all <- rbind(fit.all,dat.fit)
  }
}

fit.all$strain[fit.all$orf_name == "BF_control"] = "Reference"
fit.all$strain[fit.all$orf_name != "BF_control"] = "Query"

for (hr in sort(unique(stats.all$hours))) {
  for (orf in unique(stats.all$orf_name)) {
    stats.all$pix_mean[stats.all$hours == hr & stats.all$orf_name == orf] <-
      mean(fit.all$average[fit.all$hours == hr & fit.all$orf_name == orf], na.rm = T)
  }
}

ref.all <- NULL
i = 0
for (hr in sort(unique(fit.all$hours))) {
  for (pos in 385:768) {
    i <- i + 1
    ref.all$pix_mean[i] <-
      mean(fit.all$average[fit.all$pos %in% c(c(110000,120000,130000,140000,210000,220000,230000,240000) + pos) & fit.all$hours == hr], na.rm = T)
    ref.all$cs_mean[i] <-
      mean(fit.all$fitness[fit.all$pos %in% c(c(110000,120000,130000,140000,210000,220000,230000,240000) + pos) & fit.all$hours == hr], na.rm = T)
    ref.all$hours[i] <- hr
  }
}
ref.all <- data.frame(ref.all)

for (hr in sort(unique(ref.all$hours))) {
  pix_mean <- mean(ref.all$pix_mean[ref.all$hours == hr],na.rm = T)
  stats.all$from[stats.all$hours == hr & stats.all$pix_mean < pix_mean] = 'less'
  stats.all$from[stats.all$hours == hr & stats.all$pix_mean > pix_mean] = 'more'
}

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
  # labs(title = "LIDetector",
  labs(title = "BEAN Method",
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
  # labs(title = "LIDetector",
  labs(title = "BEAN Method",
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
oridata <- dbGetQuery(conn, "select * from 4C3_GA3_RND_6144_DATA")
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
ggsave(sprintf("%s%s_OCS_ALL.png",
               out_path,expt_name),
       width = 10, height = 10)

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
save(sig.dat, file = "figs/mess/sig_bean.RData")

##### LOAD BACK LID AND BEAN DIS DATA
load("figs/mess/sig_lid.RData")
sig.lid <- sig.dat
load("figs/mess/sig_bean.RData")
sig.bean <- sig.dat

dis.dat <- cbind(sig.lid$hours, sig.lid$con_ben, sig.lid$con_del, sig.lid$pre_ben, sig.lid$pre_del, sig.bean$pre_ben, sig.bean$pre_del)
colnames(dis.dat) <- c('hours', 'con_ben', 'con_del',
                       'lid_ben', 'lid_del',
                       'bean_ben', 'bean_del')
dis.dat <- data.frame(dis.dat)

pre.ben <- ggplot(dis.dat) +
  geom_abline() +
  geom_point(aes(x = con_ben, y = lid_ben, col = 'LID'), size = 5) +
  geom_point(aes(x = con_ben, y = bean_ben, col = 'BEAN'), size = 5) +
  labs(title = 'How good are the predictions?',
       subtitle = 'Beneficial Predictions',
       x = 'Condition Positive',
       y = 'Predicted Positive') +
  scale_color_discrete(name = 'Method') +
  theme_linedraw() +
  coord_cartesian(xlim = c(0, 100),
                  ylim = c(0, 100))

pre.del <- ggplot(dis.dat) +
  geom_abline() +
  geom_point(aes(x = con_del, y = lid_del, col = 'LID'), size = 5) +
  geom_point(aes(x = con_del, y = bean_del, col = 'BEAN'), size = 5) +
  labs(title = '',
       subtitle = 'Deleterious Predictions',
       x = 'Condition Positive',
       y = '') +
  scale_color_discrete(name = 'Method') +
  theme_linedraw() +
  coord_cartesian(xlim = c(0, 100),
                  ylim = c(0, 100))

ggarrange(pre.ben, pre.del, nrow = 1,
          common.legend = T, legend = 'bottom')
ggsave(sprintf("%s%s_PREDICTIONS.png",
               out_path,expt_name),
       width = 10, height = 6)


