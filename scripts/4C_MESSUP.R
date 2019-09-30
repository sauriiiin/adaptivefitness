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
out_path = 'figs/lid_paper/';
dat.dir <- "/home/sbp29/R/Projects/proto_plots/rawdata/4C3_GA3_RND_LID/"
expt_name <- '4C3_GA3_RND'
pvals = seq(0,1,0.005)

##### MAKING THE MESS
# alldat = dbGetQuery(conn, 'select b.*, a.orf_name, a.hours, a.average
#               from 4C3_GA3_CC2_6144_FITNESS a, 4C3_pos2coor6144 b
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
# dbWriteTable(conn, "4C3_GA3_CC2_RND_6144_DATA", alldat, overwrite = T)
# # 
# # # dbExecute(conn, 'update 4C3_GA3_RND_6144_DATA
# # #     set average = NULL
# # #     where pos in
# # #     (select pos from 4C3_borderpos)')
# # # 
# dbExecute(conn, 'create table 4C3_GA3_CC2_RND_6144_JPEG
#     select pos, hours, rnd_hrs, average
#     from 4C3_GA3_CC2_RND_6144_DATA')
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
  hours <- c(hours, as.numeric(s[5]))
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
save(sig.dat, file = "figs/lid_paper/sig_lid.RData")

##### LOAD BACK LID AND BEAN DIS DATA
load("figs/lid_paper/sig_lid.RData")
sig.lid <- sig.dat
load("figs/lid_paper/sig_bean.RData")
sig.bean <- sig.dat
load("figs/lid_paper/rnd_pred.RData")

dis.dat <- cbind(sig.lid$hours, sig.lid$con_ben, sig.lid$con_del, sig.lid$pre_ben, sig.lid$pre_del, 100 -  sig.lid$pre_ben - sig.lid$pre_del,
                 sig.bean$pre_ben, sig.bean$pre_del, 100 - sig.bean$pre_ben - sig.bean$pre_del)
colnames(dis.dat) <- c('hours', 'con_ben', 'con_del',
                       'lid_ben', 'lid_del', 'lid_neu',
                       'bean_ben', 'bean_del', 'bean_neu')
dis.dat <- data.frame(dis.dat)

rmse.lid.ben <- sqrt(mean((dis.dat$con_ben - dis.dat$lid_ben)^2, na.rm = T))
rmse.lid.del <- sqrt(mean((dis.dat$con_del - dis.dat$lid_del)^2, na.rm = T))
rmse.bean.ben <- sqrt(mean((dis.dat$con_ben - dis.dat$bean_ben)^2, na.rm = T))
rmse.bean.del <- sqrt(mean((dis.dat$con_del - dis.dat$bean_del)^2, na.rm = T))

sig.lid$method <- 'lid'
sig.bean$method <- 'mcat'

sig.dat <- rbind(sig.lid, sig.bean)
rnd.pred <- sig.dat[c(1:5,8)]
# save(rnd.pred, file = "figs/lid_paper/rnd_pred.RData")

R2.lid.ben <- summary(lm(con_ben ~ lid_ben, data = dis.dat))[8]
R2.bean.ben <- summary(lm(con_ben ~ bean_ben, data = dis.dat))[8]
R2.lid.del <- summary(lm(con_del ~ lid_del, data = dis.dat))[8]
R2.bean.del <- summary(lm(con_del ~ bean_del, data = dis.dat))[8]

pre.ben <- ggplot(dis.dat) +
  geom_abline() +
  geom_point(aes(x = con_ben, y = lid_ben, col = 'LID'), size = 5) +
  geom_point(aes(x = con_ben, y = bean_ben, col = 'MCAT'), size = 5) +
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
  geom_point(aes(x = con_del, y = lid_del, col = 'LID'), size = 5) +
  geom_point(aes(x = con_del, y = bean_del, col = 'MCAT'), size = 5) +
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

ggplot(rnd.pred) +
  geom_point(aes(x = con_ben, y = pre_ben, col = method))

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



