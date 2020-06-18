##### POWER DYNAMICS
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 05/23/2019

##### INITIALIZE
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyverse)
# library(egg)
library(ggpubr)
library(stringr)
out_path = 'figs/lid_paper/4C4/';
dat.dir <- "/home/sbp29/R/Projects/adaptivefitness/rawdata/4C4_TR_CC/"
expt_name <- '4C4_FS_NONORM'
pvals = seq(0,1,0.005)

# getmode <- function(v) {
#   uniqv <- unique(v[!is.na(v)])
#   uniqv[which.max(tabulate(match(v, uniqv)))]
# }

##### GET DATA
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
  reps <- c(reps, as.numeric(s[4]))
  # reps <- 8
  hours <- c(hours, as.numeric(s[5]))
  atmpt <- c(atmpt, as.numeric(s[6]))
}
# refs <- unique(refs)
reps <- sort(unique(reps))
hours <- sort(unique(hours))
atmpt <- sort(unique(atmpt))

hhours <- c(0,1.02,1.38,2.9,4.02,4.89,6.14,6.88,7.85,8.96,9.97,11.04)
# hhours <- c(1.02,1.38,2.9,4.02,4.89,6.14,6.88,7.85,8.96,9.97,11.04)

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
      # dat.stats <- read.csv(paste0(dat.dir,
      #                              sprintf('%s_%d_%d_STATS_P.csv',expt_name,rep,hr)),
      #                       na.strings = "NaN")

      # dat.stats <- dat.stats[dat.stats$hours != 0,]
      dat.stats$cont_hrs <- hhours[i]
      # dat.stats$cont_hrs <- hr
      dat.stats$rep <- rep
      dat.fit <- read.csv(paste0(dat.dir,
                                 sprintf('%s_%d_%d_%d_FITNESS.csv',expt_name,rep,hr,iii)),
                          na.strings = "NaN")
      # dat.fit <- read.csv(paste0(dat.dir,
      #                            sprintf('%s_%d_%d_FITNESS.csv',expt_name,rep,hr)),
      #                     na.strings = "NaN")
      dat.fit$cont_hrs <- hhours[i]
      # dat.fit$cont_hrs <- hr
      dat.fit$attempt <- iii
      dat.fit$rep <- rep
      dat.fit$se <- dat.fit$average - dat.fit$bg
      
      dat.stats$attempt <- iii
      dat.stats$replicates <- rep
      dat.stats$pthresh <- quantile(sort(dat.stats$p[dat.stats$hours == hhours[i]]),.05)
      # dat.stats$pthresh <- quantile(sort(dat.stats$p[dat.stats$hours == hr]),.05)
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

stats.all <- stats.all[!stats.all$orf_name %in% sprintf('MASK%d',1:1000),]
fit.all <- fit.all[!fit.all$orf_name %in% sprintf('MASK%d',1:1000),]
# save(stats.all, file = sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_STATS.RData", expt_name))
# save(fit.all, file = sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_FITNESS.RData", expt_name))

load(sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_STATS.RData", expt_name))
load(sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_FITNESS.RData", expt_name))
atmpt <- sort(unique(stats.all$attempt))
reps <- sort(unique(stats.all$replicates))

##### THE STATS DATA ANALYSIS
# stats.all$es <- round(stats.all$es,4)
dat.cnt.all <- NULL
spe.all <- NULL
sen.all <- NULL
for (a in atmpt) {
  for (rep in unique(reps)) {
    # effect_size <- sort(unique(round(stats.all$es[stats.all$rep == rep],2)))
    # dat.pow <- NULL
    # 
    # ##### POWER CALCULATIONS
    # for (es in effect_size) {
    #   N <- sum(stats.all$rep == rep & 
    #              stats.all$es > es-5e-03 & stats.all$es < es+5e-03)
    #   if (N > 0) {
    #     TP <- sum(stats.all$rep == rep &
    #                 stats.all$es > es-5e-03 & stats.all$es < es+5e-03 &
    #                 stats.all$p <= stats.all$pthresh & stats.all$hours != stats.all$cont_hrs)
    #     FP <- sum(stats.all$rep == rep &
    #                 stats.all$es > es-5e-03 & stats.all$es < es+5e-03 &
    #                 stats.all$p <= stats.all$pthresh & stats.all$hours == stats.all$cont_hrs)
    #     FPR <- FP/sum(stats.all$rep == rep &
    #                     stats.all$es > es-5e-03 & stats.all$es < es+5e-03 &
    #                     stats.all$hours == stats.all$cont_hrs) * 100
    #     TN <- sum(stats.all$rep == rep &
    #                 stats.all$es > es-5e-03 & stats.all$es < es+5e-03 &
    #                 stats.all$p > stats.all$pthresh & stats.all$hours == stats.all$cont_hrs)
    #     FN <- sum(stats.all$rep == rep &
    #                 stats.all$es > es-5e-03 & stats.all$es < es+5e-03 &
    #                 stats.all$p > stats.all$pthresh & stats.all$hours != stats.all$cont_hrs)
    #     POW <- TP/N * 100
    #     SEN <- TP/sum(stats.all$rep == rep &
    #                     stats.all$es > es-5e-03 & stats.all$es < es+5e-03 &
    #                     stats.all$hours != stats.all$cont_hrs) * 100
    #     SPE <- TN/sum(stats.all$rep == rep &
    #                     stats.all$es > es-5e-03 & stats.all$es < es+5e-03 &
    #                     stats.all$hours == stats.all$cont_hrs) * 100
    #     ACC <- (TP + TN)/N * 100
    #     dat.pow <- rbind(dat.pow, c(es,N,TP,FP,FPR,TN,FN,POW,SEN,SPE,ACC))
    #   }
    # }
    # dat.pow <- data.frame(dat.pow)
    # colnames(dat.pow) <- c("es","N","TP","FP","FPR","TN","FN","pow","sen","spe","acc")
    
    ##### BOX PLOTS OF EFFECT DISTRIBUTION
    stats.tmp <- stats.all[stats.all$rep == rep & stats.all$attempt == a &
                             stats.all$cont_hrs > 2 & stats.all$hours > 2,]
    
    for (ii in unique(stats.tmp$cont_hrs)) {
      for (pp in sort(unique(stats.tmp$p[stats.tmp$hours == ii & stats.tmp$cont_hrs == stats.tmp$hours]))) {
        stats.tmp$fpr[stats.tmp$hours == ii & stats.tmp$cont_hrs == stats.tmp$hours & stats.tmp$p == pp] <-
          dim(stats.tmp[stats.tmp$hours == ii & stats.tmp$cont_hrs == stats.tmp$hours & stats.tmp$p <= pp,])[1]/
          dim(stats.tmp[stats.tmp$hours == ii & stats.tmp$cont_hrs == stats.tmp$hours,])[1]
      }
    }
    
    srt <- sort(unique(stats.tmp$cen))
    dat.srt <- NULL
    t <- 1
    for (cen in srt) {
      stats.tmp$pos <- t
      dat.srt <- rbind(dat.srt, stats.tmp[stats.tmp$cen == cen,])
      t <- t + 1
    }
    
    ##### LINE PLOT OF EFFECTS
    temp <- NULL
    dat.cnt <- data.frame()
    for (pos in unique(dat.srt$pos)) {
      temp$hours <- dat.srt$hours[dat.srt$pos == pos][1]
      temp$cont_hrs <- dat.srt$cont_hrs[dat.srt$pos == pos][1]
      temp$cen <- mean(dat.srt$cen[dat.srt$pos == pos])
      temp$Deleterious <- sum(dat.srt$pos == pos & dat.srt$effect == 'Deleterious', na.rm = T)
      temp$Neutral <- sum(dat.srt$pos == pos & dat.srt$effect == 'Neutral', na.rm = T)
      temp$Beneficial <- sum(dat.srt$pos == pos & dat.srt$effect == 'Beneficial', na.rm = T)
      temp$Deleterious_p <- sum(dat.srt$pos == pos & dat.srt$effect_p == 'Deleterious', na.rm = T)
      temp$Neutral_p <- sum(dat.srt$pos == pos & dat.srt$effect_p == 'Neutral', na.rm = T)
      temp$Beneficial_p <- sum(dat.srt$pos == pos & dat.srt$effect_p == 'Beneficial', na.rm = T)
      dat.cnt <- rbind(dat.cnt,temp)
    }
    
    dat.cnt$rep <- rep
    dat.cnt$atmpt <- a
    dat.cnt.all <- rbind(dat.cnt.all, dat.cnt)
    
    dat.cnt2 <- NULL
    dat.cnt2$hours <- dat.cnt$hours
    dat.cnt2$cont_hrs <- dat.cnt$cont_hrs
    dat.cnt2$cen <- dat.cnt$cen
    dat.cnt2$Deleterious <- dat.cnt$Deleterious
    dat.cnt2$Beneficial <- dat.cnt2$Deleterious + dat.cnt$Beneficial
    dat.cnt2$Neutral <- dat.cnt2$Beneficial + dat.cnt$Neutral
    dat.cnt2$Deleterious_p <- dat.cnt$Deleterious_p
    dat.cnt2$Beneficial_p <- dat.cnt2$Deleterious_p + dat.cnt$Beneficial_p
    dat.cnt2$Neutral_p <- dat.cnt2$Beneficial_p + dat.cnt$Neutral_p
    dat.cnt2 <- data.frame(dat.cnt2)
    
    dat.cnt2$rep <- rep
    dat.cnt2$atmpt <- a
    sen.all <- rbind(sen.all,  dat.cnt2)
    
    stats.tmp$rep <- rep
    stats.tmp$atmpt <- a
    spe.all <- rbind(spe.all, stats.tmp)
    
    # save(dat.cnt2,
    #      file = sprintf('%s%sSPEDATA.RData',out_path,expt_name))
    # save(stats.tmp,
    #      file = sprintf('%s%sSENDATA.RData',out_path,expt_name))
    
    # t <- dat.cnt2$Beneficial[1]
    # plt.pow.fdr <- ggplot(dat.cnt2) +
    #   geom_area(aes(x = cen, y = Neutral, fill = 'Neutral'), alpha = 1) +
    #   geom_area(aes(x = cen, y = Beneficial, fill = 'Beneficial'), alpha = 1) +
    #   geom_area(aes(x = cen, y = Deleterious, fill = 'Deleterious'), alpha = 1) +
    #   geom_vline(xintercept = seq(0,2,0.025), col = '#757575', lwd = 0.5, alpha =0.5) +
    #   geom_hline(yintercept = c(t * seq(0,1,0.05)), col = '#757575', lwd = 0.5, alpha =0.5) +
    #   labs(title = "Which effects are detected?",
    #        subtitle = sprintf('with %d technical replicate (FDR = 0.05)',rep),
    #        x = 'Effect Size',
    #        y = 'Sensitivity') +
    #   scale_y_continuous(breaks = c(t * seq(0,1,0.1)),
    #                      minor_breaks = c(t * seq(0,1,0.05)),
    #                      labels = c('0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')) +
    #   scale_x_continuous(breaks = seq(0,2,0.05),
    #                      minor_breaks = seq(0,2,0.025)) +
    #   scale_color_manual(name = 'Effects',
    #                      breaks = c('Beneficial','Neutral','Deleterious'),
    #                      values = c('Deleterious'='#3F51B5',
    #                                 'Neutral'='#212121',
    #                                 'Beneficial'='#FFC107')) +
    #   scale_fill_manual(name = 'Effects',
    #                     breaks = c('Beneficial','Neutral','Deleterious'),
    #                     values = c('Deleterious'='#3F51B5',
    #                                'Neutral'='#212121',
    #                                'Beneficial'='#FFC107')) +
    #   theme_linedraw() +
    #   theme(axis.text.x = element_text(size=8),
    #         axis.title.x = element_text(size=10),
    #         axis.text.y = element_text(size=8),
    #         axis.title.y = element_text(size=10),
    #         # axis.title = element_blank(),
    #         legend.text = element_text(size=8),
    #         legend.title = element_text(size=10),
    #         legend.position = "bottom",
    #         plot.title = element_text(size=12),
    #         plot.subtitle = element_text(size=10)) +
    #   guides(color = guide_legend(override.aes = list(size=2)),
    #          shape = guide_legend(override.aes = list(size=2))) +
    #   coord_cartesian(xlim = c(0.8,1.2),
    #                   ylim = c(0,t))
    # 
    # plt.pow.p <- ggplot(dat.cnt2) +
    #   geom_area(aes(x = cen, y = Neutral_p, fill = 'Neutral'), alpha = 1) +
    #   geom_area(aes(x = cen, y = Beneficial_p, fill = 'Beneficial'), alpha = 1) +
    #   geom_area(aes(x = cen, y = Deleterious_p, fill = 'Deleterious'), alpha = 1) +
    #   geom_vline(xintercept = seq(0,2,0.025), col = '#757575', lwd = 0.5, alpha =0.5) +
    #   geom_hline(yintercept = c(t * seq(0,1,0.05)), col = '#757575', lwd = 0.5, alpha =0.5) +
    #   labs(title = "",
    #        subtitle = sprintf('with %d technical replicate (p <= 0.05)',rep),
    #        x = 'Effect Size',
    #        y = 'Sensitivity') +
    #   scale_y_continuous(breaks = c(t * seq(0,1,0.1)),
    #                      minor_breaks = c(t * seq(0,1,0.05)),
    #                      labels = c('0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')) +
    #   scale_x_continuous(breaks = seq(0,2,0.05),
    #                      minor_breaks = seq(0,2,0.025)) +
    #   scale_color_manual(name = 'Effects',
    #                      breaks = c('Beneficial','Neutral','Deleterious'),
    #                      values = c('Deleterious'='#3F51B5',
    #                                 'Neutral'='#212121',
    #                                 'Beneficial'='#FFC107')) +
    #   scale_fill_manual(name = 'Effects',
    #                     breaks = c('Beneficial','Neutral','Deleterious'),
    #                     values = c('Deleterious'='#3F51B5',
    #                                'Neutral'='#212121',
    #                                'Beneficial'='#FFC107')) +
    #   theme_linedraw() +
    #   theme(axis.text.x = element_text(size=8),
    #         axis.title.x = element_text(size=10),
    #         axis.text.y = element_text(size=8),
    #         axis.title.y = element_blank(),
    #         # axis.title = element_blank(),
    #         legend.text = element_text(size=8),
    #         legend.title = element_text(size=10),
    #         legend.position = "bottom",
    #         plot.title = element_text(size=12),
    #         plot.subtitle = element_text(size=10)) +
    #   guides(color = guide_legend(override.aes = list(size=2)),
    #          shape = guide_legend(override.aes = list(size=2))) +
    #   coord_cartesian(xlim = c(0.8,1.2),
    #                   ylim = c(0,t))
    # 
    # stats.tmp$hours <- as.character(stats.tmp$hours)
    # stats.tmp$cont_hrs <- as.character(stats.tmp$cont_hrs)
    # 
    # plt.fpr <- ggplot(stats.tmp[stats.tmp$hours == stats.tmp$cont_hrs,]) +
    #   geom_abline(col = 'red', linetype = 'dashed', lwd = 1) +
    #   geom_line(aes(x = p, y = fpr, col = hours), lwd = 1.2) +
    #   geom_segment(aes(x=0,xend=0.2,y=0,yend=0), col = 'black', lwd = 0.8) +
    #   geom_segment(aes(x=0.2,xend=0.2,y=0,yend=0.2), col = 'black', lwd = 0.8) +
    #   geom_segment(aes(x=0.2,xend=0,y=0.2,yend=0.2), col = 'black', lwd = 0.8) +
    #   geom_segment(aes(x=0,xend=0,y=0.2,yend=0), col = 'black', lwd = 0.8) +
    #   labs(title = "Does FPR follow random expectation?",
    #        subtitle = "for p between 0 and 1",
    #        x = "p-value cut-off",
    #        y = "False Positive Rate") +
    #   scale_x_continuous(breaks = seq(0,1,0.2),
    #                      minor_breaks = seq(0,1,0.05),
    #                      limits = c(0,1)) +
    #   scale_y_continuous(breaks = seq(0,1,0.2),
    #                      minor_breaks = seq(0,1,0.05),
    #                      limits = c(0,1)) +
    #   scale_color_discrete(name = "Hours",
    #                      breaks=as.character(c(2.9,4.02,4.89,6.14,6.88,7.85,8.96,9.97,11.04))) +
    #   theme_linedraw() +
    #   theme(axis.text.x = element_text(size=8),
    #         axis.title.x = element_text(size=10),
    #         axis.text.y = element_text(size=8),
    #         axis.title.y = element_text(size=10),
    #         # axis.title = element_blank(),
    #         legend.text = element_text(size=8),
    #         legend.title = element_text(size=10),
    #         legend.position = "bottom",
    #         plot.title = element_text(size=12),
    #         plot.subtitle = element_text(size=10)) +
    #   guides(color = guide_legend(override.aes = list(size=4)),
    #          shape = guide_legend(override.aes = list(size=4))) +
    #   coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
    # 
    # plt.fpr.z <- ggplot(stats.tmp[stats.tmp$hours == stats.tmp$cont_hrs,]) +
    #   geom_abline(col = 'red', linetype = 'dashed', lwd = 1) +
    #   geom_line(aes(x = p, y = fpr, col = hours), lwd = 1.2) +
    #   labs(title = "",
    #        subtitle = "for p < 0.2",
    #        x = "p-value cut-off",
    #        y = "False Positive Rate") +
    #   scale_x_continuous(breaks = seq(0,1,0.05),
    #                      minor_breaks = seq(0,1,0.01),
    #                      limits = c(0,1)) +
    #   scale_y_continuous(breaks = seq(0,1,0.05),
    #                      minor_breaks = seq(0,1,0.01),
    #                      limits = c(0,1)) +
    #   scale_color_discrete(name = "Hours",
    #                      breaks=as.character(c(2.9,4.02,4.89,6.14,6.88,7.85,8.96,9.97,11.04))) +
    #   theme_linedraw() +
    #   theme(axis.text.x = element_text(size=8),
    #         axis.title.x = element_text(size=10),
    #         axis.text.y = element_text(size=8),
    #         axis.title.y = element_blank(),
    #         # axis.title = element_blank(),
    #         legend.text = element_text(size=8),
    #         legend.title = element_text(size=10),
    #         legend.position = "bottom",
    #         plot.title = element_text(size=12),
    #         plot.subtitle = element_text(size=10)) +
    #   guides(color = guide_legend(override.aes = list(size=4)),
    #          shape = guide_legend(override.aes = list(size=4))) +
    #   coord_cartesian(xlim = c(0, 0.2), ylim = c(0, 0.2))
    # 
    # c.fpr <- ggpubr::ggarrange(plt.fpr, plt.fpr.z,
    #                    common.legend = TRUE, legend = "bottom")
    # c.pow <- ggpubr::ggarrange(plt.pow.fdr, plt.pow.p,
    #                    common.legend = T, legend = 'bottom')
    # edis <- ggpubr::ggarrange(c.fpr, c.pow,
    #                   nrow = 2, ncol = 1)
    # annotate_figure(edis,
    #                 top = text_grob(expt_name))
    # 
    # # save(edis, file = sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_POWDY_FIG.Rdata", expt_name))
    # 
    # ggsave(sprintf("%s%s_EDIS_%d_%d.jpg",out_path,expt_name,rep,a),
    #        height = 22, width = 20, units = "cm",
    #        dpi = 300)
  } 
}

# save(dat.cnt.all, file = sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_DAT_CNT.RData", expt_name))
# load(sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_DAT_CNT.RData", expt_name))

# spe.all$ref <- 1/4
# sen.all$ref <- 1/4

save(spe.all, file = sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_SPECIFICITY.RData", expt_name))
save(sen.all, file = sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_SENSITIVITY.RData", expt_name))

# dat.cnt.all$power[dat.cnt.all$hours < dat.cnt.all$cont_hrs] <- dat.cnt.all$Deleterious[dat.cnt.all$hours < dat.cnt.all$cont_hrs]/910 * 100
# dat.cnt.all$power[dat.cnt.all$hours > dat.cnt.all$cont_hrs] <- dat.cnt.all$Beneficial[dat.cnt.all$hours > dat.cnt.all$cont_hrs]/910 * 100
# dat.cnt.all$abs_cen <- abs(1-dat.cnt.all$cen) * 100
# 
# sens5_rep <- ggplot(dat.cnt.all[round(dat.cnt.all$abs_cen) <= 5,]) +
#   geom_boxplot(aes(x = rep, y = power, group = rep), fill = 'grey90') +
#   geom_smooth(aes(x = rep, y = power), method = 'loess', se = F, lwd = 1.2) +
#   labs(title = 'Sensitivity in detecting 5% fitness effects',
#        subtitle = 'Increases with number of replicates') +
#   scale_x_continuous(name = 'No. of Replicates',
#                      breaks = seq(0,16,2),
#                      minor_breaks = seq(0,16,1)) +
#   scale_y_continuous(name = 'Sensitivity',
#                      breaks = seq(0,100,10),
#                      minor_breaks = seq(0,105,5)) +
#   coord_cartesian(ylim = c(0,100)) +
#   theme_linedraw()
# # ggsave(sprintf("%s%s_SENS5.jpg",out_path,expt_name),sens5_rep,
# #        height = 10, width = 10, units = "cm",
# #        dpi = 300)
# 
# save(sens5_rep, file = sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_SENS5_FIG.Rdata", expt_name))


###### VIRTUAL PLATE COLONY SIZE DENSITY
# load(sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_FITNESS.RData", expt_name))
# 
# fit.all$colony[fit.all$orf_name == 'BF_control'] <- 'Reference'
# fit.all$colony[fit.all$orf_name != 'BF_control'] <- 'Query'
# 
# virP1 <- NULL
# i <- 1
# for (ref.hr in unique(fit.all$cont_hrs)) {
#   for (que.hr in unique(fit.all$hours[fit.all$cont_hrs == ref.hr])) {
#     fit.all$vir_plate[fit.all$cont_hrs == ref.hr & fit.all$hours == que.hr] <- i
#     virP1$vir_plate[i] <- i
#     virP1$ref_hrs[i] <- ref.hr
#     virP1$que_hrs[i] <- que.hr
#     virP1$ref_med[i] <- median(fit.all$average[fit.all$cont_hrs == ref.hr &
#                                                  fit.all$hours == que.hr &
#                                                  fit.all$colony == 'Reference'],
#                                na.rm = T)
#     virP1$que_med[i] <- median(fit.all$average[fit.all$cont_hrs == ref.hr &
#                                                  fit.all$hours == que.hr &
#                                                  fit.all$colony == 'Query'],
#                                na.rm = T)
#     i <- i + 1
#   }
# }
# virP1 <- data.frame(virP1)
# virP1$es <- virP1$que_med/virP1$ref_med
# 
# # save(fit.all, file = sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_VIR_PLATE1.Rdata", expt_name))
# # save(virP1, file = sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_VIR_PLATE1_ES.Rdata", expt_name))
# 
# vir_plates1 <- ggplot(fit.all) +
#   geom_line(aes(x = average, col = colony),
#             stat = 'density', trim = T, lwd = 1.2) +
#   labs(title = 'Virtual Plates',
#        subtitle = 'With Bimodal Colony Size Distribution',
#        x = 'Colony Size (pix. count)',
#        y = 'Density') +
#   scale_color_manual(name = 'Colony\nType',
#                      breaks = c("Reference","Query"),
#                      values = c("Reference" = "#3F51B5",
#                                 "Query" = "#FFC107")) +
#   theme_linedraw() +
#   facet_wrap(.~vir_plate)
# ggsave(sprintf("%s%s_VIR_PLATE1.jpg",out_path,expt_name),vir_plates1,
#        height = 15, width = 15,
#        dpi = 300)
# 
# # ##### VIRTUAL PLATE EFFECT SIZES
# load(sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_VIR_PLATE1_ES.Rdata", expt_name))
# hhours <- c(0,1.02,1.38,2.9,4.02,4.89,6.14,6.88,7.85,8.96,9.97,11.04)
# 
# es.heatmap <- ggplot(virP1) +
#   geom_tile(aes(x = sprintf('%0.1f',que_hrs), y = sprintf('%0.1f',ref_hrs), fill = round((es-1)*100,0)),
#             col = 'black') +
#   geom_text(aes(x = sprintf('%0.1f',que_hrs), y = sprintf('%0.1f',ref_hrs), label = round((es-1)*100,0)),
#             col = 'black') +
#   scale_x_discrete(limits = sprintf('%0.1f',hhours)) +
#   scale_y_discrete(limits = sprintf('%0.1f',hhours)) +
#   scale_fill_gradient2(name = 'Fitness\nEffect',
#                       low = "#303F9F", high = "#FFC107", mid = "white",
#                       trans = 'pseudo_log',
#                       midpoint = 0,
#                       breaks = c(-80,-20,-5,0,5,20,80,400),
#                       guide = F) +
#   labs(x = 'Query Colony Size Time (hour)',
#        y = 'Reference Colony Size Time (hour)') +
#   theme_linedraw()
# ggsave(sprintf("%s%s_VIR_PLATE1_ES.jpg",out_path,expt_name), es.heatmap,
#        height = 5, width = 5,
#        dpi = 300)
# 
# que.growth <- ggplot(virP1[virP1$ref_hrs == 0,],
#        aes(x = seq(1,12,1), y = que_med)) +
#   geom_line(lwd = 1.2) +
#   scale_x_discrete(limits = sprintf('%0.1f',hhours)) +
#   scale_y_continuous(breaks = seq(0,500,250)) +
#   labs(x = 'Time (hour)',
#        y = 'Query Colony Size\n(pixels)') +
#   coord_cartesian(ylim = c(0,500)) +
#   theme_linedraw() +
#   theme(axis.title = element_blank())
# ggsave(sprintf("%s%s_VIR_PLATE1_QUE.jpg",out_path,expt_name), que.growth,
#        height = 1, width = 5,
#        dpi = 300)
# 
# ref.growth <- ggplot(virP1[virP1$que_hrs == 0,],
#                      aes(x = seq(1,12,1), y = ref_med)) +
#   geom_line(lwd = 1.2) +
#   scale_x_discrete(limits = sprintf('%0.1f',hhours)) +
#   scale_y_continuous(breaks = seq(0,500,250)) +
#   labs(x = 'Time (hour)',
#        y = 'Reference Colony Size (pixels)') +
#   coord_cartesian(ylim = c(0,500)) +
#   theme_linedraw() +
#   theme(axis.title = element_blank())
# ggsave(sprintf("%s%s_VIR_PLATE1_REF.jpg",out_path,expt_name), ref.growth,
#        height = 1, width = 5,
#        dpi = 300)
# 
# # 
# # ##### VIRTUAL PLATE ILLUSTRATION
# load("/home/sbp29/R/Projects/adaptivefitness/output/4C4_FS_NONORM_FITNESS.RData")
# fit.all$colony[fit.all$orf_name == 'BF_control'] <- 'Reference'
# fit.all$colony[fit.all$orf_name != 'BF_control'] <- 'Query'
# 
# vir5 <- fit.all[fit.all$cont_hrs == 2.9 & fit.all$hours == 2.9 & fit.all$plate == 1,]
# vir5$kind <- '1. Time = 2.9 hr'
# vir11 <- fit.all[fit.all$cont_hrs == 11.04 & fit.all$hours == 11.04 & fit.all$plate == 1,]
# vir11$kind <- '2. Time = 11.04 hr'
# vir115 <- fit.all[fit.all$cont_hrs == 11.04 & fit.all$hours == 2.9 & fit.all$plate == 1,]
# vir115$kind <- '3. Virtual Plate'
# 
# vir.all <- rbind(vir5,vir11,vir115)
# 
# ggplot(vir.all) +
#   geom_point(aes(x = col, y = row, size = average, col = colony)) +
#   theme_linedraw() +
#   scale_color_manual(name = 'Colony\nType',
#                      breaks = c("Reference","Query"),
#                      values = c("Reference" = "#3F51B5",
#                                 "Query" = "#FFC107")) +
#   scale_size_continuous(range = c(0,4), guide = F) +
#   scale_x_continuous(breaks = seq(0,96,4)) +
#   scale_y_continuous(breaks = seq(0,64,4),trans = 'reverse') +
#   labs(x = 'Column', y = 'Row') +
#   coord_cartesian(xlim = c(10,33), ylim = c(10,25)) +
#   facet_wrap(.~kind, nrow = 1)
# ggsave(sprintf("%s%s_VIR_PLATE1_EXAMPLE.jpg",out_path,expt_name),
#        height = 3, width = 12,
#        dpi = 300)
