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
dat.dir <- "/home/sbp29/R/Projects/adaptivefitness/rawdata/4C4_FS_UP1/"
expt_name <- '4C4_FS_UP1'
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
# atmpt <- 1:4

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
      # dat.stats <- read.csv(paste0(dat.dir,
      #                              sprintf('%s_%d_STATS_P.csv',expt_name,hr)),
      #                       na.strings = "NaN")
      # dat.stats <- dat.stats[!dat.stats$orf_name %in% sprintf('MASK%d',1:1000),]
      dat.stats <- dat.stats[dat.stats$hours != 0,]
      dat.stats$cont_hrs <- hhours[i]
      dat.stats$rep <- rep
      dat.fit <- read.csv(paste0(dat.dir,
                                 sprintf('%s_%d_%d_%d_FITNESS.csv',expt_name,rep,hr,iii)),
                          na.strings = "NaN")
      # dat.fit <- read.csv(paste0(dat.dir,
      #                            sprintf('%s_%d_FITNESS.csv',expt_name,hr)),
      #                     na.strings = "NaN")
      # dat.fit <- dat.fit[!dat.fit$orf_name %in% sprintf('MASK%d',1:1000),]
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

stats.all <- stats.all[!stats.all$orf_name %in% sprintf('MASK%d',1:1000),]
fit.all <- fit.all[!fit.all$orf_name %in% sprintf('MASK%d',1:1000),]
# save(stats.all, file = sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_STATS.RData", expt_name))
# save(fit.all, file = sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_FITNESS.RData", expt_name))

load(sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_STATS.RData", expt_name))
load(sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_FITNESS.RData", expt_name))

##### THE STATS DATA ANALYSIS
# stats.all$es <- round(stats.all$es,4)
dat.cnt.all <- NULL
for (a in atmpt) {
  for (rep in unique(reps)) {
    effect_size <- sort(unique(round(stats.all$es[stats.all$rep == rep],2)))
    dat.pow <- NULL
    
    ##### POWER CALCULATIONS
    for (es in effect_size) {
      N <- sum(stats.all$rep == rep & 
                 stats.all$es > es-5e-03 & stats.all$es < es+5e-03)
      if (N > 0) {
        TP <- sum(stats.all$rep == rep &
                    stats.all$es > es-5e-03 & stats.all$es < es+5e-03 &
                    stats.all$p <= stats.all$pthresh & stats.all$hours != stats.all$cont_hrs)
        FP <- sum(stats.all$rep == rep &
                    stats.all$es > es-5e-03 & stats.all$es < es+5e-03 &
                    stats.all$p <= stats.all$pthresh & stats.all$hours == stats.all$cont_hrs)
        FPR <- FP/sum(stats.all$rep == rep &
                        stats.all$es > es-5e-03 & stats.all$es < es+5e-03 &
                        stats.all$hours == stats.all$cont_hrs) * 100
        TN <- sum(stats.all$rep == rep &
                    stats.all$es > es-5e-03 & stats.all$es < es+5e-03 &
                    stats.all$p > stats.all$pthresh & stats.all$hours == stats.all$cont_hrs)
        FN <- sum(stats.all$rep == rep &
                    stats.all$es > es-5e-03 & stats.all$es < es+5e-03 &
                    stats.all$p > stats.all$pthresh & stats.all$hours != stats.all$cont_hrs)
        POW <- TP/N * 100
        SEN <- TP/sum(stats.all$rep == rep &
                        stats.all$es > es-5e-03 & stats.all$es < es+5e-03 &
                        stats.all$hours != stats.all$cont_hrs) * 100
        SPE <- TN/sum(stats.all$rep == rep &
                        stats.all$es > es-5e-03 & stats.all$es < es+5e-03 &
                        stats.all$hours == stats.all$cont_hrs) * 100
        ACC <- (TP + TN)/N * 100
        dat.pow <- rbind(dat.pow, c(es,N,TP,FP,FPR,TN,FN,POW,SEN,SPE,ACC))
      }
    }
    dat.pow <- data.frame(dat.pow)
    colnames(dat.pow) <- c("es","N","TP","FP","FPR","TN","FN","pow","sen","spe","acc")
    
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
    t <- dat.cnt2$Beneficial[1]
    
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
    #                      values = c('Deleterious'='#F44336',
    #                                 'Neutral'='#303F9F',
    #                                 'Beneficial'='#4CAF50')) +
    #   scale_fill_manual(name = 'Effects',
    #                     breaks = c('Beneficial','Neutral','Deleterious'),
    #                     values = c('Deleterious'='#F44336',
    #                                'Neutral'='#303F9F',
    #                                'Beneficial'='#4CAF50')) +
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
    #                      values = c('Deleterious'='#F44336',
    #                                 'Neutral'='#303F9F',
    #                                 'Beneficial'='#4CAF50')) +
    #   scale_fill_manual(name = 'Effects',
    #                     breaks = c('Beneficial','Neutral','Deleterious'),
    #                     values = c('Deleterious'='#F44336',
    #                                'Neutral'='#303F9F',
    #                                'Beneficial'='#4CAF50')) +
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
    
    stats.tmp$hours <- as.character(stats.tmp$hours)
    stats.tmp$cont_hrs <- as.character(stats.tmp$cont_hrs)
    
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
    #   # scale_color_manual(name = "Hours",
    #   #                    breaks=c("13","14","16","17","18"),
    #   #                    values=c("13"="#D32F2F","14"="#536DFE","16"="#388E3C","17"="#795548","18"="#00BCD4",
    #   #                             "0"="transparent","8"="transparent","9"="transparent","10"="transparent","11"="transparent")) +
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
    #   # scale_color_manual(name = "Hours",
    #   #                    breaks=c("13","14","16","17","18"),
    #   #                    values=c("13"="#D32F2F","14"="#536DFE","16"="#388E3C","17"="#795548","18"="#00BCD4",
    #   #                             "0"="transparent","8"="transparent","9"="transparent","10"="transparent","11"="transparent")) +
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
    # c.fpr <- ggarrange(plt.fpr, plt.fpr.z,
    #                    common.legend = TRUE, legend = "bottom")
    # c.pow <- ggarrange(plt.pow.fdr, plt.pow.p,
    #                    common.legend = T, legend = 'bottom')
    # edis <- ggarrange(c.fpr, c.pow,
    #                   nrow = 2, ncol = 1)
    # annotate_figure(edis,
    #                 top = text_grob(expt_name))
    
    # save(edis, file = sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_POWDY_FIG.Rdata", expt_name))
    
    # ggsave(sprintf("%s%s_EDIS_%d_%d.jpg",out_path,expt_name,rep,a),
    #        height = 22, width = 20, units = "cm",
    #        dpi = 300)
  } 
}

# save(dat.cnt.all, file = sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_DAT_CNT.RData", expt_name))
load(sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_DAT_CNT.RData", expt_name))

dat.cnt.all$power[dat.cnt.all$hours < dat.cnt.all$cont_hrs] <- dat.cnt.all$Deleterious_p[dat.cnt.all$hours < dat.cnt.all$cont_hrs]/910 * 100
dat.cnt.all$power[dat.cnt.all$hours > dat.cnt.all$cont_hrs] <- dat.cnt.all$Beneficial_p[dat.cnt.all$hours > dat.cnt.all$cont_hrs]/910 * 100

dat.cnt.all$abs_cen <- abs(1-dat.cnt.all$cen) * 100

sens5_rep <- ggplot(dat.cnt.all[round(dat.cnt.all$abs_cen) == 5,]) +
  geom_boxplot(aes(x = rep, y = power, group = rep), fill = 'grey90') +
  geom_smooth(aes(x = rep, y = power), method = 'loess', se = F, lwd = 1.2) +
  labs(title = 'Sensitivity in detecting 5% fitness effects',
       subtitle = 'Increases with number of replicates') +
  scale_x_continuous(name = 'No. of Replicates',
                     breaks = seq(0,16,2),
                     minor_breaks = seq(0,16,1)) +
  scale_y_continuous(name = 'Sensitivity',
                     breaks = seq(0,100,10),
                     minor_breaks = seq(0,105,5)) +
  coord_cartesian(ylim = c(0,100)) +
  theme_linedraw()
ggsave(sprintf("%s%s_SENS5.jpg",out_path,expt_name),sens5_rep,
       height = 10, width = 10, units = "cm",
       dpi = 300)

save(sens5_rep, file = sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_SENS5_FIG.Rdata", expt_name))

# dat.cnt.all[4:9] <- dat.cnt.all[4:9]/912 * 100
# dat.cnt.all$Sensitivity <- dat.cnt.all$Beneficial + dat.cnt.all$Deleterious
# dat.cnt.all$Sensitivity_p <- dat.cnt.all$Beneficial_p + dat.cnt.all$Deleterious_p
# 
# # sen <- data.frame(x = NULL, y = NULL)
# for (rep in unique(dat.cnt.all$rep)) {
#   temp <- approx(dat.cnt.all$cen[dat.cnt.all$rep == rep], dat.cnt.all$Sensitivity_p[dat.cnt.all$rep == rep],
#                  xout = seq(0.7,1.3,0.05),
#                  method = 'linear')
#   sen <- rbind(sen, cbind(data.frame(temp), rep, data.frame(ref = 25*0.25)))
# }
# 
# ##### NON TECH REP ANALYSIS
# expt_name = '4C3_GA1_TRBLBR'
# out_path = "figs/lid_paper/EDIS/"
# norep.dat <- read.csv("/home/sbp29/R/Projects/proto_plots/rawdata/4C3_GA1_TRBLBR_POWANA.csv",
#                       na.strings = "NaN")
# colnames(norep.dat) <- c("rep","cont_hrs","hours","cen","pow",
#                          "Neutral","Deleterious","Beneficial",
#                          "Neutral_p","Deleterious_p","Beneficial_p")
# 
# 
# for (rep in unique(norep.dat$rep)) {
#   temp <- norep.dat[norep.dat$rep == rep,]
#   t <- temp$Beneficial[1]
#   
#   plt.pow.fdr <- ggplot(temp) +
#     geom_area(aes(x = cen, y = Beneficial, fill = 'Beneficial'), alpha = 0.6) +
#     # geom_point(aes(x = cen, y = Beneficial, fill = 'Beneficial'), shape = 21, col = 'black') +
#     geom_area(aes(x = cen, y = Deleterious, fill = 'Deleterious'), alpha = 0.6) +
#     # geom_point(aes(x = cen, y = Deleterious, fill = 'Deleterious'), shape = 21, col = 'black') +
#     geom_area(aes(x = cen, y = Neutral, fill = 'Neutral'), alpha = 0.6) +
#     # geom_point(aes(x = cen, y = Neutral, fill = 'Neutral'), shape = 21, col = 'black') +
#     labs(title = "Which effects are detected?",
#          subtitle = sprintf('with %d replicate (FDR = 0.05)',rep),
#          x = 'Mean Relative Fitness',
#          y = 'Effect Distribution') +
#     scale_y_continuous(breaks = c(t * seq(0,1,0.1)),
#                        minor_breaks = c(t * seq(0,1,0.05)),
#                        labels = c('0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')) +
#     scale_x_continuous(breaks = seq(0,2,0.05),
#                        minor_breaks = seq(0,2,0.025)) +
#     scale_color_manual(name = 'Effects',
#                        breaks = c('Beneficial','Neutral','Deleterious'),
#                        values = c('Deleterious'='#D32F2F',
#                                   'Neutral'='#303F9F',
#                                   'Beneficial'='#4CAF50')) +
#     scale_fill_manual(name = 'Effects',
#                       breaks = c('Beneficial','Neutral','Deleterious'),
#                       values = c('Deleterious'='#D32F2F',
#                                  'Neutral'='#303F9F',
#                                  'Beneficial'='#4CAF50')) +
#     theme_linedraw() +
#     theme(axis.text.x = element_text(size=8),
#           axis.title.x = element_text(size=10),
#           axis.text.y = element_text(size=8),
#           axis.title.y = element_text(size=10),
#           # axis.title = element_blank(),
#           legend.text = element_text(size=8),
#           legend.title = element_text(size=10),
#           legend.position = "bottom",
#           plot.title = element_text(size=12),
#           plot.subtitle = element_text(size=10)) +
#     guides(color = guide_legend(override.aes = list(size=2)),
#            shape = guide_legend(override.aes = list(size=2))) +
#     coord_cartesian(xlim = c(0.8,1.2),
#                     ylim = c(0,t))
#   # ggsave(sprintf("%s%s_EffectDis%d.png",out_path,expt_name,rep),
#   #        width = 10,height = 10) 
#   
#   plt.pow.p <- ggplot(temp) +
#     geom_area(aes(x = cen, y = Beneficial_p, fill = 'Beneficial'), alpha = 0.6) +
#     # geom_point(aes(x = cen, y = Beneficial, fill = 'Beneficial'), shape = 21, col = 'black') +
#     geom_area(aes(x = cen, y = Deleterious_p, fill = 'Deleterious'), alpha = 0.6) +
#     # geom_point(aes(x = cen, y = Deleterious, fill = 'Deleterious'), shape = 21, col = 'black') +
#     geom_area(aes(x = cen, y = Neutral_p, fill = 'Neutral'), alpha = 0.6) +
#     # geom_point(aes(x = cen, y = Neutral, fill = 'Neutral'), shape = 21, col = 'black') +
#     labs(title = "",
#          subtitle = sprintf('with %d replicate (p = 0.05)',rep),
#          x = 'Mean Relative Fitness',
#          y = 'Effect Distribution') +
#     scale_y_continuous(breaks = c(t * seq(0,1,0.1)),
#                        minor_breaks = c(t * seq(0,1,0.05)),
#                        labels = c('0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')) +
#     scale_x_continuous(breaks = seq(0,2,0.05),
#                        minor_breaks = seq(0,2,0.025)) +
#     scale_color_manual(name = 'Effects',
#                        breaks = c('Beneficial','Neutral','Deleterious'),
#                        values = c('Deleterious'='#D32F2F',
#                                   'Neutral'='#303F9F',
#                                   'Beneficial'='#4CAF50')) +
#     scale_fill_manual(name = 'Effects',
#                       breaks = c('Beneficial','Neutral','Deleterious'),
#                       values = c('Deleterious'='#D32F2F',
#                                  'Neutral'='#303F9F',
#                                  'Beneficial'='#4CAF50')) +
#     theme_linedraw() +
#     theme(axis.text.x = element_text(size=8),
#           axis.title.x = element_text(size=10),
#           axis.text.y = element_text(size=8),
#           axis.title.y = element_blank(),
#           # axis.title = element_blank(),
#           legend.text = element_text(size=8),
#           legend.title = element_text(size=10),
#           legend.position = "bottom",
#           plot.title = element_text(size=12),
#           plot.subtitle = element_text(size=10)) +
#     guides(color = guide_legend(override.aes = list(size=2)),
#            shape = guide_legend(override.aes = list(size=2))) +
#     coord_cartesian(xlim = c(0.8,1.2),
#                     ylim = c(0,t))
#   
#   c.pow <- ggarrange(plt.pow.fdr, plt.pow.p,
#                      common.legend = T, legend = 'bottom')
#   annotate_figure(c.pow,
#                   top = text_grob(expt_name))
#   
#   ggsave(sprintf("%s%s_EDIS2_%d.jpg",out_path,expt_name,rep),
#          height = 11, width = 20, units = "cm",
#          dpi = 300)
# }
# 
# ##### BEAN VS LID
# ## fit.all2 is LID data
# # r = 18
# # q = 18
# 
# # es <- mean(fit.all$fitness[fit.all$cont_hrs == r & fit.all$hours == q &
# #                              fit.all$orf_name != "BF_control" & fit.all$x6144plate_1 == 1],na.rm = T)/
# #   mean(fit.all$fitness[fit.all$cont_hrs == r & fit.all$hours == q &
# #                          fit.all$orf_name == "BF_control" & fit.all$x6144plate_1 == 1],na.rm = T)
# # 
# # bean <- ggplot() +
# #   geom_line(data = fit.all[fit.all$cont_hrs == r & fit.all$hours == q & 
# #                              fit.all$orf_name == "BF_control" & fit.all$x6144plate_1 == 1,],
# #             aes(x = fitness, col = "Reference"), stat = "density", lwd = 2) +
# #   geom_line(data = fit.all[fit.all$cont_hrs == r & fit.all$hours == q &
# #                              fit.all$orf_name != "BF_control" & fit.all$x6144plate_1 == 1,],
# #             aes(x = fitness, col = "Query"), stat = "density", lwd = 2) +
# #   scale_color_discrete(name = "Strain") +
# #   labs(title = sprintf("BEAN (ES = %.3f)",es),
# #        x = "Fitness", y = "Density") +
# #   theme_linedraw() +
# #   coord_cartesian(xlim = c(0.5,1.5))
# # 
# # es2 <- mean(fit.all2$fitness[fit.all2$cont_hrs == r & fit.all2$hours == q &
# #                                fit.all2$orf_name != "BF_control" & fit.all2$x6144plate_1 == 1],na.rm = T)/
# #   mean(fit.all2$fitness[fit.all2$cont_hrs == r & fit.all2$hours == q &
# #                           fit.all2$orf_name == "BF_control" & fit.all2$x6144plate_1 == 1],na.rm = T)
# # 
# # lid <- ggplot() +
# #   geom_line(data = fit.all2[fit.all2$cont_hrs == r & fit.all2$hours == q &
# #                               fit.all2$orf_name == "BF_control" & fit.all2$x6144plate_1 == 1,],
# #             aes(x = fitness, col = "Reference"), stat = "density", lwd = 2) +
# #   geom_line(data = fit.all2[fit.all2$cont_hrs == r & fit.all2$hours == q &
# #                               fit.all2$orf_name != "BF_control" & fit.all2$x6144plate_1 == 1,],
# #             aes(x = fitness, col = "Query"), stat = "density", lwd = 2) +
# #   scale_color_discrete(name = "Strain") +
# #   labs(title = sprintf("LID (ES = %.3f)",es2),
# #        x = "Fitness", y = "Density") +
# #   theme_linedraw() +
# #   theme(axis.title.y = element_blank()) +
# #   coord_cartesian(xlim = c(0.5,1.5))
# # 
# # fd <- ggarrange(bean, lid,
# #                 nrow = 1, ncol = 2,
# #                 common.legend = T, legend = "right")
# # annotate_figure(fd,
# #                 top = text_grob(sprintf("Fitness Distribution\nt(R) = %d | t(Q) = %d", r, q)))
# # ggsave(sprintf("figs/%s_FDIS_COMP_%d_%d.jpg",expt_name,r,q),
# #        height = 10, width = 20, units = "cm",
# #        dpi = 300)
# 
# ###
# cfit.dat <- read.csv("/home/sbp29/R/Projects/proto_plots/rawdata/4C3_GA3_CONTFITALL_8.csv",
#                      na.strings = "NaN")
# colnames(cfit.dat) <- c("cont_hrs","hours","fitness")
# 
# q = 8
# r = 18
# 
# # for (q in unique(stats.all$hours)) {
# #   for (r in unique(stats.all$cont_hrs)) {
#     ggplot(cfit.dat[cfit.dat$hours == q & cfit.dat$cont_hrs == r,]) +
#       geom_line(aes(x = fitness, col = "Reference"), stat = "density", lwd = 2) +
#       geom_line(data = stats.all[stats.all$hours == q & stats.all$cont_hrs == r,],
#                 aes(x = cs_mean, col = "Query"), stat = "density", lwd = 2) +
#       scale_color_discrete(name = "Strain") +
#       labs(title = sprintf('BEAN | t(R) = %d | t(Q) = %d', r, q),
#            subtitle = sprintf('D = %s | N = %s | B = %s',
#                               sprintf("%s",dat.cnt[dat.cnt$hours == q & dat.cnt$cont_hrs == r,][4:6])[1],
#                               sprintf("%s",dat.cnt[dat.cnt$hours == q & dat.cnt$cont_hrs == r,][4:6])[2],
#                               sprintf("%s",dat.cnt[dat.cnt$hours == q & dat.cnt$cont_hrs == r,][4:6])[3]),
#            x = "Fitness", y = "Density") +
#       theme_linedraw() +
#       coord_cartesian(xlim = c(0.8,1.2))
#     ggsave(sprintf("figs/%s_FDIS_%d_%d.jpg",expt_name,r,q),
#            height = 20, width = 20, units = "cm",
#            dpi = 300)
# #   }
# # }
