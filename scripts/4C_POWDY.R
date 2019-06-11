##### POWER DYNAMICS
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 05/23/2019

##### INITIALIZE
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyverse)
library(egg)
library(stringr)
out_path = 'figs/lid_paper/';
dat.dir <- "/home/sbp29/R/Projects/proto_plots/rawdata/4C3_GA1_LID/"
expt_name <- '4C3_GA1'
pvals = seq(0,1,0.005)

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

##### GET DATA
stats.files <- list.files(path = dat.dir,
                         pattern = "P.csv", recursive = TRUE)
fit.files <- list.files(path = dat.dir,
                        pattern = "S.csv", recursive = TRUE)
hours <- NULL
reps <- NULL
for (s in strsplit(stats.files,'_')) {
  # reps <- c(reps, as.numeric(s[4]))
  reps <- 8
  hours <- c(hours, as.numeric(s[3]))
}
reps <- unique(reps)
hours <- unique(hours)

stats.all <- NULL
fit.all <- NULL
fpr.all <- NULL
##### PUTTING IT TOGETHER
for (ii in 1:length(reps)) {
  rep <- reps[ii]
  for (i in 1:length(hours)) {
    hr <- hours[i]
    
    # dat.stats <- read.csv(paste0(dat.dir,
    #                              sprintf('%s_%d_%d_STATS_P.csv',expt_name,rep,hr)),
    #                       na.strings = "NaN")
    dat.stats <- read.csv(paste0(dat.dir,
                                 sprintf('%s_%d_STATS_P.csv',expt_name,hr)),
                          na.strings = "NaN")
    dat.stats <- dat.stats[dat.stats$hours != 0,]
    dat.stats$cont_hrs <- hr
    dat.stats$rep <- rep
    # dat.fit <- read.csv(paste0(dat.dir,
    #                            sprintf('%s_%d_%d_FITNESS.csv',expt_name,rep,hr)),
    #                     na.strings = "NaN")
    dat.fit <- read.csv(paste0(dat.dir,
                               sprintf('%s_%d_FITNESS.csv',expt_name,hr)),
                        na.strings = "NaN")
    dat.fit$cont_hrs <- hr
    dat.fit$rep <- rep
    dat.fit$se <- dat.fit$average - dat.fit$bg
    
    cont.mean <- mean(dat.fit$fitness[dat.fit$hours == hr & dat.fit$orf_name == 'BF_control' & !is.na(dat.fit$fitness)])
    dat.stats$es <-round(dat.stats$cs_mean/cont.mean,4)
    dat.stats$pthresh <- quantile(sort(dat.stats$p[dat.stats$hours == hr]),.05)
    for (ii in unique(dat.stats$hours)) {
      # dat.stats$cen[dat.stats$hours == ii] <- median(dat.stats$cs_mean[dat.stats$hours == ii])
      dat.stats$cen[dat.stats$hours == ii] <- mean(dat.stats$es[dat.stats$hours == ii])
    }

    dat.stats$effect[dat.stats$p <= dat.stats$pthresh & dat.stats$cs_mean > cont.mean] <- 'Beneficial'
    dat.stats$effect[dat.stats$p <= dat.stats$pthresh & dat.stats$cs_mean < cont.mean] <- 'Deleterious'
    dat.stats$effect[is.na(dat.stats$effect)] <- 'Neutral'
    stats.all <- rbind(stats.all,dat.stats)
    fit.all <- rbind(fit.all,dat.fit)
  }
}

##### THE STATS DATA ANALYSIS
# save(stats.all, file = "/home/sbp29/R/Projects/proto_plots/rawdata/4C3_GA1_STATS.RData")
# save(fit.all, file = "/home/sbp29/R/Projects/proto_plots/rawdata/4C3_GA1_FITNESS.RData")
stats.all$es <- round(stats.all$es,4)

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
  
  ##### POWER PLOT
  dat.pow <- data.frame(dat.pow)
  colnames(dat.pow) <- c("es","N","TP","FP","FPR","TN","FN","pow","sen","spe","acc")
  
  # pd.lid <- ggplot(dat.pow) +
  #   geom_line(aes(x = es, y = pow),col = '#D32F2F', linetype = 'twodash', lwd = 2) +
  #   scale_x_continuous(breaks = seq(-2,2,0.05),
  #                      minor_breaks = seq(-2,2,0.01)) +
  #   scale_y_continuous(breaks = seq(0,100,10),
  #                      minor_breaks = seq(0,100,5)) +
  #   labs(title = "Power Dynamics",
  #        subtitle = "Using Cubic Spline @ 5% FPR",
  #        x = "Effect Size",
  #        y = "Power") +
  #   theme_linedraw() +
  #   theme(axis.text.x = element_text(size=15),
  #         axis.title.x = element_text(size=20),
  #         axis.text.y = element_text(size=15),
  #         axis.title.y = element_text(size=20),
  #         legend.position = c(0.8,0.2),
  #         legend.background = element_blank(),
  #         legend.text = element_text(size=15),
  #         legend.title =  element_text(size=15),
  #         plot.title = element_text(size=25,hjust = 0),
  #         plot.subtitle = element_text(size=20,hjust = 0)) +
  #   coord_cartesian(xlim = c(0.8,1.2),
  #                   ylim = c(0,100))
  # ggsave(sprintf("%spower_cs%d.png",out_path,rep),
  #        pd.lid,
  #        width = 10,height = 10)
  
  ##### FPR PLOT
  # pdata = data.frame()
  # i = 1
  # for (hr in unique(stats.all$cont_hrs)) {
  #   temp = stats.all[stats.all$cont_hrs == hr & stats.all$hours == hr,]
  #   n = dim(temp)[1]
  #   for (p in pvals) {
  #     pdata = rbind(pdata, c(hr, p, sum(temp$p <= p)/n))
  #   }
  #   if (i == 1) {
  #     colnames(pdata) = c('hours','p','fpr')
  #     g = ggplot() + 
  #       geom_line(data = pdata[pdata$hours == hr,], aes(x = p, y = fpr, col = as.character(hours)), lwd = 1.5)
  #   } else {
  #     g = g + geom_line(data = pdata[pdata$hours == hr,], aes(x = p, y = fpr, col = as.character(hours)), lwd = 1.5)
  #   }
  #   i = i + 1
  # }
  # 
  # g + geom_line(data = pdata, aes(x = p, y = p, col = 'red'),
  #               linetype = 'dashed', lwd = 1.1, alpha = 0.7) +
  #   labs(title = "False positive rate",
  #        x = "p-value cut-off",
  #        y = "False Positive Rate") +
  #   scale_x_continuous(breaks = seq(0,1,0.1),
  #                      minor_breaks = seq(0,1,0.05),
  #                      limits = c(0,1)) +
  #   scale_y_continuous(breaks = seq(0,1,0.1),
  #                      minor_breaks = seq(0,1,0.05),
  #                      limits = c(0,1)) +
  #   scale_color_manual(name = "Hours",
  #                      breaks=c("8","10","14","16","18"),
  #                      values=c("8"="#D32F2F","10"="#536DFE","14"="#388E3C","16"="#795548","18"="#00BCD4","red"="red",
  #                               "0"="transparent","9"="transparent","11"="transparent","13"="transparent","17"="transparent")) +
  #   theme_linedraw() +
  #   theme(axis.text.x = element_text(size=15),
  #         axis.title.x = element_text(size=20),
  #         axis.text.y = element_text(size=15),
  #         axis.title.y = element_text(size=20),
  #         # axis.title = element_blank(),
  #         legend.text = element_text(size=15),
  #         legend.title = element_text(size=20),
  #         legend.position = "bottom",
  #         plot.title = element_text(size=25,hjust = 0),
  #         plot.subtitle = element_text(size=20,hjust = 0)) +
  #   coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
  # ggsave(sprintf("%sfpr%d.png",out_path,rep),
  #        g,
  #        width = 10,height = 10)
  
  
  ##### BOX PLOTS OF EFFECT DISTRIBUTION
  stats.tmp <- stats.all[stats.all$rep == rep &
                           stats.all$cont_hrs > 11 & stats.all$hours > 11,]
  
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
  
  # ef <- ggplot(dat.srt) +
  #   geom_bar(aes(pos, fill = effect),alpha = 0.9) +
  #   geom_vline(xintercept = c(11,24), linetype = 'dashed', col = 'red') +
  #   scale_fill_manual(name = 'Effects :',
  #                     breaks = c('Beneficial','Neutral','Deleterious'),
  #                     values = c('Deleterious'='#D32F2F',
  #                                'Neutral'='#303F9F',
  #                                'Beneficial'='#4CAF50')) +
  #   labs(x = "",
  #        y = str_wrap("Effect Distribution",10)) +
  #   scale_y_continuous(breaks = c(913 * seq(0,1,0.20)),
  #                      minor_breaks = c(913 * seq(0,1,0.05)),
  #                      labels = c('0%','20%','40%','60%','80%','100%')) +
  #   # scale_x_continuous(breaks = seq(unique(dat.srt$pos)[1],unique(dat.srt$pos)[length(unique(dat.srt$pos))],
  #   #                                 (unique(dat.srt$pos)[length(unique(dat.srt$pos))]- unique(dat.srt$pos)[1])/10),
  #   #                    labels = as.character(round(srt[seq(unique(dat.srt$pos)[1],unique(dat.srt$pos)[length(unique(dat.srt$pos))],
  #   #                                                  (unique(dat.srt$pos)[length(unique(dat.srt$pos))]- unique(dat.srt$pos)[1])/10)],2))) +
  #   scale_x_continuous(breaks = seq(0,50,4),
  #                      minor_breaks = seq(0,50,1)) +
  #   theme_linedraw() +
  #   theme(axis.text.x = element_blank(),
  #         axis.title.x = element_text(size=20),
  #         axis.ticks.x = element_blank(),
  #         axis.text.y = element_text(size=15),
  #         axis.title.y = element_text(size=20,
  #                                     angle = 90,
  #                                     vjust = 0.5),
  #         legend.position = "top",
  #         legend.background = element_rect(color = 'grey60'),
  #         legend.text = element_text(size=15),
  #         legend.title =  element_text(size=15),
  #         plot.title = element_text(size=25,hjust = 0.5),
  #         plot.subtitle = element_text(size=20,hjust = 0.5)) +
  #   coord_cartesian(xlim = c(1,36))
  # ggsave(sprintf("%seffect_dis%d.png",out_path,rep),
  #        width = 12,height = 8)
  
  # rf <- ggplot(dat.srt) +
  #   geom_line(aes(x = pos, y = cen, col = cen), lwd = 1.2) +
  #   geom_hline(yintercept = 0.9, linetype = 'dashed', col = 'red') +
  #   geom_text(aes(18, 0.8, label = "10% LOSS", hjust = 0.5), size = 4) +
  #   geom_hline(yintercept = 1.1, linetype = 'dashed', col = 'red') +
  #   geom_text(aes(18, 1.2, label = "10% GAIN", hjust = 0.5), size = 4) +
  #   geom_vline(xintercept = c(11,24), linetype = 'dashed', col = 'red') +
  #   labs(y = str_wrap("Mean Relative Fitness",10),
  #        x = "Virtual Plates") +
  #   # scale_x_continuous(breaks = seq(unique(dat.srt$pos)[1],unique(dat.srt$pos)[length(unique(dat.srt$pos))],
  #   #                                 (unique(dat.srt$pos)[length(unique(dat.srt$pos))]- unique(dat.srt$pos)[1])/10),
  #   #                    labels = as.character(round(srt[seq(unique(dat.srt$pos)[1],unique(dat.srt$pos)[length(unique(dat.srt$pos))],
  #   #                                                        (unique(dat.srt$pos)[length(unique(dat.srt$pos))]- unique(dat.srt$pos)[1])/10)],2))) +
  #   scale_x_continuous(breaks = seq(0,50,4),
  #                      minor_breaks = seq(0,50,1)) +
  #   scale_color_continuous(low = "#FFC107", high = "black",
  #                          guide = F) +
  #   theme_linedraw() +
  #   theme(axis.text.x = element_text(size=15),
  #         axis.title.x = element_text(size=20),
  #         axis.text.y = element_text(size=15),
  #         axis.title.y = element_text(size=20,
  #                                     angle = 90,
  #                                     vjust = 0.5),
  #         legend.position = "top",
  #         legend.background = element_rect(color = 'grey60'),
  #         legend.text = element_text(size=15),
  #         legend.title =  element_text(size=15),
  #         plot.title = element_text(size=25,hjust = 0),
  #         plot.subtitle = element_text(size=20,hjust = 0)) +
  #   coord_cartesian(xlim = c(1,36))
  # 
  # 
  # e.dis <- ggarrange(ef, rf,
  #                    nrow = 2,
  #                    heights = c(4,1))
  # ggsave(sprintf("%seffect_his%d.png",out_path,rep),
  #        e.dis,
  #        width = 10,height = 10)
  
  
  ##### LINE PLOT OF EFFECTS
  temp <- NULL
  dat.cnt <- data.frame()
  for (pos in unique(dat.srt$pos)) {
    temp$cen <- mean(dat.srt$cen[dat.srt$pos == pos])
    temp$Deleterious <- sum(dat.srt$pos == pos & dat.srt$effect == 'Deleterious')
    temp$Neutral <- sum(dat.srt$pos == pos & dat.srt$effect == 'Neutral')
    temp$Beneficial <- sum(dat.srt$pos == pos & dat.srt$effect == 'Beneficial')
    dat.cnt <- rbind(dat.cnt,temp)
  }
  
  dat.cnt2 <- NULL
  dat.cnt2$cen <- dat.cnt$cen
  dat.cnt2$Neutral <- dat.cnt$Neutral
  dat.cnt2$Deleterious <- dat.cnt2$Neutral + dat.cnt$Deleterious
  dat.cnt2$Beneficial <- dat.cnt2$Deleterious + dat.cnt$Beneficial
  dat.cnt2 <- data.frame(dat.cnt2)
  t <- dat.cnt2$Beneficial[1]
  
  plt.pow <- ggplot(dat.cnt2) +
    geom_area(aes(x = cen, y = Beneficial, fill = 'Beneficial'), alpha = 0.6) +
    # geom_point(aes(x = cen, y = Beneficial, fill = 'Beneficial'), shape = 21, col = 'black') +
    geom_area(aes(x = cen, y = Deleterious, fill = 'Deleterious'), alpha = 0.6) +
    # geom_point(aes(x = cen, y = Deleterious, fill = 'Deleterious'), shape = 21, col = 'black') +
    geom_area(aes(x = cen, y = Neutral, fill = 'Neutral'), alpha = 0.6) +
    # geom_point(aes(x = cen, y = Neutral, fill = 'Neutral'), shape = 21, col = 'black') +
    labs(title = "",
         subtitle = sprintf('With %d Technical Replicate',rep),
         x = 'Mean Relative Fitness',
         y = 'Effect Distribution') +
    scale_y_continuous(breaks = c(t * seq(0,1,0.1)),
                       minor_breaks = c(t * seq(0,1,0.05)),
                       labels = c('0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')) +
    scale_x_continuous(breaks = seq(0,2,0.05),
                       minor_breaks = seq(0,2,0.025)) +
    scale_color_manual(name = 'Effects',
                       breaks = c('Beneficial','Neutral','Deleterious'),
                       values = c('Deleterious'='#D32F2F',
                                  'Neutral'='#303F9F',
                                  'Beneficial'='#4CAF50')) +
    scale_fill_manual(name = 'Effects',
                      breaks = c('Beneficial','Neutral','Deleterious'),
                      values = c('Deleterious'='#D32F2F',
                                 'Neutral'='#303F9F',
                                 'Beneficial'='#4CAF50')) +
    theme_linedraw() +
    theme(axis.text.x = element_text(size=8),
          axis.title.x = element_text(size=10),
          axis.text.y = element_text(size=8),
          axis.title.y = element_text(size=10),
          # axis.title = element_blank(),
          legend.text = element_text(size=8),
          legend.title = element_text(size=10),
          legend.position = "bottom",
          plot.title = element_text(size=12),
          plot.subtitle = element_text(size=10)) +
    guides(color = guide_legend(override.aes = list(size=2)),
           shape = guide_legend(override.aes = list(size=2))) +
    coord_cartesian(xlim = c(0.8,1.2),
                    ylim = c(0,t))
  # ggsave(sprintf("%s%s_EffectDis%d.png",out_path,expt_name,rep),
  #        width = 10,height = 10) 
  
  stats.tmp$hours <- as.character(stats.tmp$hours)
  stats.tmp$cont_hrs <- as.character(stats.tmp$cont_hrs)
  
  plt.fpr <- ggplot(stats.tmp[stats.tmp$hours == stats.tmp$cont_hrs,]) +
    geom_abline(col = 'red', linetype = 'dashed', lwd = 1) +
    geom_line(aes(x = p, y = fpr, col = hours), lwd = 1.2) +
    labs(title = expt_name,
         subtitle = "False positive rate",
         x = "p-value cut-off",
         y = "False Positive Rate") +
    scale_x_continuous(breaks = seq(0,1,0.05),
                       minor_breaks = seq(0,1,0.01),
                       limits = c(0,1)) +
    scale_y_continuous(breaks = seq(0,1,0.05),
                       minor_breaks = seq(0,1,0.01),
                       limits = c(0,1)) +
    scale_color_manual(name = "Hours",
                       breaks=c("13","14","16","17","18"),
                       values=c("13"="#D32F2F","14"="#536DFE","16"="#388E3C","17"="#795548","18"="#00BCD4",
                                "0"="transparent","8"="transparent","9"="transparent","10"="transparent","11"="transparent")) +
    theme_linedraw() +
    theme(axis.text.x = element_text(size=8),
          axis.title.x = element_text(size=10),
          axis.text.y = element_text(size=8),
          axis.title.y = element_text(size=10),
          # axis.title = element_blank(),
          legend.text = element_text(size=8),
          legend.title = element_text(size=10),
          legend.position = "bottom",
          plot.title = element_text(size=12),
          plot.subtitle = element_text(size=10)) +
    guides(color = guide_legend(override.aes = list(size=6)),
           shape = guide_legend(override.aes = list(size=6))) +
    coord_cartesian(xlim = c(0, 0.2), ylim = c(0, 0.2))
  
  ggarrange(plt.fpr, plt.pow)
  ggsave(sprintf("%s%s_EDIS_%d.jpg",out_path,expt_name,rep),
         height = 12, width = 20, units = "cm",
         dpi = 300)
}

###### THE FITNESS DATA ANALYSIS
temp.fit <- fit.all[fit.all$cont_hrs == 18 & fit.all$hours ==18,]
rmse <- sqrt(mean((abs(temp.fit$se))^2, na.rm = T))
mean.cs <- mean(temp.fit$average, na.rm = T)
rmse/mean.cs * 100
# ss.tot <- sum((temp.fit$average - mean.cs)^2,na.rm = T)
# ss.res <- sum(temp.fit$se^2,na.rm = T)
# 
# R2.lid <- 1 - ss.res/ss.tot

temp.fit$rand_se <- temp.fit$average - temp.fit$average[sample(1:length(temp.fit$average))]
rand_rmse <- sqrt(mean((abs(temp.fit$rand_se))^2, na.rm = T))
rand_rmse/mean.cs * 100
# R2.rnd <- 1 - sum(temp.fit$rand_se^2,na.rm = T)/ss.tot




