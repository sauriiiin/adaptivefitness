dat.dir <- "/home/sbp29/R/Projects/adaptivefitness/rawdata/4C4_FS_CC/"
expt_name <- '4C4_FS_CC'

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
  reps <- c(reps, as.numeric(s[4]))
  hours <- c(hours, as.numeric(s[5]))
  atmpt <- c(atmpt, as.numeric(s[6]))
}
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
      dat.stats <- dat.stats[dat.stats$hours != 0,]
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
      for (h in unique(dat.stats$hours)) {
        cont.mean <- mean(dat.fit$fitness[dat.fit$hours == h & dat.fit$orf_name == 'BF_control' & !is.na(dat.fit$fitness)])
        
        ref.pix.m <- mean(dat.fit$average[dat.fit$hours == h & dat.fit$orf_name == 'BF_control' & !is.na(dat.fit$fitness)])
        que.pix.m <- mean(dat.fit$average[dat.fit$hours == h & dat.fit$orf_name != 'BF_control' & !is.na(dat.fit$fitness)])
        dat.stats$cen[dat.stats$hours == h] <- que.pix.m/ref.pix.m
        
        dat.stats$effect_p[dat.stats$hours == h & dat.stats$p <= 0.05 & dat.stats$cs_mean > cont.mean & !is.na(dat.stats$cs_mean)] <- 'Beneficial'
        dat.stats$effect_p[dat.stats$hours == h & dat.stats$p <= 0.05 & dat.stats$cs_mean < cont.mean & !is.na(dat.stats$cs_mean)] <- 'Deleterious'
        dat.stats$effect_p[dat.stats$hours == h & is.na(dat.stats$effect_p) & !is.na(dat.stats$cs_mean)] <- 'Neutral'
      }
      stats.all <- rbind(stats.all,dat.stats)
      fit.all <- rbind(fit.all,dat.fit)
    }
  }
}

dat.cnt.all <- NULL
for (a in atmpt) {
  for (rep in unique(reps)) {
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
    
    stats.tmp$fp <- stats.tmp$p * 910
    for (cont_hr in unique(stats.tmp$cont_hrs)) {
      for (hr in unique(stats.tmp$hours[stats.tmp$cont_hrs == cont_hr])) {
        for (pp in unique(stats.tmp$p[stats.tmp$cont_hrs == cont_hr & stats.tmp$hours == hr])) {
          stats.tmp$fpr[stats.tmp$hours == ii & stats.tmp$cont_hrs == stats.tmp$hours & stats.tmp$p == pp] <-
            dim(stats.tmp[stats.tmp$hours == ii & stats.tmp$cont_hrs == stats.tmp$hours & stats.tmp$p <= pp,])[1]/
            dim(stats.tmp[stats.tmp$hours == ii & stats.tmp$cont_hrs == stats.tmp$hours,])[1]
        }
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
    
    # stats.tmp %>% count(cont_hrs, hours, cen, effect_p)
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
    
    save(dat.cnt2,
         file = sprintf('%sSPEDATA.RData',out_path))
    save(stats.tmp,
         file = sprintf('%sSENDATA.RData',out_path))
    
    t <- dat.cnt2$Beneficial[1]
    plt.pow.fdr <- ggplot(dat.cnt2) +
      geom_area(aes(x = cen, y = Neutral, fill = 'Neutral'), alpha = 1) +
      geom_area(aes(x = cen, y = Beneficial, fill = 'Beneficial'), alpha = 1) +
      geom_area(aes(x = cen, y = Deleterious, fill = 'Deleterious'), alpha = 1) +
      geom_vline(xintercept = seq(0,2,0.025), col = '#757575', lwd = 0.5, alpha =0.5) +
      geom_hline(yintercept = c(t * seq(0,1,0.05)), col = '#757575', lwd = 0.5, alpha =0.5) +
      labs(title = "Which effects are detected?",
           subtitle = sprintf('with %d technical replicate (FDR = 0.05)',rep),
           x = 'Effect Size',
           y = 'Sensitivity') +
      scale_y_continuous(breaks = c(t * seq(0,1,0.1)),
                         minor_breaks = c(t * seq(0,1,0.05)),
                         labels = c('0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')) +
      scale_x_continuous(breaks = seq(0,2,0.05),
                         minor_breaks = seq(0,2,0.025)) +
      scale_color_manual(name = 'Effects',
                         breaks = c('Beneficial','Neutral','Deleterious'),
                         values = c('Deleterious'='#3F51B5',
                                    'Neutral'='#212121',
                                    'Beneficial'='#FFC107')) +
      scale_fill_manual(name = 'Effects',
                        breaks = c('Beneficial','Neutral','Deleterious'),
                        values = c('Deleterious'='#3F51B5',
                                   'Neutral'='#212121',
                                   'Beneficial'='#FFC107')) +
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
    
    plt.pow.p <- ggplot(dat.cnt2) +
      geom_area(aes(x = cen, y = Neutral_p, fill = 'Neutral'), alpha = 1) +
      geom_area(aes(x = cen, y = Beneficial_p, fill = 'Beneficial'), alpha = 1) +
      geom_area(aes(x = cen, y = Deleterious_p, fill = 'Deleterious'), alpha = 1) +
      geom_vline(xintercept = seq(0,2,0.025), col = '#757575', lwd = 0.5, alpha =0.5) +
      geom_hline(yintercept = c(t * seq(0,1,0.05)), col = '#757575', lwd = 0.5, alpha =0.5) +
      labs(title = "",
           subtitle = sprintf('with %d technical replicate (p <= 0.05)',rep),
           x = 'Effect Size',
           y = 'Sensitivity') +
      scale_y_continuous(breaks = c(t * seq(0,1,0.1)),
                         minor_breaks = c(t * seq(0,1,0.05)),
                         labels = c('0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')) +
      scale_x_continuous(breaks = seq(0,2,0.05),
                         minor_breaks = seq(0,2,0.025)) +
      scale_color_manual(name = 'Effects',
                         breaks = c('Beneficial','Neutral','Deleterious'),
                         values = c('Deleterious'='#3F51B5',
                                    'Neutral'='#212121',
                                    'Beneficial'='#FFC107')) +
      scale_fill_manual(name = 'Effects',
                        breaks = c('Beneficial','Neutral','Deleterious'),
                        values = c('Deleterious'='#3F51B5',
                                   'Neutral'='#212121',
                                   'Beneficial'='#FFC107')) +
      theme_linedraw() +
      theme(axis.text.x = element_text(size=8),
            axis.title.x = element_text(size=10),
            axis.text.y = element_text(size=8),
            axis.title.y = element_blank(),
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
    
    stats.tmp$hours <- as.character(stats.tmp$hours)
    stats.tmp$cont_hrs <- as.character(stats.tmp$cont_hrs)
    
    plt.fpr <- ggplot(stats.tmp[stats.tmp$hours == stats.tmp$cont_hrs,]) +
      geom_abline(col = 'red', linetype = 'dashed', lwd = 1) +
      geom_line(aes(x = p, y = fpr, col = hours), lwd = 1.2) +
      geom_segment(aes(x=0,xend=0.2,y=0,yend=0), col = 'black', lwd = 0.8) +
      geom_segment(aes(x=0.2,xend=0.2,y=0,yend=0.2), col = 'black', lwd = 0.8) +
      geom_segment(aes(x=0.2,xend=0,y=0.2,yend=0.2), col = 'black', lwd = 0.8) +
      geom_segment(aes(x=0,xend=0,y=0.2,yend=0), col = 'black', lwd = 0.8) +
      labs(title = "Does FPR follow random expectation?",
           subtitle = "for p between 0 and 1",
           x = "p-value cut-off",
           y = "False Positive Rate") +
      scale_x_continuous(breaks = seq(0,1,0.2),
                         minor_breaks = seq(0,1,0.05),
                         limits = c(0,1)) +
      scale_y_continuous(breaks = seq(0,1,0.2),
                         minor_breaks = seq(0,1,0.05),
                         limits = c(0,1)) +
      scale_color_discrete(name = "Hours",
                           breaks=as.character(c(2.9,4.02,4.89,6.14,6.88,7.85,8.96,9.97,11.04))) +
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
      guides(color = guide_legend(override.aes = list(size=4)),
             shape = guide_legend(override.aes = list(size=4))) +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
    
    plt.fpr.z <- ggplot(stats.tmp[stats.tmp$hours == stats.tmp$cont_hrs,]) +
      geom_abline(col = 'red', linetype = 'dashed', lwd = 1) +
      geom_line(aes(x = p, y = fpr, col = hours), lwd = 1.2) +
      labs(title = "",
           subtitle = "for p < 0.2",
           x = "p-value cut-off",
           y = "False Positive Rate") +
      scale_x_continuous(breaks = seq(0,1,0.05),
                         minor_breaks = seq(0,1,0.01),
                         limits = c(0,1)) +
      scale_y_continuous(breaks = seq(0,1,0.05),
                         minor_breaks = seq(0,1,0.01),
                         limits = c(0,1)) +
      scale_color_discrete(name = "Hours",
                           breaks=as.character(c(2.9,4.02,4.89,6.14,6.88,7.85,8.96,9.97,11.04))) +
      theme_linedraw() +
      theme(axis.text.x = element_text(size=8),
            axis.title.x = element_text(size=10),
            axis.text.y = element_text(size=8),
            axis.title.y = element_blank(),
            # axis.title = element_blank(),
            legend.text = element_text(size=8),
            legend.title = element_text(size=10),
            legend.position = "bottom",
            plot.title = element_text(size=12),
            plot.subtitle = element_text(size=10)) +
      guides(color = guide_legend(override.aes = list(size=4)),
             shape = guide_legend(override.aes = list(size=4))) +
      coord_cartesian(xlim = c(0, 0.2), ylim = c(0, 0.2))
    
    c.fpr <- ggarrange(plt.fpr, plt.fpr.z,
                       common.legend = TRUE, legend = "bottom")
    c.pow <- ggarrange(plt.pow.fdr, plt.pow.p,
                       common.legend = T, legend = 'bottom')
    edis <- ggarrange(c.fpr, c.pow,
                      nrow = 2, ncol = 1)
    annotate_figure(edis,
                    top = text_grob(expt_name))
    
    # save(edis, file = sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_POWDY_FIG.Rdata", expt_name))
    
    # ggsave(sprintf("%s%s_EDIS_%d_%d.jpg",out_path,expt_name,rep,a),
    #        height = 22, width = 20, units = "cm",
    #        dpi = 300)
  } 
}