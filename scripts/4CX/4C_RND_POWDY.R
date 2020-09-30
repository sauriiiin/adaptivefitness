##### INITIALIZE
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyverse)
# library(egg)
library(ggpubr)
library(stringr)
out_path = 'figs/paper/';
dat.dir <- "/home/sbp29/R/Projects/adaptivefitness/rawdata/4C4_RND2_RR/"
pvals = seq(0,1,0.005)

source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")
out_path = 'figs/paper/';

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 7
txt <- 5
lbls <- 9

# ##### GET DATA
stats.files <- list.files(path = dat.dir,
                          pattern = "P.csv", recursive = TRUE)
fit.files <- list.files(path = dat.dir,
                        pattern = "S.csv", recursive = TRUE)
hours <- NULL
reps <- NULL
refs <- NULL
atmpt <- NULL
expt_name <- NULL
for (s in strsplit(stats.files,'_')) {
  expt_name <- c(expt_name, str_flatten(s[1:3],collapse = '_'))
  # refs <- c(refs, as.numeric(s[4]))
  reps <- c(reps, as.numeric(s[4]))
  # reps <- 8
  hours <- c(hours, as.numeric(s[5]))
  atmpt <- c(atmpt, as.numeric(s[6]))
}
# refs <- unique(refs)
expt_name <- unique(expt_name)
reps <- sort(unique(reps))
hours <- sort(unique(hours))
atmpt <- sort(unique(atmpt))

stats.all <- NULL
fit.all <- NULL
fpr.all <- NULL
##### PUTTING IT TOGETHER
for (exp in expt_name) {
  for (iii in 1:length(atmpt)) {
    for (ii in 1:length(reps)) {
      rep <- reps[ii]
      for (i in 1:length(hours)) {
        hr <- hours[i]

        dat.stats <- read.csv(paste0(dat.dir,
                                     sprintf('%s_%d_%d_%d_STATS_P.csv',exp,rep,hr,iii)),
                              na.strings = "NaN")
        dat.stats$rep <- rep
        dat.stats$attempt <- iii
        dat.stats$replicates <- rep
        dat.stats$expt_name <- exp
        dat.stats$pthresh <- quantile(sort(dat.stats$p),.05)

        dat.fit <- read.csv(paste0(dat.dir,
                                   sprintf('%s_%d_%d_%d_FITNESS.csv',exp,rep,hr,iii)),
                            na.strings = "NaN")
        dat.fit$attempt <- iii
        dat.fit$rep <- rep
        dat.fit$se <- dat.fit$average - dat.fit$bg
        dat.fit$expt_name <- exp


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
}

expt_name <- '4C4_RND2_RR'
# save(stats.all, file = sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_STATS.RData", expt_name))
# save(fit.all, file = sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_FITNESS.RData", expt_name))
load(sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_STATS.RData", expt_name))
load(sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_FITNESS.RData", expt_name))
# dbWriteTable(conn, '4C4_RR_RND2_STATS_ALL', stats.all, overwrite = T)

##### CORRECTING RND HRS in RND2
rnd_hrs <- dbGetQuery(conn, 'select orf_name, hours, max(rnd_hrs) rnd_hrs
                      from 4C4_FS_RND2_6144_DATA
                      where orf_name is not NULL
                      group by hours, orf_name')

for (hr in unique(stats.all$hours)) {
  for (orf in unique(stats.all$orf_name[stats.all$hours == hr & !stats.all$orf_name %in% sprintf('MASK%d',1:1000)])) {
    stats.all$rnd_hrs[stats.all$hours == hr & stats.all$orf_name == orf] <- 
      rnd_hrs$rnd_hrs[rnd_hrs$hours == hr & rnd_hrs$orf_name == orf]
  }
}


##### REF FITNESS DISTRIBUTION
# ref.all <- NULL
# i = 0
# for (exp in unique(fit.all$expt_name)) {
#   for (rep in unique(fit.all$rep[fit.all$expt_name == exp])) {
#     for (a in unique(fit.all$attempt[fit.all$expt_name == exp & fit.all$rep == rep])) {
#       for (hr in sort(unique(fit.all$hours[fit.all$expt_name == exp & fit.all$rep == rep & fit.all$attempt == a]))) {
#         temp <- fit.all[fit.all$expt_name == exp & fit.all$rep == rep & fit.all$attempt == a & fit.all$hours == hr & fit.all$orf_name == 'BF_control',]
#         for (pos in 1153:1536) {
#           i <- i + 1
#           ref.all$pix_mean[i] <-
#             mean(temp$average[temp$pos %in% c(c(110000,120000,130000,140000,210000,220000,230000,240000,
#                                                 310000,320000,330000,340000,410000,420000,430000,440000) + pos)], na.rm = T)
#           ref.all$cs_mean[i] <-
#             mean(temp$fitness[temp$pos %in% c(c(110000,120000,130000,140000,210000,220000,230000,240000,
#                                                 310000,320000,330000,340000,410000,420000,430000,440000) + pos)], na.rm = T)
#           ref.all$hours[i] <- hr
#           ref.all$expt_name[i] <- exp
#           ref.all$rep[i] <- rep
#           ref.all$attempt[i] <- a
#         }
#       }
#     }
#   }
# }
# ref.all <- data.frame(ref.all)
# save(ref.all, file = sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_REFFIT.RData", expt_name))
# load(sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_REFFIT.RData", expt_name))

##### FITNESS EFFECTS
load("/home/sbp29/R/Projects/adaptivefitness/output/4C4_FS_CC_VIR_PLATE1_ES.Rdata")
stats.all <- stats.all[!stats.all$orf_name %in% sprintf('MASK%d',1:1000),]
fit.all <- fit.all[!fit.all$orf_name %in% sprintf('MASK%d',1:1000),]

for (ref.h in unique(stats.all$hours)) {
  for (mut.h in unique(stats.all$rnd_hrs[stats.all$hours == ref.h])) {
    stats.all$es[stats.all$hours == ref.h & stats.all$rnd_hrs == mut.h] <-
      virP1$es[virP1$ref_hrs == ref.h & virP1$que_hrs == mut.h]
  }
}

#####
unique(stats.all$expt_name)
stats.all$ref[stats.all$expt_name == "4C4_FS_RND2"] <- 0.25
stats.all$ref[stats.all$expt_name == "4C4_TR_RND2"] <- 0.25*3/4
stats.all$ref[stats.all$expt_name == "4C4_TRBL_RND2"] <- 0.25*1/2
stats.all$ref[stats.all$expt_name == "4C4_TRBLBR_RND2"] <- 0.25*1/4


####
load(sprintf('%sSPECIFICITY.RData',out_path))
hi <- plyr::count(spe.data, vars = c('hours','cont_hrs','rep','ref','attempt','pthresh'))
hi <- hi[hi$cont_hrs == hi$hours,]

ggplot(hi) +
  geom_boxplot(aes(x = rep, y = pthresh)) +
  facet_wrap(.~ref)

#####
for (ref in sort(unique(stats.all$ref))) {
  for (rep in sort(unique(stats.all$rep[stats.all$ref == ref]))) {
    for (h in sort(unique(stats.all$hours[stats.all$ref == ref & stats.all$rep == rep]))) {
      stats.all$pthresh[stats.all$ref == ref & stats.all$rep == rep & stats.all$hours == h] <- 
        mean(hi$pthresh[hi$ref == ref & hi$rep == rep & hi$hours == h], na.rm = T)
    }
  }
}

stats.all <- stats.all[stats.all$hours > 0,]
stats.all$pthresh[is.na(stats.all$pthresh)] <- 0.05


stats.all$effect <- NA
stats.all$effect[stats.all$p <= stats.all$pthresh & stats.all$es > 1] <- 'Beneficial'
stats.all$effect[stats.all$p <= stats.all$pthresh & stats.all$es < 1] <- 'Deleterious'
stats.all$effect[is.na(stats.all$effect)] <- 'Neutral'

stats.all$truth[stats.all$hours > stats.all$rnd_hrs] <- 'Deleterious'
stats.all$truth[stats.all$hours < stats.all$rnd_hrs] <- 'Beneficial'
stats.all$truth[stats.all$hours == stats.all$rnd_hrs] <- 'Neutral'

stats.all$phenotype <- paste(strtrim(stats.all$effect,3), strtrim(stats.all$truth,3), sep = '/')
# stats.all$phenotype <- factor(stats.all$phenotype, levels = c('Del/Ben',
                                                            # 'Neu/Ben',
                                                            # 'Ben/Ben',
                                                            # 'Ben/Del',
                                                            # 'Neu/Del',
                                                            # 'Del/Del'))


ggplot(stats.all[stats.all$hours != stats.all$rnd_hrs,]) +
  geom_histogram(aes(x = (es - 1) * 100, y=..count../5, fill = effect), binwidth = 10) +
  labs(title = 'RND REF REP',
       x = 'Fitness Effect',
       y = '') +
  scale_x_continuous(breaks = seq(-1000,1000,100),
                     minor_breaks = seq(-1000,1000,10),
                     labels = paste(seq(-1000,1000,100),'%', sep = '')) +
  # scale_fill_manual(name = 'Phenotype',
  #                   breaks = c('Ben/Ben',
  #                              'Neu/Ben',
  #                              'Del/Ben',
  #                              'Ben/Del',
  #                              'Neu/Del',
  #                              'Del/Del'),
  #                   values = c('Neu/Del'='#536DFE',
  #                              'Del/Del'='#303F9F',
  #                              'Ben/Del'='#C5CAE9',
  #                              'Del/Ben'='#FFECB3',
  #                              'Ben/Ben'='#FFA000',
  #                              'Neu/Ben'='#FFC107'),
  #                   labels = c('Neu/Del'='False-Neutral Deleterious',
  #                              'Del/Del'='True Deleterious',
  #                              'Ben/Del'='False-Beneficial Deleterious',
  #                              'Del/Ben'='False-Deleterious Beneficial',
  #                              'Ben/Ben'='True Beneficial',
  #                              'Neu/Ben'='False-Neutral Beneficial')) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107')) +
  facet_wrap(.~ref*rep, ncol = 8) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm"),
        legend.position = 'bottom',
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
# ggsave(sprintf("%sFIGURE6_effect2.jpg",out_path),
#        height = 120, width = two.c, units = 'mm',
#        dpi = 300)


#####
head(stats.all) 
# rnd.v.data <- plyr::count(stats.all, vars = c('expt_name', 'ref', 'rep', 'es', 'effect'))
# head(rnd.v.data)
rnd.v.data <- NULL

i <- 1
for (ref in unique(stats.all$ref)) {
  for (rep in unique(stats.all$rep[stats.all$ref == ref])) {
    for (es in unique(stats.all$es[stats.all$rep == rep & stats.all$ref == ref])) {
      rnd.v.data$ref[i] <- ref
      rnd.v.data$rep[i] <- rep
      rnd.v.data$es[i] <- es
      # rnd.v.data$hours[i] <- unique(stats.all$hours[stats.all$rep == rep & stats.all$ref == ref & stats.all$es == es])
      # rnd.v.data$rnd_hrs[i] <- unique(stats.all$rnd_hrs[stats.all$rep == rep & stats.all$ref == ref & stats.all$es == es])
      rnd.v.data$Beneficial[i] <- sum(stats.all$effect[stats.all$rep == rep &
                                                         stats.all$ref == ref &
                                                         stats.all$es == es] == 'Beneficial')
      rnd.v.data$Neutral[i] <- sum(stats.all$effect[stats.all$rep == rep &
                                                         stats.all$ref == ref &
                                                         stats.all$es == es] == 'Neutral')
      rnd.v.data$Deleterious[i] <- sum(stats.all$effect[stats.all$rep == rep &
                                                         stats.all$ref == ref &
                                                         stats.all$es == es] == 'Deleterious')
      i <- i + 1
    }
  }
}
rnd.v.data <- data.frame(rnd.v.data)
head(rnd.v.data)
rnd.v.data$total <- rowSums(rnd.v.data[,4:6])
rnd.v.data[,4:6] <- rnd.v.data[,4:6]/rnd.v.data$total * 100

rnd.v.data$Beneficial <- rnd.v.data$Deleterious + rnd.v.data$Beneficial
rnd.v.data$Neutral <- rnd.v.data$Beneficial + rnd.v.data$Neutral

rnd.v.data$ref_title <- sprintf('References/Plate=%0.2f%%', rnd.v.data$ref*100)
rnd.v.data$ref_title <- factor(rnd.v.data$ref_title, levels = unique(rnd.v.data$ref_title)[order(unique(rnd.v.data$ref))])

rnd.v.data$rep_title <- sprintf('Replicates/Strain=%d', rnd.v.data$rep)
rnd.v.data$rep_title <- factor(rnd.v.data$rep_title, levels = unique(rnd.v.data$rep_title)[order(unique(rnd.v.data$rep))])

load(file = sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_RNDVDATA.RData", expt_name))

t <- 100
rnd.v.plot <- ggplot(rnd.v.data) +
  geom_area(aes(x = (es-1)*100, y = Neutral, fill = 'Neutral'), alpha = 1) +
  geom_area(aes(x = (es-1)*100, y = Beneficial, fill = 'Beneficial'), alpha = 1) +
  geom_area(aes(x = (es-1)*100, y = Deleterious, fill = 'Deleterious'), alpha = 1) +
  geom_vline(xintercept = seq(-2,2,0.025)*100, col = '#757575', lwd = 0.25, alpha =0.5) +
  geom_hline(yintercept = c(t * seq(0,1,0.05)), col = '#757575', lwd = 0.25, alpha =0.5) +
  geom_vline(xintercept = c(-5,5), col = 'red', lwd = 0.25, linetype = 'dashed') +
  # geom_text(x=0, y=t*0.9, label="LID", col = 'white', size = 3) +
  labs(x = 'Fitness Effect',
       y = 'Sensitivity') +
  scale_y_continuous(breaks = c(t * seq(0,1,0.2)),
                     minor_breaks = c(t * seq(0,1,0.05)),
                     labels = paste(t * seq(0,1,0.2), '%', sep = "")) +
  scale_x_continuous(breaks = seq(-2,2,0.05)*100,
                     minor_breaks = seq(-2,2,0.025)*100,
                     labels = paste0(round(seq(-2,2,0.05)*100), '%', sep = '')) +
  scale_color_manual(name = 'Phenotype',
                     breaks = c('Beneficial','Neutral','Deleterious'),
                     values = c('Deleterious'='#3F51B5',
                                'Neutral'='#212121',
                                'Beneficial'='#FFC107'),
                     guide = F) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107')) +
  facet_wrap(.~ref_title*rep_title, ncol = 8) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(4, "mm"),
        legend.position = 'bottom',
        legend.margin = margin(0.5,0.5,0.5,0.5, "mm"),
        strip.text = element_text(size = txt-1,
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(xlim = c(-10,10),
                  ylim = c(0,t))
# ggsave(sprintf("%sFIGURE6_effect3.jpg",out_path), rnd.v.plot,
#        height = 110, width = two.c, units = 'mm',
#        dpi = 300)


##### REF REP
load(sprintf('%sREFREPDATA.RData',out_path))
refrep$power[refrep$hours < refrep$cont_hrs] <- refrep$Deleterious[refrep$hours < refrep$cont_hrs]/910 * 100
refrep$power[refrep$hours > refrep$cont_hrs] <- refrep$Beneficial[refrep$hours > refrep$cont_hrs]/910 * 100
refrep$abs_cen <- abs(1-refrep$cen) * 100
refrep$rep <- as.factor(refrep$rep)
# refrep$ref_prop <- as.factor(refrep$ref_prop)
refrep$ref <- refrep$ref_prop
refrep$ref_title <- sprintf('References/Plate=%0.2f%%', refrep$ref*100)
refrep$ref_title <- factor(refrep$ref_title, levels = unique(refrep$ref_title)[order(unique(refrep$ref))])

rnd.v.data$power <- NULL
rnd.v.data$power[rnd.v.data$es < 1] <- rnd.v.data$Deleterious[rnd.v.data$es < 1]
rnd.v.data$power[rnd.v.data$es > 1] <- rnd.v.data$Beneficial[rnd.v.data$es > 1]
rnd.v.data$abs_cen <- abs(1-rnd.v.data$es) * 100

rnd.v.data$rep <- as.factor(rnd.v.data$rep)
rnd.v.data$ref <- as.factor(rnd.v.data$ref)
# rnd.v.data$abs_cen <- as.factor(rnd.v.data$abs_cen)

# save(rnd.v.data, file = sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_RNDVDATA.RData", expt_name))

head(rnd.v.data)

# lid.rnd2.rr.sen <- ggplot(rnd.v.data[round(rnd.v.data$abs_cen) == 5,],
#        aes(x = rep, y = power, fill = ref)) +
#   geom_boxplot(outlier.shape = NA) +
#   # geom_smooth(aes(x = rep, y = power), method = 'loess', se = F, lwd = 1.2) +
#   labs(title = 'Results from random colony-size distribution',
#        x = 'No. of Replicates',
#        y = 'Sensitivity') +
#   scale_x_discrete(breaks = seq(0,16,2)) +
#   scale_y_continuous(breaks = seq(0,100,10),
#                      minor_breaks = seq(0,105,5),
#                      labels = paste(seq(0,100,10),'%',sep='')) +
#   scale_fill_manual(name = 'Reference Proportion',
#                     breaks = c(0.0625,0.125,0.1875,0.25),
#                     values = c('#607D8B','#757575','#BDBDBD','#FFFFFF'),
#                     labels = c(sprintf('%0.2f%%',c(0.0625,0.125,0.1875,0.25)*100))) +
#   coord_cartesian(ylim = c(0,100)) +
#   theme_linedraw() +
#   theme( plot.title = element_text(size = titles),
#          axis.title = element_text(size = titles),
#         axis.text = element_text(size = txt),
#         legend.title = element_text(size = titles),
#         legend.text = element_text(size = txt),
#         legend.key = element_rect(size = 0.5),
#         legend.key.size = unit(.8, 'lines'),
#         legend.position = c(0.5,0.1),
#         legend.margin = margin(0.5,0.5,0.5,0.5, "mm")) +
#   guides(fill = guide_legend(nrow=1))
# 
# 
# lid.rr.sen <- ggplot(refrep[round(refrep$abs_cen) == 5,],
#                       aes(x = rep,
#                           y = power,
#                           fill = ref_prop)) +
#   geom_boxplot(outlier.shape = NA) +
#   # geom_smooth(aes(x = rep, y = power), method = 'loess', se = F, lwd = 1.2) +
#   labs(title = 'Results from bimodal colony-size distribution',
#        x = 'No. of Replicates',
#        y = 'Sensitivity') +
#   scale_x_discrete(breaks = seq(0,16,2)) +
#   scale_y_continuous(breaks = seq(0,100,10),
#                      minor_breaks = seq(0,105,5),
#                      labels = paste(seq(0,100,10),'%',sep='')) +
#   scale_fill_manual(name = 'Reference Proportion',
#                     breaks = c(0.0625,0.125,0.1875,0.25),
#                     values = c('#607D8B','#757575','#BDBDBD','#FFFFFF'),
#                     labels = c(sprintf('%0.2f%%',c(0.0625,0.125,0.1875,0.25)*100))) +
#   coord_cartesian(ylim = c(0,100)) +
#   theme_linedraw() +
#   theme(plot.title = element_text(size = titles),
#         axis.title = element_text(size = titles),
#         axis.text = element_text(size = txt),
#         legend.title = element_text(size = titles),
#         legend.text = element_text(size = txt),
#         legend.key = element_rect(size = 0.5),
#         legend.key.size = unit(.8, 'lines'),
#         legend.position = c(0.5,0.1),
#         legend.margin = margin(0.5,0.5,0.5,0.5, "mm")) +
#   guides(fill = guide_legend(nrow=1))
# 
# fig6 <- ggpubr::ggarrange(lid.rr.sen, lid.rnd2.rr.sen,
#           nrow = 1,
#           common.legend = T,
#           legend = 'bottom',
#           labels = c('a','b'),
#           font.label = list(face = 'bold', size = lbls),
#           hjust=-1)

es = 5
fig6 <- ggplot(rnd.v.data[round(rnd.v.data$abs_cen) == es,],
       aes(x = rep,
           y = power)) +
  geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1),
                size = 0.3, width = 0.5) +
  geom_line(stat = 'summary', fun = 'mean', group = 1,
            size = 0.3) +
  geom_point(stat = 'summary', fun = 'mean', group = 1,
             aes(col = 'Random CS Distribution'),
             size = 1) +
  geom_errorbar(data = refrep[round(refrep$abs_cen) == es,],
                aes(x = rep,
                    y = power),
                stat="summary", fun.data="mean_se", fun.args = list(mult = 1),
                size = 0.3, width = 0.5) +
  geom_line(data = refrep[round(refrep$abs_cen) == es,],
            aes(x = rep,
                y = power),
            stat = 'summary', fun = 'mean', group = 1,
            size = 0.3) +
  geom_point(data = refrep[round(refrep$abs_cen) == es,],
             stat = 'summary', fun = 'mean', group = 1,
             aes(x = rep,
                 y = power,
                 col = 'Bimodal CS Distribution'),
             size = 1) +
  scale_color_manual(name = 'Virtual Plates',
                     breaks = c('Bimodal CS Distribution', 'Random CS Distribution'),
                     values = c('#7C4DFF', '#E64A19')) +
  scale_y_continuous(breaks = seq(0,100,10),
                     labels = paste0(seq(0,100,10),'%',sep = '')) +
  labs(x = 'Replicates per strain',
       y = 'Sensitivity') +
  facet_wrap(.~ref_title, nrow = 1) +
  theme_linedraw() +
  theme( plot.title = element_text(size = titles),
         axis.title = element_text(size = titles),
         axis.text = element_text(size = txt),
         legend.title = element_text(size = titles),
         legend.text = element_text(size = txt+1),
         legend.key = element_rect(size = 0.5),
         legend.key.size = unit(.8, 'lines'),
         legend.position = 'bottom',
         legend.margin = margin(0.5,0.5,0.5,0.5, "mm"),
         strip.text = element_text(size = txt+2,
                                   margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%sFIGURE6.jpg",out_path), fig6,
       height = 70, width = two.c, units = 'mm',
       dpi = 300)


es = 7
fig6 <- ggplot(rnd.v.data[round(rnd.v.data$abs_cen) == es,],
               aes(x = rep,
                   y = power)) +
  geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1),
                size = 0.3, width = 0.5) +
  geom_line(stat = 'summary', fun = 'mean', group = 1,
            size = 0.3) +
  geom_point(stat = 'summary', fun = 'mean', group = 1,
             aes(col = 'Random CS Distribution'),
             size = 1) +
  geom_errorbar(data = refrep[round(refrep$abs_cen) == es,],
                aes(x = rep,
                    y = power),
                stat="summary", fun.data="mean_se", fun.args = list(mult = 1),
                size = 0.3, width = 0.5) +
  geom_line(data = refrep[round(refrep$abs_cen) == es,],
            aes(x = rep,
                y = power),
            stat = 'summary', fun = 'mean', group = 1,
            size = 0.3) +
  geom_point(data = refrep[round(refrep$abs_cen) == es,],
             stat = 'summary', fun = 'mean', group = 1,
             aes(x = rep,
                 y = power,
                 col = 'Bimodal CS Distribution'),
             size = 1) +
  scale_color_manual(name = 'Virtual Plates',
                     breaks = c('Bimodal CS Distribution', 'Random CS Distribution'),
                     values = c('#7C4DFF', '#E64A19')) +
  scale_y_continuous(breaks = seq(0,100,10),
                     labels = paste0(seq(0,100,10),'%',sep = '')) +
  labs(x = 'Replicates per strain',
       y = 'Sensitivity') +
  facet_wrap(.~ref_title, nrow = 1) +
  theme_linedraw() +
  theme( plot.title = element_text(size = titles),
         axis.title = element_text(size = titles),
         axis.text = element_text(size = txt),
         legend.title = element_text(size = titles),
         legend.text = element_text(size = txt+1),
         legend.key = element_rect(size = 0.5),
         legend.key.size = unit(.8, 'lines'),
         legend.position = 'bottom',
         legend.margin = margin(0.5,0.5,0.5,0.5, "mm"),
         strip.text = element_text(size = txt+2,
                                   margin = margin(0.1,0,0.1,0, "mm")))

fig6.2 <- ggarrange(fig6, rnd.v.plot,
                    labels = c('a','b'),
                    label.args = list(gp=gpar(font = 2, fontsize = lbls),
                                      hjust=-1),
                    heights = c(1,2))
ggsave(sprintf("%sFIGURE_S7.jpg",out_path), fig6.2,
       height = two.c, width = two.c, units = 'mm',
       dpi = 300)


##### VP2 Fitness Distribution
load("/home/sbp29/R/Projects/adaptivefitness/output/4C4_FS_CC_VIR_PLATE1_ES.Rdata")
hhours <- c(0,1.02,1.38,2.9,4.02,4.89,6.14,6.88,7.85,8.96,9.97,11.04)

virP1$text <- '11x11 Fitness Effect Matrix'
es.heatmap <- ggplot(virP1[virP1$ref_hrs > 0 & virP1$que_hrs > 0,]) +
  geom_tile(aes(x = sprintf('%0.1f',que_hrs), y = sprintf('%0.1f',ref_hrs), fill = round((es-1)*100,0)),
            col = 'black') +
  geom_text(aes(x = sprintf('%0.1f',que_hrs), y = sprintf('%0.1f',ref_hrs), label = round((es-1)*100,0)),
            col = 'black', size = 2) +
  scale_x_discrete(limits = sprintf('%0.1f',hhours[2:12])) +
  scale_y_discrete(limits = sprintf('%0.1f',sort(hhours[2:12], decreasing = T))) +
  scale_fill_gradient2(name = 'Fitness\nEffect %',
                       low = "#3F51B5", high = "#FFC107", mid = "white",
                       trans = 'pseudo_log',
                       midpoint = 0,
                       breaks = c(-50,-5,0,5,50,500)) +
  labs(title = 'Bimodal Colony-Size Distribution',
       x = 'Mutant Colony-Size Time (hour)',
       y = 'Reference Colony-Size Time (hour)') +
  facet_wrap(.~text) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'right',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt+1,
                                  margin = margin(0.1,0,0.1,0, "mm")))

library(ggridges)
# density.data <- dbGetQuery(conn, 'select a.pos, c.orf_name, a.hours, b.rnd_hrs, a.average, b.average rnd_avg
#                            from 4C4_FS_6144_JPEG a, 4C4_FS_RND2_6144_JPEG b, 4C4_pos2orf_name c
#                            where a.hours = b.hours and a.pos = b.pos and a.pos = c.pos')
# save(density.data, file = sprintf('%sCS_DENSITY_DATA.RData', out_path))
load(file = sprintf('%sCS_DENSITY_DATA.RData', out_path))

density.data$colony[density.data$orf_name == 'BF_control'] <- 'Reference'
density.data$colony[density.data$orf_name != 'BF_control' & !is.na(density.data$orf_name)] <- 'Mutant'

density.data$title <- sprintf('Reference Colony\nSize Time = %0.1f hr', density.data$hours)
density.data$title <- factor(density.data$title, levels = sprintf('Reference Colony\nSize Time = %0.1f hr',
                                                                  c(unique(sort(density.data$hours)))))

# density.data$title <- sprintf('t[R] = %0.1f hr', density.data$hours)
# for (i in 1:length(density.data$title)) {
#   density.data$title[i] <- expression(density.data$title[i])
# }
# density.data$title <- factor(density.data$title, levels = sprintf('Reference Colony\nSize Time = %0.1f hr',
#                                                                   c(unique(sort(density.data$hours)))))

den.plot.cn <- ggplot(density.data[!is.na(density.data$colony) &
                      density.data$hours > 0,],
       aes(x = average, y = colony, group = colony, fill = colony)) +
  geom_density_ridges(quantile_lines = TRUE,
                      scale = 3, alpha = 0.8, size = 0.3,
                      vline_size = 0.2, vline_color = "black") +
  labs(title = '',
       x = 'Colony Size (pixel counts)',
       y = 'Colony Type') +
  scale_fill_manual(name = 'Colony Type',
                    breaks = c("Reference","Mutant"),
                    labels = c("Reference","Mutant"),
                    values = c("Reference" = "#9E9E9E",
                               "Mutant" = "#7B1FA2"),
                    guide = F) +
  coord_cartesian(xlim = c(0,650)) +
  facet_wrap(.~title) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt+1,
                                  margin = margin(0.1,0,0.1,0, "mm")))
# ggsave(sprintf("%sFIGURE2A.jpg",out_path), den.plot.cn,
#        height = 70, width = one.c, units = 'mm',
#        dpi = 300)


den.plot.cp <- ggplot(density.data[!is.na(density.data$colony) &
                                     density.data$hours > 0,],
                      aes(x = rnd_avg, y = colony, group = colony, fill = colony)) +
  geom_density_ridges(quantile_lines = TRUE,
                      scale = 3, alpha = 0.8, size = 0.3,
                      vline_size = 0.2, vline_color = "black") +
  labs(title = '',
       x = 'Colony Size (pixel counts)',
       y = 'Colony Type') +
  scale_fill_manual(name = 'Colony Type',
                    breaks = c("Reference","Mutant"),
                    labels = c("Reference","Mutant"),
                    values = c("Reference" = "#9E9E9E",
                               "Mutant" = "#7B1FA2"),
                    guide = F) +
  coord_cartesian(xlim = c(0,650)) +
  facet_wrap(.~title) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt+1,
                                  margin = margin(0.1,0,0.1,0, "mm")))
# ggsave(sprintf("%sFIGURE2D.jpg",out_path), den.plot.cp,
#        height = 70, width = one.c, units = 'mm',
#        dpi = 300)


fig2 <- ggpubr::ggarrange(NULL, NULL,
                          den.plot.cn, es.heatmap,
                          den.plot.cp, NULL, 
                  nrow = 3, ncol = 2,
                  labels = c('a','b','c','d','e',''),
                  heights = c(0.9,1,1),
                  font.label = list(face = 'bold', size = lbls),
                  hjust=-1)
ggsave(sprintf("%sFIGURE2.jpg",out_path), fig2,
       height = two.c*1.2, width = two.c, units = 'mm',
       dpi = 300)

##### MCAT RND2

##### RND_HRS
mcat.pval <- dbGetQuery(conn, 'select * from 4C4_FS_RND2_BEAN_6144_PVALUE
                        where orf_name != "BFC100"')
for (hr in unique(mcat.pval$hours)) {
  for (orf in unique(mcat.pval$orf_name[mcat.pval$hours == hr & !mcat.pval$orf_name %in% sprintf('MASK%d',1:1000)])) {
    mcat.pval$rnd_hrs[mcat.pval$hours == hr & mcat.pval$orf_name == orf] <-
      rnd_hrs$rnd_hrs[rnd_hrs$hours == hr & rnd_hrs$orf_name == orf]
  }
}

for (ref.h in unique(mcat.pval$hours)) {
  for (mut.h in unique(mcat.pval$rnd_hrs[mcat.pval$hours == ref.h])) {
    mcat.pval$es[mcat.pval$hours == ref.h & mcat.pval$rnd_hrs == mut.h] <-
      virP1$es[virP1$ref_hrs == ref.h & virP1$que_hrs == mut.h]
  }
}

load(sprintf('%sMCAT_SPECIFICITY_DATA.RData',out_path))
hi <- plyr::count(stats.tmp, vars = c('hours','cont_hrs','rep','attempt','pthresh'))
hi <- hi[hi$cont_hrs == hi$hours,]

ggplot(hi) +
  geom_point(aes(x = hours, y = pthresh))

mcat.pval$ref <- 0.25
mcat.pval$rep <- 16
hi$ref <- 0.25

##### PTHRESH
for (ref in sort(unique(mcat.pval$ref))) {
  for (rep in sort(unique(mcat.pval$rep[mcat.pval$ref == ref]))) {
    for (h in sort(unique(mcat.pval$hours[mcat.pval$ref == ref & mcat.pval$rep == rep]))) {
      mcat.pval$pthresh[mcat.pval$ref == ref & mcat.pval$rep == rep & mcat.pval$hours == h] <- 
        mean(hi$pthresh[hi$ref == ref & hi$rep == rep & hi$hours == h], na.rm = T)
    }
  }
}

mcat.pval <- mcat.pval[mcat.pval$hours > 0,]
mcat.pval$pthresh[is.na(mcat.pval$pthresh)] <- 0.05

mcat.pval$effect <- NA
mcat.pval$effect[mcat.pval$p <= mcat.pval$pthresh & mcat.pval$es > 1] <- 'Beneficial'
mcat.pval$effect[mcat.pval$p <= mcat.pval$pthresh & mcat.pval$es < 1] <- 'Deleterious'
mcat.pval$effect[is.na(mcat.pval$effect)] <- 'Neutral'

mcat.pval$truth[mcat.pval$hours > mcat.pval$rnd_hrs] <- 'Deleterious'
mcat.pval$truth[mcat.pval$hours < mcat.pval$rnd_hrs] <- 'Beneficial'
mcat.pval$truth[mcat.pval$hours == mcat.pval$rnd_hrs] <- 'Neutral'

mcat.pval$phenotype <- paste(strtrim(mcat.pval$effect,3), strtrim(mcat.pval$truth,3), sep = '/')
# mcat.pval$phenotype <- factor(mcat.pval$phenotype, levels = c('Del/Ben',
#                                                               'Neu/Ben',
#                                                               'Ben/Ben',
#                                                               'Ben/Del',
#                                                               'Neu/Del',
#                                                               'Del/Del'))



mcat.rnd2.res <- mcat.pval
lid.rnd2.res <- stats.all

save(mcat.rnd2.res, lid.rnd2.res, file = sprintf('%sRND2_RESULTS.RData',out_path))
load(file = sprintf('%sRND2_RESULTS.RData',out_path))

mcat.rnd2.res$effect <- factor(mcat.rnd2.res$effect, levels = c('Beneficial','Neutral','Deleterious'))
lid.rnd2.res$effect <- factor(lid.rnd2.res$effect, levels = c('Beneficial','Neutral','Deleterious'))
mcat.rnd2.res$truth <- factor(mcat.rnd2.res$truth, levels = c('Beneficial','Neutral','Deleterious'))
lid.rnd2.res$truth <- factor(lid.rnd2.res$truth, levels = c('Beneficial','Neutral','Deleterious'))

mcat.pie <- plyr::count(mcat.rnd2.res, vars = 'effect')
mcat.pie$freq <- mcat.pie$freq/sum(mcat.pie$freq) * 100
truth.pie <- plyr::count(mcat.rnd2.res, vars = 'truth')
truth.pie$freq <- truth.pie$freq/sum(truth.pie$freq) * 100
lid.pie <- plyr::count(lid.rnd2.res[lid.rnd2.res$ref == 0.25 & lid.rnd2.res$rep == 16 & lid.rnd2.res$attempt == 1,], vars = 'effect')
lid.pie$freq <- lid.pie$freq/sum(lid.pie$freq) * 100

fig5a.l <- ggplot(lid.pie, aes(x = "", y = freq, fill = effect)) +
  geom_bar(stat = 'identity') +
  coord_polar("y", start = 0) +
  geom_label_repel(aes(label = sprintf("%0.2f%%",freq)),
                   position = position_stack(vjust = 0.5),
                   colour ='white',
                   label.size = 0.15,
                   show.legend = F) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107'),
                    guide = F) +
  labs(title = 'LID CLASSIFICATION') +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(size = titles, hjust = 0.5),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = 'bottom',
        strip.text = element_text(size = titles),
        legend.box.spacing = unit(0.5,"mm"))

fig5a.m <- ggplot(truth.pie, aes(x = "", y = freq, fill = truth)) +
  geom_bar(stat = 'identity') +
  coord_polar("y", start = 0) +
  geom_label_repel(aes(label = sprintf("%0.2f%%",freq)),
                   position = position_stack(vjust = 0.5),
                   colour ='white',
                   label.size = 0.15,
                   show.legend = F) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107'),
                    guide = F) +
  labs(title = 'ACTUAL CLASSIFICATION') +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(size = titles, hjust = 0.5),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = 'bottom',
        strip.text = element_text(size = titles),
        legend.box.spacing = unit(0.5,"mm"))

fig5a.r <- ggplot(mcat.pie, aes(x = "", y = freq, fill = effect)) +
  geom_bar(stat = 'identity') +
  coord_polar("y", start = 0) +
  geom_label_repel(aes(label = sprintf("%0.2f%%",freq)),
                   position = position_stack(vjust = 0.5),
                   colour ='white',
                   label.size = 0.15,
                   show.legend = F) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107'),
                    guide = F) +
  labs(title = 'MCAT CLASSIFICATION') +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(size = titles, hjust = 0.5),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = 'bottom',
        strip.text = element_text(size = titles),
        legend.box.spacing = unit(0.5,"mm"))

mcatbar <- ggplot(mcat.rnd2.res,
                  aes(x = as.factor(hours), fill = effect)) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..),
            col = '#F5F5F5', size = 2,
            position = position_stack(vjust = 0.5)) +
  labs(title = 'MCAT',
       x = 'Reference Population Time Point (hours)',
       y = '') +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107'),
                    guide = F) +
  scale_x_discrete(labels = sprintf('%0.1f',sort(unique(mcat.rnd2.res$hours)))) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm"))

lidbar <- ggplot(lid.rnd2.res[lid.rnd2.res$ref == 0.25 & lid.rnd2.res$rep == 16 & lid.rnd2.res$attempt == 1,],
                 aes(x = as.factor(hours), fill = effect)) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..),
            col = '#F5F5F5', size = 2,
            position = position_stack(vjust = 0.5)) +
  labs(title = 'LID',
       x = 'Reference Population Time Point (hours)',
       y = '') +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107'),
                    guide = F) +
  scale_x_discrete(labels = sprintf('%0.1f',sort(unique(lid.rnd2.res$hours)))) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm"))

truthbar <- ggplot(lid.rnd2.res[lid.rnd2.res$ref == 0.25 & lid.rnd2.res$rep == 16 & lid.rnd2.res$attempt == 1,],
                   aes(x = as.factor(hours), fill = truth)) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..),
            col = '#F5F5F5', size = 2,
            position = position_stack(vjust = 0.5)) +
  labs(title = 'TRUTH',
       x = 'Reference Population Time Point (hours)',
       y = 'Number of Mutants') +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107'),
                    guide = F) +
  scale_x_discrete(labels = sprintf('%0.1f',sort(unique(lid.rnd2.res$hours)))) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm"))

mcat.rnd2.res$effect <- factor(mcat.rnd2.res$effect, levels = c('Beneficial','Deleterious','Neutral'))
lid.rnd2.res$effect <- factor(lid.rnd2.res$effect, levels = c('Beneficial','Deleterious','Neutral'))
lid.rnd2.res$truth <- factor(lid.rnd2.res$truth, levels = c('Beneficial','Deleterious','Neutral'))

mcat.es.histo <- ggplot(mcat.rnd2.res) +
  geom_histogram(aes(x = (es - 1) * 100, y=..count../5, fill = effect), binwidth = 10) +
  labs(title = 'MCAT',
       x = 'Fitness Effect',
       y = '') +
  scale_x_continuous(breaks = seq(-1000,1000,100),
                     minor_breaks = seq(-1000,1000,10),
                     labels = paste(seq(-1000,1000,100),'%', sep = '')) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107')) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm"),
        legend.position = 'bottom',
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))

lid.es.histo <- ggplot(lid.rnd2.res[lid.rnd2.res$ref == 0.25 & lid.rnd2.res$rep == 16 & lid.rnd2.res$attempt == 1,]) +
  geom_histogram(aes(x = (es - 1) * 100, y=..count../5, fill = effect), binwidth = 10) +
  labs(title = 'LID',
       x = 'Fitness Effect',
       y = '') +
  scale_x_continuous(breaks = seq(-1000,1000,100),
                     minor_breaks = seq(-1000,1000,10),
                     labels = paste(seq(-1000,1000,100),'%', sep = '')) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107')) +
  # facet_wrap(.~ref*rep, ncol = 8) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm"),
        legend.position = 'bottom',
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))

truth.es.histo <- ggplot(lid.rnd2.res[lid.rnd2.res$ref == 0.25 & lid.rnd2.res$rep == 16 & lid.rnd2.res$attempt == 1,]) +
  geom_histogram(aes(x = (es - 1) * 100, y=..count../5, fill = truth), binwidth = 10) +
  labs(title = 'TRUTH',
       x = 'Fitness Effect',
       y = 'Number of Mutants') +
  scale_x_continuous(breaks = seq(-1000,1000,100),
                     minor_breaks = seq(-1000,1000,10),
                     labels = paste(seq(-1000,1000,100),'%', sep = '')) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107')) +
  # facet_wrap(.~ref*rep, ncol = 8) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm"),
        legend.position = 'bottom',
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))

# fig5.a <- ggpubr::ggarrange(fig5a.l, fig5a.m, fig5a.r,
#                             nrow = 1)
# fig5.b <- ggpubr::ggarrange(lidbar, truthbar, mcatbar,
#                               nrow = 1)
# fig5.c <- ggpubr::ggarrange(lid.es.histo, mcat.es.histo,
#                             nrow = 1,
#                            common.legend = T,
#                            legend = 'bottom')
# 
# fig5 <- ggarrange(fig5.a, fig5.b, fig5.c,
#                   nrow = 3,
#                   ncol = 1,
#                   labels = c('a','b','c'),
#                   label.args = list(gp=gpar(font = 2, fontsize = lbls),
#                                     hjust=-1))
# ggsave(sprintf("%sFIGURE5.jpg",out_path), fig5,
#        height = two.c, width = two.c, units = 'mm',
#        dpi = 300)


fig5 <- ggpubr::ggarrange(fig5a.m, fig5a.l, fig5a.r,
                          truthbar, lidbar, mcatbar,
                          truth.es.histo, lid.es.histo, mcat.es.histo,
                  nrow = 3,
                  ncol = 3,
                  common.legend = T,
                  legend = 'bottom',
                  labels = c('a','b','c','d','e','f','g','h','i'),
                  font.label = list(face = 'bold', size = lbls),
                  hjust=-1)
ggsave(sprintf("%sFIGURE5.jpg",out_path), fig5,
       height = two.c, width = two.c, units = 'mm',
       dpi = 300)

mcat.pie$power <- mcat.pie$freq/truth.pie$freq * 100
lid.pie$power <- lid.pie$freq/truth.pie$freq * 100
mcat.pie
mean(mcat.pie$power[c(1,3)])
lid.pie
mean(lid.pie$power[c(1,3)])

sum(lid.rnd2.res$phenotype[lid.rnd2.res$ref == 0.25 & lid.rnd2.res$rep == 16] == 'Neu/Neu')
sum(lid.rnd2.res$truth[lid.rnd2.res$ref == 0.25 & lid.rnd2.res$rep == 16] == 'Neutral')

plyr::count(lid.rnd2.res$phenotype[lid.rnd2.res$ref == 0.25 & lid.rnd2.res$rep == 16 & lid.rnd2.res$attempt == 1])
plyr::count(mcat.rnd2.res$phenotype[mcat.rnd2.res$ref == 0.25 & mcat.rnd2.res$rep == 16])
plyr::count(lid.rnd2.res$truth[lid.rnd2.res$ref == 0.25 & lid.rnd2.res$rep == 16 & lid.rnd2.res$attempt == 1])


