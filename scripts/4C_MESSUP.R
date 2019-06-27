##### 4C_MESSUP
##### MESSING UP THE FINAL PLATES
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 06/26/2019

##### INITIALIZE
library(RMariaDB)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyverse)
# library(egg)
library(ggpubr)
library(stringr)
source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")
out_path = 'figs/';
dat.dir <- "/home/sbp29/R/Projects/proto_plots/rawdata/4C3_GA1_RND_LID/"
expt_name <- '4C3_GA1_RND'
pvals = seq(0,1,0.005)

##### MAKING UP THE MESS
alldat = dbGetQuery(conn, 'select b.*, a.orf_name, a.hours, a.bg, a.average, a.fitness
              from 4C3_GA1_RAW_6144_FITNESS a, 4C3_pos2coor6144 b
              where a.pos = b.pos
              order by hours, 6144plate, 6144col, 6144row')
alldat$rnd_hrs <- alldat$hours

# for (i in 1:dim(alldat)[1]) {
#   if (alldat$orf_name[i] != "BF_control" & !is.na(alldat$orf_name[i])) {
#     temp <- sample_n(alldat[alldat$pos == alldat$pos[i] & alldat$hours != alldat$hours[i],], 1)
#     alldat$rnd_hrs[i] <- temp$hours
#     alldat$average[i] <- temp$average
#   }
# }

for (orf in unique(alldat$orf_name[!is.na(alldat$orf_name) & alldat$orf_name != "BF_control"])) {
  for (hr in unique(alldat$hours)) {
    hr.tmp <- sample(unique(alldat$hours[alldat$hours != hr]),1)
    alldat$average[!is.na(alldat$orf_name) & alldat$orf_name == orf & alldat$hours == hr] <-
      alldat$average[!is.na(alldat$orf_name) & alldat$orf_name == orf & alldat$hours == hr.tmp]
    alldat$rnd_hrs[!is.na(alldat$orf_name) & alldat$orf_name == orf & alldat$hours == hr] <-
      alldat$hours[!is.na(alldat$orf_name) & alldat$orf_name == orf & alldat$hours == hr.tmp]
  }
}

dbWriteTable(conn, "4C3_GA1_RND_6144_DATA", alldat, overwrite = T)
dbWriteTable(conn, "4C3_GA1_RND_BEAN_6144_DATA", alldat, overwrite = T)

dbExecute(conn, 'update 4C3_GA1_RND_6144_DATA
    set average = NULL
    where pos in
    (select pos from 4C3_borderpos)')

dbExecute(conn, 'create table 4C3_GA1_RND_6144_JPEG
    select pos, hours, rnd_hrs, average
    from 4C3_GA1_RND_6144_DATA')

dbExecute(conn, 'create table 4C3_GA1_RND_BEAN_6144_JPEG
    select pos, hours, rnd_hrs, average
    from 4C3_GA1_RND_BEAN_6144_DATA')

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
  for (i in 1:length(hours)) {
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
fit.all[fit.all$orf_name == "BF_control",]

ref.all <- NULL
i = 0
for (hr in sort(unique(fit.all$hours))) {
  for (pos in 1152:1536) {
    i <- i + 1
    ref.all$pix_mean[i] <-
      mean(fit.all$average[fit.all$pos %in% c(c(110000,120000,130000,140000,210000,220000,230000,240000) + pos) & fit.all$hours == hr], na.rm = T)
    ref.all$cs_mean[i] <-
      mean(fit.all$fitness[fit.all$pos %in% c(c(110000,120000,130000,140000,210000,220000,230000,240000) + pos) & fit.all$hours == hr], na.rm = T)
    ref.all$hours[i] <- hr
  }
}
ref.all <- data.frame(ref.all)

hr = 11

cs <- ggplot() +
  geom_point(data = fit.all[fit.all$hours == hr,],
             aes(x = x6144col_1, y = x6144row_1, col = average), shape = 15) +
  scale_x_continuous(breaks = seq(1,96,1)) +
  scale_y_continuous(breaks = seq(1,64,1),trans = 'reverse') +
  scale_color_distiller(name = "PIX",
                        limits = c(quantile(fit.all$average[fit.all$hours == hr],0.05,na.rm = T)[[1]],
                                   quantile(fit.all$average[fit.all$hours == hr],0.95,na.rm = T)[[1]]),
                        palette = "Spectral") +
  labs(title = sprintf("LID at %d hours",hr),
       subtitle = "Colony Sizes",
       x = "Columns",
       y = "Rows") +
  theme_linedraw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank())

f <- ggplot() +
  geom_point(data = fit.all[fit.all$hours == hr,],
             aes(x = x6144col_1, y = x6144row_1, col = fitness), shape = 15) +
  scale_x_continuous(breaks = seq(1,96,1)) +
  scale_y_continuous(breaks = seq(1,64,1),trans = 'reverse') +
  scale_color_distiller(name = "FIT",
                        limits = c(quantile(fit.all$fitness[fit.all$hours == hr],0.05,na.rm = T)[[1]],
                                   quantile(fit.all$fitness[fit.all$hours == hr],0.95,na.rm = T)[[1]]),
                        palette = "PiYG") +
  labs(title = "",
       subtitle = "Fitness",
       x = "Columns",
       y = "Rows") +
  theme_linedraw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank())

csVf <- ggplot() +
  geom_point(data = ref.all[ref.all$hours == hr,],
             aes(x = pix_mean, y = cs_mean, col = "Reference"), alpha = 0.5) +
  # geom_smooth(data = fit.all[fit.all$cont_hrs == hr & fit.all$hours == hr & fit.all$orf_name == "BF_control",],
  #            aes(x = average, y = fitness), method = lm, se = F) +
  geom_point(data = stats.all[stats.all$hours == hr,],
             aes(x = pix_mean, y = cs_mean, col = effect_p)) +
  # geom_smooth(data = stats.all[stats.all$hours == hr,],
  #             aes(x = pix_mean, y = cs_mean), method = lm, se = F) +
  scale_color_manual(name = "Strain Kind",
                     breaks = c("Reference", "Beneficial", "Neutral", "Deleterious"),
                     values = c("Reference" = "#607D8B",
                                "Beneficial" = "#388E3C",
                                "Neutral" = "#303F9F",
                                "Deleterious" = "#D32F2F")) +
  scale_x_continuous(breaks = seq(0,1000,100), minor_breaks = seq(0,1000,25)) +
  scale_y_continuous(breaks = seq(-10,10,0.2), minor_breaks = seq(-10,10,0.05)) +
  labs(title = "",
       subtitle = "CS and Fitness Correlation",
       x = "Colony Size (pix)",
       y = "Fitness") +
  theme_linedraw() +
  coord_cartesian(xlim = c(100,500),
                  ylim = c(0.2,2))

ggarrange(cs,f,csVf,
          nrow = 1, ncol = 3, widths = c(1.8,1.8,1.4))

ggsave(sprintf("%s%s_CSFIT_%d.png",
               out_path,expt_name,hr),
       width = 18,height = 4.5)
