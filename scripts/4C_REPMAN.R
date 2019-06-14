##### REPLICATE MANAGEMENT
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 06/12/2019

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
reps <- seq(1,20,1)
for (s in strsplit(stats.files,'_')) {
  hours <- c(hours, as.numeric(s[3]))
}
hours <- unique(hours)

stats.all <- NULL
fit.all <- NULL
fpr.all <- NULL
##### PUTTING IT TOGETHER
for (i in 1:length(hours)) {
  hr <- hours[i]
  dat.stats <- read.csv(paste0(dat.dir,
                               sprintf('%s_%d_STATS_P.csv',expt_name,hr)),
                        na.strings = "NaN")
  dat.stats$cont_hrs <- hr
  
  dat.fit <- read.csv(paste0(dat.dir,
                             sprintf('%s_%d_FITNESS.csv',expt_name,hr)),
                      na.strings = "NaN")
  dat.fit$cont_hrs <- hr
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
  
  dat.stats$effect_p[dat.stats$p <= 0.05 & dat.stats$cs_mean > cont.mean] <- 'Beneficial'
  dat.stats$effect_p[dat.stats$p <= 0.05 & dat.stats$cs_mean < cont.mean] <- 'Deleterious'
  dat.stats$effect_p[is.na(dat.stats$effect_p)] <- 'Neutral'
  
  stats.all <- rbind(stats.all,dat.stats)
  fit.all <- rbind(fit.all,dat.fit)
}