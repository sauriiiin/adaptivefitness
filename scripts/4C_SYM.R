##### SYMMETRY ISSUE
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 07/30/2019

##### INITIALIZE
library(RMariaDB)
library(ggplot2)
library(ggExtra)
library(gridExtra)
library(ggpubr)
library(stringr)
source("R/functions/initialize.sql.R")

##### GET/SET DATA
conn <- initialize.sql("saurin_test")

expt_name = '4C3_GA3'
expt = 'FS1-GA3'
out_path = 'figs/sym/';
density = 6144;

dat.dir <- "/home/sbp29/R/Projects/proto_plots/rawdata/4C3_GA3_LID/"
dat.dir_cc <- "/home/sbp29/R/Projects/proto_plots/rawdata/4C3_GA3_CC_LID/"

cont.name = 'BF_control'
tablename_p2o  = '4C3_pos2orf_name3'
tablename_bpos = '4C3_borderpos'
tablename_fit = sprintf('%s_%d_FITNESS',expt_name,density)
tablename_fits = sprintf('%s_%d_FITNESS_STATS',expt_name,density)
tablename_fit_cc = sprintf('%s_CC_%d_FITNESS',expt_name,density)
tablename_fits_cc = sprintf('%s_CC_%d_FITNESS_STATS',expt_name,density)

##### VIRTUAL PLATE DATA
pvals = seq(0,1,0.005)
stats.files <- list.files(path = dat.dir,
                          pattern = "P.csv", recursive = TRUE)
stats.files_cc <- list.files(path = dat.dir,
                             pattern = "P.csv", recursive = TRUE)
hours <- NULL
reps <- NULL
for (s in strsplit(stats.files,'_')) {
  # reps <- c(reps, as.numeric(s[4]))
  reps <- 8
  hours <- c(hours, as.numeric(s[4]))
}
reps <- unique(reps)
hours <- sort(unique(hours))

stats.all <- NULL
stats.all_cc <- NULL

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
    stats.all <- rbind(stats.all,dat.stats)
    
    dat.stats_cc <- read.csv(paste0(dat.dir_cc,
                                    sprintf('%s_CC_%d_%d_STATS_P.csv',expt_name,rep,hr)),
                             na.strings = "NaN")
    dat.stats_cc <- dat.stats_cc[dat.stats_cc$hours != 0,]
    dat.stats_cc$cont_hrs <- hr
    dat.stats_cc$rep <- rep
    stats.all_cc <- rbind(stats.all_cc,dat.stats_cc)
  }
}

##### REF DATA
contpos = dbGetQuery(conn, sprintf('select pos from %s
                               where orf_name = "%s" and pos < 10000
                               and pos not in
                               (select pos from %s)',
                              tablename_p2o,cont.name,tablename_bpos))
contpos <- contpos$pos
restpos = dbGetQuery(conn, sprintf('select pos from %s
                               where orf_name != "%s" and pos < 10000
                               and pos not in
                               (select pos from %s)',
                              tablename_p2o,cont.name,tablename_bpos))
restpos <- restpos$pos

refpos <- NULL
quepos <- NULL
for (i in c(110000,120000,130000,140000,210000,220000,230000,240000)) {
  refpos <- cbind(refpos, contpos + i)
  quepos <- cbind(quepos, restpos + i)
}

hr.r <- 14
hr.q <- 16

contfit <- NULL
contfit_cc <- NULL
contavg <- NULL
contavg_cc <- NULL
for (ii in 1:dim(refpos)[[1]]) {
  temp <- dbGetQuery(conn, sprintf('select * from %s
                                   where hours = %d and pos in (%s)',
                                   tablename_fit, hr.r,
                                   str_c(refpos[ii,], collapse=",")))
  temp_cc <- dbGetQuery(conn, sprintf('select * from %s
                                   where hours = %d and pos in (%s)',
                                   tablename_fit_cc, hr.r,
                                   str_c(refpos[ii,], collapse=",")))
  if (sum(temp$fitness, na.rm = T) > 0) {
    md <- mad(temp$fitness, na.rm = T)
    m <- median(temp$fitness, na.rm = T)
    contfit <- c(contfit, mean(temp$fitness[temp$fitness > m - 3*md & temp$fitness < m + 3*md], na.rm = T))
    md <- mad(temp_cc$fitness, na.rm = T)
    m <- median(temp_cc$fitness, na.rm = T)
    contfit_cc <- c(contfit_cc, mean(temp_cc$fitness[temp_cc$fitness > m - 3*md & temp_cc$fitness < m + 3*md], na.rm = T))

    md <- mad(temp$average, na.rm = T)
    m <- median(temp$average, na.rm = T)
    contavg <- c(contavg, mean(temp$average[temp$average > m - 3*md & temp$average < m + 3*md], na.rm = T))
    md <- mad(temp_cc$average, na.rm = T)
    m <- median(temp_cc$average, na.rm = T)
    contavg_cc <- c(contavg_cc, mean(temp_cc$average[temp_cc$average > m - 3*md & temp_cc$average < m + 3*md], na.rm = T))
  }
}
contfit <- data.frame(contfit)
contfit_cc <- data.frame(contfit_cc)
contavg <- data.frame(contavg)
contavg_cc <- data.frame(contavg_cc)
colnames(contfit) <- 'cs_mean'
colnames(contfit_cc) <- 'cs_mean'
colnames(contavg) <- 'pix_mean'
colnames(contavg_cc) <- 'pix_mean'

qfit <- NULL
qfit_cc <- NULL
qavg <- NULL
qavg_cc <- NULL
for (ii in 1:dim(quepos)[[1]]) {
  temp <- dbGetQuery(conn, sprintf('select * from %s
                                   where hours = %d and pos in (%s)',
                                   tablename_fit, hr.q,
                                   str_c(quepos[ii,], collapse=",")))
  temp_cc <- dbGetQuery(conn, sprintf('select * from %s
                                   where hours = %d and pos in (%s)',
                                   tablename_fit_cc, hr.q,
                                   str_c(quepos[ii,], collapse=",")))
  if (sum(temp$fitness, na.rm = T) > 0) {
    md <- mad(temp$fitness, na.rm = T)
    m <- median(temp$fitness, na.rm = T)
    qfit <- c(qfit, mean(temp$fitness[temp$fitness > m - 3*md & temp$fitness < m + 3*md], na.rm = T))
    md <- mad(temp_cc$fitness, na.rm = T)
    m <- median(temp_cc$fitness, na.rm = T)
    qfit_cc <- c(qfit_cc, mean(temp_cc$fitness[temp_cc$fitness > m - 3*md & temp_cc$fitness < m + 3*md], na.rm = T))

    md <- mad(temp$average, na.rm = T)
    m <- median(temp$average, na.rm = T)
    qavg <- c(qavg, mean(temp$average[temp$average > m - 3*md & temp$average < m + 3*md], na.rm = T))
    md <- mad(temp_cc$average, na.rm = T)
    m <- median(temp_cc$average, na.rm = T)
    qavg_cc <- c(qavg_cc, mean(temp_cc$average[temp_cc$average > m - 3*md & temp_cc$average < m + 3*md], na.rm = T))
  }
}
qfit <- data.frame(qfit)
qfit_cc <- data.frame(qfit_cc)
qavg <- data.frame(qavg)
qavg_cc <- data.frame(qavg_cc)
colnames(qfit) <- 'cs_mean'
colnames(qfit_cc) <- 'cs_mean'
colnames(qavg) <- 'pix_mean'
colnames(qavg_cc) <- 'pix_mean'

##### PLOTTING FITNESS
pix.raw <- ggplot() +
  geom_line(data = qavg, aes(x = pix_mean, col = 'Query'),
                 stat = 'density', lwd = 1.2) +
  geom_line(data = contavg, aes(x = pix_mean, col = 'Ref'),
            stat = 'density', lwd = 1.2) +
  scale_color_discrete(name = 'Colony Type') +
  # scale_fill_discrete(name = 'Colony Type') +
  labs(title = 'Colony Size Distribution',
       subtitle = sprintf('%s | R = %d Hrs | Q = %d Hrs | ES = %0.3f',
                          expt, hr.r, hr.q, mean(qavg$pix_mean)/mean(contavg$pix_mean)),
       x = "Pixel Count",
       y = "Count") +
  coord_cartesian(xlim = c(200, 600)) +
  theme_linedraw() +
  theme(axis.title.x = element_blank())

pix.cc <- ggplot() +
  geom_line(data = qavg_cc, aes(x = pix_mean, col = 'Query'),
                 stat = 'density', lwd = 1.2) +
  geom_line(data = contavg_cc, aes(x = pix_mean, col = 'Ref'),
                 stat = 'density', lwd = 1.2) +
  scale_color_discrete(name = 'Colony Type') +
  labs(title = '',
       subtitle = sprintf('%s-CC | R = %d Hrs | Q = %d Hrs | ES = %0.3f',
                          expt, hr.r, hr.q, mean(qavg_cc$pix_mean)/mean(contavg_cc$pix_mean)),
       x = "Pixel Count",
       y = "Count") +
  coord_cartesian(xlim = c(200, 600)) +
  theme_linedraw()

ggarrange(pix.raw, pix.cc,
          ncol = 1,
          common.legend = T,
          legend = 'bottom')
# ggsave(sprintf("%s/%s_PIXDIS_%d_%d.jpg",out_path,expt_name,hr.r,hr.q),
#        height = 12, width = 7,
#        dpi = 300)


fit.raw <- ggplot() +
  geom_histogram(data = stats.all[stats.all$hours == hr.q & stats.all$cont_hrs == hr.r,],
                 aes(x = cs_mean, fill = 'Query'),
                 col = 'black', binwidth = 0.01, alpha = 0.7) +
  geom_histogram(data = contfit, aes(x = cs_mean, fill = 'Ref'),
                 col = 'black', binwidth = 0.01, alpha = 0.7) +
  # geom_line(data = stats.all[stats.all$hours == hr.q & stats.all$cont_hrs == hr.r,],
  #                aes(x = cs_mean, col = 'Query'),
  #                stat = 'density', lwd = 1.2) +
  # geom_line(data = contfit,
  #                aes(x = cs_mean, col = 'Ref'),
  #                stat = 'density', lwd = 1.2) +
  geom_vline(xintercept = quantile(contfit, c(0.025, 0.975), na.rm = T),
             col = 'red', linetype = 'dashed') +
  scale_fill_discrete(name = 'Colony Type') +
  labs(title = 'Fitness Distribution',
       subtitle = sprintf('%s | R = %d Hrs | Q = %d Hrs | ES = %0.3f',
                          expt, hr.r, hr.q, mean(stats.all$cs_mean[stats.all$hours == hr.q & stats.all$cont_hrs == hr.r])/mean(contfit$cs_mean)),
       x = "Fitness",
       y = "Count") +
  coord_cartesian(xlim = c(0.85, 1.15)) +
  theme_linedraw() +
  theme(axis.title.x = element_blank())

fit.cc <- ggplot() +
  geom_histogram(data = stats.all_cc[stats.all_cc$hours == hr.q & stats.all_cc$cont_hrs == hr.r,],
                 aes(x = cs_mean, fill = 'Query'),
                 col = 'black', binwidth = 0.01, alpha = 0.7) +
  geom_histogram(data = contfit_cc, aes(x = cs_mean, fill = 'Ref'),
                 col = 'black', binwidth = 0.01, alpha = 0.7) +
  # geom_line(data = stats.all_cc[stats.all_cc$hours == hr.q & stats.all_cc$cont_hrs == hr.r,],
  #           aes(x = cs_mean, col = 'Query'),
  #           stat = 'density', lwd = 1.2) +
  # geom_line(data = contfit_cc,
  #           aes(x = cs_mean, col = 'Ref'),
  #           stat = 'density', lwd = 1.2) +
  geom_vline(xintercept = quantile(contfit_cc, c(0.025, 0.975), na.rm = T),
             col = 'red', linetype = 'dashed') +
  scale_fill_discrete(name = 'Colony Type') +
  labs(title = '',
       subtitle = sprintf('%s-CC | R = %d Hrs | Q = %d Hrs| ES = %0.3f',
                          expt, hr.r, hr.q, mean(stats.all$cs_mean[stats.all_cc$hours == hr.q & stats.all$cont_hrs == hr.r])/mean(contfit_cc$cs_mean)),
       x = "Fitness",
       y = "Count") +
  coord_cartesian(xlim = c(0.85, 1.15)) +
  theme_linedraw()

ggarrange(fit.raw, fit.cc,
          ncol = 1,
          common.legend = T,
          legend = 'bottom')
ggsave(sprintf("%s/%s_FDIS_%d_%d.jpg",out_path,expt_name,hr.r,hr.q),
       height = 12, width = 7,
       dpi = 300)



