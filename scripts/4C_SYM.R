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

expt_name = '4C3_GA1'
expt = 'FS1-GA1'
out_path = 'figs/sym/';
density = 6144;

dat.dir <- "/home/sbp29/R/Projects/proto_plots/rawdata/4C3_GA1_LID/"
dat.dir_cc <- "/home/sbp29/R/Projects/proto_plots/rawdata/4C3_GA1_CC_LID/"

cont.name = 'BF_control'
tablename_p2o  = '4C3_pos2orf_name1'
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
  hours <- c(hours, as.numeric(s[3]))
}
reps <- unique(reps)
hours <- sort(unique(hours))

stats.all <- NULL
stats.all_cc <- NULL

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
reppos <- NULL
for (i in c(110000,120000,130000,140000,210000,220000,230000,240000)) {
  reppos <- cbind(reppos, contpos + i)
}

hr.r <- 17
hr.q <- 18

contfit <- NULL
contfit_cc <- NULL
for (ii in 1:dim(reppos)[[1]]) {
  temp <- dbGetQuery(conn, sprintf('select fitness from %s 
                                   where hours = %d and pos in (%s)',
                                   tablename_fit, hr.r,
                                   str_c(reppos[ii,], collapse=",")))
  temp_cc <- dbGetQuery(conn, sprintf('select fitness from %s 
                                   where hours = %d and pos in (%s)',
                                   tablename_fit_cc, hr.r,
                                   str_c(reppos[ii,], collapse=",")))
  if (sum(temp$fitness, na.rm = T) > 0) {
    md <- mad(temp$fitness, na.rm = T)
    m <- median(temp$fitness, na.rm = T)
    contfit <- c(contfit, mean(temp$fitness[temp$fitness > m - 3*md & temp$fitness < m + 3*md], na.rm = T))
    md <- mad(temp_cc$fitness, na.rm = T)
    m <- median(temp_cc$fitness, na.rm = T)
    contfit_cc <- c(contfit_cc, mean(temp_cc$fitness[temp_cc$fitness > m - 3*md & temp_cc$fitness < m + 3*md], na.rm = T))
  }
}
contfit <- data.frame(contfit)
contfit_cc <- data.frame(contfit_cc)
colnames(contfit) <- 'cs_mean'
colnames(contfit_cc) <- 'cs_mean'

# qfit <- dbGetQuery(conn, sprintf('select * from %s
#                                  where hours = %d',
#                                  tablename_fits, hr))
# 
# qfit_cc <- dbGetQuery(conn, sprintf('select * from %s
#                                  where hours = %d',
#                                  tablename_fits_cc, hr))

##### PLOTTING FITNESS
fit.raw <- ggplot() +
  geom_histogram(data = stats.all[stats.all$hours == hr.q & stats.all$cont_hrs == hr.r,],
                 aes(x = cs_mean, fill = 'Query'),
                 col = 'black', binwidth = 0.01, alpha = 0.7) +
  geom_histogram(data = contfit, aes(x = cs_mean, fill = 'Ref'),
                 col = 'black', binwidth = 0.01, alpha = 0.7) +
  geom_vline(xintercept = quantile(contfit, c(0.025, 0.975), na.rm = T),
             col = 'red', linetype = 'dashed') +
  scale_fill_discrete(name = 'Colony Type') +
  labs(title = 'Fitness Distribution',
       subtitle = sprintf('%s | R = %d Hrs | Q = %d Hrs', expt, hr.r, hr.q),
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
  geom_vline(xintercept = quantile(contfit_cc, c(0.025, 0.975), na.rm = T),
             col = 'red', linetype = 'dashed') +
  scale_fill_discrete(name = 'Colony Type') +
  labs(title = '',
       subtitle = sprintf('%s-CC | R = %d Hrs | Q = %d Hrs', expt, hr.r, hr.q),
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

