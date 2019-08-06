##### TALK FIGURES
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 08/03/2019

##### INITIALIZE
library(RMariaDB)
library(readxl)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(grid)
library(tidyverse)
# library(egg)
source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

##### INITIALIZE
expt_name = '4C3_GA3'
expt = 'FS1-3'
out_path = 'figs/talk/';
density = 6144;

##### GET DATA
conn <- initialize.sql("saurin_test")
# tablename_fit = sprintf('%s_%d_FITNESS',expt_name,density);

alldat = dbGetQuery(conn, 'select c.*, a.*
                          from 4C3_GA3_CC_6144_FITNESS a, 4C3_pos2coor6144 c
                          where a.pos = c.pos
                          order by a.hours, 6144plate, 6144col, 6144row')

alldat.b = dbGetQuery(conn, 'select c.*, b.*
                          from 4C3_GA3_BEAN_6144_FITNESS b, 4C3_pos2coor6144 c
                          where b.pos = c.pos
                          order by b.hours, 6144plate, 6144col, 6144row')
dbDisconnect(conn)

##### ERROR
rmse.dat <- NULL
for (hr in sort(unique(alldat$hours))) {
  errdat <- alldat[alldat$hours == hr & !is.na(alldat$average),]
  errdat.b <- alldat.b[alldat.b$hours == hr & !is.na(alldat.b$average),]
  rmse.dat <- rbind(rmse.dat, cbind(hr,
                                    sqrt(mean((errdat$average - errdat$bg)^2, na.rm = T))/mean(errdat$average, na.rm = T) * 100,
                                    sqrt(mean((errdat.b$average - errdat.b$bg)^2, na.rm = T))/mean(errdat.b$average, na.rm = T) * 100,
                                    sqrt(mean((errdat$average - sample(errdat$average))^2, na.rm = T))/mean(errdat$average, na.rm = T)* 100))
}
rmse.dat <- data.frame(rmse.dat)
colnames(rmse.dat) <- c('hours','LID','BEAN','RND')

# my_comparisons <- list( c("BEAN", "LID"), c("BEAN", "RND"), c("LID", "BEAN") )

ggplot(rmse.dat[rmse.dat$hours > 10,]) +
  geom_boxplot(aes(x = 'LID', y = LID, fill = 'LID')) +
  geom_boxplot(aes(x = 'BEAN', y = BEAN, fill = 'BEAN')) +
  geom_boxplot(aes(x = 'RND', y = RND, fill = 'RND')) +
  # stat_compare_means(comparisons = my_comparisons, method = "t.test") +
  coord_cartesian(ylim = c(0,30)) +
  scale_fill_discrete(name = 'Method', guide = F) +
  labs(title = 'Error in Background Prediction',
       x = 'Method',
       y = 'RMSE %') +
  theme_linedraw()
# ggsave(sprintf("%s%s_ERROR_BOX.png",
#                out_path,expt_name),
#        width = 6, height = 5)

# t.test(rmse.dat$LID, rmse.dat$RND)
