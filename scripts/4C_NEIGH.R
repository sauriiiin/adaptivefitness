##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 05/16/2019

##### INITIALIZE
library(RMariaDB)
library(ggplot2)
library(ggExtra)
library(gridExtra)
source("R/functions/initialize.sql.R")

##### GET/SET DATA
expt_name = '4C3_GA1'
expt = 'FS1-1'
out_path = 'figs/neigh/';
density = 6144;

##### CHECK POSITION WISE VARIABILITY
conn <- initialize.sql("saurin_test")

expt = 'FS1-1'
density = 6144;

tablename_fit = sprintf('%s_%d_FITNESS',expt_name,density);
tablename_p2o = '4C3_pos2orf_name1';
tablename_bpos = '4C3_borderpos';

p2c_info = NULL
p2c_info[1] = '4C3_pos2coor6144'
p2c_info[2] = '6144plate'
p2c_info[3] = '6144col'
p2c_info[4] = '6144row'

hours = dbGetQuery(conn, sprintf('select distinct hours from %s order by hours asc', tablename_fit))
p2c = dbGetQuery(conn, sprintf('select * from %s a order by a.%s, a.%s, a.%s',
                               p2c_info[1],
                               p2c_info[2],
                               p2c_info[3],
                               p2c_info[4]))
n_plates = dbGetQuery(conn, sprintf('select distinct %s from %s a order by %s asc',
                                    p2c_info[2],
                                    p2c_info[1],
                                    p2c_info[2]))

for (hr in hours[[1]][8:length(hours[[1]])]) {
  for (pl in n_plates[[1]]) {
    alldat = dbGetQuery(conn, sprintf('select a.*, b.*
                                      from %s a, %s b
                                      where a.hours = %d
                                      and a.pos = b.pos
                                      and b.%s = %d
                                      order by b.%s, b.%s',
                                      tablename_fit,
                                      p2c_info[1],hr,p2c_info[2],
                                      pl,
                                      p2c_info[3],p2c_info[4]))
    
    alldat$source[alldat$`6144row`%%2==1 & alldat$`6144col`%%2==1] = 'TL'
    alldat$source[alldat$`6144row`%%2==0 & alldat$`6144col`%%2==1] = 'BL'
    alldat$source[alldat$`6144row`%%2==1 & alldat$`6144col`%%2==0] = 'TR'
    alldat$source[alldat$`6144row`%%2==0 & alldat$`6144col`%%2==0] = 'BR'
    
    alldat$colony[alldat$orf_name == 'BF_control'] = 'Reference'
    alldat$colony[alldat$orf_name != 'BF_control'] = 'Query'
    alldat$colony[is.na(alldat$orf_name)] = 'Gap'
    
    mn_avg <- NULL
    dif_avg <- NULL
    var_avg <- NULL
    source <- NULL
    cnt = 1
    for (sr in unique(alldat$source)) {
      temp <- alldat[alldat$source == sr,]
      for (i in seq(2,length(unique(temp$`6144col`)))) {
        col <- unique(temp$`6144col`)[i]
        lf <- unique(temp$`6144col`)[i-1]
        rt <- unique(temp$`6144col`)[i+1]
        for (ii in seq(2,length(unique(temp$`6144row`[temp$`6144col` == col])))) {
          row <- unique(temp$`6144row`[temp$`6144col` == col])[ii]
          up <- unique(temp$`6144row`[temp$`6144col` == col])[ii-1]
          dw <- unique(temp$`6144row`[temp$`6144col` == col])[ii+1]
          if (alldat$orf_name[alldat$`6144col` == col & alldat$`6144row` == row] == 'BF_control') {
            if (!is.na(alldat$average[alldat$`6144col` == col & alldat$`6144row` == row])) {
              a <- alldat$average[alldat$`6144col` == col & alldat$`6144row` == row]
              u <- alldat$average[alldat$`6144col` == col & alldat$`6144row` == up]
              d <- alldat$average[alldat$`6144col` == col & alldat$`6144row` == dw]
              l <- alldat$average[alldat$`6144col` == lf & alldat$`6144row` == row]
              r <- alldat$average[alldat$`6144col` == rt & alldat$`6144row` == row]
              alldat$var[alldat$`6144col` == col & alldat$`6144row` == row] <- sd(c(a,u,d,l,r),na.rm = T)/mean(c(a,u,d,l,r),na.rm = T)
              alldat$neigh[alldat$`6144col` == col & alldat$`6144row` == row] <-  mean(c(u,d,l,r),na.rm = T)
              alldat$diff[alldat$`6144col` == col & alldat$`6144row` == row] <- a - mean(c(u,d,l,r),na.rm = T)
              mn_avg <- c(mn_avg, mean(c(u,d,l,r),na.rm = T))
              dif_avg <- c(dif_avg, (a - mean(c(u,d,l,r),na.rm = T))/a*100)
              var_avg <- c(var_avg, sd(c(a,u,d,l,r),na.rm = T)/mean(c(a,u,d,l,r),na.rm = T))
              source <- c(source, alldat$source[alldat$`6144col` == col & alldat$`6144row` == row])
              cnt <- cnt + 1
            }
          }
        }
      }
      diff_std <- sd(alldat$diff[alldat$source == sr],na.rm = T)
      diff_mean <- mean(alldat$diff[alldat$source == sr],na.rm = T)
      alldat$outlier[alldat$source == sr & !is.na(alldat$average) & alldat$diff > (diff_mean + 2*diff_std)] = 'Right'
      alldat$outlier[alldat$source == sr & !is.na(alldat$average) & alldat$diff < (diff_mean - 2*diff_std)] = 'Left'
    }