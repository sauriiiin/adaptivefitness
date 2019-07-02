##### COMPETITION
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 07/02/2019
##### Correcting for competion on the plate

##### INITIALIZE
library(RMariaDB)
library(ggplot2)
library(ggExtra)
library(gridExtra)
source("R/functions/initialize.sql.R")

##### GET/SET DATA
expt_name = '4C3_GA1'
expt = 'FS1-1'
out_path = 'figs/comp/';
density = 6144;

##### CHECK POSITION WISE VARIABILITY
conn <- initialize.sql("saurin_test")

tablename_jpeg = sprintf('%s_%d_JPEG',expt_name,density);
tablename_fit = sprintf('%s_%d_FITNESS',expt_name,density);
tablename_p2o = '4C3_pos2orf_name1';
tablename_bpos = '4C3_borderpos';

p2c_info = NULL
p2c_info[1] = '4C3_pos2coor6144'
p2c_info[2] = '6144plate'
p2c_info[3] = '6144col'
p2c_info[4] = '6144row'

hours = dbGetQuery(conn, sprintf('select distinct hours from %s order by hours asc', tablename_fit))

jpegdat <- data.frame()

for (hr in hours$hours) {
  for (pl in n_plates$`6144plate`) {
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
    
    alldat$source[alldat$`6144row`%%2==1 & alldat$`6144col`%%2==1] = '1TL'
    alldat$source[alldat$`6144row`%%2==0 & alldat$`6144col`%%2==1] = '3BL'
    alldat$source[alldat$`6144row`%%2==1 & alldat$`6144col`%%2==0] = '2TR'
    alldat$source[alldat$`6144row`%%2==0 & alldat$`6144col`%%2==0] = '4BR'
    
    alldat$colony[alldat$orf_name == 'BF_control'] = 'Reference'
    alldat$colony[alldat$orf_name != 'BF_control'] = 'Query'
    alldat$colony[is.na(alldat$orf_name)] = 'Gap'
    
    for (sr in unique(alldat$source)) {
      temp <- alldat[alldat$source == sr,]
      for (i in seq(2,length(unique(temp$`6144col`)))) {
        col <- unique(temp$`6144col`)[i]
        lf <- tail(temp$`6144col`[temp$`6144col` < col & temp$orf_name == 'BF_control' & !is.na(temp$orf_name)],1)
        rt <- temp$`6144col`[temp$`6144col` > col & temp$orf_name == 'BF_control' & !is.na(temp$orf_name)][1]
        for (ii in seq(2,length(unique(temp$`6144row`[temp$`6144col` == col])))) {
          row <- unique(temp$`6144row`[temp$`6144col` == col])[ii]
          up <- tail(temp$`6144row`[temp$`6144row` < row & temp$orf_name == 'BF_control' & !is.na(temp$orf_name)],1)
          dw <- temp$`6144row`[temp$`6144row` > row & temp$orf_name == 'BF_control' & !is.na(temp$orf_name)][1]
          if (!is.na(alldat$average[alldat$`6144col` == col & alldat$`6144row` == row])) {
            a <- alldat$average[alldat$`6144col` == col & alldat$`6144row` == row]
            u <- alldat$average[alldat$`6144col` == col & alldat$`6144row` == up]
            d <- alldat$average[alldat$`6144col` == col & alldat$`6144row` == dw]
            l <- alldat$average[alldat$`6144col` == lf & alldat$`6144row` == row]
            r <- alldat$average[alldat$`6144col` == rt & alldat$`6144row` == row]
            alldat$var[alldat$`6144col` == col & alldat$`6144row` == row] <- sd(c(a,u,d,l,r),na.rm = T)/mean(c(a,u,d,l,r),na.rm = T)
            alldat$neigh[alldat$`6144col` == col & alldat$`6144row` == row] <-  mean(c(u,d,l,r),na.rm = T)
            alldat$diff[alldat$`6144col` == col & alldat$`6144row` == row] <- a - mean(c(u,d,l,r),na.rm = T)
          }
          # cnt <- cnt + 1
        }
      }
      diff_std <- sd(alldat$diff[alldat$source == sr],na.rm = T)
      diff_mean <- mean(alldat$diff[alldat$source == sr],na.rm = T)
      alldat$outlier[alldat$source == sr & !is.na(alldat$average) & alldat$diff > (diff_mean + 3*diff_std)] = 'Bigger'
      alldat$outlier[alldat$source == sr & !is.na(alldat$average) & alldat$diff < (diff_mean - 3*diff_std)] = 'Smaller'
      alldat$outlier[is.na(alldat$outlier)] = 'Normal'
      
    }
    
    ggplot(alldat[alldat$average > 0,]) +
      geom_point(aes(x = `6144col`, y = `6144row`, shape = colony, col = gap)) +
      scale_x_continuous(breaks = seq(1,96,1),limits = c(1,96)) +
      scale_y_continuous(breaks = seq(1,64,1),limits = c(64,1),trans = 'reverse') +
      # scale_color_discrete(guide = F) +
      scale_shape_discrete(guide = F) +
      theme_linedraw() +
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank())
    
    
    jpegdat <- rbind(jpegdat, alldat)
    
  }
}



