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
expt_name = '4C3_96R'
expt = 'FS1-96R'
out_path = 'figs/comp/';
density = 6144;

##### CHECK POSITION WISE VARIABILITY
conn <- initialize.sql("saurin_test")

tablename_jpeg = sprintf('%s_%d_JPEG',expt_name,density);
tablename_jpeg_mca = sprintf('%s_MCA_%d_JPEG',expt_name,density);
tablename_fit = sprintf('%s_%d_FITNESS',expt_name,density);
tablename_p2o = '4C3_96R_pos2orf_name';
tablename_bpos = '4C3_borderpos';

p2c_info = NULL
p2c_info[1] = '4C3_96R_pos2coor6144'
p2c_info[2] = '6144plate'
p2c_info[3] = '6144col'
p2c_info[4] = '6144row'

hours = dbGetQuery(conn, sprintf('select distinct hours from %s order by hours asc', tablename_fit))
n_plates = dbGetQuery(conn, sprintf('select distinct %s from %s a order by %s asc',
                                    p2c_info[2],
                                    p2c_info[1],
                                    p2c_info[2]))

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
    
    # alldat$average[alldat$colony == 'Gap'] = 0

    alldat$outlier <- NULL
    alldat$var <- NULL
    alldat$neigh <- NULL
    alldat$diff <- NULL

    # temp <- alldat[alldat$orf_name == 'BF_control' & !is.na(alldat$orf_name),]
    temp <- alldat
    # cnt <- 0
    for (i in seq(1,dim(temp)[1])) {
      col <- temp$`6144col`[i]
      row <- temp$`6144row`[i]
      
      lf <- tail(temp$`6144col`[temp$`6144row` == row & temp$`6144col` < col & temp$orf_name == 'BF_control' & !is.na(temp$orf_name)],1)
      rt <- temp$`6144col`[temp$`6144row` == row & temp$`6144col` > col & temp$orf_name == 'BF_control' & !is.na(temp$orf_name)][1]
      if (length(lf) == 0) {
        up = NaN
      }
      if (length(rt) == 0) {
        up = NaN
      }
      
      up <- tail(temp$`6144row`[temp$`6144col` == col & temp$`6144row` < row & temp$orf_name == 'BF_control' & !is.na(temp$orf_name)],1)
      dw <- temp$`6144row`[temp$`6144col` == col & temp$`6144row` > row & temp$orf_name == 'BF_control' & !is.na(temp$orf_name)][1]
      if (length(up) == 0) {
        up = NaN
      }
      if (length(dw) == 0) {
        up = NaN
      }
      
      a <- alldat$average[alldat$`6144col` == col & alldat$`6144row` == row]
      u <- mean(alldat$average[alldat$`6144col` == col & alldat$`6144row` == up], na.rm = T)
      d <- mean(alldat$average[alldat$`6144col` == col & alldat$`6144row` == dw], na.rm = T)
      l <- mean(alldat$average[alldat$`6144col` == lf & alldat$`6144row` == row], na.rm = T)
      r <- mean(alldat$average[alldat$`6144col` == rt & alldat$`6144row` == row], na.rm = T)
      alldat$var[alldat$`6144col` == col & alldat$`6144row` == row] <- sd(c(a,u,d,l,r),na.rm = T)/median(c(a,u,d,l,r),na.rm = T)
      alldat$neigh[alldat$`6144col` == col & alldat$`6144row` == row] <-  median(c(u,d,l,r),na.rm = T)
      # alldat$wt_neigh[alldat$`6144col` == col & alldat$`6144row` == row] <-
      #   weighted.mean(c(u,d,l,r), c(1/abs(row-up),1/abs(row-dw),1/abs(col-lf),1/abs(col-rt)),na.rm = T)
      alldat$diff[alldat$`6144col` == col & alldat$`6144row` == row] <- a - median(c(u,d,l,r),na.rm = T)
      
    }
    
    diff_std <- sd(alldat$diff,na.rm = T)
    diff_mean <- mean(alldat$diff,na.rm = T)
    alldat$coef <- alldat$average/alldat$neigh
    
    med <- median(alldat$coef, na.rm = T)
    lim.low <- median(alldat$coef, na.rm = T) - 3*sd(alldat$coef, na.rm = T)
    lim.hig <- median(alldat$coef, na.rm = T) + 3*sd(alldat$coef, na.rm = T)
    
    alldat$outlier[alldat$coef > lim.hig & !is.na(alldat$coef)] = 'Bigger'
    alldat$outlier[alldat$coef < lim.low & !is.na(alldat$coef)] = 'Smaller'
    alldat$outlier[is.na(alldat$outlier)] = 'Normal'
    
    ggplot() +
      # geom_histogram(data = alldat, aes(x = coef))
      geom_point(data = alldat[alldat$colony == 'Gap',], aes(x = `6144col`, y = `6144row`), col = 'Black') +
      geom_point(data = alldat[alldat$coef > lim.hig,], aes(x = `6144col`, y = `6144row`), col = 'Red') +
      geom_point(data = alldat[alldat$coef < lim.low,], aes(x = `6144col`, y = `6144row`), col = 'Green')

    
    alldat$nearBig <- 'N'
    for (o in alldat$pos) {
      c = alldat$`6144col`[alldat$pos == o]
      r = alldat$`6144row`[alldat$pos == o]
      if (alldat$outlier[alldat$`6144col` == c & alldat$`6144row` == r] == 'Bigger') {
        alldat$nearBig[alldat$`6144col` == c - 2 & alldat$`6144row` == r - 2 |
                     alldat$`6144col` == c - 1 & alldat$`6144row` == r - 2 |
                     alldat$`6144col` == c & alldat$`6144row` == r - 2 |
                     alldat$`6144col` == c + 1 & alldat$`6144row` == r - 2 |
                     alldat$`6144col` == c + 2 & alldat$`6144row` == r - 2 |
                     alldat$`6144col` == c - 2 & alldat$`6144row` == r - 1 |
                     alldat$`6144col` == c + 2 & alldat$`6144row` == r - 1 |
                     alldat$`6144col` == c - 2 & alldat$`6144row` == r |
                     alldat$`6144col` == c + 2 & alldat$`6144row` == r |
                     alldat$`6144col` == c - 2 & alldat$`6144row` == r + 1 |
                     alldat$`6144col` == c + 2 & alldat$`6144row` == r + 1 |
                     alldat$`6144col` == c - 2 & alldat$`6144row` == r + 2 |
                     alldat$`6144col` == c - 1 & alldat$`6144row` == r + 2 |
                     alldat$`6144col` == c & alldat$`6144row` == r + 2 |
                     alldat$`6144col` == c + 1 & alldat$`6144row` == r + 2 |
                     alldat$`6144col` == c + 2 & alldat$`6144row` == r + 2] = 'B2'
        # alldat$cor[alldat$`6144col` == c - 2 & alldat$`6144row` == r - 2 |
        #              alldat$`6144col` == c - 1 & alldat$`6144row` == r - 2 |
        #              alldat$`6144col` == c & alldat$`6144row` == r - 2 |
        #              alldat$`6144col` == c + 1 & alldat$`6144row` == r - 2 |
        #              alldat$`6144col` == c + 2 & alldat$`6144row` == r - 2 |
        #              alldat$`6144col` == c - 2 & alldat$`6144row` == r - 1 |
        #              alldat$`6144col` == c + 2 & alldat$`6144row` == r - 1 |
        #              alldat$`6144col` == c - 2 & alldat$`6144row` == r |
        #              alldat$`6144col` == c + 2 & alldat$`6144row` == r |
        #              alldat$`6144col` == c - 2 & alldat$`6144row` == r + 1 |
        #              alldat$`6144col` == c + 2 & alldat$`6144row` == r + 1 |
        #              alldat$`6144col` == c - 2 & alldat$`6144row` == r + 2 |
        #              alldat$`6144col` == c - 1 & alldat$`6144row` == r + 2 |
        #              alldat$`6144col` == c & alldat$`6144row` == r + 2 |
        #              alldat$`6144col` == c + 1 & alldat$`6144row` == r + 2 |
        #              alldat$`6144col` == c + 2 & alldat$`6144row` == r + 2] = 
        #   alldat$coef[alldat$`6144col` == c & alldat$`6144row` == r]
      }
    }
    for (o in alldat$pos) {
      c = alldat$`6144col`[alldat$pos == o]
      r = alldat$`6144row`[alldat$pos == o]
      if (alldat$outlier[alldat$`6144col` == c & alldat$`6144row` == r] == 'Bigger') {
        alldat$nearBig[alldat$`6144col` == c - 1 & alldat$`6144row` == r - 1 |
                     alldat$`6144col` == c & alldat$`6144row` == r - 1 |
                     alldat$`6144col` == c + 1 & alldat$`6144row` == r - 1 |
                     alldat$`6144col` == c - 1 & alldat$`6144row` == r |
                     alldat$`6144col` == c + 1 & alldat$`6144row` == r |
                     alldat$`6144col` == c - 1 & alldat$`6144row` == r + 1 |
                     alldat$`6144col` == c & alldat$`6144row` == r + 1 |
                     alldat$`6144col` == c + 1 & alldat$`6144row` == r + 1] = 'B1'
        # alldat$cor[alldat$`6144col` == c - 1 & alldat$`6144row` == r - 1 |
        #              alldat$`6144col` == c & alldat$`6144row` == r - 1 |
        #              alldat$`6144col` == c + 1 & alldat$`6144row` == r - 1 |
        #              alldat$`6144col` == c - 1 & alldat$`6144row` == r |
        #              alldat$`6144col` == c + 1 & alldat$`6144row` == r |
        #              alldat$`6144col` == c - 1 & alldat$`6144row` == r + 1 |
        #              alldat$`6144col` == c & alldat$`6144row` == r + 1 |
        #              alldat$`6144col` == c + 1 & alldat$`6144row` == r + 1] = 
        #   alldat$coef[alldat$`6144col` == c & alldat$`6144row` == r]
      }
    }
    alldat$nearBig[alldat$outlier == 'Bigger'] = 'B'
    # alldat$nearBig[alldat$colony == 'Gap']
    
    ggplot(alldat) +
      geom_point(aes(x = `6144col`, y = `6144row`, col = nearBig, size = coef))
    
    alldat$nearSick <- 'N'
    for (o in alldat$pos) {
      c = alldat$`6144col`[alldat$pos == o]
      r = alldat$`6144row`[alldat$pos == o]
      if (alldat$outlier[alldat$`6144col` == c & alldat$`6144row` == r] == 'Smaller') {
        alldat$nearSick[alldat$`6144col` == c - 2 & alldat$`6144row` == r - 2 |
                         alldat$`6144col` == c - 1 & alldat$`6144row` == r - 2 |
                         alldat$`6144col` == c & alldat$`6144row` == r - 2 |
                         alldat$`6144col` == c + 1 & alldat$`6144row` == r - 2 |
                         alldat$`6144col` == c + 2 & alldat$`6144row` == r - 2 |
                         alldat$`6144col` == c - 2 & alldat$`6144row` == r - 1 |
                         alldat$`6144col` == c + 2 & alldat$`6144row` == r - 1 |
                         alldat$`6144col` == c - 2 & alldat$`6144row` == r |
                         alldat$`6144col` == c + 2 & alldat$`6144row` == r |
                         alldat$`6144col` == c - 2 & alldat$`6144row` == r + 1 |
                         alldat$`6144col` == c + 2 & alldat$`6144row` == r + 1 |
                         alldat$`6144col` == c - 2 & alldat$`6144row` == r + 2 |
                         alldat$`6144col` == c - 1 & alldat$`6144row` == r + 2 |
                         alldat$`6144col` == c & alldat$`6144row` == r + 2 |
                         alldat$`6144col` == c + 1 & alldat$`6144row` == r + 2 |
                         alldat$`6144col` == c + 2 & alldat$`6144row` == r + 2] = 'S2'
      }
    }
    for (o in alldat$pos) {
      c = alldat$`6144col`[alldat$pos == o]
      r = alldat$`6144row`[alldat$pos == o]
      if (alldat$outlier[alldat$`6144col` == c & alldat$`6144row` == r] == 'Smaller') {
        alldat$nearSick[alldat$`6144col` == c - 1 & alldat$`6144row` == r - 1 |
                         alldat$`6144col` == c & alldat$`6144row` == r - 1 |
                         alldat$`6144col` == c + 1 & alldat$`6144row` == r - 1 |
                         alldat$`6144col` == c - 1 & alldat$`6144row` == r |
                         alldat$`6144col` == c + 1 & alldat$`6144row` == r |
                         alldat$`6144col` == c - 1 & alldat$`6144row` == r + 1 |
                         alldat$`6144col` == c & alldat$`6144row` == r + 1 |
                         alldat$`6144col` == c + 1 & alldat$`6144row` == r + 1] = 'S1'
      }
    }
    alldat$nearSick[alldat$outlier == 'Smaller'] = 'S'
    
    ggplot(alldat) +
      geom_point(aes(x = `6144col`, y = `6144row`, col = nearSick, size = coef))
  
    alldat$mca <- alldat$average
    ### NORM 1
    # alldat$mca[alldat$nearBig != 'N'] <- alldat$average[alldat$nearBig != 'N']/alldat$coef[alldat$nearBig != 'N']
    ### NORM 2
    # alldat$mca[alldat$nearBig != 'N' & alldat$nearBig != 'B2' &
    #                  alldat$nearSick != 'N' & alldat$nearSick != 'S2'] <- alldat$average[alldat$nearBig != 'N' & alldat$nearBig != 'B2' &
    #                                                                                        alldat$nearSick != 'N' & alldat$nearSick != 'S2']/
    #   alldat$coef[alldat$nearBig != 'N' & alldat$nearBig != 'B2' & alldat$nearSick != 'N' & alldat$nearSick != 'S2']
    ### NORM 3
    # alldat$mca[alldat$nearBig == 'B' | alldat$nearBig == 'B1' |
    #              alldat$nearSick == 'S' | alldat$nearSick == 'S1'] <- alldat$average[alldat$nearBig == 'B' | alldat$nearBig == 'B1' |
    #                                                                                    alldat$nearSick == 'S' | alldat$nearSick == 'S1']/
    #   alldat$coef[alldat$nearBig == 'B' | alldat$nearBig == 'B1' | alldat$nearSick == 'S' | alldat$nearSick == 'S1']
    ### NORM 4
    alldat$mca[alldat$nearBig == 'B'] <- alldat$average[alldat$nearBig == 'B']/alldat$coef[alldat$nearBig == 'B']
    
    ggplot() +
      # geom_point(aes(x = `6144col`, y = `6144row`, col = outlier))
      geom_line(data = alldat, aes(x = average, col = nearBig), stat = 'density') +
      geom_line(data = alldat, aes(x = mca, col = nearBig), stat = 'density', lwd = 1.2)

    # ggplot(alldat) +
    #   geom_line(aes(x = mca), stat = 'density', col = 'Red') +
    #   geom_line(aes(x = average), stat = 'density', col = 'Blue') +
    #   coord_cartesian(xlim = c(200, 600))
    #   
    
    jpegdat <- rbind(jpegdat, alldat)
  }
}

jpegdat <- data.frame(jpegdat$pos, jpegdat$hours, jpegdat$mca)
colnames(jpegdat) <- c('pos','hours','average')
dbWriteTable(conn, tablename_jpeg_mca, jpegdat, overwrite = T)

ggplot(alldat[!is.na(alldat$average),]) +
  # geom_point(aes(x = average, y =  neigh, col = nearSick)) +
  geom_line(aes(x = average, col = nearSick), stat = 'density', lwd = 1.2) +
  # facet_wrap(.~source) +
  # labs(title = 'Small Colonies and Neighbors',
  #      x = 'Colony Size (pix)',
  #      y = 'Density') +
  scale_color_manual(name = '',
                     breaks = c('N','S','S1','S2'),
                     values = c('N'='#BDBDBD',
                                'S' = '#FFA000',
                                'S1' = '#009688',
                                'S2' = '#673AB7'),
                     labels = c('Rest','Small','Deg1','Deg2')) +
  theme_linedraw() +
  theme(legend.position = 'bottom')

ggplot(alldat[!is.na(alldat$average),]) +
  geom_line(aes(x = average, col = nearBig), stat = 'density', lwd = 1.2) +
  facet_wrap(.~source) +
  labs(title = 'Big Colonies and Neighbors',
       x = 'Colony Size (pix)',
       y = 'Density') +
  scale_color_manual(name = '',
                     breaks = c('N','B','B1','B2'),
                     values = c('N'='#BDBDBD',
                                'B' = '#FFA000',
                                'B1' = '#009688',
                                'B2' = '#673AB7'),
                     labels = c('Rest','Big','Deg1','Deg2')) +
  theme_linedraw() +
  theme(legend.position = 'bottom')

