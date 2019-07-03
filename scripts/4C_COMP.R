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
    
    alldat$average[alldat$colony == 'Gap'] = 0
    
    alldat$outlier <- NULL
    alldat$var <- NULL
    alldat$neigh <- NULL
    alldat$diff <- NULL
    
    for (sr in unique(alldat$source)) {
      temp <- alldat[alldat$source == sr,]
      for (i in seq(1,length(unique(temp$`6144col`)))) {
        col <- unique(temp$`6144col`)[i]
        lf <- tail(temp$`6144col`[temp$`6144col` < col & temp$orf_name == 'BF_control' & !is.na(temp$orf_name)],1)
        rt <- temp$`6144col`[temp$`6144col` > col & temp$orf_name == 'BF_control' & !is.na(temp$orf_name)][1]
        for (ii in seq(1,length(unique(temp$`6144row`[temp$`6144col` == col])))) {
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
      # alldat$outlier[alldat$source == sr & !is.na(alldat$average) & alldat$diff > (diff_mean + 3*diff_std)] = 'Bigger'
      # alldat$outlier[alldat$source == sr & !is.na(alldat$average) & alldat$diff < (diff_mean - 3*diff_std)] = 'Smaller'
      alldat$outlier[alldat$source == sr & !is.na(alldat$average) &
                       alldat$diff > quantile(alldat$diff[alldat$source == sr & alldat$orf_name == 'BF_control'], 0.995, na.rm = T)[[1]]] = 'Bigger'
      alldat$outlier[alldat$source == sr & !is.na(alldat$average) &
                       alldat$diff < quantile(alldat$diff[alldat$source == sr & alldat$orf_name == 'BF_control'], 0.005, na.rm = T)[[1]]] = 'Smaller'
      alldat$outlier[is.na(alldat$outlier)] = 'Normal'
      
    }
    
    ggplot(alldat) +
      geom_point(aes(x = `6144col`, y = `6144row`, shape = colony, size = outlier, col = outlier)) +
      scale_x_continuous(breaks = seq(1,96,1),limits = c(1,96)) +
      scale_y_continuous(breaks = seq(1,64,1),limits = c(64,1),trans = 'reverse') +
      scale_size_manual(values = c('Smaller' = 2, 'Normal' = 1, 'Bigger' = 2)) +
      # scale_color_discrete(guide = F) +
      scale_shape_discrete(guide = F) +
      theme_linedraw() +
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank())
    
    ggplot(alldat) +
      geom_line(aes(x = average, col = outlier), stat = 'density')
    
    
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
      }
    }
    
    alldat$nearBig[alldat$outlier == 'Bigger'] = 'B'
    
    alldat$nearSmall <- 'N'
    for (o in alldat$pos) {
      c = alldat$`6144col`[alldat$pos == o]
      r = alldat$`6144row`[alldat$pos == o]
      if (alldat$outlier[alldat$`6144col` == c & alldat$`6144row` == r] == 'Smaller') {
        alldat$nearSmall[alldat$`6144col` == c - 2 & alldat$`6144row` == r - 2 |
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
        alldat$nearSmall[alldat$`6144col` == c - 1 & alldat$`6144row` == r - 1 |
                         alldat$`6144col` == c & alldat$`6144row` == r - 1 |
                         alldat$`6144col` == c + 1 & alldat$`6144row` == r - 1 |
                         alldat$`6144col` == c - 1 & alldat$`6144row` == r |
                         alldat$`6144col` == c + 1 & alldat$`6144row` == r |
                         alldat$`6144col` == c - 1 & alldat$`6144row` == r + 1 |
                         alldat$`6144col` == c & alldat$`6144row` == r + 1 |
                         alldat$`6144col` == c + 1 & alldat$`6144row` == r + 1] = 'S1'
      }
    }
    
    alldat$nearSmall[alldat$outlier == 'Smaller'] = 'S'
    
    neigh.small <- ggplot(alldat) +
      # geom_point(aes(x = `6144col`, y = `6144row`, shape = colony, col = nearSmall)) +
      geom_point(aes(x = `6144col`, y = `6144row`, shape = colony, col = nearSmall)) +
      scale_x_continuous(breaks = seq(1,96,1),limits = c(1,96)) +
      scale_y_continuous(breaks = seq(1,64,1),limits = c(64,1),trans = 'reverse') +
      labs(title = 'Small Colonies and Neighbors') +
      scale_color_manual(name = '',
                         breaks = c('S','S1','S2'),
                         values = c('N'='#BDBDBD',
                                    'S' = '#FFA000',
                                    'S1' = '#009688',
                                    'S2' = '#673AB7'),
                         labels = c('Small','Deg1','Deg2')) +
      scale_shape_manual(breaks = c('Reference','Query','Gap'),
                         values = c('Reference' = 19,
                                    'Query' = 19,
                                    'Gap' = 0),
                         guide = F) +
      theme_linedraw() +
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            legend.position = 'bottom')
    
    neigh.big <- ggplot(alldat) +
      geom_point(aes(x = `6144col`, y = `6144row`, shape = colony, col = nearBig)) +
      scale_x_continuous(breaks = seq(1,96,1),limits = c(1,96)) +
      scale_y_continuous(breaks = seq(1,64,1),limits = c(64,1),trans = 'reverse') +
      labs(title = 'Big Colonies and Neighbors') +
      scale_color_manual(name = '',
                         breaks = c('B','B1','B2'),
                         values = c('N'='#BDBDBD',
                                    'B' = '#FFA000',
                                    'B1' = '#009688',
                                    'B2' = '#673AB7'),
                         labels = c('Big','Deg1','Deg2')) +
      scale_shape_manual(breaks = c('Reference','Query','Gap'),
                         values = c('Reference' = 19,
                                    'Query' = 19,
                                    'Gap' = 0),
                         guide = F) +
      theme_linedraw() +
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            legend.position = 'bottom')
    
    ggplot(alldat[alldat$orf_name == 'BF_control' & alldat$outlier != 'Bigger' & !is.na(alldat$average),]) +
      geom_line(aes(x = fitness, col = nearBig), stat = 'density') +
      geom_line(data = alldat[alldat$outlier == 'Bigger' & !is.na(alldat$average),],
               aes(x = fitness, col = 'B'), stat = 'density') +
      facet_wrap(.~source)
    
    ggplot(alldat[alldat$orf_name == 'BF_control' & !is.na(alldat$average),]) +
      geom_line(aes(x = average, col = nearSmall), stat = 'density') +
      facet_wrap(.~source) +
      coord_cartesian(ylim = c(0,0.1))
    
    ggplot(alldat) +
      geom_point(aes(x = average, y = fitness, col = nearSmall))
    
    jpegdat <- rbind(jpegdat, alldat)
    
  }
}

median(alldat$average[alldat$nearSmall == 'N'],na.rm = T)


