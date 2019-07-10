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
    
    alldat$average[alldat$colony == 'Gap'] = 0
    
    alldat$health <- NULL
    for (sr in unique(alldat$source)) {
      alldat$health[alldat$colony != 'Refernce' &
                      alldat$source == sr &
                      alldat$average < quantile(alldat$average[alldat$colony == 'Reference' &
                                                                 alldat$source == sr], 0.005, na.rm = T)[[1]] &
                      !is.na(alldat$average)] = 'Sick'
      alldat$health[alldat$colony != 'Reference' &
                      alldat$source == sr &
                      alldat$average > quantile(alldat$average[alldat$colony == 'Reference' &
                                                                 alldat$source == sr], 0.995, na.rm = T)[[1]] &
                      !is.na(alldat$average)] = 'Healthy'
      # alldat$health[is.na(alldat$health[alldat$source == sr & alldat$colony != 'Reference'])] = 'Normal'
    }
    
    alldat$health[is.na(alldat$health) & alldat$colony != 'Reference'] = 'Normal'
    # alldat$health[alldat$colony == 'Gap'] = 'Gap'
    alldat$health[alldat$colony == 'Reference'] = 'Reference'
    
    # ggplot(alldat) +
    #   geom_point(aes(x = `6144col`, y = `6144row`, col = health))

    alldat$outlier <- NULL
    alldat$var <- NULL
    alldat$neigh <- NULL
    alldat$diff <- NULL

    temp <- alldat[alldat$orf_name == 'BF_control' & !is.na(alldat$orf_name),]
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
          alldat$var[alldat$`6144col` == col & alldat$`6144row` == row] <- sd(c(a,u,d,l,r),na.rm = T)/median(c(a,u,d,l,r),na.rm = T)
          alldat$neigh[alldat$`6144col` == col & alldat$`6144row` == row] <-  median(c(u,d,l,r),na.rm = T)
          alldat$diff[alldat$`6144col` == col & alldat$`6144row` == row] <- a - median(c(u,d,l,r),na.rm = T)
        }
        # cnt <- cnt + 1
      }
    }
    diff_std <- sd(alldat$diff,na.rm = T)
    diff_mean <- mean(alldat$diff,na.rm = T)
    alldat$outlier[!is.na(alldat$average) &
                     alldat$var > quantile(alldat$var[alldat$orf_name == 'BF_control'], 0.95, na.rm = T)[[1]] &
                     alldat$diff > 0] = 'Bigger'
    alldat$outlier[!is.na(alldat$average) &
                     alldat$var > quantile(alldat$var[alldat$orf_name == 'BF_control'], 0.95, na.rm = T)[[1]] &
                     alldat$diff < 0] = 'Smaller'
    alldat$outlier[is.na(alldat$outlier)] = 'Normal'
    
    alldat$coef[alldat$outlier == 'Bigger'] <-
      mean(alldat$average[alldat$outlier == 'Bigger']/alldat$neigh[alldat$outlier == 'Bigger'], na.rm = T)
    alldat$coef[alldat$outlier == 'Smaller'] <-
      mean(alldat$average[alldat$outlier == 'Smaller']/alldat$neigh[alldat$outlier == 'Smaller'], na.rm = T)
    alldat$coef[alldat$outlier == 'Normal'] <-
      mean(alldat$average[alldat$outlier == 'Normal']/alldat$neigh[alldat$outlier == 'Normal'], na.rm = T)

    # ggplot(alldat) +
    #   geom_point(aes(x = `6144col`, y = `6144row`, shape = colony, size = outlier, col = outlier)) +
    #   scale_x_continuous(breaks = seq(1,96,1),limits = c(1,96)) +
    #   scale_y_continuous(breaks = seq(1,64,1),limits = c(64,1),trans = 'reverse') +
    #   scale_size_manual(values = c('Smaller' = 2, 'Normal' = 2, 'Bigger' = 2)) +
    #   scale_shape_manual(breaks = c('Reference','Query','Gap'),
    #                      values = c('Reference' = 19,
    #                                 'Query' = 19,
    #                                 'Gap' = 0),
    #                      guide = F) +
    #   # scale_color_discrete(guide = F) +
    #   theme_linedraw() +
    #   theme(axis.title = element_blank(),
    #         axis.text = element_blank(),
    #         axis.ticks = element_blank())
    # 
    # ggplot(alldat) +
    #   geom_line(aes(x = average, col = outlier), stat = 'density') +
    #   geom_line(aes(x = neigh, col = outlier), stat = 'density', linetype = 'dotted')
    # 
    # ggplot(alldat) +
    #   # geom_point(aes(x = `6144col`, y = `6144row`, col = coef))
    #   geom_line(aes(x = coef), stat = 'density') +
    #   labs(title = sprintf('Median = %0.3f',median(alldat$coef, na.rm = T)))
    
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
    
    ggplot(alldat) +
      geom_point(aes(x = `6144col`, y = `6144row`, col = nearBig))
    
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
      geom_point(aes(x = `6144col`, y = `6144row`, col = outlier))
  
    # alldat$mca <- alldat$average
    # alldat$mca[alldat$nearBig == 'B'] <- alldat$average[alldat$nearBig == 'B'] * median(alldat$average[alldat$nearBig == 'N'], na.rm = T)/
    #   median(alldat$average[alldat$nearBig == 'B'], na.rm = T)
    # alldat$mca[alldat$nearBig == 'B1'] <- alldat$average[alldat$nearBig == 'B1'] * median(alldat$average[alldat$nearBig == 'N'], na.rm = T)/
    #   median(alldat$average[alldat$nearBig == 'B1'], na.rm = T)
    # alldat$mca[alldat$nearBig == 'B2'] <- alldat$average[alldat$nearBig == 'B2'] * median(alldat$average[alldat$nearBig == 'N'], na.rm = T)/
    #   median(alldat$average[alldat$nearBig == 'B2'], na.rm = T)
    
    # for (sr in unique(alldat$source)) {
    #   alldat$mca[alldat$nearSick == 'S1' & alldat$source == sr] <-
    #     alldat$average[alldat$nearSick == 'S1' & alldat$source == sr] * median(alldat$average[alldat$nearSick == 'N' & alldat$source == sr], na.rm = T)/
    #     median(alldat$average[alldat$nearSick == 'S1' & alldat$source == sr], na.rm = T)
    #   alldat$mca[alldat$nearSick == 'S2' & alldat$source == sr] <-
    #     alldat$average[alldat$nearSick == 'S2' & alldat$source == sr] * median(alldat$average[alldat$nearSick == 'N' & alldat$source == sr], na.rm = T)/
    #     median(alldat$average[alldat$nearSick == 'S2' & alldat$source == sr], na.rm = T)
    # }
    alldat$mca <- alldat$average/comp_coef #change this
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

