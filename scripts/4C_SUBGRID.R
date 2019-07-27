##### FITNESS IN SUBGRIDS
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 07/19/2019
##### Dividing the plates in subgrids and comparing fitness and raw pixel distributions

##### INITIALIZE
library(RMariaDB)
library(ggplot2)
library(ggExtra)
library(gridExtra)
library(ggpubr)
source("R/functions/initialize.sql.R")

##### GET/SET DATA
expt_name = '4C3_GA3'
expt = 'FS1-GA3'
out_path = 'figs/comp/';
density = 6144;

##### CHECK POSITION WISE VARIABILITY
conn <- initialize.sql("saurin_test")

tablename_jpeg = sprintf('%s_%d_JPEG',expt_name,density);
tablename_jpeg_cc = sprintf('%s_CC_%d_JPEG',expt_name,density);
tablename_fit = sprintf('%s_%d_FITNESS',expt_name,density);
tablename_fit_cc = sprintf('%s_CC_%d_FITNESS',expt_name,density);
tablename_mdfr = sprintf('%s_CC_%d_MDFR',expt_name,density);

p2c_info = NULL
p2c_info[1] = '4C3_pos2coor6144'
p2c_info[2] = '6144plate'
p2c_info[3] = '6144col'
p2c_info[4] = '6144row'

p2c <- dbGetQuery(conn, sprintf('select * from %s
                                order by %s, %s, %s',
                  p2c_info[1], p2c_info[2], p2c_info[3], p2c_info[4]))

# p2c$source[p2c$`6144row`%%2==1 & p2c$`6144col`%%2==1] = '1TL'
# p2c$source[p2c$`6144row`%%2==0 & p2c$`6144col`%%2==1] = '3BL'
# p2c$source[p2c$`6144row`%%2==1 & p2c$`6144col`%%2==0] = '2TR'
# p2c$source[p2c$`6144row`%%2==0 & p2c$`6144col`%%2==0] = '4BR'

p2c$source = 'ALL'

pos <- NULL
neigh <- matrix(ncol = 8, nrow = 6144*2)
neigh_sr <- matrix(ncol = 8, nrow = 6144*2)
i <- 1
for (pl in sort(unique(p2c$`6144plate`))) {
  for (sr in sort(unique(p2c$source))) {
    temp <- p2c[p2c$source == sr & p2c$`6144plate` == pl,]
    for (c in sort(unique(temp$`6144col`))) {
      for (r in sort(unique(temp$`6144row`[temp$`6144col` == c]))) {
        pos <- rbind(pos, temp$pos[temp$`6144col` == c & temp$`6144row` == r])

        len <- length(temp$pos[temp$`6144col` == c - 1 & temp$`6144row` %in% c(r - 1, r, r + 1) |
                                 temp$`6144col` == c & temp$`6144row` %in% c(r - 1, r + 1) |
                                 temp$`6144col` == c + 1 & temp$`6144row` %in% c(r - 1, r, r + 1)])
        neigh[i,1:len] <- temp$pos[temp$`6144col` == c - 1 & temp$`6144row` %in% c(r - 1, r, r + 1) |
                                     temp$`6144col` == c & temp$`6144row` %in% c(r - 1, r + 1) |
                                     temp$`6144col` == c + 1 & temp$`6144row` %in% c(r - 1, r, r + 1)]
        len_sr <- length(temp$pos[temp$`6144col` == c - 2 & temp$`6144row` %in% c(r - 2, r, r + 2) |
                                 temp$`6144col` == c & temp$`6144row` %in% c(r - 2, r + 2) |
                                 temp$`6144col` == c + 2 & temp$`6144row` %in% c(r - 2, r, r + 2)])
        neigh_sr[i,1:len_sr] <- temp$pos[temp$`6144col` == c - 2 & temp$`6144row` %in% c(r - 2, r, r + 2) |
                                           temp$`6144col` == c & temp$`6144row` %in% c(r - 2, r + 2) |
                                           temp$`6144col` == c + 2 & temp$`6144row` %in% c(r - 2, r, r + 2)]
        i <- i + 1
      }
    } 
  }
}
grids <- cbind(pos, neigh)
grids_sr <- cbind(pos, neigh_sr)

alldat <- dbGetQuery(conn, sprintf('select a.*, b.*
                                  from %s a, %s b
                                  where a.pos = b.pos
                                  order by a.hours, b.%s, b.%s, b.%s',
                                  tablename_fit,
                                  p2c_info[1],p2c_info[2],
                                  p2c_info[3],p2c_info[4]))

alldat$colony[alldat$orf_name == 'BF_control'] = 'Reference'
alldat$colony[alldat$orf_name != 'BF_control'] = 'Query'
alldat$colony[is.na(alldat$orf_name)] = 'Gap'

alldat$source[alldat$`6144row`%%2==1 & alldat$`6144col`%%2==1] = '1TL'
alldat$source[alldat$`6144row`%%2==0 & alldat$`6144col`%%2==1] = '3BL'
alldat$source[alldat$`6144row`%%2==1 & alldat$`6144col`%%2==0] = '2TR'
alldat$source[alldat$`6144row`%%2==0 & alldat$`6144col`%%2==0] = '4BR'

hr = 18
pl = 1
tempdat <- alldat[alldat$hours == hr & alldat$`6144plate` == pl,]
# tempdat$bg[tempdat$`6144col` < 5 | tempdat$`6144col` > 92 | tempdat$`6144row` < 5 | tempdat$`6144row` > 60] <- NA
tempdat$average[is.na(tempdat$orf_name)] <- 0
for (i in seq(1,dim(grids)[1])) {
  tempdat$neigh[tempdat$hours == hr & tempdat$pos == grids[i]] <-
    mean(tempdat$average[tempdat$hours == hr & tempdat$pos %in% grids[i,2:9]], na.rm = T)
  tempdat$neigh_sr[tempdat$hours == hr & tempdat$pos == grids_sr[i]] <-
    mean(tempdat$average[tempdat$hours == hr & tempdat$pos %in% grids_sr[i,2:9]], na.rm = T)
}

tempdat$neigh[is.na(tempdat$average) & !is.na(tempdat$orf_name)] <- NA
tempdat$neigh_sr[is.na(tempdat$average) & !is.na(tempdat$orf_name)] <- NA

tempdat$score <- tempdat$average/(tempdat$neigh + tempdat$neigh_sr)

# for (i in seq(1,dim(grids)[1])) {
#   tempdat$score_neigh[tempdat$hours == hr & tempdat$pos == grids[i]] <-
#     (mean(tempdat$score[tempdat$hours == hr & tempdat$pos %in% grids[i,2:9]], na.rm = T) +
#        mean(tempdat$score[tempdat$hours == hr & tempdat$pos %in% grids_sr[i,2:9]], na.rm = T))/2
# }

ggplot(tempdat[tempdat$colony != 'Gap',]) +
  # geom_point(aes(x = `6144col`, y = `6144row`, col = score_neigh/score))
  # geom_point(aes(x = score, y = score_neigh))
  geom_line(aes(x = score - score_neigh, col = colony), stat = 'density')

md <- mad(tempdat$score, na.rm =T)
ll <- median(tempdat$score, na.rm =T) - 2*md
ul <- median(tempdat$score, na.rm =T) + 2*md

tempdat$outlier <- NULL
tempdat$outlier[tempdat$score < ll & !is.na(tempdat$score)] <- 'Sick'
tempdat$outlier[tempdat$score > ul & !is.na(tempdat$score)] <- 'Healthy'
tempdat$outlier[is.na(tempdat$outlier)] <- 'Normal'

tempdat$nearsick <- NULL
tempdat$nearsick_sr <- NULL
for (p in tempdat$pos[tempdat$outlier == 'Sick']) {
  tempdat$nearsick[tempdat$pos %in% grids[grids[,1] == p, 2:9]] <- 'Y'
  tempdat$nearsick_sr[tempdat$pos %in% grids_sr[grids_sr[,1] == p, 2:9]] <- 'Y'
}
tempdat$nearsick[is.na(tempdat$orf_name)] <- 'N'
tempdat$nearsick[is.na(tempdat$nearsick)] <- 'N'
tempdat$nearsick_sr[is.na(tempdat$orf_name)] <- 'N'
tempdat$nearsick_sr[is.na(tempdat$nearsick_sr)] <- 'N'

tempdat$average_cc <- tempdat$average
tempdat$average_cc[tempdat$nearsick == 'Y'] <- tempdat$average_cc[tempdat$nearsick == 'Y'] *
  median(tempdat$score, na.rm =T)/tempdat$score[tempdat$nearsick == 'Y']
tempdat$average_cc2 <- tempdat$average*median(tempdat$score[tempdat$colony == 'Reference'], na.rm =T)/tempdat$score
# tempdat$fitness_cc2 <- tempdat$fitness*tempdat$score_neigh/tempdat$score

plt.scr <- ggplot(tempdat[tempdat$colony != 'Gap',]) +
  geom_line(aes(x = average/bg, col = colony), stat = 'density', lwd = 1.2) +
  geom_vline(xintercept = quantile(tempdat$average[tempdat$colony == 'Reference']/tempdat$bg[tempdat$colony == 'Reference'], c(0.05,0.95), na.rm = T), col = 'blue') +
  geom_vline(xintercept = quantile(tempdat$average[tempdat$colony == 'Query']/tempdat$bg[tempdat$colony == 'Query'], c(0.05,0.95), na.rm = T), col = 'red') +
  labs(title = sprintf('ES = %0.3f', median(tempdat$average[tempdat$colony == 'Query']/tempdat$bg[tempdat$colony == 'Query'], na.rm = T)/
                         median(tempdat$average[tempdat$colony == 'Reference']/tempdat$bg[tempdat$colony == 'Reference'], na.rm = T)),
       x = 'Raw Pixel Counts',
       y = 'Frequency') +
  theme_linedraw()

plt.pix <- ggplot(tempdat[tempdat$colony != 'Gap',]) +
  geom_line(aes(x = average_cc2/bg, col = colony), stat = 'density', lwd = 1.2) +
  geom_vline(xintercept = quantile(tempdat$average_cc2[tempdat$colony == 'Reference']/tempdat$bg[tempdat$colony == 'Reference'], c(0.05,0.95), na.rm = T), col = 'blue') +
  geom_vline(xintercept = quantile(tempdat$average_cc2[tempdat$colony == 'Query']/tempdat$bg[tempdat$colony == 'Query'], c(0.05,0.95), na.rm = T), col = 'red') +
  labs(title = sprintf('ES = %0.3f', median(tempdat$average_cc2[tempdat$colony == 'Query']/tempdat$bg[tempdat$colony == 'Query'], na.rm = T)/
                         median(tempdat$average_cc2[tempdat$colony == 'Reference']/tempdat$bg[tempdat$colony == 'Reference'], na.rm = T)),
       x = 'CC Pixel Counts',
       y = 'Frequency') +
  theme_linedraw()

annotate_figure(ggarrange(plt.scr, plt.pix, align = "hv",
                          common.legend = T, legend = 'right',
                          ncol = 2, nrow = 1),
                top = sprintf('%s: Competition', expt))
# ggsave(sprintf('%s%s_COMP_SCORE.png',out_path,expt_name),
#        height = 7, width = 7)

ggplot(tempdat[tempdat$`6144plate` == pl,]) +
  geom_point(aes(x = `6144col`, y = `6144row`, shape = colony), col ='lightblue') +
  geom_point(data = tempdat[tempdat$`6144plate` == pl &
                              tempdat$colony == 'Gap',],
             aes(x = `6144col`, y = `6144row`, shape = colony), col = 'Black') +
  # geom_point(data = tempdat[tempdat$`6144plate` == pl &
  #                           tempdat$score > ul,],
  #          aes(x = `6144col`, y = `6144row`, shape = colony), col = 'Red') +
  geom_point(data = tempdat[tempdat$`6144plate` == pl &
                              tempdat$score < ll,],
             aes(x = `6144col`, y = `6144row`, shape = colony), col = 'darkgreen')


##### PUTTING IT ALL TOGETHER
compdat <- data.frame()

for (hr in sort(unique(alldat$hours))) {
  for (pl in sort(unique(alldat$`6144plate`[alldat$hours == hr]))) {
    tempdat <- alldat[alldat$hours == hr & alldat$`6144plate` == pl,]
    tempdat$average[is.na(tempdat$orf_name)] <- 0
    # tempdat$bg[tempdat$`6144col` < 5 | tempdat$`6144col` > 92 | tempdat$`6144row` < 5 | tempdat$`6144row` > 60] <- NA

    for (i in seq(1,dim(grids)[1])) {
      tempdat$neigh[tempdat$pos == grids[i]] <- mean(tempdat$average[tempdat$pos %in% grids[i,2:9]], na.rm = T)
      tempdat$neigh_sr[tempdat$pos == grids_sr[i]] <- mean(tempdat$average[tempdat$pos %in% grids_sr[i,2:9]], na.rm = T)
    }
    
    tempdat$neigh[is.na(tempdat$average) & !is.na(tempdat$orf_name)] <- NA
    tempdat$neigh_sr[is.na(tempdat$average) & !is.na(tempdat$orf_name)] <- NA
    tempdat$score <- tempdat$average/((tempdat$neigh + tempdat$neigh_sr)/2)
    
    for (i in seq(1,dim(grids)[1])) {
      tempdat$score_neigh[tempdat$pos == grids[i]] <-
        (mean(tempdat$score[tempdat$pos %in% grids[i,2:9]], na.rm = T) +
          mean(tempdat$score[tempdat$pos %in% grids_sr[i,2:9]], na.rm = T))/2
    }
    
    md <- mad(tempdat$score[tempdat$orf_name == 'BF_control'], na.rm =T)
    ll <- median(tempdat$score[tempdat$orf_name == 'BF_control'], na.rm =T) - 3*md
    ul <- median(tempdat$score[tempdat$orf_name == 'BF_control'], na.rm =T) + 3*md
    
    ggplot(tempdat[tempdat$`6144plate` == pl,]) +
      geom_point(aes(x = `6144col`, y = `6144row`, shape = colony), col ='yellow') +
      geom_point(data = tempdat[tempdat$`6144plate` == pl &
                                  tempdat$colony == 'Gap',],
                 aes(x = `6144col`, y = `6144row`, shape = colony), col = 'Blue') +
      geom_point(data = tempdat[tempdat$`6144plate` == pl &
                                  tempdat$score > ul,],
                 aes(x = `6144col`, y = `6144row`, shape = colony), col = 'Green') +
      geom_point(data = tempdat[tempdat$`6144plate` == pl &
                                  tempdat$score < ll,],
                 aes(x = `6144col`, y = `6144row`, shape = colony), col = 'red')
    
    tempdat$outlier <- NULL
    tempdat$outlier[tempdat$score < ll & !is.na(tempdat$score)] <- 'Sick'
    tempdat$outlier[tempdat$score > ul & !is.na(tempdat$score)] <- 'Healthy'
    tempdat$outlier[is.na(tempdat$outlier)] <- 'Normal'
    
    tempdat$reallysick <- NULL
    for (p in tempdat$pos[tempdat$outlier == 'Sick']) {
      if (sum(tempdat$outlier[tempdat$pos %in% grids[grids[,1] == p, 2:9] | tempdat$pos %in% grids_sr[grids_sr[,1] == p, 2:9]] == 'Healthy', na.rm = T) > 0) {
        tempdat$reallysick[tempdat$pos == p] <- 'Y'
      }
    }
    tempdat$reallysick[is.na(tempdat$reallysick) & tempdat$outlier == 'Sick'] <- 'N'
    tempdat$reallysick[is.na(tempdat$reallysick)] <- 'Normal'
    
    tempdat$nearsick <- NULL
    # tempdat$nearsick_sr <- NULL
    for (p in tempdat$pos[tempdat$reallysick == 'Y']) {
      tempdat$nearsick[tempdat$pos %in% grids[grids[,1] == p, 2:9]] <- 'Y'
      # tempdat$nearsick_sr[tempdat$pos %in% grids_sr[grids_sr[,1] == p, 2:9]] <- 'Y'
    }
    tempdat$nearsick[is.na(tempdat$orf_name)] <- 'N'
    tempdat$nearsick[is.na(tempdat$nearsick)] <- 'N'
    # # tempdat$nearsick_sr[is.na(tempdat$orf_name)] <- 'N'
    # # tempdat$nearsick_sr[is.na(tempdat$nearsick_sr)] <- 'N'
    
    tempdat$average_cc2 <- tempdat$average
    tempdat$average_cc2[tempdat$nearsick == 'Y'] <- tempdat$average_cc2[tempdat$nearsick == 'Y'] * 
      tempdat$score_neigh[tempdat$nearsick == 'Y']/tempdat$score[tempdat$nearsick == 'Y']
    
    # tempdat$fitness_cc2 <- tempdat$fitness * median(tempdat$score, na.rm = T)/tempdat$score
    
    # tempdat$modifier <- median(tempdat$score, na.rm = T)/tempdat$score
    
    compdat <- rbind(compdat, tempdat)
  }
}

compdat$average[is.na(compdat$orf_name)] <- NA
compdat$average_cc2[is.na(compdat$orf_name)] <- NA
compdat$modifier[is.na(compdat$orf_name)] <- NA

jpegdat <- data.frame(compdat$pos, compdat$hours, compdat$average, compdat$average_cc2)
colnames(jpegdat) <- c('pos','hours','average_raw', 'average')
dbWriteTable(conn, tablename_jpeg_cc, jpegdat, overwrite = T)

# fitdat <- data.frame(compdat$orf_name, compdat$pos, compdat$hours, compdat$bg, compdat$average, compdat$average_cc2, compdat$fitness_cc2)
# colnames(fitdat) <- c('orf_name','pos','hours','bg','average','average_cc2','fitness')
# dbWriteTable(conn, tablename_fit_cc, fitdat, overwrite = T)
# 
# mdfrdat <- data.frame(compdat$pos, compdat$hours, compdat$modifier)
# colnames(mdfrdat) <- c('pos','hours', 'mdfr')
# dbWriteTable(conn, tablename_mdfr, mdfrdat, overwrite = T)

