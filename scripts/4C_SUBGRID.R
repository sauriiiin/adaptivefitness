##### FITNESS IN SUBGRIDS
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 07/19/2019
##### Dividing the plates in subgrids and comparing fitness and raw pixel distributions

##### INITIALIZE
library(RMariaDB)
library(ggplot2)
library(ggExtra)
library(gridExtra)
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
# tablename_p2o = '4C3_pos2orf_name1';
# tablename_bpos = '4C3_borderpos';

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
neigh <- matrix(ncol = 24, nrow = 6144*2)
neigh_sr <- matrix(ncol = 24, nrow = 6144*2)
i <- 1
for (pl in sort(unique(p2c$`6144plate`))) {
  for (sr in sort(unique(p2c$source))) {
    temp <- p2c[p2c$source == sr & p2c$`6144plate` == pl,]
    for (c in sort(unique(temp$`6144col`))) {
      for (r in sort(unique(temp$`6144row`[temp$`6144col` == c]))) {
        pos <- rbind(pos, temp$pos[temp$`6144col` == c & temp$`6144row` == r])
        # len <- length(temp$pos[temp$`6144col` == c - 2 & temp$`6144row` %in% c(r - 2, r - 1, r, r + 1, r + 2) |
        #                          temp$`6144col` == c - 1 & temp$`6144row` %in% c(r - 2, r - 1, r, r + 1, r + 2) |
        #                          temp$`6144col` == c & temp$`6144row` %in% c(r - 2, r - 1, r + 1, r + 2) |
        #                          temp$`6144col` == c + 1 & temp$`6144row` %in% c(r - 2, r - 1, r, r + 1, r + 2) |
        #                          temp$`6144col` == c + 2 & temp$`6144row` %in% c(r - 2, r - 1, r, r + 1, r + 2)])
        # neigh[i,1:len] <- temp$pos[temp$`6144col` == c - 2 & temp$`6144row` %in% c(r - 2, r - 1, r, r + 1, r + 2) |
        #                              temp$`6144col` == c - 1 & temp$`6144row` %in% c(r - 2, r - 1, r, r + 1, r + 2) |
        #                              temp$`6144col` == c & temp$`6144row` %in% c(r - 2, r - 1, r + 1, r + 2) |
        #                              temp$`6144col` == c + 1 & temp$`6144row` %in% c(r - 2, r - 1, r, r + 1, r + 2) |
        #                              temp$`6144col` == c + 2 & temp$`6144row` %in% c(r - 2, r - 1, r, r + 1, r + 2)]
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
tempdat$average[is.na(tempdat$orf_name)] <- 0
# for (hr in sort(unique(alldat$hours))) {
  for (i in seq(1,dim(grids)[1])) {
   tempdat$neigh[tempdat$hours == hr & tempdat$pos == grids[i]] <-
     mean(tempdat$average[tempdat$hours == hr & tempdat$pos %in% grids[i,2:25]], na.rm = T)
   tempdat$neigh_sr[tempdat$hours == hr & tempdat$pos == grids_sr[i]] <-
     mean(tempdat$average[tempdat$hours == hr & tempdat$pos %in% grids_sr[i,2:25]], na.rm = T)
  }
# }

tempdat$neigh[is.na(tempdat$average) & !is.na(tempdat$orf_name)] <- NA
tempdat$neigh_sr[is.na(tempdat$average) & !is.na(tempdat$orf_name)] <- NA

tempdat$score <- tempdat$average/(tempdat$neigh + tempdat$neigh_sr)

for (i in seq(1,dim(grids)[1])) {
  tempdat$score_neigh[tempdat$hours == hr & tempdat$pos == grids[i]] <-
    (quantile(tempdat$score[tempdat$hours == hr & tempdat$pos %in% grids[i,2:25]], 0.5, na.rm = T)[[1]] +
       quantile(tempdat$score[tempdat$hours == hr & tempdat$pos %in% grids_sr[i,2:25]], 0.5, na.rm = T)[[1]])/2
}

ggplot(tempdat) +
  geom_histogram(aes(x = score), binwidth = .02) +
  labs(title = sprintf('%s: Competition', expt),
       x = 'Competition Score',
       y = 'Frequency') +
  theme_linedraw()
ggsave(sprintf('%s%s_COMP_SCORE.png',out_path,expt_name),
       height = 7, width = 7)

md <- mad(tempdat$score, na.rm =T)
ll <- median(tempdat$score, na.rm =T) - 2*md
ul <- median(tempdat$score, na.rm =T) + 2*md

# score.med <- median(tempdat$score, na.rm = T)
ggplot(tempdat) +
  # geom_line(aes(x = score), stat = 'density', binwidth = 30) +
  # geom_vline(xintercept = c(ll, ul)) +
  # geom_point(data = tempdat[tempdat$`6144plate` == 1 & tempdat$score > ul,],
  #            aes(x = score, y = 0), col = 'Red')
  geom_line(aes(x = average), stat = 'density', lwd = 1.2) +
  geom_point(data = tempdat[tempdat$`6144plate` == 1 & tempdat$score > ul,],
             aes(x = average, y = 0), col = 'Red') +
  geom_point(data = tempdat[tempdat$`6144plate` == 1 & tempdat$score > ul,],
             aes(x = average * score_neigh/score, y = 0.0001), col = 'Blue') +
  geom_point(data = tempdat[tempdat$`6144plate` == 1 & tempdat$score < ll,],
             aes(x = average, y = 0.0002), col = 'darkgreen') +
  geom_point(data = tempdat[tempdat$`6144plate` == 1 & tempdat$score < ll,],
             aes(x = average * score_neigh/score, y = 0.0003), col = 'Orange')
  # geom_point(data = tempdat[tempdat$`6144plate` == 1 & tempdat$score < ul,],
  #            aes(x = average, y = score), col = 'Red') + #average * score_neigh/score
  # geom_point(data = tempdat[tempdat$`6144plate` == 1 & tempdat$score > ul,],
  #            aes(x = average, y = score), col = 'Blue') #average * score_neigh/score

ggplot(tempdat[tempdat$`6144plate` == 1,]) +
  geom_point(aes(x = `6144col`, y = `6144row`, shape = colony), col ='lightblue') +
  geom_point(data = tempdat[tempdat$`6144plate` == 1 &
                            tempdat$score > ul,],
           aes(x = `6144col`, y = `6144row`, shape = colony), col = 'Red') +
  geom_point(data = tempdat[tempdat$`6144plate` == 1 &
                              tempdat$score < ll,],
             aes(x = `6144col`, y = `6144row`, shape = colony), col = 'darkgreen') +
  geom_point(data = tempdat[tempdat$`6144plate` == 1 &
                              tempdat$colony == 'Gap',],
             aes(x = `6144col`, y = `6144row`, shape = colony), col = 'Black')

tempdat$average_cc <- tempdat$average
tempdat$average_cc[tempdat$`6144plate` == 1 & !is.na(tempdat$score) & tempdat$score > ul] <-
  tempdat$average[tempdat$`6144plate` == 1 & !is.na(tempdat$score) & tempdat$score > ul] *
  tempdat$score_neigh[tempdat$`6144plate` == 1 & !is.na(tempdat$score) & tempdat$score > ul]/
  tempdat$score[tempdat$`6144plate` == 1 & !is.na(tempdat$score) & tempdat$score > ul]
# tempdat$average_cc[tempdat$`6144plate` == 1 & !is.na(tempdat$score) & tempdat$score < ll] <-
#   tempdat$average[tempdat$`6144plate` == 1 & !is.na(tempdat$score) & tempdat$score < ll] *
#   tempdat$score_neigh[tempdat$`6144plate` == 1 & !is.na(tempdat$score) & tempdat$score < ll]/
#   tempdat$score[tempdat$`6144plate` == 1 & !is.na(tempdat$score) & tempdat$score < ll]

ggplot(tempdat) +
  # geom_line(aes(x = average, col = colony), stat = 'density', lwd = 1.2) +
  geom_line(aes(x = average_cc, col = colony), stat = 'density', lwd = 1.2) +
  coord_cartesian(xlim = c(200,600),
                  ylim = c(0, 0.01))


##### PUTTING IT ALL TOGETHER
compdat <- data.frame()

for (hr in sort(unique(alldat$hours))) {
  for (pl in sort(unique(alldat$`6144plate`[alldat$hours == hr]))) {
    tempdat <- alldat[alldat$hours == hr & alldat$`6144plate` == pl,]
    tempdat$average[is.na(tempdat$orf_name)] <- 0

    for (i in seq(1,dim(grids)[1])) {
      tempdat$neigh[tempdat$pos == grids[i]] <- mean(tempdat$average[tempdat$pos %in% grids[i,2:25]], na.rm = T)
      tempdat$neigh_sr[tempdat$pos == grids_sr[i]] <- mean(tempdat$average[tempdat$pos %in% grids_sr[i,2:25]], na.rm = T)
    }
    
    tempdat$neigh[is.na(tempdat$average) & !is.na(tempdat$orf_name)] <- NA
    tempdat$neigh_sr[is.na(tempdat$average) & !is.na(tempdat$orf_name)] <- NA
    tempdat$score <- tempdat$average/(tempdat$neigh + tempdat$neigh_sr)
    
    md <- mad(tempdat$score, na.rm =T)
    ll <- median(tempdat$score, na.rm =T) - 3*md
    ul <- median(tempdat$score, na.rm =T) + 3*md
    
    for (i in seq(1,dim(grids)[1])) {
      tempdat$score_neigh[tempdat$pos == grids[i]] <-
        quantile(tempdat$score[tempdat$pos %in% grids[i,2:25]], 0.5, na.rm = T)[[1]]
    }
    
    # ggplot(tempdat) +
    #   # geom_line(aes(x = score), stat = 'density', binwidth = 30) +
    #   # geom_vline(xintercept = c(ll, ul)) +
    #   # geom_point(data = tempdat[tempdat$`6144plate` == 1 & tempdat$score > ul,],
    #   #            aes(x = score, y = 0), col = 'Red')
    #   geom_line(aes(x = average), stat = 'density', lwd = 1.2) +
    #   geom_point(data = tempdat[tempdat$`6144plate` == 1 & tempdat$score > ul,],
    #              aes(x = average, y = 0), col = 'Red') +
    #   geom_point(data = tempdat[tempdat$`6144plate` == 1 & tempdat$score > ul,],
    #              aes(x = average * score_neigh/score, y = 0.0001), col = 'Blue')
    # # geom_point(data = tempdat[tempdat$`6144plate` == 1 & tempdat$score > ul,],
    # #          aes(x = average, y = average * score_neigh/score), col = 'Blue')
    # 
    # ggplot(tempdat[tempdat$`6144plate` == 1,]) +
    #   geom_point(aes(x = `6144col`, y = `6144row`, shape = colony), col ='Red') +
    #   geom_point(data = tempdat[tempdat$`6144plate` == 1 &
    #                               tempdat$colony == 'Gap',],
    #              aes(x = `6144col`, y = `6144row`, shape = colony), col = 'Blue') +
    #   geom_point(data = tempdat[tempdat$`6144plate` == 1 &
    #                               tempdat$score > ul,],
    #              aes(x = `6144col`, y = `6144row`, shape = colony))
    
    tempdat$average[tempdat$`6144plate` == 1 & !is.na(tempdat$score) & tempdat$score > ul] <-
      tempdat$average[tempdat$`6144plate` == 1 & !is.na(tempdat$score) & tempdat$score > ul] *
      tempdat$score_neigh[tempdat$`6144plate` == 1 & !is.na(tempdat$score) & tempdat$score > ul]/
      tempdat$score[tempdat$`6144plate` == 1 & !is.na(tempdat$score) & tempdat$score > ul]
    
    # ggplot(tempdat) +
    #   geom_line(aes(x = average), stat = 'density')
    
    compdat <- rbind(compdat, tempdat)
  }
}

compdat$average[is.na(compdat$orf_name)] <- NA

newdat <- data.frame(compdat$pos, compdat$hours, compdat$average)
colnames(newdat) <- c('pos','hours','average')
dbWriteTable(conn, tablename_jpeg_cc, newdat, overwrite = T)



