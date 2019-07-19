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
tablename_fit = sprintf('%s_%d_FITNESS',expt_name,density);
tablename_p2o = '4C3_pos2orf_name3';
tablename_bpos = '4C3_borderpos';

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
i <- 1
for (pl in sort(unique(p2c$`6144plate`))) {
  for (sr in sort(unique(p2c$source))) {
    temp <- p2c[p2c$source == sr & p2c$`6144plate` == pl,]
    for (c in sort(unique(temp$`6144col`))) {
      for (r in sort(unique(temp$`6144row`[temp$`6144col` == c]))) {
        pos <- rbind(pos, temp$pos[temp$`6144col` == c & temp$`6144row` == r])
        len <- length(temp$pos[temp$`6144col` == c - 2 & temp$`6144row` %in% c(r - 2, r - 1, r, r + 1, r + 2) |
                                 temp$`6144col` == c - 1 & temp$`6144row` %in% c(r - 2, r - 1, r, r + 1, r + 2) |
                                 temp$`6144col` == c & temp$`6144row` %in% c(r - 2, r - 1, r + 1, r + 2) |
                                 temp$`6144col` == c + 1 & temp$`6144row` %in% c(r - 2, r - 1, r, r + 1, r + 2) |
                                 temp$`6144col` == c + 2 & temp$`6144row` %in% c(r - 2, r - 1, r, r + 1, r + 2)])
        neigh[i,1:len] <- temp$pos[temp$`6144col` == c - 2 & temp$`6144row` %in% c(r - 2, r - 1, r, r + 1, r + 2) |
                                     temp$`6144col` == c - 1 & temp$`6144row` %in% c(r - 2, r - 1, r, r + 1, r + 2) |
                                     temp$`6144col` == c & temp$`6144row` %in% c(r - 2, r - 1, r + 1, r + 2) |
                                     temp$`6144col` == c + 1 & temp$`6144row` %in% c(r - 2, r - 1, r, r + 1, r + 2) |
                                     temp$`6144col` == c + 2 & temp$`6144row` %in% c(r - 2, r - 1, r, r + 1, r + 2)]
        i <- i + 1
      }
    } 
  }
}
grids <- cbind(pos, neigh)

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

alldat <- alldat[alldat$hours == hr,]
# for (hr in sort(unique(alldat$hours))) {
  for (i in seq(1,dim(grids)[1])) {
   alldat$neigh[alldat$hours == hr & alldat$pos == grids[i]] <- median(alldat$average[alldat$hours == hr & alldat$pos %in% grids[i,2:25]], na.rm = T)
  }
# }

alldat$neigh[is.na(alldat$average) & !is.na(alldat$orf_name)] <- NA
alldat$diff <- alldat$average - alldat$neigh

ggplot(alldat) +
  # geom_line(aes(x = average, col = 'Pix'), stat = 'density') +
  # geom_line(aes(x = neigh, col = 'Neigh'), stat = 'density')
  geom_abline() +
  geom_point(aes(x = average, y = neigh, col = colony, shape = outlier), alpha =0.5)
  # geom_histogram(aes(x = diff), fill = 'Blue') +
  # coord_cartesian(ylim = c(0, 10))

quantile(alldat$diff, c(0.05,0.95), na.rm = T)

m <- mean(alldat$diff, na.rm = T)
s <- sd(alldat$diff, na.rm = T)


alldat$outlier[alldat$diff > m + 2*s] = 'Big'
alldat$outlier[alldat$diff < m - 2*s] = 'Small'

ggplot(alldat[alldat$`6144plate` == 1,]) +
  geom_point(aes(x = `6144col`, y = `6144row`, col = outlier, shape = colony))
