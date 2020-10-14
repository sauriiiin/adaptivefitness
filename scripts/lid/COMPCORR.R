##### COMPCORR
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 02/11/2020
##### Dividing the plates in subgrids and comparing fitness and raw pixel distributions
##### For competition correction

##### INITIALIZE
library(RMariaDB)
library(ggplot2)
library(ggExtra)
library(gridExtra)
library(ggpubr)
source("R/functions/initialize.sql.R")

##### GET/SET DATA
expt_name = 'SDPG_GLU_FS_R1'
expt = 'SDPG GLU FS R1'
density = 1536;

##### CHECK POSITION WISE VARIABILITY
conn <- initialize.sql("saurin_test")

tablename_p2o = 'SDPG_pos2orf_name'
tablename_jpeg = sprintf('%s_%d_RAW',expt_name,density);
tablename_jpeg_cc = sprintf('%s_CC_%d_RAW',expt_name,density);

p2c_info = NULL
p2c_info[1] = 'SDPG_pos2coor'
p2c_info[2] = 'plate'
p2c_info[3] = 'col'
p2c_info[4] = 'row'

p2c <- dbGetQuery(conn, sprintf('select * from %s
                                where density = %d
                                order by %s, %s, %s',
                                p2c_info[1],
                                density,
                                p2c_info[2], p2c_info[3], p2c_info[4]))
p2c$source = 'ALL'

pos <- NULL
neigh <- matrix(ncol = 8, nrow = density*length(unique(p2c$plate)))
neigh_sr <- matrix(ncol = 8, nrow = density*length(unique(p2c$plate)))
i <- 1
for (pl in sort(unique(p2c$plate))) {
  for (sr in sort(unique(p2c$source))) {
    temp <- p2c[p2c$source == sr & p2c$plate == pl,]
    for (c in sort(unique(temp$col))) {
      for (r in sort(unique(temp$row[temp$col == c]))) {
        pos <- rbind(pos, temp$pos[temp$col == c & temp$row == r])
        
        len <- length(temp$pos[temp$col == c - 1 & temp$row %in% c(r - 1, r, r + 1) |
                                 temp$col == c & temp$row %in% c(r - 1, r + 1) |
                                 temp$col == c + 1 & temp$row %in% c(r - 1, r, r + 1)])
        neigh[i,1:len] <- temp$pos[temp$col == c - 1 & temp$row %in% c(r - 1, r, r + 1) |
                                     temp$col == c & temp$row %in% c(r - 1, r + 1) |
                                     temp$col == c + 1 & temp$row %in% c(r - 1, r, r + 1)]
        len_sr <- length(temp$pos[temp$col == c - 2 & temp$row %in% c(r - 2, r, r + 2) |
                                    temp$col == c & temp$row %in% c(r - 2, r + 2) |
                                    temp$col == c + 2 & temp$row %in% c(r - 2, r, r + 2)])
        neigh_sr[i,1:len_sr] <- temp$pos[temp$col == c - 2 & temp$row %in% c(r - 2, r, r + 2) |
                                           temp$col == c & temp$row %in% c(r - 2, r + 2) |
                                           temp$col == c + 2 & temp$row %in% c(r - 2, r, r + 2)]
        i <- i + 1
      }
    } 
  }
}
grids <- cbind(data.frame(pos), data.frame(neigh))
grids_sr <- cbind(data.frame(pos), data.frame(neigh_sr))

grids <- data.frame(lapply(grids, as.numeric))
grids_sr <- data.frame(lapply(grids_sr, as.numeric))

alldat <- dbGetQuery(conn, sprintf('select b.*, c.orf_name, a.hours, a.average
                                   from %s a, %s b, %s c
                                   where a.pos = b.pos and b.pos = c.pos
                                   and a.hours > 0
                                   order by a.hours, b.%s, b.%s, b.%s',
                                   tablename_jpeg,
                                   p2c_info[1],tablename_p2o,
                                   p2c_info[2],
                                   p2c_info[3],p2c_info[4]))

alldat$colony[alldat$orf_name == 'BF_control'] = 'Reference'
alldat$colony[alldat$orf_name != 'BF_control'] = 'Query'

alldat$orf_name[alldat$orf_name == ''] = NA
alldat$colony[is.na(alldat$orf_name)] = 'Gap'

alldat$source[alldat$row%%2==1 & alldat$col%%2==1] = '1TL'
alldat$source[alldat$row%%2==0 & alldat$col%%2==1] = '3BL'
alldat$source[alldat$row%%2==1 & alldat$col%%2==0] = '2TR'
alldat$source[alldat$row%%2==0 & alldat$col%%2==0] = '4BR'

##### COMPETITION CORRECTION
compdat <- data.frame()
sick_N <- NULL
healthy_N <- NULL

for (hr in sort(unique(alldat$hours))) {
  for (pl in sort(unique(alldat$plate[alldat$hours == hr]))) {
    tempdat <- alldat[alldat$hours == hr & alldat$plate == pl,]
    tempdat$average[is.na(tempdat$orf_name)] <- 0
    
    for (i in seq(1,dim(grids)[1])) {
      tempdat$neigh[tempdat$pos == grids[i,1]] <- mean(tempdat$average[tempdat$pos %in% grids[i,2:9]], na.rm = T)
      tempdat$neigh_sr[tempdat$pos == grids_sr[i,1]] <- mean(tempdat$average[tempdat$pos %in% grids_sr[i,2:9]], na.rm = T)
    }
    
    tempdat$neigh[is.na(tempdat$average) & !is.na(tempdat$orf_name)] <- NA
    tempdat$neigh_sr[is.na(tempdat$average) & !is.na(tempdat$orf_name)] <- NA
    tempdat$score <- tempdat$average/((tempdat$neigh + tempdat$neigh_sr)/2)
    
    md <- mad(tempdat$score[tempdat$orf_name == 'BF_control'], na.rm =T)
    ll <- median(tempdat$score[tempdat$orf_name == 'BF_control'], na.rm =T) - 3*md
    ul <- median(tempdat$score[tempdat$orf_name == 'BF_control'], na.rm =T) + 3*md
    
    tempdat$size <- NULL
    tempdat$size[tempdat$score < ll & !is.na(tempdat$score)] <- 'Small'
    tempdat$size[tempdat$score > ul & !is.na(tempdat$score)] <- 'Big'
    tempdat$size[is.na(tempdat$size)] <- 'Normal'
    
    tempdat$sick <- NULL
    tempdat$healthy_neigh <- NULL
    for (p in tempdat$pos[tempdat$size == 'Small']) {
      # N <- sum(tempdat$size[tempdat$pos %in% grids[grids[,1] == p, 2:9] | tempdat$pos %in% grids_sr[grids_sr[,1] == p, 2:9]] == 'Big', na.rm = T)
      N <- sum(tempdat$size[tempdat$pos %in% grids[grids[,1] == p, 2:9]] == 'Big', na.rm = T)
      if (N > 0) {
        tempdat$sick[tempdat$pos == p] <- 'Y'
        tempdat$healthy_neigh[tempdat$pos == p] <- N
        healthy_N <- rbind(healthy_N, N)
      }
    }
    tempdat$sick[is.na(tempdat$sick) & tempdat$size == 'Small'] <- 'N'
    tempdat$sick[is.na(tempdat$sick)] <- 'Normal'
    
    tempdat$healthy <- NULL
    tempdat$sick_neigh <- NULL
    for (p in tempdat$pos[tempdat$size == 'Big']) {
      # N <- sum(tempdat$size[tempdat$pos %in% grids[grids[,1] == p, 2:9] | tempdat$pos %in% grids_sr[grids_sr[,1] == p, 2:9]] == 'Small', na.rm = T)
      N <- sum(tempdat$size[tempdat$pos %in% grids[grids[,1] == p, 2:9]] == 'Small', na.rm = T)
      if (N > 0) {
        tempdat$healthy[tempdat$pos == p] <- 'Y'
        tempdat$sick_neigh[tempdat$pos == p] <- N
        sick_N <- rbind(sick_N, N)
      }
    }
    tempdat$healthy[is.na(tempdat$healthy) & tempdat$size == 'Big'] <- 'N'
    tempdat$healthy[is.na(tempdat$healthy)] <- 'Normal'
    
    tempdat$comp.b <- NULL
    tempdat$comp.s <- NULL
    tempdat$comp <- NULL
    tempdat$driver <- NULL
    
    b2s <- NULL
    s2b <- NULL
    for (p in tempdat$pos[tempdat$sick == 'Y' | tempdat$healthy == 'Y']){
      if (tempdat$healthy[tempdat$pos == p] == 'Y') {
        # N <- sum(tempdat$healthy_neigh[tempdat$pos %in% grids[grids[,1] == p, 2:9] | tempdat$pos %in% grids_sr[grids_sr[,1] == p, 2:9]] + 1 >
        #            tempdat$sick_neigh[tempdat$pos == p], na.rm = T)
        N <- sum(tempdat$healthy_neigh[tempdat$pos %in% grids[grids[,1] == p, 2:9]] + 1 > tempdat$sick_neigh[tempdat$pos == p], na.rm = T)
        # a healthy colony should have no sick neighbors which have more healthy neighbors than it has sick ones
        b2s <- rbind(b2s, c(tempdat$sick_neigh[tempdat$pos == p], tempdat$healthy_neigh[tempdat$pos %in% grids[grids[,1] == p, 2:9]]))
        if (N == 0) {
          tempdat$comp[tempdat$pos %in% grids[grids[,1] == p, 2:9]] <- 'CB'
          tempdat$comp.b[tempdat$pos %in% grids[grids[,1] == p, 2:9]] <- 'CB'
          tempdat$driver[tempdat$pos == p] <- 'Big'
        }
      } else {
        # N <- sum(tempdat$sick_neigh[tempdat$pos %in% grids[grids[,1] == p, 2:9] | tempdat$pos %in% grids_sr[grids_sr[,1] == p, 2:9]] <
        #            tempdat$healthy_neigh[tempdat$pos == p], na.rm = T)
        N <- sum(tempdat$sick_neigh[tempdat$pos %in% grids[grids[,1] == p, 2:9]] < tempdat$healthy_neigh[tempdat$pos == p], na.rm = T)
        # a sick colony should have atleast one healthy neighbor that has less sick nieghbors than it has healthy ones
        s2b <- rbind(s2b, c(tempdat$healthy_neigh[tempdat$pos == p], tempdat$sick_neigh[tempdat$pos %in% grids[grids[,1] == p, 2:9]]))
        if (N > 0) {
          tempdat$comp[tempdat$pos %in% grids[grids[,1] == p, 2:9]] <- 'CS'
          tempdat$comp.s[tempdat$pos %in% grids[grids[,1] == p, 2:9]] <- 'CS'
          tempdat$driver[tempdat$pos == p] <- 'Small'
        }
      }
    }
    # tempdat$comp[is.na(tempdat$orf_name)] = NA
    if (sum(!is.na(tempdat$comp.b)) > 0) {
      tempdat$comp[!is.na(tempdat$comp.b) & !is.na(tempdat$comp.s)] = NA
    } else {
      tempdat$comp.b <- 'NA'
    }
    tempdat$comp[tempdat$sick == 'Y' & tempdat$comp == 'CS'] = NA
    tempdat$comp[tempdat$healthy == 'Y' & tempdat$comp == 'CB'] = NA
    tempdat$comp[is.na(tempdat$comp)] = 'No'

    tempdat$average_cc <- tempdat$average
    tempdat$average_cc[tempdat$comp == 'CS'] <- tempdat$average_cc[tempdat$comp == 'CS'] *
      median(tempdat$average_cc[tempdat$orf_name == 'BF_control'], na.rm = T)/median(tempdat$average_cc[tempdat$comp == 'CS'], na.rm = T)
    tempdat$average_cc[tempdat$comp == 'CB'] <- tempdat$average_cc[tempdat$comp == 'CB'] *
      median(tempdat$average_cc[tempdat$orf_name == 'BF_control'], na.rm = T)/median(tempdat$average_cc[tempdat$comp == 'CB'], na.rm = T)
    
    compdat <- rbind(compdat, tempdat)
  }
}
compdat$average[is.na(compdat$orf_name)] <- NA
compdat$average_cc[is.na(compdat$orf_name)] <- NA

jpegdat <- data.frame(compdat$pos, compdat$hours, compdat$average, compdat$average_cc)
colnames(jpegdat) <- c('pos','hours','average_raw', 'average')
dbWriteTable(conn, tablename_jpeg_cc, jpegdat, overwrite = T)


