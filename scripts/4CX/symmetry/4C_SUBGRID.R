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
expt_name = 'OESP1_FS'
expt = 'OESP1 FS'
out_path = 'figs/comp/';
density = 6144;

##### CHECK POSITION WISE VARIABILITY
conn <- initialize.sql("saurin_test")

tablename_p2o = 'OESP1_pos2orf_name'
tablename_jpeg = sprintf('%s_%d_RAW',expt_name,density);
tablename_jpeg_cc = sprintf('%s_CC_%d_RAW',expt_name,density);
tablename_fit = sprintf('%s_%d_FITNESS',expt_name,density);
tablename_fit_cc = sprintf('%s_CC_%d_FITNESS',expt_name,density);
# tablename_mdfr = sprintf('%s_CC_%d_MDFR',expt_name,density);

p2c_info = NULL
p2c_info[1] = 'OESP1_pos2coor'
p2c_info[2] = 'plate'
p2c_info[3] = 'col'
p2c_info[4] = 'row'

p2c <- dbGetQuery(conn, sprintf('select * from %s
                                where density = %d
                                order by %s, %s, %s',
                  p2c_info[1],
                  density,
                  p2c_info[2], p2c_info[3], p2c_info[4]))

# p2c$source[p2c$row%%2==1 & p2c$col%%2==1] = '1TL'
# p2c$source[p2c$row%%2==0 & p2c$col%%2==1] = '3BL'
# p2c$source[p2c$row%%2==1 & p2c$col%%2==0] = '2TR'
# p2c$source[p2c$row%%2==0 & p2c$col%%2==0] = '4BR'

p2c$source = 'ALL'

pos <- NULL
neigh <- matrix(ncol = 8, nrow = 6144*8)
neigh_sr <- matrix(ncol = 8, nrow = 6144*8)
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
grids <- cbind(pos, neigh)
grids_sr <- cbind(pos, neigh_sr)

alldat <- dbGetQuery(conn, sprintf('select b.*, c.orf_name, a.hours, a.average
                                  from %s a, %s b, %s c
                                  where a.pos = b.pos and b.pos = c.pos
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

# alldat.cc <- dbGetQuery(conn, sprintf('select a.*, b.*
#                                   from %s a, %s b
#                                   where a.pos = b.pos
#                                   order by a.hours, b.%s, b.%s, b.%s',
#                                   tablename_fit_cc,
#                                   p2c_info[1],p2c_info[2],
#                                   p2c_info[3],p2c_info[4]))
# 
# alldat.cc$colony[alldat.cc$orf_name == 'BF_control'] = 'Reference'
# alldat.cc$colony[alldat.cc$orf_name != 'BF_control'] = 'Query'
# alldat.cc$colony[is.na(alldat.cc$orf_name)] = 'Gap'
# 
# alldat.cc$source[alldat.cc$row%%2==1 & alldat.cc$col%%2==1] = '1TL'
# alldat.cc$source[alldat.cc$row%%2==0 & alldat.cc$col%%2==1] = '3BL'
# alldat.cc$source[alldat.cc$row%%2==1 & alldat.cc$col%%2==0] = '2TR'
# alldat.cc$source[alldat.cc$row%%2==0 & alldat.cc$col%%2==0] = '4BR'

# alldat.b <- dbGetQuery(conn, sprintf('select a.*, b.*
#                                   from %s a, %s b
#                                   where a.pos = b.pos and a.hours = 17 and b.6144plate = 1
#                                   order by a.hours, b.%s, b.%s, b.%s',
#                                   '4C3_GA3_BEAN_6144_FITNESS',
#                                   p2c_info[1],p2c_info[2],
#                                   p2c_info[3],p2c_info[4]))
# 
# alldat.b$colony[alldat.b$orf_name == 'BF_control'] = 'Reference'
# alldat.b$colony[alldat.b$orf_name != 'BF_control'] = 'Query'
# alldat.b$colony[is.na(alldat.b$orf_name)] = 'Gap'
# 
# alldat.b$source[alldat.b$row%%2==1 & alldat.b$col%%2==1] = '1TL'
# alldat.b$source[alldat.b$row%%2==0 & alldat.b$col%%2==1] = '3BL'
# alldat.b$source[alldat.b$row%%2==1 & alldat.b$col%%2==0] = '2TR'
# alldat.b$source[alldat.b$row%%2==0 & alldat.b$col%%2==0] = '4BR'

##### COMPETITION CORRECTION
compdat <- data.frame()
sick_N <- NULL
healthy_N <- NULL

for (hr in sort(unique(alldat$hours))) {
  for (pl in sort(unique(alldat$plate[alldat$hours == hr]))) {
    tempdat <- alldat[alldat$hours == hr & alldat$plate == pl,]
    tempdat$average[is.na(tempdat$orf_name)] <- 0

    for (i in seq(1,dim(grids)[1])) {
      tempdat$neigh[tempdat$pos == grids[i]] <- mean(tempdat$average[tempdat$pos %in% grids[i,2:9]], na.rm = T)
      tempdat$neigh_sr[tempdat$pos == grids_sr[i]] <- mean(tempdat$average[tempdat$pos %in% grids_sr[i,2:9]], na.rm = T)
    }
    
    tempdat$neigh[is.na(tempdat$average) & !is.na(tempdat$orf_name)] <- NA
    tempdat$neigh_sr[is.na(tempdat$average) & !is.na(tempdat$orf_name)] <- NA
    tempdat$score <- tempdat$average/((tempdat$neigh + tempdat$neigh_sr)/2)
    
    md <- mad(tempdat$score[tempdat$orf_name == 'BF_control'], na.rm =T)
    ll <- median(tempdat$score[tempdat$orf_name == 'BF_control'], na.rm =T) - 3*md
    ul <- median(tempdat$score[tempdat$orf_name == 'BF_control'], na.rm =T) + 3*md
    
    # ggplot(tempdat) +
    #   geom_line(aes(x = score, col = colony), stat = 'density', lwd = 1.2) +
    #   # geom_histogram(aes(x = score, fill = colony)) +
    #   geom_vline(xintercept = c(ll,ul), col = 'Red', linetype = 'dashed') +
    #   labs(title = expt,
    #        subtitle = sprintf('Competition Scores of Plate #%d at %0.2f Hrs', pl, hr),
    #        x = 'Competition Score',
    #        y = 'Density') +
    #   scale_color_discrete(name = 'Colony Type') +
    #   theme_linedraw()
    # ggsave(sprintf("%scomp_score_%d_%d.jpg",out_path,round(hr),pl),
    #        width = 8, height = 6,
    #        dpi = 300)

    # ggplot(tempdat) +
    #   geom_point(aes(x = col, y = row, shape = colony, col = 'Normal'), size = 1) +
    #   geom_point(data = tempdat[tempdat$score > ul,],
    #              aes(x = col, y = row, shape = colony, col = 'Big')) +
    #   geom_point(data = tempdat[tempdat$score < ll,],
    #              aes(x = col, y = row, shape = colony, col = 'Small')) +
    #   # geom_point(data = tempdat[tempdat$colony == 'Gap',],
    #   #            aes(x = col, y = row, shape = colony, col = 'Gap')) +
    #   geom_text(data = tempdat, aes(x = col, y = row, label = average), check_overlap = TRUE, size = 0.83) +
    #   scale_x_continuous(breaks = seq(0,96,4)) +
    #   scale_y_continuous(breaks = seq(0,64,4),trans = 'reverse') +
    #   labs(title = sprintf('%s: Potential Drivers of Competition',expt),
    #        subtitle = sprintf('Plate #%d at %0.2f Hrs', pl, hr),
    #        x = 'Columns',
    #        y = 'Rows') +
    #   scale_color_manual(name = 'Size',
    #                      breaks = c('Normal','Small','Big','Gap'),
    #                      values = c('Normal'='#9E9E9E','Small'='#673AB7','Big'='#FFC107','Gap'='#FF5252')) +
    #   scale_shape_discrete(name = 'Colony Type') +
    #   theme_linedraw()
    # ggsave(sprintf("%spotentialdrivers_%d_%d.jpg",out_path,round(hr),pl),
    #        width = 8, height = 5,
    #        dpi = 300)
      
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
    
    # ggplot() +
    #   geom_point(data = tempdat,
    #              aes(x = col, y = row, shape = colony, col = 'Normal'), size = 1) +
    #   # geom_point(data = tempdat[tempdat$colony == 'Gap',],
    #   #             aes(x = col, y = row, shape = colony, col = 'Gap')) +
    #   geom_point(data = tempdat[tempdat$healthy_neigh > 0 & !is.na(tempdat$healthy_neigh),],
    #              aes(x = col, y = row, shape = colony, col = 'Small')) +
    #   geom_point(data = tempdat[tempdat$sick_neigh > 0 & !is.na(tempdat$sick_neigh),],
    #              aes(x = col, y = row, shape = colony, col = 'Big')) +
    #   geom_text(data = tempdat, aes(x = col, y = row, label = average), check_overlap = TRUE, size = 0.83) +
    #   scale_x_continuous(breaks = seq(0,96,4)) +
    #   scale_y_continuous(breaks = seq(0,64,4),trans = 'reverse') +
    #   labs(title = sprintf('%s: Drivers of Competition',expt),
    #        subtitle = sprintf('Plate #%d at %0.2f Hrs', pl, hr),
    #        x = 'Columns',
    #        y = 'Rows') +
    #   scale_color_manual(name = 'Size',
    #                      breaks = c('Normal','Small','Big','Gap'),
    #                      values = c('Normal'='#9E9E9E','Small'='#673AB7','Big'='#FFC107','Gap'='#FF5252')) +
    #   scale_shape_discrete(name = 'Colony Type') +
    #   theme_linedraw()
    # ggsave(sprintf("%s%s_drivers_%d_%d.jpg",out_path,expt_name,
    #                round(hr),pl),
    #        width = 8, height = 5,
    #        dpi = 300)
    
    # ggplot() +
    #   geom_point(data = tempdat,
    #              aes(x = col, y = row, shape = colony, col = 'Normal'), size = 1) +
    #   # geom_point(data = tempdat[tempdat$colony == 'Gap',],
    #   #             aes(x = col, y = row, shape = colony, col = 'Gap')) +
    #   geom_point(data = tempdat[tempdat$healthy_neigh > 0 & !is.na(tempdat$healthy_neigh),],
    #              aes(x = col, y = row, shape = colony, col = 'Small')) +
    #   geom_point(data = tempdat[tempdat$sick_neigh > 0 & !is.na(tempdat$sick_neigh),],
    #              aes(x = col, y = row, shape = colony, col = 'Big')) +
    #   geom_text(data = tempdat, aes(x = col, y = row, label = round(score, 2)), check_overlap = TRUE, size = 0.83) +
    #   scale_x_continuous(breaks = seq(0,96,4)) +
    #   scale_y_continuous(breaks = seq(0,64,4),trans = 'reverse') +
    #   labs(title = sprintf('%s: Drivers of Competition',expt),
    #        subtitle = sprintf('Plate #%d at %0.2f Hrs', pl, hr),
    #        x = 'Columns',
    #        y = 'Rows') +
    #   scale_color_manual(name = 'Size',
    #                      breaks = c('Normal','Small','Big','Gap'),
    #                      values = c('Normal'='#9E9E9E','Small'='#673AB7','Big'='#FFC107','Gap'='#FF5252')) +
    #   scale_shape_discrete(name = 'Colony Type') +
    #   theme_linedraw()
    # ggsave(sprintf("%sdrivers_score_%d_%d.jpg",out_path,round(hr),pl),
    #        width = 8, height = 5,
    #        dpi = 300)
    
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
    
    # ggplot() +
    #   geom_point(data = tempdat[tempdat$comp == 'No',],aes(x = col, y = row, col = comp, shape = colony), size = 1) +
    #   geom_point(data = tempdat[tempdat$comp != 'No',],aes(x = col, y = row, col = comp, shape = colony)) +
    #   geom_point(data = tempdat[!is.na(tempdat$driver),],
    #              aes(x = col, y = row, col = driver, shape = colony)) +
    #   geom_text(data = tempdat, aes(x = col, y = row, label = average), check_overlap = TRUE, size = 0.83) +
    #   scale_x_continuous(breaks = seq(0,96,4)) +
    #   scale_y_continuous(breaks = seq(0,64,4),trans = 'reverse') +
    #   labs(title = sprintf('%s: Competition Correction',expt),
    #        subtitle = sprintf('Plate #%d at %0.2f Hrs', pl, hr),
    #        x = 'Columns',
    #        y = 'Rows') +
    #   scale_color_manual(name = 'Comp Pos',
    #                      breaks = c('No','Small','Big','CS','CB'),
    #                      values = c('No'='#9E9E9E','Small'='#673AB7','Big'='#FFC107','CS'='#8BC34A','CB'='#536DFE'),
    #                      labels = c('Normal','Small','Big','Next2S','Next2B')) +
    #   scale_shape_discrete(name = 'Colony Type') +
    #   theme_linedraw()
    # ggsave(sprintf("%s%scompcorr_%d_%d.jpg",out_path,expt_name,
    #                round(hr),pl),
    #        width = 8, height = 5,
    #        dpi = 300)
    
    # ggplot(tempdat[tempdat$colony != 'Gap',]) +
    #   geom_histogram(aes(x = average, fill = comp), position = 'stack')
    
    # tempdat$nearsick <- NULL
    # tempdat$nearsick_sr <- NULL
    # for (p in tempdat$pos[tempdat$sick == 'Y']) {
    #   tempdat$nearsick[tempdat$pos %in% grids[grids[,1] == p, 2:9]] <- 'N1'
    #   tempdat$nearsick_sr[tempdat$pos %in% grids_sr[grids_sr[,1] == p, 2:9]] <- 'N2'
    # }
    # tempdat$nearsick[is.na(tempdat$orf_name)] <- 'N'
    # tempdat$nearsick[is.na(tempdat$nearsick)] <- 'N'
    # tempdat$nearsick[tempdat$sick == 'Y'] <- 'N'
    # tempdat$nearsick_sr[is.na(tempdat$orf_name)] <- 'N'
    # tempdat$nearsick_sr[is.na(tempdat$nearsick_sr)] <- 'N'
    # tempdat$nearsick_sr[tempdat$nearsick == 'N1'] <- 'N'
    # tempdat$nearsick_sr[tempdat$sick == 'Y'] <- 'N'
    
    # ggplot(tempdat) +
    #   geom_point(aes(x = col, y = row, shape = colony, col = 'Normal'), size = 1) +
    #   geom_point(data = tempdat[tempdat$colony == 'Gap',],
    #              aes(x = col, y = row, shape = colony, col = 'Gap')) +
    #   geom_point(data = tempdat[tempdat$nearsick_sr == 'N2',],
    #              aes(x = col, y = row, shape = colony, col = nearsick_sr)) +
    #   geom_point(data = tempdat[tempdat$nearsick == 'N1',],
    #              aes(x = col, y = row, shape = colony, col = nearsick)) +
    #   geom_point(data = tempdat[tempdat$sick == 'Y',],
    #              aes(x = col, y = row, shape = colony, col = sick)) +
    #   # geom_point(data = tempdat[tempdat$colony == 'Gap',],
    #   #            aes(x = col, y = row, shape = colony, col = 'Gap')) +
    #   scale_x_continuous(breaks = seq(0,96,4)) +
    #   scale_y_continuous(breaks = seq(0,64,4),trans = 'reverse') +
    #   labs(title = expt,
    #        subtitle = sprintf('Corrected Positions of Plate #%d at %d Hrs', pl, hr),
    #        x = 'Columns',
    #        y = 'Rows') +
    #   scale_color_manual(name = 'Health',
    #                      breaks = c('Normal','Y','N1','N2','Gap'),
    #                      values = c('Normal'='#9E9E9E','Y'='#2196F3','N1'='#FFC107','N2'='#673AB7','Gap'='#FF5252'),
    #                      labels = c('Normal','Small','Pres Neigh','Past Neigh','Gap')) +
    #   scale_shape_discrete(name = 'Colony Type') +
    #   theme_linedraw()
    # ggsave(sprintf("%snearsick_%d_%d.jpg",out_path,hr,pl),
    #        width = 8, height = 5,
    #        dpi = 300)
    
    tempdat$average_cc <- tempdat$average
    # tempdat$average_cc[tempdat$nearsick == 'N1'] <- tempdat$average_cc[tempdat$nearsick == 'N1'] *
    #   median(tempdat$average_cc[tempdat$orf_name == 'BF_control'], na.rm = T)/median(tempdat$average_cc[tempdat$nearsick == 'N1'], na.rm = T)
    # tempdat$average_cc[tempdat$nearsick_sr == 'N2'] <- tempdat$average_cc[tempdat$nearsick_sr == 'N2'] *
    #   median(tempdat$average_cc[tempdat$orf_name == 'BF_control'], na.rm = T)/median(tempdat$average_cc[tempdat$nearsick_sr == 'N2'], na.rm = T)
    tempdat$average_cc[tempdat$comp == 'CS'] <- tempdat$average_cc[tempdat$comp == 'CS'] *
      median(tempdat$average_cc[tempdat$orf_name == 'BF_control'], na.rm = T)/median(tempdat$average_cc[tempdat$comp == 'CS'], na.rm = T)
    tempdat$average_cc[tempdat$comp == 'CB'] <- tempdat$average_cc[tempdat$comp == 'CB'] *
      median(tempdat$average_cc[tempdat$orf_name == 'BF_control'], na.rm = T)/median(tempdat$average_cc[tempdat$comp == 'CB'], na.rm = T)

    
    # ggplot() +
    #   geom_line(data = tempdat[tempdat$colony != 'Gap',],
    #             aes(x = average, col = 'RAW'), stat = 'density') +
    #   geom_line(data = tempdat[tempdat$colony != 'Gap',],
    #             aes(x = average_cc, col = 'CC'), stat = 'density') +
    #   geom_vline(xintercept = quantile(tempdat$average, c(0.025, 0.5, 0.975), na.rm = T), col = 'Blue') +
    #   geom_vline(xintercept = quantile(tempdat$average_cc, c(0.025, 0.5,  0.975), na.rm = T), col = 'Red') +
    #   facet_grid(.~colony)
    #   geom_line(data = tempdat[tempdat$nearsick == 'N1',],
    #             aes(x = average_cc, col = 'Pres Neigh'), stat = 'density') +
    #   geom_line(data = tempdat[tempdat$nearsick_sr == 'N2',],
    #             aes(x = average_cc, col = 'Past Neigh'), stat = 'density')
    # 
    # ggplot() +
    #   geom_point(data = tempdat[tempdat$nearsick == 'N1',],
    #             aes(x = average, y = average/score * median(score, na.rm = T), col = 'Pres Neigh')) +
    #   geom_point(data = tempdat[tempdat$nearsick_sr == 'N2',],
    #             aes(x = average, y = average/score * median(score, na.rm = T), col = 'Past Neigh')) +
    #   geom_abline()
    
    compdat <- rbind(compdat, tempdat)
  }
}

compdat$average[is.na(compdat$orf_name)] <- NA
compdat$average_cc[is.na(compdat$orf_name)] <- NA

jpegdat <- data.frame(compdat$pos, compdat$hours, compdat$average, compdat$average_cc)
colnames(jpegdat) <- c('pos','hours','average_raw', 'average')
dbWriteTable(conn, tablename_jpeg_cc, jpegdat, overwrite = T)

#####

plt.scr <- ggplot(tempdat[tempdat$colony != 'Gap',]) +
  geom_line(aes(x = average, col = colony), stat = 'density', lwd = 1.2) +
  geom_vline(xintercept = quantile(tempdat$average[tempdat$colony == 'Reference'], c(0.05,0.95), na.rm = T), col = 'blue') +
  geom_vline(xintercept = quantile(tempdat$average[tempdat$colony == 'Query'], c(0.05,0.95), na.rm = T), col = 'red') +
  labs(title = sprintf('ES = %0.3f', median(tempdat$average[tempdat$colony == 'Query'], na.rm = T)/
                         median(tempdat$average[tempdat$colony == 'Reference'], na.rm = T)),
       x = 'Raw Pixel Counts',
       y = 'Frequency') +
  theme_linedraw() +
  coord_cartesian(xlim = c(200,600))

plt.pix <- ggplot(tempdat[tempdat$colony != 'Gap',]) +
  geom_line(aes(x = average_cc, col = colony), stat = 'density', lwd = 1.2) +
  geom_vline(xintercept = quantile(tempdat$average_cc[tempdat$colony == 'Reference'], c(0.05,0.95), na.rm = T), col = 'blue') +
  geom_vline(xintercept = quantile(tempdat$average_cc[tempdat$colony == 'Query'], c(0.05,0.95), na.rm = T), col = 'red') +
  labs(title = sprintf('ES = %0.3f', median(tempdat$average_cc[tempdat$colony == 'Query'], na.rm = T)/
                         median(tempdat$average_cc[tempdat$colony == 'Reference'], na.rm = T)),
       x = 'CC Pixel Counts',
       y = 'Frequency') +
  theme_linedraw() +
  coord_cartesian(xlim = c(200,600))

annotate_figure(ggarrange(plt.scr, plt.pix, align = "hv",
                          common.legend = T, legend = 'right',
                          ncol = 2, nrow = 1),
                top = text_grob(sprintf('%s: Competition Correction of Plate #%d at %0.2f Hrs',
                                        expt, pl, hr)))

ggsave(sprintf("%scorrected_%d_%d.jpg",out_path,hr,pl),
       width = 10, height = 5,
       dpi = 300)

ggplot() +
  geom_line(data = alldat[alldat$hours == hr & alldat$plate == pl & alldat$pos %in%
                            compdat$pos[compdat$hours == hr & compdat$plate == pl & compdat$comp == 'CS'],],
            aes(x = fitness, col = 'LID W/O CC'), stat = 'density', lwd = 1.2) +
  geom_line(data = alldat.cc[alldat.cc$hours == hr & alldat.cc$plate == pl & alldat.cc$pos %in%
                               compdat$pos[compdat$hours == hr & compdat$plate == pl & compdat$comp == 'CS'],],
            aes(x = fitness, col = 'LID W CC'), stat = 'density', lwd = 1.2) +
  labs(title = 'Impact of Competition Correction (CC)',
       subtitle = 'On fitness of colonies growing near small colonies or gaps',
       x = 'Fitness',
       y = 'Density') +
  scale_color_discrete(name = '',
                       breaks = c('LID W CC',
                                  'LID W/O CC',
                                  'BEAN')) +
  coord_cartesian(xlim = c(0.7,1.3)) +
  theme_linedraw()
ggsave(sprintf("%s%s_COMP_CORR.png",
               out_path,expt_name),
       width = 6, height = 5)


ggplot() +
  geom_abline(col = 'red') +
  geom_point(data = alldat[alldat$hours == hr & alldat$plate == pl & alldat$pos %in%
                            compdat$pos[compdat$hours == hr & compdat$plate == pl & compdat$comp == 'CS'],],
            aes(x = fitness, 
                y = alldat.cc$fitness[alldat.cc$hours == hr & alldat.cc$plate == pl & alldat.cc$pos %in%
                                compdat$pos[compdat$hours == hr & compdat$plate == pl & compdat$comp == 'CS']],
                col = 'Next2S')) +
  geom_point(data = alldat[alldat$hours == hr & alldat$plate == pl & alldat$pos %in%
                             compdat$pos[compdat$hours == hr & compdat$plate == pl & compdat$comp == 'CB'],],
             aes(x = fitness, 
                 y = alldat.cc$fitness[alldat.cc$hours == hr & alldat.cc$plate == pl & alldat.cc$pos %in%
                                         compdat$pos[compdat$hours == hr & compdat$plate == pl & compdat$comp == 'CB']],
                 col = 'Next2B')) +
  coord_cartesian(xlim = c(0.7,1.3),
                  ylim = c(0.7,1.3))

ggplot() +
  geom_histogram(data = data.frame(sick_N),
                 aes(x = sick_N, fill = 'Small Neigh'), alpha = 0.8) +
  geom_histogram(data = data.frame(healthy_N),
                 aes(x = healthy_N, fill = 'Big Neigh'), alpha = 0.8)
  
  
  