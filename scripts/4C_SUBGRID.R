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
# tablename_mdfr = sprintf('%s_CC_%d_MDFR',expt_name,density);

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

# alldat.cc <- dbGetQuery(conn, sprintf('select a.*, b.*
#                                   from %s a, %s b
#                                   where a.pos = b.pos and a.hours = 18
#                                   order by a.hours, b.%s, b.%s, b.%s',
#                                   tablename_fit_cc,
#                                   p2c_info[1],p2c_info[2],
#                                   p2c_info[3],p2c_info[4]))
# 
# alldat.cc$colony[alldat.cc$orf_name == 'BF_control'] = 'Reference'
# alldat.cc$colony[alldat.cc$orf_name != 'BF_control'] = 'Query'
# alldat.cc$colony[is.na(alldat.cc$orf_name)] = 'Gap'
# 
# alldat.cc$source[alldat.cc$`6144row`%%2==1 & alldat.cc$`6144col`%%2==1] = '1TL'
# alldat.cc$source[alldat.cc$`6144row`%%2==0 & alldat.cc$`6144col`%%2==1] = '3BL'
# alldat.cc$source[alldat.cc$`6144row`%%2==1 & alldat.cc$`6144col`%%2==0] = '2TR'
# alldat.cc$source[alldat.cc$`6144row`%%2==0 & alldat.cc$`6144col`%%2==0] = '4BR'
# 
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
# alldat.b$source[alldat.b$`6144row`%%2==1 & alldat.b$`6144col`%%2==1] = '1TL'
# alldat.b$source[alldat.b$`6144row`%%2==0 & alldat.b$`6144col`%%2==1] = '3BL'
# alldat.b$source[alldat.b$`6144row`%%2==1 & alldat.b$`6144col`%%2==0] = '2TR'
# alldat.b$source[alldat.b$`6144row`%%2==0 & alldat.b$`6144col`%%2==0] = '4BR'

##### COMPETITION CORRECTION
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
    
    # ggplot(tempdat) +
    #   geom_line(aes(x = score, col = colony), stat = 'density', lwd = 1.2) +
    #   # geom_histogram(aes(x = score, fill = colony)) +
    #   geom_vline(xintercept = c(ll,ul), col = 'Red', linetype = 'dashed') +
    #   labs(title = expt,
    #        subtitle = sprintf('Competition Scores of Plate #%d at %d Hrs', pl, hr),
    #        x = 'Competition Score',
    #        y = 'Density') +
    #   scale_color_discrete(name = 'Colony Type') +
    #   theme_linedraw()
    # ggsave(sprintf("%scomp_score_%d_%d.jpg",out_path,hr,pl),
    #        width = 8, height = 6,
    #        dpi = 300)
    
    # ggplot(tempdat) +
    #   geom_point(aes(x = `6144col`, y = `6144row`, shape = colony, col = 'Normal'), size = 1) +
    #   geom_point(data = tempdat[tempdat$colony == 'Gap',],
    #              aes(x = `6144col`, y = `6144row`, shape = colony, col = 'Gap')) +
    #   geom_point(data = tempdat[tempdat$score > ul,],
    #              aes(x = `6144col`, y = `6144row`, shape = colony, col = 'Healthy')) +
    #   geom_point(data = tempdat[tempdat$score < ll,],
    #              aes(x = `6144col`, y = `6144row`, shape = colony, col = 'Sick')) +
    #   scale_x_continuous(breaks = seq(0,96,4)) +
    #   scale_y_continuous(breaks = seq(0,64,4),trans = 'reverse') +
    #   labs(title = expt,
    #        subtitle = sprintf('Colony Health Layout of Plate #%d at %d Hrs', pl, hr),
    #        x = 'Columns',
    #        y = 'Rows') +
    #   scale_color_manual(name = 'Health',
    #                      breaks = c('Normal','Sick','Healthy','Gap'),
    #                      values = c('Normal'='#9E9E9E','Sick'='#673AB7','Healthy'='#FFC107','Gap'='#FF5252')) +
    #   scale_shape_discrete(name = 'Colony Type') +
    #   theme_linedraw()
    # ggsave(sprintf("%shealth_dis_%d_%d.jpg",out_path,hr,pl),
    #        width = 8, height = 5,
    #        dpi = 300)
      
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
    tempdat$nearsick_sr <- NULL
    for (p in tempdat$pos[tempdat$reallysick == 'Y']) {
      tempdat$nearsick[tempdat$pos %in% grids[grids[,1] == p, 2:9]] <- 'N1'
      tempdat$nearsick_sr[tempdat$pos %in% grids_sr[grids_sr[,1] == p, 2:9]] <- 'N2'
    }
    tempdat$nearsick[is.na(tempdat$orf_name)] <- 'N'
    tempdat$nearsick[is.na(tempdat$nearsick)] <- 'N'
    tempdat$nearsick[tempdat$reallysick == 'Y'] <- 'N'
    tempdat$nearsick_sr[is.na(tempdat$orf_name)] <- 'N'
    tempdat$nearsick_sr[is.na(tempdat$nearsick_sr)] <- 'N'
    tempdat$nearsick_sr[tempdat$nearsick == 'N1'] <- 'N'
    tempdat$nearsick_sr[tempdat$reallysick == 'Y'] <- 'N'
    
    # ggplot(tempdat) +
    #   geom_point(aes(x = `6144col`, y = `6144row`, shape = colony, col = 'Normal'), size = 1) +
    #   geom_point(data = tempdat[tempdat$colony == 'Gap',],
    #              aes(x = `6144col`, y = `6144row`, shape = colony, col = 'Gap')) +
    #   geom_point(data = tempdat[tempdat$nearsick_sr == 'N2',],
    #              aes(x = `6144col`, y = `6144row`, shape = colony, col = nearsick_sr)) +
    #   geom_point(data = tempdat[tempdat$nearsick == 'N1',],
    #              aes(x = `6144col`, y = `6144row`, shape = colony, col = nearsick)) +
    #   geom_point(data = tempdat[tempdat$reallysick == 'Y',],
    #              aes(x = `6144col`, y = `6144row`, shape = colony, col = reallysick)) +
    #   # geom_point(data = tempdat[tempdat$colony == 'Gap',],
    #   #            aes(x = `6144col`, y = `6144row`, shape = colony, col = 'Gap')) +
    #   scale_x_continuous(breaks = seq(0,96,4)) +
    #   scale_y_continuous(breaks = seq(0,64,4),trans = 'reverse') +
    #   labs(title = expt,
    #        subtitle = sprintf('Corrected Positions of Plate #%d at %d Hrs', pl, hr),
    #        x = 'Columns',
    #        y = 'Rows') +
    #   scale_color_manual(name = 'Health',
    #                      breaks = c('Normal','Y','N1','N2','Gap'),
    #                      values = c('Normal'='#9E9E9E','Y'='#2196F3','N1'='#FFC107','N2'='#673AB7','Gap'='#FF5252'),
    #                      labels = c('Normal','Sick','Pres Neigh','Past Neigh','Gap')) +
    #   scale_shape_discrete(name = 'Colony Type') +
    #   theme_linedraw()
    # ggsave(sprintf("%snearsick_%d_%d.jpg",out_path,hr,pl),
    #        width = 8, height = 5,
    #        dpi = 300)
    
    tempdat$average_cc2 <- tempdat$average
    tempdat$average_cc2[tempdat$nearsick == 'N1'] <- tempdat$average_cc2[tempdat$nearsick == 'N1'] *
      median(tempdat$average_cc2[tempdat$orf_name == 'BF_control'], na.rm = T)/median(tempdat$average_cc2[tempdat$nearsick == 'N1'], na.rm = T)
    tempdat$average_cc2[tempdat$nearsick_sr == 'N2'] <- tempdat$average_cc2[tempdat$nearsick_sr == 'N2'] *
      median(tempdat$average_cc2[tempdat$orf_name == 'BF_control'], na.rm = T)/median(tempdat$average_cc2[tempdat$nearsick_sr == 'N2'], na.rm = T)
    # tempdat$average_cc2[tempdat$nearsick == 'N1' & !is.na(tempdat$score)] <- tempdat$average_cc2[tempdat$nearsick == 'N1' & !is.na(tempdat$score)] *
    #   median(tempdat$score[tempdat$orf_name == 'BF_control'], na.rm = T)/tempdat$score[tempdat$nearsick == 'N1' & !is.na(tempdat$score)]
    # tempdat$average_cc2[tempdat$nearsick_sr == 'N2' & !is.na(tempdat$score)] <- tempdat$average_cc2[tempdat$nearsick_sr == 'N2' & !is.na(tempdat$score)] *
    #   median(tempdat$score[tempdat$orf_name == 'BF_control'], na.rm = T)/tempdat$score[tempdat$nearsick_sr == 'N2' & !is.na(tempdat$score)]
    
    # ggplot() +
    #   geom_line(data = tempdat[tempdat$colony != 'Gap',],
    #             aes(x = average, col = 'RAW'), stat = 'density') +
    #   geom_line(data = tempdat[tempdat$colony != 'Gap',],
    #             aes(x = average_cc2, col = 'CC'), stat = 'density') +
    #   geom_vline(xintercept = quantile(tempdat$average, c(0.025, 0.5, 0.975), na.rm = T), col = 'Blue') +
    #   geom_vline(xintercept = quantile(tempdat$average_cc2, c(0.025, 0.5,  0.975), na.rm = T), col = 'Red') +
    #   facet_grid(.~colony)
    #   geom_line(data = tempdat[tempdat$nearsick == 'N1',],
    #             aes(x = average_cc2, col = 'Pres Neigh'), stat = 'density') +
    #   geom_line(data = tempdat[tempdat$nearsick_sr == 'N2',],
    #             aes(x = average_cc2, col = 'Past Neigh'), stat = 'density')
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
compdat$average_cc2[is.na(compdat$orf_name)] <- NA
compdat$modifier[is.na(compdat$orf_name)] <- NA

jpegdat <- data.frame(compdat$pos, compdat$hours, compdat$average, compdat$average_cc2)
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
  geom_line(aes(x = average_cc2, col = colony), stat = 'density', lwd = 1.2) +
  geom_vline(xintercept = quantile(tempdat$average_cc2[tempdat$colony == 'Reference'], c(0.05,0.95), na.rm = T), col = 'blue') +
  geom_vline(xintercept = quantile(tempdat$average_cc2[tempdat$colony == 'Query'], c(0.05,0.95), na.rm = T), col = 'red') +
  labs(title = sprintf('ES = %0.3f', median(tempdat$average_cc2[tempdat$colony == 'Query'], na.rm = T)/
                         median(tempdat$average_cc2[tempdat$colony == 'Reference'], na.rm = T)),
       x = 'CC Pixel Counts',
       y = 'Frequency') +
  theme_linedraw() +
  coord_cartesian(xlim = c(200,600))

annotate_figure(ggarrange(plt.scr, plt.pix, align = "hv",
                          common.legend = T, legend = 'right',
                          ncol = 2, nrow = 1),
                top = text_grob(sprintf('%s: Competition Correction of Plate #%d at %d Hrs',
                                        expt, pl, hr)))

ggsave(sprintf("%scorrected_%d_%d.jpg",out_path,hr,pl),
       width = 10, height = 5,
       dpi = 300)

ggplot() +
  geom_line(data = alldat.b[alldat.b$hours == 18 & alldat.b$`6144plate` == 1 & alldat.b$pos %in% tempdat$pos[tempdat$nearsick == 'N1'],],
            aes(x = fitness, col = 'BEAN'), stat = 'density', lwd = 1.2) +
  geom_line(data = alldat[alldat$hours == 18 & alldat$`6144plate` == 1 & alldat$pos %in% tempdat$pos[tempdat$nearsick == 'N1'],],
            aes(x = fitness, col = 'LID W/O CC'), stat = 'density', lwd = 1.2) +
  geom_line(data = alldat.cc[alldat.cc$hours == 18 & alldat.cc$`6144plate` == 1 & alldat.cc$pos %in% tempdat$pos[tempdat$nearsick == 'N1'],],
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

