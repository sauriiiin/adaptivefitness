##### BORDERS AND NEIGHBORS
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 06/24/2019

##### INITIALIZE
library(RMariaDB)
library(ggplot2)
library(ggExtra)
library(gridExtra)
source("R/functions/initialize.sql.R")

##### GET/SET DATA
expt_name = '4C3_GA1_RAW'
expt = 'FS1-1-RAW'
out_path = 'figs/borderandneigh/';
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
p2c = dbGetQuery(conn, sprintf('select * from %s a order by a.%s, a.%s, a.%s',
                               p2c_info[1],
                               p2c_info[2],
                               p2c_info[3],
                               p2c_info[4]))
n_plates = dbGetQuery(conn, sprintf('select distinct %s from %s a order by %s asc',
                                    p2c_info[2],
                                    p2c_info[1],
                                    p2c_info[2]))

# hr = 18
# pl = 2
jpegdat <- data.frame()

for (hr in 17:18) { #hours$hours
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
    
    alldat$border[alldat$`6144row` %in% 1 | alldat$`6144row` %in% 64 |
                    alldat$`6144col` %in% 1 | alldat$`6144col` %in% 96] = 'B1'
    alldat$border[alldat$`6144row` %in% 2 | alldat$`6144row` %in% 63 |
                    alldat$`6144col` %in% 2 | alldat$`6144col` %in% 95] = 'B2'
    alldat$border[alldat$`6144row` %in% 3 | alldat$`6144row` %in% 62 |
                    alldat$`6144col` %in% 3 | alldat$`6144col` %in% 94] = 'B3'
    alldat$border[alldat$`6144row` %in% 4 | alldat$`6144row` %in% 61 |
                    alldat$`6144col` %in% 4 | alldat$`6144col` %in% 93] = 'B4'
    alldat$border[is.na(alldat$border)] = 'N'
    
    ##### PLOTTING THE PIXEL COUNTS OF BORDERS WRT SOURCE INFORMATION
    ggplot(alldat) +
      geom_point(aes(x = average, col = border), stat = 'density') +
      facet_wrap(.~source, ncol = 2) +
      labs(x = 'Pixel Count',
           y = 'Density') +
      scale_color_discrete(name = 'Border',
                           breaks = c('B1','B2','B3','B4','N'),
                           labels = c('One','Two','Three','Four','Interior')) +
      scale_y_continuous(breaks = seq(0,0.03,0.005)) +
      theme_linedraw() +
      theme(legend.position = 'bottom') +
      coord_cartesian(xlim = c(200,700),
                      ylim = c(0,0.03))
    ggsave(sprintf("%s%s_BORDERS_%d_%d.png",
                   out_path,expt_name,hr,pl),
           width = 7,height = 7.5)
    
    ##### MEDIAN CORRECTING THE BORDERS
    alldat$average_mbc <- alldat$average
    
    for (sr in unique(alldat$source)) {
      alldat$average_mbc[alldat$border == 'B1' & alldat$source == sr]  <- alldat$average_mbc[alldat$border ==  "B1" & alldat$source == sr] *
        median(alldat$average_mbc[alldat$border ==  "N" & alldat$source == sr], na.rm = T)/median(alldat$average_mbc[alldat$border ==  "B1" & alldat$source == sr], na.rm = T)
      alldat$average_mbc[alldat$border == 'B2' & alldat$source == sr]  <- alldat$average_mbc[alldat$border ==  "B2" & alldat$source == sr] *
        median(alldat$average_mbc[alldat$border ==  "N" & alldat$source == sr], na.rm = T)/median(alldat$average_mbc[alldat$border ==  "B2" & alldat$source == sr], na.rm = T)
      alldat$average_mbc[alldat$border == 'B3' & alldat$source == sr]  <- alldat$average_mbc[alldat$border ==  "B3" & alldat$source == sr] *
        median(alldat$average_mbc[alldat$border ==  "N" & alldat$source == sr], na.rm = T)/median(alldat$average_mbc[alldat$border ==  "B3" & alldat$source == sr], na.rm = T)
      alldat$average_mbc[alldat$border == 'B4' & alldat$source == sr]  <- alldat$average_mbc[alldat$border ==  "B4" & alldat$source == sr] *
        median(alldat$average_mbc[alldat$border ==  "N" & alldat$source == sr], na.rm = T)/median(alldat$average_mbc[alldat$border ==  "B4" & alldat$source == sr], na.rm = T)
    }
    
    ggplot(alldat) +
      geom_point(aes(x = average_mbc, col = border), stat = 'density') +
      facet_wrap(.~source, ncol = 2) +
      labs(x = 'Pixel Count',
           y = 'Density') +
      scale_color_discrete(name = 'Border',
                           breaks = c('B1','B2','B3','B4','N'),
                           labels = c('One','Two','Three','Four','Interior')) +
      scale_y_continuous(breaks = seq(0,0.03,0.005)) +
      theme_linedraw() +
      theme(legend.position = 'bottom') +
      coord_cartesian(xlim = c(200,700),
                      ylim = c(0,0.03))
    ggsave(sprintf("%s%s_BORDERS_MBC_%d_%d.png",
                   out_path,expt_name,hr,pl),
           width = 7,height = 7.5)
    
    ##### WHAT HAPPENS NEAR GAPS
    alldat$gap <- 'N'
    for (o in alldat$pos) {
      c = alldat$`6144col`[alldat$pos == o]
      r = alldat$`6144row`[alldat$pos == o]
      if (alldat$colony[alldat$`6144col` == c & alldat$`6144row` == r] == 'Gap') {
        alldat$gap[alldat$`6144col` == c - 1 & alldat$`6144row` == r - 1 |
                     alldat$`6144col` == c & alldat$`6144row` == r - 1 |
                     alldat$`6144col` == c + 1 & alldat$`6144row` == r - 1 |
                     alldat$`6144col` == c - 1 & alldat$`6144row` == r |
                     alldat$`6144col` == c + 1 & alldat$`6144row` == r |
                     alldat$`6144col` == c - 1 & alldat$`6144row` == r + 1 |
                     alldat$`6144col` == c & alldat$`6144row` == r + 1 |
                     alldat$`6144col` == c + 1 & alldat$`6144row` == r + 1] = 'G1'
        alldat$gap[alldat$`6144col` == c - 2 & alldat$`6144row` == r - 2 |
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
                     alldat$`6144col` == c + 2 & alldat$`6144row` == r + 2] = 'G2'
      }
    }
    
    ggplot(alldat) +
      geom_point(aes(x = average_mbc, col = gap), stat = 'density') +
      facet_wrap(.~source, ncol = 2) +
      labs(x = 'Pixel Count',
           y = 'Density') +
      scale_color_discrete(name = 'Gaps',
                           breaks = c('G1','G2','N'),
                           labels = c('One','Two','Away')) +
      scale_y_continuous(breaks = seq(0,0.03,0.005)) +
      theme_linedraw() +
      theme(legend.position = 'bottom') +
      coord_cartesian(xlim = c(200,700),
                      ylim = c(0,0.02))
    ggsave(sprintf("%s%s_NEIGHS_%d_%d.png",
                   out_path,expt_name,hr,pl),
           width = 7,height = 7.5)
    
    ##### CORRECTING THE NEAR GAP POSITIONS
    alldat$average_mgc <- alldat$average_mbc
    
    for (sr in unique(alldat$source)) {
      alldat$average_mgc[alldat$gap == 'G1' & alldat$source == sr]  <- alldat$average_mgc[alldat$gap ==  "G1" & alldat$source == sr] *
        median(alldat$average_mgc[alldat$gap ==  "N" & alldat$source == sr], na.rm = T)/median(alldat$average_mgc[alldat$gap ==  "G1" & alldat$source == sr], na.rm = T)
      alldat$average_mgc[alldat$gap == 'G2' & alldat$source == sr]  <- alldat$average_mgc[alldat$gap ==  "G2" & alldat$source == sr] *
        median(alldat$average_mgc[alldat$gap ==  "N" & alldat$source == sr], na.rm = T)/median(alldat$average_mgc[alldat$gap ==  "G2" & alldat$source == sr], na.rm = T)
    }
    
    ggplot(alldat) +
      geom_point(aes(x = average_mgc, col = gap), stat = 'density') +
      facet_wrap(.~source, ncol = 2) +
      labs(x = 'Pixel Count',
           y = 'Density') +
      scale_color_discrete(name = 'Gaps',
                           breaks = c('G1','G2','N'),
                           labels = c('One','Two','Away')) +
      scale_y_continuous(breaks = seq(0,0.03,0.005)) +
      theme_linedraw() +
      theme(legend.position = 'bottom') +
      coord_cartesian(xlim = c(200,700),
                      ylim = c(0,0.02))
    ggsave(sprintf("%s%s_NEIGHS_MGC_%d_%d.png",
                   out_path,expt_name,hr,pl),
           width = 7,height = 7.5)
    
    jpegdat <- rbind(jpegdat, alldat)
  }
}

jpegdat <- data.frame(jpegdat$pos, jpegdat$hours, jpegdat$average_mgc)
dbWriteTable(conn, "4C3_GA1_MC_6144_JPEG", jpegdat, overwrite = T)

##### WHERE ARE THE OTHER EXTREMES
for (sr in unique(alldat$source)) {
  avg.ul <- mean(alldat$average_mgc[alldat$source == sr], na.rm = T) + 2*sd(alldat$average_mgc[alldat$source == sr], na.rm = T)
  avg.ll <- mean(alldat$average_mgc[alldat$source == sr], na.rm = T) - 2*sd(alldat$average_mgc[alldat$source == sr], na.rm = T)
  
  alldat$extreme[alldat$average_mgc > avg.ul & alldat$source == sr] <- 'Big'
  alldat$extreme[alldat$average_mgc < avg.ll & alldat$source == sr] <- 'Small' 
}

ggplot(alldat) +
  geom_point(data = alldat,
             aes(x=`6144col`, y=`6144row`,
                 col = extreme,
                 shape = colony),
             size = 3)

##### FITNESS FROM CORRECTED DATA
expt_name = '4C3_GA1_MC_BOR'
expt = 'FS1-1-MC-BOR'

##### CHECK POSITION WISE VARIABILITY
tablename_jpeg = sprintf('%s_%d_JPEG',expt_name,density);
tablename_fit = sprintf('%s_%d_FITNESS',expt_name,density);

# hr = 18
# pl = 1
for (hr in 17:18) {
  for (pl in 1:2) {
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
    
    alldat$border[alldat$`6144row` %in% 1 | alldat$`6144row` %in% 64 |
                    alldat$`6144col` %in% 1 | alldat$`6144col` %in% 96] = 'B1'
    alldat$border[alldat$`6144row` %in% 2 | alldat$`6144row` %in% 63 |
                    alldat$`6144col` %in% 2 | alldat$`6144col` %in% 95] = 'B2'
    alldat$border[alldat$`6144row` %in% 3 | alldat$`6144row` %in% 62 |
                    alldat$`6144col` %in% 3 | alldat$`6144col` %in% 94] = 'B3'
    alldat$border[alldat$`6144row` %in% 4 | alldat$`6144row` %in% 61 |
                    alldat$`6144col` %in% 4 | alldat$`6144col` %in% 93] = 'B4'
    alldat$border[is.na(alldat$border)] = 'N'
    
    alldat$gap <- 'N'
    for (o in alldat$pos) {
      c = alldat$`6144col`[alldat$pos == o]
      r = alldat$`6144row`[alldat$pos == o]
      if (alldat$colony[alldat$`6144col` == c & alldat$`6144row` == r] == 'Gap') {
        alldat$gap[alldat$`6144col` == c - 1 & alldat$`6144row` == r - 1 |
                     alldat$`6144col` == c & alldat$`6144row` == r - 1 |
                     alldat$`6144col` == c + 1 & alldat$`6144row` == r - 1 |
                     alldat$`6144col` == c - 1 & alldat$`6144row` == r |
                     alldat$`6144col` == c + 1 & alldat$`6144row` == r |
                     alldat$`6144col` == c - 1 & alldat$`6144row` == r + 1 |
                     alldat$`6144col` == c & alldat$`6144row` == r + 1 |
                     alldat$`6144col` == c + 1 & alldat$`6144row` == r + 1] = 'G1'
        alldat$gap[alldat$`6144col` == c - 2 & alldat$`6144row` == r - 2 |
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
                     alldat$`6144col` == c + 2 & alldat$`6144row` == r + 2] = 'G2'
      }
    }
    
    ggplot(alldat) +
      geom_point(aes(x = fitness, col = gap), stat = 'density') +
      facet_wrap(.~source, ncol = 2) +
      labs(x = 'Relative Fitness',
           y = 'Density') +
      scale_color_discrete(name = 'Gaps',
                           breaks = c('G1','G2','N'),
                           labels = c('One','Two','Away')) +
      theme_linedraw() +
      theme(legend.position = 'bottom') +
      coord_cartesian(xlim = c(0.5,1.5))
    ggsave(sprintf("%s%s_GAPFIT_MGC_%d_%d.png",
                   out_path,expt_name,hr,pl),
           width = 7,height = 7.5)
    
    ggplot(alldat) +
      geom_point(aes(x = fitness, col = border), stat = 'density') +
      facet_wrap(.~source, ncol = 2) +
      labs(x = 'Relative Fitness',
           y = 'Density') +
      scale_color_discrete(name = 'Border',
                           breaks = c('B1','B2','B3','B4','N'),
                           labels = c('One','Two','Three','Four','Interior')) +
      theme_linedraw() +
      theme(legend.position = 'bottom') +
      coord_cartesian(xlim = c(0.5,1.5))
    ggsave(sprintf("%s%s_BORFIT_MGC_%d_%d.png",
                   out_path,expt_name,hr,pl),
           width = 7,height = 7.5)
  }
}


for (sr in unique(alldat$source)) {
  fit.ul <- mean(alldat$fitness[alldat$source == sr], na.rm = T) + 2*sd(alldat$fitness[alldat$source == sr], na.rm = T)
  fit.ll <- mean(alldat$fitness[alldat$source == sr], na.rm = T) - 2*sd(alldat$fitness[alldat$source == sr], na.rm = T)
  
  alldat$extreme[alldat$fitness > fit.ul & alldat$source == sr] <- 'Big'
  alldat$extreme[alldat$fitness < fit.ll & alldat$source == sr] <- 'Small' 
}

ggplot(alldat[alldat$fitness > .7 & alldat$fitness < 1.3,]) +
  geom_point(aes(x=fitness),stat = 'density')

ggplot() +
  geom_point(data = alldat,
             aes(x=`6144col`, y=`6144row`),
             size = 0.2) +
  geom_point(data = alldat[alldat$fitness < 0.8 | alldat$fitness > 1.2,],
             aes(x=`6144col`, y=`6144row`, 
                 shape = colony),
             size = 3)
