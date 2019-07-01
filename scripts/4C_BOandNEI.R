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
expt_name = '4C3_GA1'
expt = 'FS1-1'
out_path = 'figs/borderandneigh/';
density = 6144;

##### CHECK POSITION WISE VARIABILITY
conn <- initialize.sql("saurin_test")

tablename_jpeg = sprintf('%s_%d_JPEG',expt_name,density);
tablename_fit = sprintf('%s_%d_FITNESS',expt_name,density);
tablename_fit_mcg = sprintf('%s_MCG_%d_FITNESS',expt_name,density);
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
    
    # alldat$border[alldat$`6144row` %in% 1 | alldat$`6144row` %in% 64 |
    #                 alldat$`6144col` %in% 1 | alldat$`6144col` %in% 96] = 'B1'
    # alldat$border[alldat$`6144row` %in% 2 | alldat$`6144row` %in% 63 |
    #                 alldat$`6144col` %in% 2 | alldat$`6144col` %in% 95] = 'B2'
    # alldat$border[alldat$`6144row` %in% 3 | alldat$`6144row` %in% 62 |
    #                 alldat$`6144col` %in% 3 | alldat$`6144col` %in% 94] = 'B3'
    # alldat$border[alldat$`6144row` %in% 4 | alldat$`6144row` %in% 61 |
    #                 alldat$`6144col` %in% 4 | alldat$`6144col` %in% 93] = 'B4'
    # alldat$border[is.na(alldat$border)] = 'N'
    # 
    # # ##### MEDIAN CORRECTING THE BORDERS
    # alldat$average_mbc <- alldat$average
    # 
    # for (sr in unique(alldat$source)) {
    #   alldat$average_mbc[alldat$border == 'B1' & alldat$source == sr]  <- alldat$average_mbc[alldat$border ==  "B1" & alldat$source == sr] *
    #     median(alldat$average_mbc[alldat$border ==  "N" & alldat$source == sr], na.rm = T)/median(alldat$average_mbc[alldat$border ==  "B1" & alldat$source == sr], na.rm = T)
    #   alldat$average_mbc[alldat$border == 'B2' & alldat$source == sr]  <- alldat$average_mbc[alldat$border ==  "B2" & alldat$source == sr] *
    #     median(alldat$average_mbc[alldat$border ==  "N" & alldat$source == sr], na.rm = T)/median(alldat$average_mbc[alldat$border ==  "B2" & alldat$source == sr], na.rm = T)
    #   alldat$average_mbc[alldat$border == 'B3' & alldat$source == sr]  <- alldat$average_mbc[alldat$border ==  "B3" & alldat$source == sr] *
    #     median(alldat$average_mbc[alldat$border ==  "N" & alldat$source == sr], na.rm = T)/median(alldat$average_mbc[alldat$border ==  "B3" & alldat$source == sr], na.rm = T)
    #   alldat$average_mbc[alldat$border == 'B4' & alldat$source == sr]  <- alldat$average_mbc[alldat$border ==  "B4" & alldat$source == sr] *
    #     median(alldat$average_mbc[alldat$border ==  "N" & alldat$source == sr], na.rm = T)/median(alldat$average_mbc[alldat$border ==  "B4" & alldat$source == sr], na.rm = T)
    # }
    # 
    # plot.raw <- ggplot(alldat) +
    #   geom_point(aes(x = average, col = border), stat = 'density') +
    #   facet_wrap(.~source, ncol = 2) +
    #   labs(title = 'Raw Data',
    #        x = 'Pixel Count',
    #        y = 'Density') +
    #   scale_color_discrete(name = 'Borders',
    #                        breaks = c('B1','B2','B3','B4','N'),
    #                        labels = c('One','Two','Three','Four','Away')) +
    #   scale_y_continuous(breaks = seq(0,0.03,0.005)) +
    #   theme_linedraw() +
    #   theme(legend.position = 'bottom') +
    #   coord_cartesian(xlim = c(200,600),
    #                   ylim = c(0,0.022))
    # 
    # plot.mcg <- ggplot(alldat) +
    #   geom_point(aes(x = average_mbc, col = border), stat = 'density') +
    #   facet_wrap(.~source, ncol = 2) +
    #   labs(title = 'Median Corrected Gaps',
    #        x = 'Pixel Count',
    #        y = 'Density') +
    #   scale_color_discrete(name = 'Gaps',
    #                        breaks = c('G1','G2','N'),
    #                        labels = c('One','Two','Away')) +
    #   scale_y_continuous(breaks = seq(0,0.03,0.005)) +
    #   theme_linedraw() +
    #   theme(legend.position = 'bottom') +
    #   coord_cartesian(xlim = c(200,600),
    #                   ylim = c(0,0.022))
    # ggarrange(plot.raw, plot.mcg,
    #           common.legend = T, legend = 'bottom')
    # ggsave(sprintf("%s%s_NEIGHS_MBC_%d_%d.png",
    #                out_path,expt_name,hr,pl),
    #        width = 14,height = 7.5)
    
    ##### WHAT HAPPENS NEAR GAPS
    alldat$average[alldat$colony == "Gap"] <- 0
    alldat$gap <- 'N'
    for (o in alldat$pos) {
      c = alldat$`6144col`[alldat$pos == o]
      r = alldat$`6144row`[alldat$pos == o]
      if (alldat$colony[alldat$`6144col` == c & alldat$`6144row` == r] == 'Gap') {
      # sr = alldat$source[alldat$`6144col` == c & alldat$`6144row` == r]
      # l = quantile(alldat$average[alldat$source == sr & alldat$orf_name == 'BF_control'], 0.1, na.rm = T)[[1]]
      # if (!is.na(alldat$average[alldat$`6144col` == c & alldat$`6144row` == r])) {
      #   if (alldat$colony[alldat$`6144col` == c & alldat$`6144row` == r] == 'Gap' |
      #       alldat$average[alldat$`6144col` == c & alldat$`6144row` == r] < l) {
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
        # }
      }
    }
    
    for (o in alldat$pos) {
      c = alldat$`6144col`[alldat$pos == o]
      r = alldat$`6144row`[alldat$pos == o]
      if (alldat$colony[alldat$`6144col` == c & alldat$`6144row` == r] == 'Gap') {
      # sr = alldat$source[alldat$`6144col` == c & alldat$`6144row` == r]
      # l = quantile(alldat$average[alldat$source == sr & alldat$orf_name == 'BF_control'], 0.1, na.rm = T)[[1]]
      # if (!is.na(alldat$average[alldat$`6144col` == c & alldat$`6144row` == r])) {
      #   if (alldat$colony[alldat$`6144col` == c & alldat$`6144row` == r] == 'Gap' |
      #       alldat$average[alldat$`6144col` == c & alldat$`6144row` == r] < l) {
          alldat$gap[alldat$`6144col` == c - 1 & alldat$`6144row` == r - 1 |
                       alldat$`6144col` == c & alldat$`6144row` == r - 1 |
                       alldat$`6144col` == c + 1 & alldat$`6144row` == r - 1 |
                       alldat$`6144col` == c - 1 & alldat$`6144row` == r |
                       alldat$`6144col` == c + 1 & alldat$`6144row` == r |
                       alldat$`6144col` == c - 1 & alldat$`6144row` == r + 1 |
                       alldat$`6144col` == c & alldat$`6144row` == r + 1 |
                       alldat$`6144col` == c + 1 & alldat$`6144row` == r + 1] = 'G1'
        # }
      }
    }
    
    ggplot(alldat[alldat$average > 0,]) +
      geom_point(aes(x = `6144col`, y = `6144row`, shape = colony, col = gap)) +
      scale_x_continuous(breaks = seq(1,96,1),limits = c(1,96)) +
      scale_y_continuous(breaks = seq(1,64,1),limits = c(64,1),trans = 'reverse') +
      # scale_color_discrete(guide = F) +
      scale_shape_discrete(guide = F) +
      theme_linedraw() +
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank())
    
    ##### CORRECTING THE NEAR GAP POSITIONS
    # alldat$average_mgc <- alldat$average_mbc
    alldat$average_mgc <- alldat$average
    
    for (sr in unique(alldat$source)) {
      alldat$average_mgc[alldat$gap == 'G1' & alldat$source == sr]  <- alldat$average_mgc[alldat$gap ==  "G1" & alldat$source == sr] *
        median(alldat$average_mgc[alldat$gap ==  "N" & alldat$source == sr & alldat$orf_name == "BF_control"], na.rm = T)/
        median(alldat$average_mgc[alldat$gap ==  "G1" & alldat$source == sr & alldat$orf_name == "BF_control"], na.rm = T)
      
      alldat$average_mgc[alldat$gap == 'G2' & alldat$source == sr]  <- alldat$average_mgc[alldat$gap ==  "G2" & alldat$source == sr] *
        median(alldat$average_mgc[alldat$gap ==  "N" & alldat$source == sr & alldat$orf_name == "BF_control"], na.rm = T)/
        median(alldat$average_mgc[alldat$gap ==  "G2" & alldat$source == sr & alldat$orf_name == "BF_control"], na.rm = T)
    }
    
    plot.raw <- ggplot(alldat) +
      geom_line(aes(x = average, col = gap), stat = 'density', lwd = 1.2) +
      facet_wrap(.~source, ncol = 2) +
      labs(title = 'Raw Data',
           x = 'Pixel Count',
           y = 'Density') +
      scale_color_discrete(name = 'Gaps',
                           breaks = c('G1','G2','N'),
                           labels = c('One','Two','Away')) +
      scale_y_continuous(breaks = seq(0,0.03,0.005)) +
      theme_linedraw() +
      theme(legend.position = 'bottom') +
      coord_cartesian(xlim = c(200,600),
                      ylim = c(0,0.022))
    # 
    # plot.mcg <- ggplot(alldat) +
    #   geom_line(aes(x = average_mgc, col = gap), stat = 'density', lwd = 1.2) +
    #   facet_wrap(.~source, ncol = 2) +
    #   labs(title = 'Median Corrected Gaps',
    #        x = 'Pixel Count',
    #        y = 'Density') +
    #   scale_color_discrete(name = 'Gaps',
    #                        breaks = c('G1','G2','N'),
    #                        labels = c('One','Two','Away')) +
    #   scale_y_continuous(breaks = seq(0,0.03,0.005)) +
    #   theme_linedraw() +
    #   theme(legend.position = 'bottom') +
    #   coord_cartesian(xlim = c(200,600),
    #                   ylim = c(0,0.022))
    # ggarrange(plot.raw, plot.mcg,
    #           common.legend = T, legend = 'bottom')
    # ggsave(sprintf("%s%s_NEIGHS_MGC_%d_%d.png",
    #                out_path,expt_name,hr,pl),
    #        width = 14,height = 7.5)
    
    jpegdat <- rbind(jpegdat, alldat)
  }
}

jpegdat <- data.frame(jpegdat$pos, jpegdat$hours, jpegdat$average_mgc)
colnames(jpegdat) <- c('pos','hours','average')
dbWriteTable(conn, "4C3_GA1_MCG_6144_JPEG", jpegdat, overwrite = T)

##### BEFORE AND AFTER
alldat = dbGetQuery(conn, sprintf('select a.*, b.*, c.*
                                      from %s a, %s b, %s c
                                      where a.pos = b.pos
                                      and b.pos = c.pos
                                      order by a.hours, c.%s, c.%s, c.%s',
                                  tablename_fit, tablename_fit_mcg,
                                  p2c_info[1],p2c_info[2],
                                  p2c_info[3],p2c_info[4]))

alldat$source[alldat$`6144row`%%2==1 & alldat$`6144col`%%2==1] = '1TL'
alldat$source[alldat$`6144row`%%2==0 & alldat$`6144col`%%2==1] = '3BL'
alldat$source[alldat$`6144row`%%2==1 & alldat$`6144col`%%2==0] = '2TR'
alldat$source[alldat$`6144row`%%2==0 & alldat$`6144col`%%2==0] = '4BR'

alldat$colony[alldat$orf_name == 'BF_control'] = 'Reference'
alldat$colony[alldat$orf_name != 'BF_control'] = 'Query'
alldat$colony[is.na(alldat$orf_name)] = 'Gap'

