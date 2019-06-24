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

hr = 18
pl = 1

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

# alldat$border[alldat$`6144row` %in% 1:4 | alldat$`6144row` %in% 61:64 |
#                 alldat$`6144col` %in% 1:4 | alldat$`6144col` %in% 93:96] = 'Y'
# alldat$border[is.na(alldat$border)] = 'N'

alldat$border[alldat$`6144row` %in% 1 | alldat$`6144row` %in% 64 |
                alldat$`6144col` %in% 1 | alldat$`6144col` %in% 96] = 'B1'
alldat$border[alldat$`6144row` %in% 2 | alldat$`6144row` %in% 63 |
                alldat$`6144col` %in% 2 | alldat$`6144col` %in% 95] = 'B2'
alldat$border[alldat$`6144row` %in% 3 | alldat$`6144row` %in% 62 |
                alldat$`6144col` %in% 3 | alldat$`6144col` %in% 94] = 'B3'
alldat$border[alldat$`6144row` %in% 4 | alldat$`6144row` %in% 61 |
                alldat$`6144col` %in% 4 | alldat$`6144col` %in% 93] = 'B4'
alldat$border[is.na(alldat$border)] = 'N'

ggplot(alldat) +
  geom_point(aes(x = average, col = border), stat = 'density') +
  facet_wrap(.~source, ncol = 2) +
  labs(x = 'Pixel Count',
       y = 'Density') +
  scale_color_discrete(name = 'Border',
                     breaks = c('B1','B2','B3','B4','N'),
                     labels = c('One','Two','Three','Four','Interior')) +
  theme_linedraw() +
  theme(legend.position = 'bottom')

ggsave(sprintf("%s%s_BORDERS_%d_%d.png",
               out_path,expt_name,hr,pl),
       width = 10,height = 11)

alldat$average_mbc <- alldat$average
alldat$average_mbc[alldat$border == 'B1']  <- alldat$average_mbc[alldat$border ==  "B1"] *
  median(alldat$average_mbc[alldat$border ==  "N"], na.rm = T)/median(alldat$average_mbc[alldat$border ==  "B1"], na.rm = T)
alldat$average_mbc[alldat$border == 'B2']  <- alldat$average_mbc[alldat$border ==  "B2"] *
  median(alldat$average_mbc[alldat$border ==  "N"], na.rm = T)/median(alldat$average_mbc[alldat$border ==  "B2"], na.rm = T)
alldat$average_mbc[alldat$border == 'B3']  <- alldat$average_mbc[alldat$border ==  "B3"] *
  median(alldat$average_mbc[alldat$border ==  "N"], na.rm = T)/median(alldat$average_mbc[alldat$border ==  "B3"], na.rm = T)
alldat$average_mbc[alldat$border == 'B4']  <- alldat$average_mbc[alldat$border ==  "B4"] *
  median(alldat$average_mbc[alldat$border ==  "N"], na.rm = T)/median(alldat$average_mbc[alldat$border ==  "B4"], na.rm = T)

ggplot(alldat) +
  geom_point(aes(x = average_mbc, col = border), stat = 'density') +
  facet_grid(.~source)
