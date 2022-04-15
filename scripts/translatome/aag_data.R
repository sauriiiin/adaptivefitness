##### AAG-ATG MUTANT ANALYSIS
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 04/11/2022

source('/home/sbp29/R/Projects/adaptivefitness/scripts/translatome/initialize.R')
source('/home/sbp29/R/Projects/adaptivefitness/scripts/translatome/fit_sum.R')

data.fit.aag <- dbGetQuery(conn, 'select a.*, b.density, b.plate_no, b.plate_col, b.plate_row
                                         from TR_AAG_FS_ALL_1536_FITNESS a, TR_AAG_pos2coor b
                                         where a.pos = b.pos
                                         order by a.hours, b.plate_no, b.plate_col, b.plate_row')
data.fit.aag$hours <- data.fit.aag$hours + 2
data.fit.aag$condition[data.fit.aag$plate_no %in% c(1,8)]  <- 'YPDA'
data.fit.aag$condition[data.fit.aag$plate_no %in% c(2,9)]  <- 'DM'
data.fit.aag$condition[data.fit.aag$plate_no %in% c(3,10)] <- 'SA'
data.fit.aag$condition[data.fit.aag$plate_no %in% c(4,11)] <- 'FL'
data.fit.aag$condition[data.fit.aag$plate_no %in% c(5,12)] <- 'HU'
data.fit.aag$condition[data.fit.aag$plate_no %in% c(6,13)] <- 'TN'
data.fit.aag$condition[data.fit.aag$plate_no %in% c(7,14)] <- 'HO'

data.fit.aag$attempt[data.fit.aag$plate_no %in% c(1:7)]   <- 'ONE'
data.fit.aag$attempt[data.fit.aag$plate_no %in% c(8:14)]  <- 'TWO'

head(data.fit.aag)

temp <- fit_sum(data.fit.aag, orf_types, 4)
data.fit.aag <- temp[[1]]
data.fit.aag.sum <- temp[[2]]

data.fit.aag$category[data.fit.aag$strain_id > 100000 & data.fit.aag$strain_id < 200000] <- 'KanMX'
data.fit.aag$category[data.fit.aag$strain_id > 200000 & data.fit.aag$strain_id < 300000] <- 'AAG'
data.fit.aag$category[data.fit.aag$strain_id > 300000] <- 'WT'

data.fit.aag.sum$category[data.fit.aag.sum$strain_id > 100000 & data.fit.aag.sum$strain_id < 200000] <- 'KanMX'
data.fit.aag.sum$category[data.fit.aag.sum$strain_id > 200000 & data.fit.aag.sum$strain_id < 300000] <- 'AAG'
data.fit.aag.sum$category[data.fit.aag.sum$strain_id > 300000] <- 'WT'

dbWriteTable(conn, 'TR_AAG_CS_FIT', data.fit.aag, overwrite = T)
dbWriteTable(conn, 'TR_AAG_CS_FIT_SUMMARY', data.fit.aag.sum, overwrite = T)


data.fit.aag %>%
  filter(hours == 90, strain_id > 100000) %>%
  ggplot(aes(x = category, y = fitness)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(condition ~ orf_name)

