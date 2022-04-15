##### DELETIONS SCREEN DATA
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 04/05/2022

source('/home/sbp29/R/Projects/adaptivefitness/scripts/translatome/initialize.R')
source('/home/sbp29/R/Projects/adaptivefitness/scripts/translatome/fit_sum.R')

########################## INIT DATA ##################################
##### COLONY SIZE TO FITNESS DATA
tr.conds <- data.frame(arms = c('R1','R1','R1','R1','R1','R2','R1'),
                       conds = c('YPDA','DM','HU','HO','TN','FL','SA'))

data.fit.init <- NULL
for (c in unique(tr.conds$conds)) {
  for (a in tr.conds$arms[tr.conds$conds == c]) {
    temp.fit <- dbGetQuery(conn, sprintf('select a.*, b.density, b.plate_no, b.plate_col, b.plate_row
                                         from TR_DEL_FS_%s_%s_1536_FITNESS a, TR_DEL_pos2coor b
                                         where a.pos = b.pos
                                         order by a.hours, b.plate_no, b.plate_col, b.plate_row', c, a))
    temp.fit$condition <- c
    data.fit.init <- rbind(data.fit.init, temp.fit)
  }
}
data.fit.init$rep <- as.numeric(str_trunc(as.character(data.fit.init$pos), 4, side = 'left', ellipsis = ''))

data.fit.init <- data.fit.init[!(data.fit.init$condition == 'FL' & data.fit.init$hours == 192) &
                                 !(data.fit.init$condition == 'GA' & data.fit.init$hours == 67) &
                                 !(data.fit.init$condition == 'SA' & data.fit.init$hours %in% c(51,63,72)),]

temp <- fit_sum(data.fit.init, orf_types, 4)
data.fit.init <- temp[[1]]
data.fit.init.sum <- temp[[2]]

data.fit.init$category[data.fit.init$orf_name %in% c('BOR','NULL')] <- NA
data.fit.init$category[data.fit.init$orf_name %in% c('HO')] <- 'Reference'
data.fit.init$category[data.fit.init$orf_name %notin% c('BOR','NULL','HO')] <- 'KanMX'
data.fit.init.sum$category[data.fit.init.sum$orf_name %in% c('BOR','NULL')] <- NA
data.fit.init.sum$category[data.fit.init.sum$orf_name %in% c('HO')] <- 'Reference'
data.fit.init.sum$category[data.fit.init.sum$orf_name %notin% c('BOR','NULL','HO')] <- 'KanMX'

dbWriteTable(conn, 'TR_DEL_INIT_CS_FIT', data.fit.init, overwrite = T)
dbWriteTable(conn, 'TR_DEL_INIT_CS_FIT_SUMMARY', data.fit.init.sum, overwrite = T)

###### COLONY SIZE TO GROWTH TO FITNESS
data.growth.init <- dbGetQuery(conn, 'select a.*, b.density, b.plate_no, b.plate_col, b.plate_row
                                         from TR_DEL_GROWTH_FS_ALL_1536_FITNESS a, TR_DEL_pos2coor b
                                         where a.pos = b.pos
                                         order by a.hours, b.plate_no, b.plate_col, b.plate_row')
data.growth.init$condition[data.growth.init$hours == 1] <- 'DM'
data.growth.init$condition[data.growth.init$hours == 2] <- 'FL'
data.growth.init$condition[data.growth.init$hours == 3] <- 'HO'
data.growth.init$condition[data.growth.init$hours == 4] <- 'HU'
data.growth.init$condition[data.growth.init$hours == 5] <- 'SA'
data.growth.init$condition[data.growth.init$hours == 6] <- 'TN'
data.growth.init$condition[data.growth.init$hours == 7] <- 'YPDA'

temp <- fit_sum(data.growth.init, orf_types, 4)
data.growth.init <- temp[[1]]
data.growth.init.sum <- temp[[2]]

data.growth.init$category[data.growth.init$orf_name %in% c('BOR','NULL')] <- NA
data.growth.init$category[data.growth.init$orf_name %in% c('HO')] <- 'Reference'
data.growth.init$category[data.growth.init$orf_name %notin% c('BOR','NULL','HO')] <- 'KanMX'
data.growth.init.sum$category[data.growth.init.sum$orf_name %in% c('BOR','NULL')] <- NA
data.growth.init.sum$category[data.growth.init.sum$orf_name %in% c('HO')] <- 'Reference'
data.growth.init.sum$category[data.growth.init.sum$orf_name %notin% c('BOR','NULL','HO')] <- 'KanMX'

dbWriteTable(conn, 'TR_DEL_INIT_GR_FIT', data.growth.init, overwrite = T)
dbWriteTable(conn, 'TR_DEL_INIT_GR_FIT_SUMMARY', data.growth.init.sum, overwrite = T)


########################## VALID DATA ##################################
##### COLONY SIZE TO FITNESS
tr.conds <- read_xlsx(path = '/home/sbp29/RAW_Data/TranslatomeDEL_VAL/TR_DEL_VAL_MS_INFO.xlsx') %>% data.frame()
tr.conds <- tr.conds %>% filter(stage_id == 'FS', density == 1536)

data.fit.valid <- NULL
for (a in unique(tr.conds$arm)) {
  temp.fit <- dbGetQuery(conn, sprintf('select a.*, b.density, b.plate_no, b.plate_col, b.plate_row
                                         from TR_DEL_VAL_MS1536_FS_%s_1536_FITNESS a, TR_DEL_VAL_MS1536_pos2coor b
                                         where a.pos = b.pos
                                         order by a.hours, b.plate_no, b.plate_col, b.plate_row', a))
  temp.fit$condition <- a
  data.fit.valid <- rbind(data.fit.valid, temp.fit)
}

temp <- fit_sum(data.fit.valid, orf_types, 4)
data.fit.valid <- temp[[1]]
data.fit.valid.sum <- temp[[2]]

data.fit.valid$category[data.fit.valid$orf_name %in% c('BOR','NULL')] <- NA
data.fit.valid$category[data.fit.valid$orf_name %in% c('HO')] <- 'Reference'
data.fit.valid$category[data.fit.valid$orf_name %notin% c('BOR','NULL','HO')] <- 'KanMX'
data.fit.valid.sum$category[data.fit.valid.sum$orf_name %in% c('BOR','NULL')] <- NA
data.fit.valid.sum$category[data.fit.valid.sum$orf_name %in% c('HO')] <- 'Reference'
data.fit.valid.sum$category[data.fit.valid.sum$orf_name %notin% c('BOR','NULL','HO')] <- 'KanMX'

dbWriteTable(conn, 'TR_DEL_VALID_CS_FIT', data.fit.valid, overwrite = T)
dbWriteTable(conn, 'TR_DEL_VALID_CS_FIT_SUMMARY', data.fit.valid.sum, overwrite = T)

###### COLONY SIZE TO GROWTH TO FITNESS
data.growth.valid <- dbGetQuery(conn, 'select a.*, b.density, b.plate_no, b.plate_col, b.plate_row
                                         from TR_DEL_VAL_MS_GROWTH_FS_ALL_1536_FITNESS a, TR_DEL_VAL_MS1536_pos2coor b
                                         where a.pos = b.pos
                                         order by a.hours, b.plate_no, b.plate_col, b.plate_row')
data.growth.valid$condition[data.growth.valid$hours == 1] <- 'DM'
data.growth.valid$condition[data.growth.valid$hours == 2] <- 'FL'
data.growth.valid$condition[data.growth.valid$hours == 3] <- 'HO'
data.growth.valid$condition[data.growth.valid$hours == 4] <- 'HU'
data.growth.valid$condition[data.growth.valid$hours == 5] <- 'SA'
data.growth.valid$condition[data.growth.valid$hours == 6] <- 'TN'
data.growth.valid$condition[data.growth.valid$hours == 7] <- 'YPDA'

temp <- fit_sum(data.growth.valid, orf_types, 4)
data.growth.valid <- temp[[1]]
data.growth.valid.sum <- temp[[2]]

data.growth.valid$category[data.growth.valid$orf_name %in% c('BOR','NULL')] <- NA
data.growth.valid$category[data.growth.valid$orf_name %in% c('HO')] <- 'Reference'
data.growth.valid$category[data.growth.valid$orf_name %notin% c('BOR','NULL','HO')] <- 'KanMX'
data.growth.valid.sum$category[data.growth.valid.sum$orf_name %in% c('BOR','NULL')] <- NA
data.growth.valid.sum$category[data.growth.valid.sum$orf_name %in% c('HO')] <- 'Reference'
data.growth.valid.sum$category[data.growth.valid.sum$orf_name %notin% c('BOR','NULL','HO')] <- 'KanMX'

dbWriteTable(conn, 'TR_DEL_VALID_GR_FIT', data.growth.valid, overwrite = T)
dbWriteTable(conn, 'TR_DEL_VALID_GR_FIT_SUMMARY', data.growth.valid.sum, overwrite = T)


##################### AAG/ATG MUTANT DATA ##########################
##### COLONY SIZE TO FITNESS
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
# data.fit.aag$attempt[data.fit.aag$plate_no %in% c(1:7)]   <- 'ONE'
# data.fit.aag$attempt[data.fit.aag$plate_no %in% c(8:14)]  <- 'TWO'

temp <- fit_sum(data.fit.aag, orf_types, 4)
data.fit.aag <- temp[[1]]
data.fit.aag.sum <- temp[[2]]

data.fit.aag$category[data.fit.aag$strain_id > 100000 & data.fit.aag$strain_id < 200000] <- 'KanMX'
data.fit.aag$category[data.fit.aag$strain_id > 200000 & data.fit.aag$strain_id < 300000] <- 'AAG'
data.fit.aag$category[data.fit.aag$strain_id > 300000] <- 'WT'
data.fit.aag$category[data.fit.aag$orf_name == 'BY4741'] <- 'Reference'
data.fit.aag$orf_name[data.fit.aag$orf_name == 'BY4741_HO'] <- 'HO'
data.fit.aag$category[data.fit.aag$orf_name == 'HO'] <- 'Reference2'

data.fit.aag.sum$category[data.fit.aag.sum$strain_id > 100000 & data.fit.aag.sum$strain_id < 200000] <- 'KanMX'
data.fit.aag.sum$category[data.fit.aag.sum$strain_id > 200000 & data.fit.aag.sum$strain_id < 300000] <- 'AAG'
data.fit.aag.sum$category[data.fit.aag.sum$strain_id > 300000] <- 'WT'
data.fit.aag.sum$category[data.fit.aag.sum$orf_name == 'BY4741'] <- 'Reference'
data.fit.aag.sum$orf_name[data.fit.aag.sum$orf_name == 'BY4741_HO'] <- 'HO'
data.fit.aag.sum$category[data.fit.aag.sum$orf_name == 'HO'] <- 'Reference2'

##################### COMBINING ALL DATA ###########################
##### FITNESS DATA
data.fit <- rbind(cbind(data.fit.init, attempt = 'init', data = 'cs'),
                  cbind(data.fit.valid, attempt = 'valid', data = 'cs'),
                  cbind(data.growth.init, attempt = 'init', data = 'growth'),
                  cbind(data.growth.valid, attempt = 'valid', data = 'growth'),
                  cbind(data.fit.aag, attempt = 'aag', data = 'cs'))

data.fit$pos <- as.numeric(data.fit$pos)
data.fit$rep <- as.numeric(data.fit$rep)

data.fit <- arrange(data.fit,data,attempt,condition,pos,hours) %>% 
  mutate(avg_ratio=(lag(average,0)/lag(average,1)))
data.fit$avg_ratio[data.fit$data != 'cs'] <- NA
head(data.fit)

dbWriteTable(conn, 'TRANSLATOME_DEL_FITNESS_DATA', data.fit, overwrite = T)

##### FITNESS SUMMARY DATA
data.fit.sum <- rbind(cbind(data.fit.init.sum, attempt = 'init', data = 'cs'),
                      cbind(data.fit.valid.sum, attempt = 'valid', data = 'cs'),
                      cbind(data.growth.init.sum, attempt = 'init', data = 'growth'),
                      cbind(data.growth.valid.sum, attempt = 'valid', data = 'growth'),
                      cbind(data.fit.aag.sum, attempt = 'aag', data = 'cs'))
data.fit.sum$rep <- as.numeric(data.fit.sum$rep)

# data.fit.lim <- data.fit.sum %>%
#   filter(orf_name == 'HO') %>%
#   group_by(data, attempt, condition, hours) %>%
#   summarise(fit.ll = quantile(fit.median, 0.025, na.rm = T),
#             fit.m = median(fit.median, na.rm = T),
#             fit.ul = quantile(fit.median, 0.975, na.rm = T),
#             .groups = 'keep') %>%
#   data.frame()
# data.fit.sum <- merge(data.fit.sum, data.fit.lim, 
#                       by = c('data','attempt','condition','hours'))
# 
# data.fit.sum$phenotype[data.fit.sum$fit.median > data.fit.sum$fit.ul] <- 'Beneficial'
# data.fit.sum$phenotype[data.fit.sum$fit.median < data.fit.sum$fit.ll] <- 'Deleterious'
# data.fit.sum$phenotype[is.na(data.fit.sum$phenotype) & !is.na(data.fit.sum$fit.median)] <- 'Neutral'
# head(data.fit.sum)

dbWriteTable(conn, 'TRANSLATOME_DEL_FITNESS_SUMMARY_DATA', data.fit.sum, overwrite = T)

