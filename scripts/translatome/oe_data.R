##### INITIAL VS VALIDATION DATA
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 03/23/2022

source('/home/sbp29/R/Projects/adaptivefitness/scripts/translatome/initialize.R')
source('/home/sbp29/R/Projects/adaptivefitness/scripts/translatome/fit_sum.R')

########################## INIT DATA ##################################

data.p2c.init <- dbGetQuery(conn, 'select * from TRANS_OE_pos2coor')
data.p2c.init$pos <- as.numeric(data.p2c.init$pos)
data.p2c.init <- data.p2c.init %>% filter(density == 384)
data.p2c.init$rep <- as.numeric(str_trunc(as.character(data.p2c.init$pos), 6, side = 'left', ellipsis = ''))
data.p2c.init <- data.p2c.init[,c(-1,-2)]
colnames(data.p2c.init) <- c('plate_no_384','plate_row_384','plate_col_384','rep')

##### COLONY SIZE TO FITNESS DATA
tr.conds <- data.frame(arms = c('ONE','ONE','ONE','ONE','TWO','TWO','TWO'),
                       conds = c('GA','SA','HO','HU','DM','FL','TN'))

smudge <- dbGetQuery(conn, 'select * from TR_OE_ONE_FS_GA_smudgebox')
data.fit.init <- NULL
for (a in unique(tr.conds$arms)) {
  for (c in tr.conds$conds[tr.conds$arms == a]) {
    temp.fit <- dbGetQuery(conn, sprintf('select a.*, b.density, b.plate_no, b.plate_col, b.plate_row
                                         from TR_OE_FS_%s_%s_6144_FITNESS a, TRANS_OE_pos2coor b
                                         where a.pos = b.pos
                                         order by a.hours, b.plate_no, b.plate_col, b.plate_row', a, c))
    temp.fit$condition <- c
    data.fit.init <- rbind(data.fit.init, temp.fit)
  }
}
data.fit.init$average[data.fit.init$condition == 'HO' & data.fit.init$plate %in% c(2,3,10,18,19,23,26,32)] <- NA
data.fit.init$average[data.fit.init$condition %in% c('GA','HO','HU','SA') & data.fit.init$plate == 14] <- NA
data.fit.init$average[data.fit.init$condition %in% c('GA','HO','HU','SA') & data.fit.init$pos %in% smudge$pos] <- NA
data.fit.init$fitness[data.fit.init$condition == 'HO' & data.fit.init$plate %in% c(2,3,10,18,19,23,26,32)] <- NA
data.fit.init$fitness[data.fit.init$condition %in% c('GA','HO','HU','SA') & data.fit.init$plate == 14] <- NA
data.fit.init$fitness[data.fit.init$condition %in% c('GA','HO','HU','SA') & data.fit.init$pos %in% smudge$pos] <- NA

data.fit.init <- data.fit.init[!(data.fit.init$condition == 'FL' & data.fit.init$hours == 192) &
                                 !(data.fit.init$condition == 'GA' & data.fit.init$hours == 67) &
                                 !(data.fit.init$condition == 'SA' & data.fit.init$hours %in% c(51,63,72)),]

temp <- fit_sum(data.fit.init, orf_types, 6)
data.fit.init <- temp[[1]]
data.fit.init.sum <- temp[[2]]

data.fit.init <- merge(data.fit.init, data.p2c.init, by = 'rep')

dbWriteTable(conn, 'TR_OE_INIT_CS_FIT', data.fit.init, overwrite = T)
dbWriteTable(conn, 'TR_OE_INIT_CS_FIT_SUMMARY', data.fit.init.sum, overwrite = T)

# ##### COLONY SIZE TO AUC TO FITNESS
# data.gc.fit.init <- dbGetQuery(conn, 'select a.*, b.density, b.plate_no, b.plate_col, b.plate_row
#                     from TR_OE_GC_6144_FITNESS a, TRANS_OE_pos2coor b
#                     where a.pos = b.pos
#                     order by a.hours, b.plate_no, b.plate_col, b.plate_row')
# 
# data.gc.fit.init$condition[data.gc.fit.init$hours == 1] <- 'GA'
# data.gc.fit.init$condition[data.gc.fit.init$hours == 2] <- 'SA'
# data.gc.fit.init$condition[data.gc.fit.init$hours == 3] <- 'HO'
# data.gc.fit.init$condition[data.gc.fit.init$hours == 4] <- 'HU'
# data.gc.fit.init$condition[data.gc.fit.init$hours == 5] <- 'DM'
# data.gc.fit.init$condition[data.gc.fit.init$hours == 6] <- 'FL'
# data.gc.fit.init$condition[data.gc.fit.init$hours == 7] <- 'TN'
# 
# temp <- fit_sum(data.gc.fit.init, orf_types, 6)
# data.gc.fit.init <- temp[[1]]
# data.gc.fit.init.sum <- temp[[2]]
# 
# data.gc.fit.init <- merge(data.gc.fit.init, data.p2c.lin, by = 'rep')
# 
# dbWriteTable(conn, 'TR_OE_INIT_CS_GC_FIT', data.gc.fit.init, overwrite = T)
# dbWriteTable(conn, 'TR_OE_INIT_CS_GC_FIT_SUMMARY', data.gc.fit.init.sum, overwrite = T)

##### COLONY SIZE TO FITNESS DATA (TWO GA)
data.fit.init.two.ga <- NULL
data.fit.init.two.ga <- dbGetQuery(conn, 'select a.*, b.density, b.plate_no, b.plate_col, b.plate_row
                                         from TR_OE_FS_TWO_GA_6144_FITNESS a, TRANS_OE_pos2coor b
                                         where a.pos = b.pos
                                         order by a.hours, b.plate_no, b.plate_col, b.plate_row')
data.fit.init.two.ga$condition <- 'GA'

temp <- fit_sum(data.fit.init.two.ga, orf_types, 6)
data.fit.init.two.ga <- temp[[1]]
data.fit.init.two.ga.sum <- temp[[2]]

data.fit.init.two.ga <- merge(data.fit.init.two.ga, data.p2c.init, by = 'rep')

dbWriteTable(conn, 'TR_OE_INIT_TWO_GA_CS_FIT', data.fit.init.two.ga, overwrite = T)
dbWriteTable(conn, 'TR_OE_INIT_TWO_GA_CS_FIT_SUMMARY', data.fit.init.two.ga.sum, overwrite = T)

###### COLONY SIZE TO GROWTH TO FITNESS
data.growth.init <- dbGetQuery(conn, 'select a.*, b.density, b.plate_no, b.plate_col, b.plate_row
                                         from TR_OE_GROWTH_FS_ALL_6144_FITNESS a, TRANS_OE_pos2coor b
                                         where a.pos = b.pos
                                         order by a.hours, b.plate_no, b.plate_col, b.plate_row')
data.growth.init$condition[data.growth.init$hours == 1] <- 'DM'
data.growth.init$condition[data.growth.init$hours == 2] <- 'FL'
data.growth.init$condition[data.growth.init$hours == 3] <- 'GA'
data.growth.init$condition[data.growth.init$hours == 4] <- 'HO'
data.growth.init$condition[data.growth.init$hours == 5] <- 'HU'
data.growth.init$condition[data.growth.init$hours == 6] <- 'SA'
data.growth.init$condition[data.growth.init$hours == 7] <- 'TN'

temp <- fit_sum(data.growth.init, orf_types, 6)
data.growth.init <- temp[[1]]
data.growth.init.sum <- temp[[2]]

data.growth.init <- merge(data.growth.init, data.p2c.init, by = 'rep')

dbWriteTable(conn, 'TR_OE_INIT_GR_FIT', data.growth.init, overwrite = T)
dbWriteTable(conn, 'TR_OE_INIT_GR_FIT_SUMMARY', data.growth.init.sum, overwrite = T)


########################## VALID DATA ##################################

data.p2c.valid <- dbGetQuery(conn, 'select * from TR_OE_VAL_MS_pos2coor')
data.p2c.valid <- data.p2c.valid %>% filter(density == 384)
data.p2c.valid$rep <- as.numeric(str_trunc(as.character(data.p2c.valid$pos), 6, side = 'left', ellipsis = ''))
data.p2c.valid <- data.p2c.valid[,c(-1,-2)]
colnames(data.p2c.valid) <- c('plate_no_384','plate_row_384','plate_col_384','rep')

##### COLONY SIZE TO FITNESS
tr.conds <- readxl::read_xlsx(path = '/home/sbp29/RAW_Data/TranslatomeOE_VAL/TR_OE_VAL_MS_INFO2.xlsx') %>% data.frame()
data.fit.valid <- NULL
for (a in unique(tr.conds$arm)) {
  temp.cs <- dbGetQuery(conn, sprintf('select a.*, b.density, b.plate_no, b.plate_col, b.plate_row
                                         from TR_OE_VAL_MS_FS_%s_6144_FITNESS a, TR_OE_VAL_MS_pos2coor b
                                         where a.pos = b.pos
                                         order by a.hours, b.plate_no, b.plate_col, b.plate_row', a))
  temp.cs$condition <- a
  data.fit.valid <- rbind(data.fit.valid, temp.cs)
}
temp <- fit_sum(data.fit.valid, orf_types, 6)
data.fit.valid <- temp[[1]]
data.fit.valid.sum <- temp[[2]]

data.fit.valid <- merge(data.fit.valid, data.p2c.valid, by = 'rep')

dbWriteTable(conn, 'TR_OE_VALID_CS_FIT', data.fit.valid, overwrite = T)
dbWriteTable(conn, 'TR_OE_VALID_CS_FIT_SUMMARY', data.fit.valid.sum, overwrite = T)

##### COLONY SIZE TO AUC TO FITNESS
data.gc.fit.valid <- dbGetQuery(conn, 'select a.*, b.density, b.plate_no, b.plate_col, b.plate_row
                    from TR_OE_VAL_MS_GC_6144_FITNESS a, TR_OE_VAL_MS_pos2coor b
                    where a.pos = b.pos
                    order by a.hours, b.plate_no, b.plate_col, b.plate_row')

data.gc.fit.valid$condition[data.gc.fit.valid$hours == 1] <- 'GA'
data.gc.fit.valid$condition[data.gc.fit.valid$hours == 2] <- 'SA'
data.gc.fit.valid$condition[data.gc.fit.valid$hours == 3] <- 'HU'
data.gc.fit.valid$condition[data.gc.fit.valid$hours == 4] <- 'DM'
data.gc.fit.valid$condition[data.gc.fit.valid$hours == 5] <- 'FL'
data.gc.fit.valid$condition[data.gc.fit.valid$hours == 6] <- 'TN'

temp <- fit_sum(data.gc.fit.valid, orf_types, 6)
data.gc.fit.valid <- temp[[1]]
data.gc.fit.valid.sum <- temp[[2]]

data.gc.fit.valid <- merge(data.gc.fit.valid, data.p2c.valid, by = 'rep')

dbWriteTable(conn, 'TR_OE_VALID_CS_GC_FIT', data.gc.fit.valid, overwrite = T)
dbWriteTable(conn, 'TR_OE_VALID_CS_GC_FIT_SUMMARY', data.gc.fit.valid.sum, overwrite = T)

###### COLONY SIZE TO GROWTH TO FITNESS
data.growth.valid <- dbGetQuery(conn, 'select a.*, b.density, b.plate_no, b.plate_col, b.plate_row
                                         from TR_OE_VAL_MS_GROWTH_FS_ALL_6144_FITNESS a, TR_OE_VAL_MS_pos2coor b
                                         where a.pos = b.pos
                                         order by a.hours, b.plate_no, b.plate_col, b.plate_row')
data.growth.valid$condition[data.growth.valid$hours == 1] <- 'DM'
data.growth.valid$condition[data.growth.valid$hours == 2] <- 'FL'
data.growth.valid$condition[data.growth.valid$hours == 3] <- 'GA'
data.growth.valid$condition[data.growth.valid$hours == 4] <- 'HU'
data.growth.valid$condition[data.growth.valid$hours == 5] <- 'SA'
data.growth.valid$condition[data.growth.valid$hours == 6] <- 'TN'

temp <- fit_sum(data.growth.valid, orf_types, 6)
data.growth.valid <- temp[[1]]
data.growth.valid.sum <- temp[[2]]

data.growth.valid <- merge(data.growth.valid, data.p2c.valid, by = 'rep')

dbWriteTable(conn, 'TR_OE_VALID_GR_FIT', data.growth.valid, overwrite = T)
dbWriteTable(conn, 'TR_OE_VALID_GR_FIT_SUMMARY', data.growth.valid.sum, overwrite = T)


########################## LIN DATA ##################################

data.p2c.lin <- dbGetQuery(conn, 'select * from LIN_OE_2202_pos2coor')
data.p2c.lin <- data.p2c.lin %>% filter(density == 384)
data.p2c.lin$rep <- as.numeric(str_trunc(as.character(data.p2c.lin$pos), 6, side = 'left', ellipsis = ''))
data.p2c.lin <- data.p2c.lin[,c(-1,-2)]
colnames(data.p2c.lin) <- c('plate_no_384','plate_row_384','plate_col_384','rep')

##### COLONY SIZE TO FITNESS
tr.conds <- readxl::read_xlsx(path = '/home/sbp29/RAW_Data/LinOE2202/LIN_OE_2202_INFO.xlsx') %>% data.frame() %>% filter(stage_id == 'FS')
data.fit.lin <- NULL
for (a in unique(tr.conds$arm)) {
  temp.cs <- dbGetQuery(conn, sprintf('select a.*, b.density, b.plate_no, b.plate_col, b.plate_row
                                         from LIN_OE_2202_FS_%s_6144_FITNESS a, LIN_OE_2202_pos2coor b
                                         where a.pos = b.pos
                                         order by a.hours, b.plate_no, b.plate_col, b.plate_row', a))
  temp.cs$condition <- a
  data.fit.lin <- rbind(data.fit.lin, temp.cs)
}
data.fit.lin <- data.fit.lin[!(data.fit.lin$condition == 'GA' & data.fit.lin$plate_no == 12),]
temp <- fit_sum(data.fit.lin, orf_types, 6)
data.fit.lin <- temp[[1]]
data.fit.lin.sum <- temp[[2]]

data.fit.lin <- merge(data.fit.lin, data.p2c.lin, by = 'rep')

dbWriteTable(conn, 'TR_OE_LIN_CS_FIT', data.fit.lin, overwrite = T)
dbWriteTable(conn, 'TR_OE_LIN_CS_FIT_SUMMARY', data.fit.lin.sum, overwrite = T)

##### COLONY SIZE TO AUC TO FITNESS
data.gc.fit.lin <- dbGetQuery(conn, 'select a.*, b.density, b.plate_no, b.plate_col, b.plate_row
                    from TR_OE_LIN_GC_6144_FITNESS a, LIN_OE_2202_pos2coor b
                    where a.pos = b.pos and orf_name != "NULL" and fitness is not NULL
                    order by a.hours, b.plate_no, b.plate_col, b.plate_row')

data.gc.fit.lin$condition[data.gc.fit.lin$hours == 1] <- 'GA'
data.gc.fit.lin$condition[data.gc.fit.lin$hours == 2] <- 'SA'
data.gc.fit.lin$condition[data.gc.fit.lin$hours == 3] <- 'CF'

temp <- fit_sum(data.gc.fit.lin, orf_types, 6)
data.gc.fit.lin <- temp[[1]]
data.gc.fit.lin.sum <- temp[[2]]

data.gc.fit.lin <- merge(data.gc.fit.lin, data.p2c.lin, by = 'rep')

dbWriteTable(conn, 'TR_OE_LIN_CS_GC_FIT', data.gc.fit.lin, overwrite = T)
dbWriteTable(conn, 'TR_OE_LIN_CS_GC_FIT_SUMMARY', data.gc.fit.lin.sum, overwrite = T)

##### COLONY SIZE TO GROWTH TO FITNESS
data.growth.lin <- dbGetQuery(conn, 'select a.*, b.density, b.plate_no, b.plate_col, b.plate_row
                                         from LIN_OE_2202_GROWTH_FS_ALL_6144_FITNESS a, LIN_OE_2202_pos2coor b
                                         where a.pos = b.pos
                                         order by a.hours, b.plate_no, b.plate_col, b.plate_row')
data.growth.lin$condition[data.growth.lin$hours == 1] <- 'CF'
data.growth.lin$condition[data.growth.lin$hours == 2] <- 'GA'
data.growth.lin$condition[data.growth.lin$hours == 3] <- 'SA'

temp <- fit_sum(data.growth.lin, orf_types, 6)
data.growth.lin <- temp[[1]]
data.growth.lin.sum <- temp[[2]]

data.growth.lin <- merge(data.growth.lin, data.p2c.lin, by = 'rep')

dbWriteTable(conn, 'TR_OE_LIN_GR_FIT', data.growth.lin, overwrite = T)
dbWriteTable(conn, 'TR_OE_LIN_GR_FIT_SUMMARY', data.growth.lin.sum, overwrite = T)

##################### GETTING DATA FROM MYSQL
data.fit.init <- dbGetQuery(conn, 'select * from TR_OE_INIT_CS_FIT')
data.fit.init.two.ga <- dbGetQuery(conn, 'select * from TR_OE_INIT_TWO_GA_CS_FIT')
data.fit.valid <- dbGetQuery(conn, 'select * from TR_OE_VALID_CS_FIT')
data.fit.lin <- dbGetQuery(conn, 'select * from TR_OE_LIN_CS_FIT')

data.gc.fit.valid <- dbGetQuery(conn, 'select * from TR_OE_VALID_CS_GC_FIT')
data.gc.fit.lin <- dbGetQuery(conn, 'select * from TR_OE_LIN_CS_GC_FIT')

data.growth.init <- dbGetQuery(conn, 'select * from TR_OE_INIT_GR_FIT')
data.growth.valid <- dbGetQuery(conn, 'select * from TR_OE_VALID_GR_FIT')
data.growth.lin <- dbGetQuery(conn, 'select * from TR_OE_LIN_GR_FIT')


data.fit.init.sum <- dbGetQuery(conn, 'select * from TR_OE_INIT_CS_FIT_SUMMARY')
data.fit.init.two.ga.sum <- dbGetQuery(conn, 'select * from TR_OE_INIT_TWO_GA_CS_FIT_SUMMARY')
data.fit.valid.sum <- dbGetQuery(conn, 'select * from TR_OE_VALID_CS_FIT_SUMMARY')
data.fit.lin.sum <- dbGetQuery(conn, 'select * from TR_OE_LIN_CS_FIT_SUMMARY')

data.gc.fit.valid.sum <- dbGetQuery(conn, 'select * from TR_OE_VALID_CS_GC_FIT_SUMMARY')
data.gc.fit.lin.sum <- dbGetQuery(conn, 'select * from TR_OE_LIN_CS_GC_FIT_SUMMARY')

data.growth.init.sum <- dbGetQuery(conn, 'select * from TR_OE_INIT_GR_FIT_SUMMARY')
data.growth.valid.sum <- dbGetQuery(conn, 'select * from TR_OE_VALID_GR_FIT_SUMMARY')
data.growth.lin.sum <- dbGetQuery(conn, 'select * from TR_OE_LIN_GR_FIT_SUMMARY')

##################### COMBINING ALL DATA ###########################
##### FITNESS DATA
data.fit.init.two.ga$condition <- 'GA2'

data.fit <- rbind(cbind(data.fit.init, attempt = 'init', data = 'cs'),
                  cbind(data.fit.init.two.ga, attempt = 'init', data = 'cs'),
                  cbind(data.fit.valid, attempt = 'valid', data = 'cs'),
                  cbind(data.fit.lin, attempt = 'lin', data = 'cs'),
                  # cbind(data.gc.fit.init, attempt = 'init', data = 'auc'),
                  cbind(data.gc.fit.valid, attempt = 'valid', data = 'auc'),
                  cbind(data.gc.fit.lin, attempt = 'lin', data = 'auc'),
                  cbind(data.growth.init, attempt = 'init', data = 'growth'),
                  cbind(data.growth.valid, attempt = 'valid', data = 'growth'),
                  cbind(data.growth.lin, attempt = 'lin', data = 'growth'))

data.fit$pos <- as.numeric(data.fit$pos)
data.fit$rep <- as.numeric(data.fit$rep)

data.fit <- arrange(data.fit,data,attempt,condition,pos,hours) %>% 
  mutate(avg_ratio=(lag(average,0)/lag(average,1)))
data.fit$avg_ratio[data.fit$data != 'cs'] <- NA
head(data.fit)

dbWriteTable(conn, 'TRANSLATOME_OE_FITNESS_DATA', data.fit)

##### FITNESS SUMMARY DATA
data.fit.init.two.ga.sum$condition <- 'GA2'
data.fit.sum <- rbind(cbind(data.fit.init.sum, attempt = 'init', data = 'cs'),
                      cbind(data.fit.init.two.ga.sum, attempt = 'init', data = 'cs'),
                      cbind(data.fit.valid.sum, attempt = 'valid', data = 'cs'),
                      cbind(data.fit.lin.sum, attempt = 'lin', data = 'cs'),
                      # cbind(data.gc.fit.init.sum, attempt = 'init', data = 'auc'),
                      cbind(data.gc.fit.valid.sum, attempt = 'valid', data = 'auc'),
                      cbind(data.gc.fit.lin.sum, attempt = 'lin', data = 'auc'),
                      cbind(data.growth.init.sum, attempt = 'init', data = 'growth'),
                      cbind(data.growth.valid.sum, attempt = 'valid', data = 'growth'),
                      cbind(data.growth.lin.sum, attempt = 'lin', data = 'growth'))
data.fit.sum$rep <- as.numeric(data.fit.sum$rep)

# data.fit.lim <- data.fit.sum %>%
#   filter(orf_name == 'BF_control') %>%
#   group_by(data, attempt, condition, hours) %>%
#   summarize(fit.ll = quantile(fit.median, 0.025, na.rm = T),
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

dbWriteTable(conn, 'TRANSLATOME_OE_FITNESS_SUMMARY_DATA', data.fit.sum)

