tr.conds <- data.frame(arms = c('ONE','ONE','ONE','ONE','TWO','TWO','TWO'),
                       conds = c('GA','SA','HO','HU','DM','FL','TN'))

##### GATHER DATA
data.growth.all <- dbGetQuery(conn, 'select a.*, b.density, b.plate_no, b.plate_col, b.plate_row
                                         from TR_OE_GROWTH_FS_ALL_6144_FITNESS a, TRANS_OE_pos2coor b
                                         where a.pos = b.pos
                                         order by a.hours, b.plate_no, b.plate_col, b.plate_row')
data.growth.all$condition[data.growth.all$hours == 1] <- 'DM'
data.growth.all$condition[data.growth.all$hours == 2] <- 'FL'
data.growth.all$condition[data.growth.all$hours == 3] <- 'GA'
data.growth.all$condition[data.growth.all$hours == 4] <- 'HU'
data.growth.all$condition[data.growth.all$hours == 5] <- 'SA'
data.growth.all$condition[data.growth.all$hours == 6] <- 'TN'
data.growth.all$condition[data.growth.all$hours == 7] <- 'HO'

data.growth.all$average[data.growth.all$condition == 'HO' & data.growth.all$plate_no %in% c(2,3,10,18,19,23,26,32)] <- NA
data.growth.all$average[data.growth.all$condition %in% c('GA','HU','HO') & data.growth.all$plate_no == 14] <- NA
data.growth.all$average[data.growth.all$condition %in% c('GA','HU','HO') & data.growth.all$pos %in% smudge$pos] <- NA

data.growth.all$fitness[data.growth.all$condition == 'HO' & data.growth.all$plate_no %in% c(2,3,10,18,19,23,26,32)] <- NA
data.growth.all$fitness[data.growth.all$condition %in% c('GA','HU','HO') & data.growth.all$plate_no == 14] <- NA
data.growth.all$fitness[data.growth.all$condition %in% c('GA','HU','HO') & data.growth.all$pos %in% smudge$pos] <- NA

data.growth.all <- merge(data.growth.all, orf_types, by = 'orf_name', all.x = T)
data.growth.all$category[data.growth.all$orf_name == 'BF_control'] <- 'Reference'
data.growth.all$category[data.growth.all$strain_id >= 1000000] <- 'ConditionControls'

data.growth.all$rep <- as.numeric(str_trunc(as.character(data.growth.all$pos), 5, side = 'left', ellipsis = ''))

data.growth.sum <- data.growth.all %>%
  group_by(condition, hours, rep, strain_id, orf_name, category) %>%
  summarize(avg.median = median(average, na.rm = T), avg.mad = mad(average, na.rm = T),
            fitness.median = median(fitness, na.rm = T), fitness.mad = mad(fitness, na.rm = T),
            .groups = 'keep') %>%
  data.frame()

##### REMOVE OUTLIERS
data.growth.all <- merge(data.growth.all, data.growth.sum,
                      by = c('condition','hours','rep','strain_id','orf_name','category'), all = T)

data.growth.all$fitness[data.growth.all$fitness < (data.growth.all$fitness.median - 2*data.growth.all$fitness.mad) |
                       data.growth.all$fitness > (data.growth.all$fitness.median + 2*data.growth.all$fitness.mad)] <- NA
data.growth.all$average[is.na(data.growth.all$fitness)] <- NA

##### REPLICATE COUNT PER REP
data.growth.cnt <- data.growth.all[!is.na(data.growth.all$fitness),] %>%
  group_by(condition, rep) %>% 
  count() %>%
  data.frame()
data.growth.all <- merge(data.growth.all, data.growth.cnt, by = c('condition','rep'), all.x = T)
data.growth.all$n[is.na(data.growth.all$n)] <- 0

##### REMOVE DATA WITH LESS THAN 4 REPLICATES
data.growth.all$fitness[data.growth.all$n < 4] <- NA

##### FITNESS SUMMARY
data.growth.sum <- data.growth.all %>%
  group_by(condition, hours, rep, strain_id, orf_name, n, category) %>%
  summarize(avg.median = median(average, na.rm = T),
            fitness.median = median(fitness, na.rm = T),
            .groups = 'keep') %>%
  data.frame()
head(data.growth.sum)

data.growth.all <- merge(data.growth.all[,c(1:19,24)], data.growth.sum,
                      by = c('condition','hours','rep','strain_id','orf_name','n','category'), all = T)

write.csv(data.growth.all[,c(1:14)], file = '/home/sbp29/R/Projects/adaptivefitness/output/translatome/pinning_artifact_correction/tr_oe_fitandgrowth_all.csv')

#####
data.growth.all$pos <- as.numeric(data.growth.all$pos)
data.growth.all %>%
  ggplot(aes(x = as.factor(pos), y = fitness, col = category)) +
  geom_point() +
  facet_wrap(.~condition) +
  theme(axis.text.x = element_blank()) 
