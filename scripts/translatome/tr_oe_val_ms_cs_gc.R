
source('/home/sbp29/R/Projects/adaptivefitness/scripts/translatome/initialize.R')
load('/home/sbp29/R/Projects/adaptivefitness/output/translatome/pinning_artifact_correction/VALID_CS_GROWTH_DATA_ALL.RData')

data.growth.curve <- merge(data.gc.res[,-5], data.gr.res[,-5], 
                           by = c('condition','pos','rep','strain_id'))
data.growth.curve <- merge(data.strains[,-6], data.growth.curve, 
                           by = c('condition','rep','strain_id'))
data.growth.curve$pos <- as.numeric(data.growth.curve$pos)
data.p2c.valid <- dbGetQuery(conn, 'select * from TR_OE_VAL_MS_pos2coor')
data.p2c.valid$pos <- as.numeric(data.p2c.valid$pos)

data.growth.curve <- merge(data.growth.curve, data.p2c.valid %>% filter(density == 6144), by = 'pos') 

dbWriteTable(conn, 'TR_OE_VALID_CS_GROWTH_CURVE', data.growth.curve, overwrite = T)

data.growth.curve$hours[data.growth.curve$condition == 'GA'] <- 1
data.growth.curve$hours[data.growth.curve$condition == 'SA'] <- 2
data.growth.curve$hours[data.growth.curve$condition == 'HU'] <- 3
data.growth.curve$hours[data.growth.curve$condition == 'DM'] <- 4
data.growth.curve$hours[data.growth.curve$condition == 'FL'] <- 5
data.growth.curve$hours[data.growth.curve$condition == 'TN'] <- 6

temp.mock.raw <- data.growth.curve[,c('pos','hours','auc_l')]
colnames(temp.mock.raw) <- c('pos','hours','average')

data.mock.raw <- NULL
for (h in unique(temp.mock.raw$hours)) {
  temp <- temp.mock.raw[temp.mock.raw$hours == h,]
  temp <- merge(temp, data.p2c.valid %>% filter(density == 6144), by = 'pos', all.y = T)
  temp$hours <- h
  data.mock.raw <- rbind(data.mock.raw, temp[,c('pos','hours','average')])
  }

dbWriteTable(conn, 'TR_OE_VAL_MS_GC_6144_CLEAN', data.mock.raw, overwrite = T)

data.growth.curve %>%
  ggplot(aes(x = pos, y = auc_e)) +
  geom_point(aes(col = category)) +
  facet_wrap(.~condition) +
  theme_minimal()

data.growth.curve %>%
  ggplot(aes(x = auc_l, y = auc_e)) +
  geom_point(aes(col = category)) +
  facet_wrap(.~condition) +
  theme_minimal()

data.growth.curve %>%
  group_by(condition) %>%
  count()
