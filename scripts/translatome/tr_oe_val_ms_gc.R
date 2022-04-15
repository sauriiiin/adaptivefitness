
source('/home/sbp29/R/Projects/adaptivefitness/scripts/translatome/initialize.R')
load('/home/sbp29/R/Projects/adaptivefitness/output/translatome/pinning_artifact_correction/LIN_GROWTH_DATA_ALL.RData')

data.growth.curve <- merge(data.gc.res[,-5], data.gr.res[,-5], 
                           by = c('condition','pos','rep','strain_id'))
data.growth.curve <- merge(data.strains[,-6], data.growth.curve, 
                           by = c('condition','rep','strain_id'))
data.growth.curve$pos <- as.numeric(data.growth.curve$pos)
data.p2c.valid <- dbGetQuery(conn, 'select * from LIN_OE_2202_pos2coor')
data.p2c.valid$pos <- as.numeric(data.p2c.valid$pos)

data.growth.curve <- merge(data.growth.curve, data.p2c.valid %>% filter(density == 6144), by = 'pos') 

dbWriteTable(conn, 'TR_OE_LIN_GROWTH_CURVE', data.growth.curve, overwrite = T)

data.growth.curve %>%
  ggplot(aes(x = pos, y = auc_l)) +
  geom_point(aes(col = category)) +
  facet_wrap(.~condition) +
  theme_minimal()

data.growth.curve %>%
  group_by(condition) %>%
  count()
