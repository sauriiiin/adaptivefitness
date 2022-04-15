source('scripts/translatome/initialize.R')

data.growth <- read.table(file = '/home/acwach/overexp/colonies_del_init', sep = ' ', header = T)

p2s <- dbGetQuery(conn, 'select a.*, b.strain_id
                from TR_DEL_pos2coor a, TR_DEL_pos2strainid b
                where a.pos = b.pos')
p2s$pos <- as.numeric(p2s$pos)

data.growth %>%
  group_by(condition, arm) %>%
  count() %>% data.frame()

data.growth %>%
  filter(condition == 'HO', arm == 'R2')

temp.growth <- NULL
for (c in unique(data.growth$condition)) {
  temp <- data.growth[data.growth$condition == c,]
  temp <- merge(p2s[p2s$density == 1536,], temp[,c('condition','position','growth')], by.x = 'pos', by.y = 'position', all.x = T)
  temp$condition[is.na(temp$condition)] <- c
  temp.growth <- rbind(temp, temp.growth)
}
data.growth <- temp.growth
unique(data.growth$condition)

data.growth$hours[data.growth$condition == 'DM'] <- 1
data.growth$hours[data.growth$condition == 'FL'] <- 2
data.growth$hours[data.growth$condition == 'HO'] <- 3
data.growth$hours[data.growth$condition == 'HU'] <- 4
data.growth$hours[data.growth$condition == 'SA'] <- 5
data.growth$hours[data.growth$condition == 'TN'] <- 6
data.growth$hours[data.growth$condition == 'YP'] <- 7

data.growth <- data.growth[,c('pos','hours','growth')]
colnames(data.growth) <- c('pos','hours','average')

dbWriteTable(conn, 'TR_DEL_VAL_MS_GROWTH_FS_ALL_1536_CLEAN', data.growth, overwrite = T)

