
source('/home/sbp29/R/Projects/adaptivefitness/scripts/translatome/initialize.R')

tr.conds <- readxl::read_xlsx(path = '/home/sbp29/RAW_Data/TranslatomeOE_VAL/TR_OE_VAL_MS_INFO2.xlsx') %>% data.frame()

##### GATHER DATA
data.fit.all <- NULL
for (a in unique(tr.conds$arm)) {
  temp.cs <- dbGetQuery(conn, sprintf('select a.*, b.density, b.plate_no, b.plate_col, b.plate_row
                                         from TR_OE_VAL_MS_FS_%s_6144_FITNESS a, TR_OE_VAL_MS_pos2coor b
                                         where a.pos = b.pos
                                         order by a.hours, b.plate_no, b.plate_col, b.plate_row', a))
  temp.cs$condition <- a
  temp.cs$saturation <- max(temp.cs$hours)
  data.fit.all <- rbind(data.fit.all, temp.cs)
}
data.fit.all <- merge(data.fit.all, orf_types, by = 'orf_name', all.x = T)
data.fit.all$category[data.fit.all$orf_name == 'BF_control'] <- 'Reference'

data.fit.all$rep <- as.numeric(str_trunc(as.character(data.fit.all$pos), 5, side = 'left', ellipsis = ''))

##### REMOVE OUTLIERS
data.fit.sum <- data.fit.all %>%
  group_by(condition, hours, rep, strain_id, orf_name, category) %>%
  summarize(avg.median = median(average, na.rm = T), avg.mad = mad(average, na.rm = T),
            fitness.median = median(fitness, na.rm = T), fitness.mad = mad(fitness, na.rm = T),
            .groups = 'keep') %>%
  data.frame()

data.fit.all <- merge(data.fit.all, data.fit.sum,
                      by = c('condition','hours','rep','strain_id','orf_name','category'), all = T)

data.fit.all$fitness[data.fit.all$fitness < (data.fit.all$fitness.median - 2*data.fit.all$fitness.mad) |
                       data.fit.all$fitness > (data.fit.all$fitness.median + 2*data.fit.all$fitness.mad)] <- NA
data.fit.all$average[is.na(data.fit.all$fitness)] <- NA

# ## REPLICATE COUNT PER REP
# data.fit.cnt <- data.fit.all[!is.na(data.fit.all$fitness),] %>%
#   group_by(condition, hours, rep) %>% 
#   count() %>%
#   data.frame()
# data.fit.all <- merge(data.fit.all, data.fit.cnt, by = c('condition','hours','rep'), all.x = T)
# data.fit.all$n[is.na(data.fit.all$n)] <- 0
# 
# ## REMOVE DATA WITH LESS THAN 4 REPLICATES
# data.fit.all$fitness[data.fit.all$n < 4] <- NA

data.fit.sum <- data.fit.all %>%
  group_by(condition, hours, rep, strain_id, orf_name, category) %>%
  summarize(avg.median = median(average, na.rm = T),  avg.std = sd(average, na.rm = T),
            fitness.median = median(fitness, na.rm = T),
            .groups = 'keep') %>%
  data.frame()
# head(data.fit.sum)

data.fit.all <- merge(data.fit.all[,c(1:20)], data.fit.sum,
                      by = c('condition','hours','rep','strain_id','orf_name','category'), all = T)
# head(data.fit.all)
# 
# write.csv(data.fit.all[,c(1:14)], 
#           file = '/home/sbp29/R/Projects/adaptivefitness/output/translatome/pinning_artifact_correction/tr_oe_val_ms_fitandcs_all.csv')
# 
# ##### PLOT FITNESS CURVE
# data.fit.all %>%
#   ggplot(aes(x = hours, y = norm_cs)) +
#   geom_line(aes(group = rep)) +
#   facet_wrap(.~condition) +
#   coord_cartesian(ylim = c(0,2000))
# 
# # data.fit.sum %>%
# #   filter(fitness.median >= 0) %>%
# #   group_by(category, strain_id) %>%
# #   count() %>%
# #   group_by(category) %>%
# #   count()
# 
# data.fit.all %>%
#   group_by(condition, hours, rep) %>%
#   summarize(cs = median(norm_cs, na.rm = T)) %>%
#   group_by(condition, hours) %>%
#   summarize(max_cs = max(cs, na.rm = T))
# 
# data.fit.all %>%
#   group_by(condition, hours, rep) %>%
#   summarize(cs = median(norm_cs, na.rm = T)) %>%
#   group_by(condition, hours) %>%
#   summarize(max_cs = max(cs, na.rm = T)) %>%
#   filter(max_cs <= 1800) %>%
#   group_by(condition) %>%
#   summarize(sat_hr = max(hours))
