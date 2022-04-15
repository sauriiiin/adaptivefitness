
expt.file <- read_xlsx(path = '/home/sbp29/RAW_Data/TranslatomeDEL_VAL/TR_DEL_VAL_MS_INFO.xlsx') %>% data.frame()
head(expt.file)

data.del.val <- NULL
for (e in unique(expt.file$expt_id)) {
  # for (s in unique(expt.file$stage_id[expt.file$expt_id == e])) {
  s <- 'FS'
  for (a in unique(expt.file$arm[expt.file$expt_id == e & expt.file$stage_id == s])) {
    for (d in unique(expt.file$density[expt.file$expt_id == e & expt.file$stage_id == s & expt.file$arm == a])) {
      temp.fit <- dbGetQuery(conn, sprintf('select a.*, b.density, b.plate_no, b.plate_col, b.plate_row
                                         from %s_%s_%s_%d_FITNESS a, %s_pos2coor b
                                         where a.pos = b.pos
                                         order by a.hours, b.plate_no, b.plate_col, b.plate_row', e, s, a, d, e))
      temp.fit$expt_id <- e
      temp.fit$stage_id <- s
      temp.fit$condition <- a
      
      data.del.val <- rbind(data.del.val, temp.fit)
    }
  }
  # }
}
data.del.val <- data.del.val[!(data.del.val$condition == 'YPDA' &
                               data.del.val$plate_no == 2),]
data.del.val <- data.del.val[!(data.del.val$condition == 'FL' &
                                 data.del.val$density %in% c(1536, 6144) &
                                 data.del.val$plate_no == 2),]
data.del.val <- data.del.val[!(data.del.val$condition == 'SA' &
                                 data.del.val$density %in% c(1536) &
                                 data.del.val$plate_no == 2),]
data.del.val <- data.del.val[!(data.del.val$condition == 'DM' &
                                 data.del.val$density %in% c(6144) &
                                 data.del.val$plate_no == 1),]
data.del.val <- data.del.val[!(data.del.val$condition == 'TN' &
                                 data.del.val$density %in% c(1536, 6144) &
                                 data.del.val$plate_no == 1),]
data.del.val <- data.del.val[!(data.del.val$condition == 'HU' &
                                 data.del.val$density %in% c(1536, 6144) &
                                 data.del.val$plate_no == 1),]

data.del.val.sum <- data.del.val %>%
  filter(orf_name %notin% c('BOR','NULL')) %>%
  group_by(expt_id, stage_id, condition, density, hours, strain_id, orf_name) %>%
  summarise(avg.median = median(average, na.rm = T), avg.mad = mad(average, na.rm = T),
            fitness.median = median(fitness, na.rm = T), fitness.mad = mad(fitness, na.rm = T),
            .groups = 'keep') %>%
  data.frame()

##### REMOVE OUTLIERS
data.del.val <- merge(data.del.val, data.del.val.sum,
                            by = c('expt_id','stage_id','condition','density','hours','strain_id','orf_name'), all = T)

data.del.val$average[data.del.val$average < (data.del.val$avg.median - 2*data.del.val$avg.mad) |
                             data.del.val$average > (data.del.val$avg.median + 2*data.del.val$avg.mad)] <- NA
data.del.val$fitness[data.del.val$fitness < (data.del.val$fitness.median - 2*data.del.val$fitness.mad) |
                             data.del.val$fitness > (data.del.val$fitness.median + 2*data.del.val$fitness.mad)] <- NA

data.del.val.sum <- data.del.val %>%
  filter(orf_name %notin% c('BOR','NULL')) %>%
  group_by(expt_id, stage_id, condition, density, hours, strain_id, orf_name) %>%
  summarise(avg.median = median(average, na.rm = T), avg.mad = mad(average, na.rm = T),
            fitness.median = median(fitness, na.rm = T), fitness.mad = mad(fitness, na.rm = T),
            .groups = 'keep') %>%
  data.frame()
head(data.del.val.sum)

unique(data.del.val.sum$hours)
##### GROWTH CURVES
data.del.val %>%
  filter(orf_name %notin% c('BOR','NULL')) %>%
  ggplot(aes(x = hours, y = average)) +
  stat_summary(aes(fill = orf_name, group = orf_name),
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.4) +
  stat_summary(aes(col = orf_name, group = orf_name),
               fun=mean, geom="line", lwd =1) +
  scale_color_discrete(guide = 'none') +
  scale_fill_discrete(guide = 'none') +
  facet_grid(density~condition,
             scales = 'free_y') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) 

data.del.val.sum %>%
  filter(orf_name %notin% c('BOR','NULL')) %>%
  ggplot(aes(x = hours, y = fitness.mad)) +
  stat_summary(aes(fill = orf_name, group = orf_name),
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.4) +
  stat_summary(aes(col = orf_name, group = orf_name),
               fun=mean, geom="line", lwd =1) +
  scale_color_discrete(guide = 'none') +
  scale_fill_discrete(guide = 'none') +
  facet_grid(density~condition,
             scales = 'free_y') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) 

unique(data.del.val$hours)
data.del.val$saturation[data.del.val$density == 384] <- 60
data.del.val$saturation[data.del.val$density == 1536] <- 60 
data.del.val$saturation[data.del.val$density == 6144 &
                                data.del.val$condition %notin% c('FL','SA')] <- 6
data.del.val$saturation[data.del.val$density == 6144 &
                                data.del.val$condition %in% c('FL')] <- 12
data.del.val$saturation[data.del.val$density == 6144 &
                                data.del.val$condition %in% c('SA')] <- 60

data.del.val.sum$saturation[data.del.val.sum$density == 384] <- 60
data.del.val.sum$saturation[data.del.val.sum$density == 1536] <- 60 
data.del.val.sum$saturation[data.del.val.sum$density == 6144 &
                                data.del.val.sum$condition %notin% c('FL','SA')] <- 6
data.del.val.sum$saturation[data.del.val.sum$density == 6144 &
                                data.del.val.sum$condition %in% c('FL')] <- 12
data.del.val.sum$saturation[data.del.val.sum$density == 6144 &
                                data.del.val.sum$condition %in% c('SA')] <- 60

##### FITNESS PLOTS
plot.fit.sum <- data.del.val %>%
  filter(hours == saturation,
         orf_name %notin% c('REF','BOR','NULL')) %>%
  ggplot(aes(x = orf_name, y = fitness)) +
  # geom_jitter(aes(col = as.factor(strain_id)), size = 0.3) +
  geom_hline(yintercept = c(0.99,1.01), col = 'red', linetype = 'dashed', size = 0.2) +
  geom_boxplot(aes(fill = orf_name), outlier.shape = NA, size = 0.3) +
  scale_color_discrete(guide = 'none') +
  scale_fill_discrete(guide = 'none') +
  facet_grid(density~condition) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0.5,1.5))
ggsave(sprintf("%s/TR_DEL_VAL_SUMMARY.jpg",fig_path), plot.fit.sum,
       height = two.c, width = two.c*3, units = 'mm',
       bg = 'white',
       dpi = 600)


###### SAVING DATA
data.del.val$rep <- as.numeric(str_trunc(as.character(data.del.val$pos), 5, side = 'left', ellipsis = ''))
write.csv(data.del.val[data.del.val$orf_name %notin% c('NULL','BOR'),c(3:14,20)],
          file = '/home/sbp29/R/Projects/adaptivefitness/output/translatome/tr_del_val_fitandcs_all.csv')

