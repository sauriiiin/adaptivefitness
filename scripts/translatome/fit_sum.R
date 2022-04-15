

fit_sum <- function(data.fit.all, orf_types, n){
  
  ##### ADD CLASSIFIERS
  data.fit.all <- merge(data.fit.all, orf_types[,c(1,7)], by = 'orf_name', all.x = T)
  data.fit.all$category[data.fit.all$orf_name == 'BF_control'] <- 'Reference'
  data.fit.all$category[data.fit.all$strain_id >= 1000000] <- 'ConditionControls'
  
  data.fit.all$rep <- as.numeric(str_trunc(as.character(data.fit.all$pos), n, side = 'left', ellipsis = ''))
  
  data.fit.all$source[data.fit.all$plate_row%%2==1 & data.fit.all$plate_col%%2==1] = 'TL'
  data.fit.all$source[data.fit.all$plate_row%%2==0 & data.fit.all$plate_col%%2==1] = 'BL'
  data.fit.all$source[data.fit.all$plate_row%%2==1 & data.fit.all$plate_col%%2==0] = 'TR'
  data.fit.all$source[data.fit.all$plate_row%%2==0 & data.fit.all$plate_col%%2==0] = 'BR'
  
  ##### REMOVE OUTLIERS
  data.fit.sum.temp <- data.fit.all %>%
    group_by(condition, hours, rep, strain_id, orf_name, category) %>%
    summarise(avg.median = median(average, na.rm = T), avg.mad = mad(average, na.rm = T),
              fitness.median = median(fitness, na.rm = T), fitness.mad = mad(fitness, na.rm = T),
              .groups = 'keep') %>%
    data.frame()
  
  data.fit.all <- merge(data.fit.all, data.fit.sum.temp,
                        by = c('condition','hours','rep','strain_id','orf_name','category'), all = T)
  data.fit.all$fitness[data.fit.all$fitness < (data.fit.all$fitness.median - 2*data.fit.all$fitness.mad) |
                         data.fit.all$fitness > (data.fit.all$fitness.median + 2*data.fit.all$fitness.mad)] <- NA
  data.fit.all$average[is.na(data.fit.all$fitness)] <- NA
  
  data.fit.all <- data.fit.all[,c(1:15)]
  
  ## REPLICATE COUNT PER REP
  data.fit.cnt <- data.fit.all[!is.na(data.fit.all$fitness),] %>%
    group_by(condition, hours, rep) %>%
    count() %>%
    data.frame()
  data.fit.all <- merge(data.fit.all, data.fit.cnt, by = c('condition','hours','rep'), all.x = T)
  data.fit.all$n[is.na(data.fit.all$n)] <- 0
  data.fit.all$fitness[data.fit.all$n < 4] <- NA
  
  ##### FITNESS SUMMARY
  data.fit.sum <- data.fit.all %>%
    group_by(condition, hours, rep, strain_id, orf_name, category) %>%
    summarise(avg.median = median(average, na.rm = T), avg.std = sd(average, na.rm = T),
              fit.median = median(fitness, na.rm = T), fit.std = sd(fitness, na.rm = T),
              .groups = 'keep') %>%
    data.frame()
  
  ##### RETURN DATA
  output <- list(data.fit.all, data.fit.sum)
  return(output)
}

