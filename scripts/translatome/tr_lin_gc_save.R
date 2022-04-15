##### OVEREXPRESSION LIN SCREEN GROWTH CURVE ANALYSIS
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 03/25/2022 

source('/home/sbp29/R/Projects/adaptivefitness/scripts/translatome/initialize.R')

data.fit.all <- dbGetQuery(conn, 'select * from TR_OE_LIN_CS_FIT')
data.fit.sum <- dbGetQuery(conn, 'select * from TR_OE_LIN_CS_FIT_SUMMARY')

##### REFERENCE LIMITS
data.fit.lim <- data.fit.all %>%
  filter(orf_name == 'BF_control') %>%
  group_by(condition, hours, orf_name, rep) %>%
  summarize(average = median(average, na.rm = T),
            fitness = median(fitness, na.rm = T),
            .groups = 'keep') %>%
  group_by(condition, hours, orf_name) %>%
  summarize(avg_ll = quantile(average, 0.025, na.rm = T),
            avg_m = median(average, na.rm = T),
            avg_ul = quantile(average, 0.975, na.rm = T),
            fitness_ll = quantile(fitness, 0.025, na.rm = T),
            fitness_m = median(fitness, na.rm = T),
            fitness_ul = quantile(fitness, 0.975, na.rm = T),
            .groups = 'keep') %>%
  data.frame()

data.fit.sum <- merge(data.fit.sum, data.fit.lim[,-3], by = c('condition','hours'))

##### GROWTH CURVE ANALYSIS
data.fit.all$source[data.fit.all$plate_row%%2==1 & data.fit.all$plate_col%%2==1] = 'TL'
data.fit.all$source[data.fit.all$plate_row%%2==0 & data.fit.all$plate_col%%2==1] = 'BL'
data.fit.all$source[data.fit.all$plate_row%%2==1 & data.fit.all$plate_col%%2==0] = 'TR'
data.fit.all$source[data.fit.all$plate_row%%2==0 & data.fit.all$plate_col%%2==0] = 'BR'

data.fit.lim.src <- data.fit.all %>%
  filter(orf_name == 'BF_control') %>%
  group_by(condition, hours, source) %>%
  summarize(avg_ll = quantile(average, 0.025, na.rm = T),
            avg_m = median(average, na.rm = T),
            avg_ul = quantile(average, 0.975, na.rm = T),
            fitness_ll = quantile(fitness, 0.025, na.rm = T),
            fitness_m = median(fitness, na.rm = T),
            fitness_ul = quantile(fitness, 0.975, na.rm = T),
            .groups = 'keep') %>%
  data.frame()

data.fit.all <- merge(data.fit.all, data.fit.lim.src, by = c('condition','hours','source'))
data.fit.all$norm_cs <- data.fit.all$fitness * data.fit.all$avg_m
data.fit.all$pos <- as.numeric(data.fit.all$pos)

data.strains <- data.fit.all %>%
  filter(orf_name != 'NULL') %>%
  group_by(condition, rep, strain_id, orf_name, category) %>%
  count() %>% data.frame()

data.gr.res <- NULL
data.gc.res <- NULL
for (c in unique(data.fit.all$condition)) {
  sprintf('Now analyzing %s.', c)
  data.pred <- NULL
  col.names <- NULL
  data.temp <- data.fit.all[data.fit.all$condition == c,]
  for (p in unique(data.temp$pos)) {
    r <- unique(data.temp$rep[data.temp$pos == p])
    s <- unique(data.temp$strain_id[data.temp$pos == p])
    temp <- data.temp[data.temp$pos == p,]
    if (sum(!is.na(temp$average)) >= length(temp$average) * 0.7) {
      lo <- loess.smooth(temp$hours, log(temp$average),
                         span = 0.6, evaluation = 100, degree = 2,
                         family = 'gaussian')
      data.pred <- cbind(data.pred,exp(lo$y))
      # data.pred <- cbind(data.pred, temp$average)
      col.names <- cbind(col.names,paste(c,p,r,s,sep = ','))
    }
  }
  data.pred <- cbind(lo$x, data.pred)
  # data.pred <- cbind(temp$hours, data.pred)
  data.pred <- data.frame(data.pred)
  colnames(data.pred) <- c('Time',col.names)
  # save(data.pred, file = sprintf('/home/sbp29/R/Projects/adaptivefitness/output/translatome/pinning_artifact_correction/growth_curve_prediction/data.lin.pred.pos.%s.RData',c))
  
  ## GROWTH CURVE ANALYSIS
  temp.gr.res <- NULL
  for (i in 2:dim(data.pred)[2]) {
    fit0 <- fit_easylinear(data.pred$Time[2:dim(data.pred)[1]], data.pred[2:dim(data.pred)[1],i], h = 8, quota = 1);
    
    temp_res <- data.frame(colnames(data.pred[i]), maxgr = coef(fit0)[[3]],
                           dtime = log(2)/coef(fit0)[[3]], ltime = coef(fit0)[[4]])
    temp.gr.res <- rbind(temp.gr.res, temp_res)
  }
  temp.gr.res <- data.frame(temp.gr.res)
  colnames(temp.gr.res) <- c('sample','gr','dtime','lag')
  temp <- str_split(temp.gr.res$sample, ',', simplify = T)
  colnames(temp) <- c('condition','pos','rep','strain_id')
  temp.gr.res <- cbind(temp, temp.gr.res)
  
  temp.gc.res <- SummarizeGrowthByPlate(data.pred)
  temp <- str_split(temp.gc.res$sample, ',', simplify = T)
  colnames(temp) <- c('condition','pos','rep','strain_id')
  temp.gc.res <- cbind(temp, temp.gc.res)
  
  data.gr.res <- rbind(data.gr.res, temp.gr.res)
  data.gc.res <- rbind(data.gc.res, temp.gc.res)
}
# data.gcgr.res <- merge(data.gc.res, data.gr.res, by = c('arm','condition','pos','rep','strain_id','sample'))
# data.gcgr.res <- merge(data.strains[,-7], data.gcgr.res, by = c('arm','condition','pos','rep','strain_id'))

save.image(file = '/home/sbp29/R/Projects/adaptivefitness/output/translatome/pinning_artifact_correction/LIN_CS_GROWTH_DATA_ALL.RData')




