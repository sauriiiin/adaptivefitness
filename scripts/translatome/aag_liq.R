##### AAG-ATG MUTANT LIQUID SCREENS
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 04/17/2022


##### INITIALIZE
source('/home/sbp29/R/Projects/adaptivefitness/scripts/translatome/initialize.R')
platemap <- read_xlsx('/home/sbp29/R/Projects/adaptivefitness/rawdata/translatome/validation/AAG_LIQ_VAL_PLATEMAP.xlsx') %>% data.frame()
platemap$smudge[platemap$pos_id == 'A12'] <- 'Y'

##### GATHER DATA
data.aag.liq <- read_xlsx('/home/sbp29/R/Projects/adaptivefitness/rawdata/translatome/validation/AAG_LIQ_VAL_VALUES.xlsx') %>% data.frame()
data.aag.liq2 <- merge(platemap, 
                       melt(data.aag.liq, id.vars = 'Time', variable.name = 'pos_id', value.name = 'OD'),
                       by = 'pos_id', all = T)
data.aag.liq2$category[data.aag.liq2$strain_id > 100000 & data.aag.liq2$strain_id < 200000] <- 'KanMX'
data.aag.liq2$category[data.aag.liq2$strain_id > 200000 & data.aag.liq2$strain_id < 300000] <- 'AAG'
data.aag.liq2$category[data.aag.liq2$strain_id > 300000] <- 'WT'


##### GROWTH ANALYSIS
data.aag.gc <- SummarizeGrowthByPlate(data.aag.liq)
data.aag.gr <- NULL
for (i in 2:dim(data.aag.liq)[2]) {
  if (sum(data.aag.liq[,i] < 0.005) < 5) {
    fit0 <- fit_easylinear(data.aag.liq$Time, data.aag.liq[,i], h = 8, quota = 1);
    
    temp_res <- data.frame(sample = colnames(data.aag.liq[i]), maxgr = coef(fit0)[[3]],
                           dtime = log(2)/coef(fit0)[[3]], ltime = coef(fit0)[[4]])
    data.aag.gr <- rbind(data.aag.gr, temp_res)
  }
}

data.aag.liq.fit <- merge(platemap, merge(data.aag.gc, data.aag.gr, by = 'sample', all = T), by.x = 'pos_id', by.y = 'sample', all = T)
data.aag.liq.fit$category[data.aag.liq.fit$strain_id > 100000 & data.aag.liq.fit$strain_id < 200000] <- 'KanMX'
data.aag.liq.fit$category[data.aag.liq.fit$strain_id > 200000 & data.aag.liq.fit$strain_id < 300000] <- 'AAG'
data.aag.liq.fit$category[data.aag.liq.fit$strain_id > 300000] <- 'WT'

##### PLOT RESULTS
data.aag.liq2 %>%
  filter(is.na(smudge), strain_id > 0) %>%
  ggplot(aes(x = Time, y = OD)) +
  stat_summary(aes(col = category), geom = 'line', fun = 'median') +
  stat_summary(aes(fill = category), 
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.4) +
  facet_grid(condition ~ orf_name)

data.aag.liq.fit %>%
  filter(is.na(smudge), strain_id > 0) %>%
  ggplot(aes(x = category, y = auc_l)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  stat_compare_means(method = 't.test', size = 3) +
  facet_grid(condition ~ orf_name) +
  coord_cartesian(ylim = c(0,10000))


##### SAVING DATA
dbWriteTable(conn, 'TR_AAG_LIQ_OD', data.aag.liq2, overwrite = T)
dbWriteTable(conn, 'TR_AAG_LIQ_RESULTS', data.aag.liq.fit, overwrite = T)


  