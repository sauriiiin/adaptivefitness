
##### 
`%notin%` <- Negate(`%in%`)

#####
# data.fit.all is all the fitness data after outlier removal

data.fit.cs <- merge(data.fit.all[,c(1:21)], data.fit.lim[,-4], by = c('arm','condition','hours'))

data.temp <- data.fit.cs %>% 
  # filter(arm == 'TWO', condition == 'FL') %>%
  group_by(arm, condition, hours, rep, strain_id, orf_name, category) %>%
  summarize(fitness.median = median(fitness, na.rm = T),
            avg_m = median(avg_m, na.rm = T),
            .groups = 'keep') %>%
  data.frame() 
data.temp$norm_cs <- data.temp$fitness.median * data.temp$avg_m

data.temp <- data.temp %>%
  group_by(arm, condition, rep, strain_id, orf_name, category, hours) %>%
  summarise(norm_cs = norm_cs, .groups = 'keep') %>% data.frame()

data.strains <- data.temp %>%
  group_by(arm, condition, rep, strain_id, orf_name, category) %>%
  count() %>% data.frame()

data.strains %>%
  group_by(arm, condition, category) %>%
  count() %>% data.frame()

data.gr.res <- NULL
data.gc.res <- NULL
for (a in unique(data.temp$arm)) {
  for (c in unique(data.temp$condition[data.temp$arm == a])) {
    data.pred <- NULL
    col.names <- NULL
    for (r in unique(data.temp$rep[data.temp$arm == a & data.temp$condition == c])) {
      for (s in unique(data.temp$strain_id[data.temp$arm == a & data.temp$condition == c & data.temp$rep == r])) {
        temp <- data.temp[data.temp$arm == a &
                            data.temp$condition == c &
                            data.temp$rep == r &
                            data.temp$strain_id == s,]
        if (sum(is.na(temp$norm_cs)) <= 5) {
          lo <- loess.smooth(temp$hours, log(temp$norm_cs),
                             span = 0.6, evaluation = 100, degree = 2,
                             family = 'gaussian')
          data.pred <- cbind(data.pred,exp(lo$y))
          # data.pred <- cbind(data.pred, temp$norm_cs)
          col.names <- cbind(col.names,paste(a,c,r,s,sep = ','))
          
          # temp.plot <- ggplot() +
          #   geom_line(data = data.frame(lo), aes(x = x, y = y)) +
          #   geom_point(data = temp, aes(x = hours, y = log(norm_cs)), shape = 1) +
          #   labs(y = 'Log( Colony Size (pixels) )',
          #        x = 'Time (hours)',
          #        title = paste(a,c,r,s,sep = ' | ')) +
          #   theme_linedraw() +
          #   theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
          #         axis.title = element_text(size = titles),
          #         axis.text = element_text(size = txt),
          #         legend.title = element_text(size = titles),
          #         legend.text = element_text(size = txt),
          #         legend.position = 'bottom',
          #         legend.key.size = unit(3, "mm"),
          #         legend.box.spacing = unit(0.5,"mm"),
          #         strip.text = element_text(size = txt,
          #                                   face = 'bold',
          #                                   margin = margin(0.1,0,0.1,0, "mm"))) +
          #   coord_cartesian(ylim = c(3,9))
          # ggsave(sprintf("%s/%s/growthcurves/%s_%s_%s.jpg",fig_path, expt.name, o, c, b), temp.plot,
          #        height = 100, width = 100, units = 'mm',
          #        dpi = 600)
        }
      }
    }
    data.pred <- cbind(lo$x, data.pred)
    # data.pred <- cbind(temp$hours, data.pred)
    data.pred <- data.frame(data.pred)
    colnames(data.pred) <- c('Time',col.names)
    head(data.pred)
    
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
    colnames(temp) <- c('arm','condition','rep','strain_id')
    temp.gr.res <- cbind(temp, temp.gr.res)
    
    temp.gc.res <- SummarizeGrowthByPlate(data.pred)
    temp <- str_split(temp.gc.res$sample, ',', simplify = T)
    colnames(temp) <- c('arm','condition','rep','strain_id')
    temp.gc.res <- cbind(temp, temp.gc.res)
    
    data.gr.res <- rbind(data.gr.res, temp.gr.res)
    data.gc.res <- rbind(data.gc.res, temp.gc.res)
  }
}
data.gcgr.res <- merge(data.gc.res, data.gr.res, by = c('arm','condition','rep','strain_id','sample'))
data.gcgr.res <- merge(data.strains[,-7], data.gcgr.res, by = c('arm','condition','rep','strain_id'))


##### REFERENCE LIMITS
head(data.gcgr.res)
data.gc.lim <- data.gcgr.res[,c(-3,-4,-6,-7,-16)] %>%
  filter(orf_name == 'BF_control') %>%
  # group_by(arm, condition, orf_name) %>%
  melt(id.vars = c('arm','condition','orf_name')) %>%
  group_by(arm, condition, orf_name, variable) %>%
  summarize(ll = quantile(value, 0.005, na.rm = T),
            m = median(value, na.rm = T),
            ul = quantile(value, 0.995, na.rm = T),
            .groups = 'keep') %>%
  data.frame()

data.gcgr.sum <- melt(data.gcgr.res[!is.na(data.gcgr.res$category),c(1:14,17:19)] %>%
       filter(category != 'Reference'), 
       id.vars = c('arm','condition','rep','strain_id','orf_name','category','sample'))
data.gcgr.sum <- merge(data.gcgr.sum, data.gc.lim[,-3], by = c('arm','condition','variable'), all.x = T)

data.gcgr.sum$es <- (data.gcgr.sum$value - data.gcgr.sum$m)/data.gcgr.sum$m * 100

data.gcgr.sum$annotation[str_detect(data.gcgr.sum$orf_name, 'smo')] <- 'Unannotated'
data.gcgr.sum$annotation[is.na(data.gcgr.sum$annotation)] <- 'Annotated'


##### FALSE DISCOVERY
data.gcgr.fdr <- merge(merge(data.gcgr.sum %>%
                               filter(variable == 'gr', value > m, category == 'Transient') %>%
                               group_by(arm, condition, category, annotation) %>%
                               count() %>% data.frame(),
                             data.gcgr.sum %>%
                               filter(variable == 'gr', value > ul, category == 'Transient') %>%
                               group_by(arm, condition, category, annotation) %>%
                               count() %>% data.frame(),
                             by = c('arm','condition','category','annotation'), 
                             suffixes = c('_m','_ul')),
                       merge(data.gcgr.sum %>%
                               filter(variable == 'gr', value > ul, category == 'Transient', es > 10) %>%
                               group_by(arm, condition, category, annotation) %>%
                               count() %>% data.frame(),
                             data.gcgr.sum %>%
                               filter(variable == 'gr', value > ul, category == 'Transient', es > 30) %>%
                               group_by(arm, condition, category, annotation) %>%
                               count() %>% data.frame(),
                             by = c('arm','condition','category','annotation'),
                             suffixes = c('_es10','_es30'), all = T),
                       by = c('arm','condition','category','annotation'), all.x = T)
data.gcgr.fdr$fp <- data.gcgr.fdr$n_m * 0.005
data.gcgr.fdr$fdr <- data.gcgr.fdr$fp/data.gcgr.fdr$n_ul * 100
data.gcgr.fdr


##### PLOTTING MUTANTS WITH HIGH GR
data.gcgr.sum %>%
  filter(variable == 'gr', value > ul, category == 'Transient')

##### FIND MUTANTS WITH PHENOTYPES
data.gcgr.sum2 <- rbind(merge(data.gcgr.sum[,-8] %>%
              filter(arm == 'ONE', condition %in% c('HO','HU','SA')),
            data.gcgr.sum[,-8] %>%
              filter(arm == 'ONE', condition %in% c('GA')),
            by = c('arm','variable','rep','strain_id','orf_name','category','annotation'),
            suffixes = c('_stress','_cont')),
      merge(data.gcgr.sum[,-8] %>%
              filter(arm == 'TWO', condition %in% c('FL','TN')),
            data.gcgr.sum[,-8] %>%
              filter(arm == 'TWO', condition %in% c('DM')),
            by = c('arm','variable','rep','strain_id','orf_name','category','annotation'),
            suffixes = c('_stress','_cont')))

head(data.gcgr.sum2)

data.gcgr.sum2 %>%
  filter(variable %in% c('auc_l','r','k','gr'), category == 'Transient', 
         value_stress > ul_stress) %>%
  group_by(arm, condition_stress, category, variable) %>%
  count() %>% data.frame()

write.csv(data.gcgr.sum2 %>% filter(variable %notin% c('gr','dtime','lag')), 
          file = '/home/sbp29/R/Projects/adaptivefitness/output/translatome/tr_oe_growth_parameters.csv')
write.csv(data.gcgr.res, 
          file = '/home/sbp29/R/Projects/adaptivefitness/output/translatome/tr_oe_growth_parameters_all.csv')
