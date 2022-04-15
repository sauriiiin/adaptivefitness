##### DELETION SCREEN ANALYSIS
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 04/05/2022

source('/home/sbp29/R/Projects/adaptivefitness/scripts/translatome/initialize.R')

##### GATHER DATA
data.fit <- dbGetQuery(conn, 'select * from TRANSLATOME_DEL_FITNESS_DATA')
data.fit.sum <- dbGetQuery(conn, 'select * from TRANSLATOME_DEL_FITNESS_SUMMARY_DATA')

##### PLOT
# ### PLOTTING REF CS FIT DATA
# plot.ref.growth <- data.fit %>%
#   filter(orf_name == 'BF_control') %>%
#   ggplot(aes(x = hours, y = average)) +
#   geom_line(data = data.fit, aes(x = hours, y = average, group = pos),
#             col = 'grey50', alpha = .5, lwd = 0.2) +
#   geom_hline(yintercept = 1200, col = 'blue', linetype = 'dashed', lwd = 0.5) +
#   stat_summary(aes(group = category), fill = 'red',
#                fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.4) +
#   stat_summary(aes(group = category), col = 'red',
#                fun=median, geom="line", lwd =1) +
#   stat_summary(aes(y = avg_ratio * 1000, group = category), col = 'black',
#                fun=median, geom="line", lwd =1) +
#   # stat_summary(col = 'blue',
#   #              fun=max, geom="point", size = 1.5) +
#   scale_y_continuous("Colony Size", sec.axis = sec_axis(~ (./1000), name = "Ratio")) +
#   facet_grid(attempt ~ condition, scales = 'free_x') +
#   coord_cartesian(ylim = c(0,2500))
# ggsave(sprintf("%s/REF_GROWTH.jpg",fig_path), plot.ref.growth,
#        height = 250, width = 500, units = 'mm',
#        bg = 'white',
#        dpi = 300)
# 
# ### PLOTTING GROWTH RATIOS
# plot.all.ratio <- data.fit %>%
#   ggplot(aes(x = hours, y = avg_ratio)) +
#   stat_summary(aes(group = category, fill = category),
#                fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.4) +
#   stat_summary(aes(group = category, col = category),
#                fun=median, geom="line", lwd =1) +
#   scale_y_continuous("Ratio") +
#   facet_grid(attempt ~ condition, scales = 'free_x') +
#   theme_linedraw() +
#   theme(axis.title = element_text(size = titles),
#         axis.text = element_text(size = txt),
#         legend.title = element_text(size = titles),
#         legend.text = element_text(size = txt),
#         legend.position = 'bottom',
#         legend.key.size = unit(3, "mm"),
#         legend.box.spacing = unit(0.5,"mm"),
#         strip.text = ggtext::element_markdown(size = txt, colour = 'black',face = 'bold',
#                                               margin = margin(0.5,0,0.1,0, "mm")))
# ggsave(sprintf("%s/ALL_RATIO.jpg",fig_path), plot.all.ratio,
#        height = 250, width = 500, units = 'mm',
#        bg = 'white',
#        dpi = 300)

##### END-POINT FITNESS
data.fit.hrs <- data.fit %>%
  filter(category == 'Reference', !is.na(avg_ratio), data == 'cs') %>%
  group_by(data, attempt, condition, hours, rep) %>%
  summarise(avg_ratio = median(avg_ratio, na.rm = T), .groups = 'keep') %>%
  group_by(data, attempt, condition, hours) %>%
  summarise(avg_ratio = median(avg_ratio, na.rm = T), .groups = 'keep') %>%
  filter(avg_ratio > 1.05, avg_ratio < 1.1) %>%
  group_by(data, attempt, condition) %>%
  summarise(hours = max(hours), .groups = 'keep') %>%
  data.frame()

data.fit.sum2 <- rbind(merge(data.fit.sum, data.fit.hrs, by = c('data','attempt','condition','hours')),
                       data.fit.sum %>% filter(data != 'cs'))


##### ORF-WISE FITNESS RESULTS
data.fit.sum3 <- data.fit.sum2 %>%
  filter(strain_id > 0) %>%
  group_by(orf_name, category) %>%
  count() %>% data.frame()
data.fit.sum3 <- data.fit.sum3[,-3]

temp.sum <- data.fit.sum2 %>%
  group_by(data, attempt, condition, orf_name, category) %>%
  summarise(fitness = median(fit.median, na.rm = T),
            .groups = 'keep') %>%
  data.frame()

col.names <- c('orf_name','category')
for (d in unique(temp.sum$data)) {
  for (a in unique(temp.sum$attempt[temp.sum$data == d])) {
    for (c in unique(temp.sum$condition[temp.sum$data == d & temp.sum$attempt == a])) {
      temp <- temp.sum[temp.sum$data == d & temp.sum$attempt == a & temp.sum$condition == c,
                            c('orf_name','category','fitness')]
      data.fit.sum3 <- merge(data.fit.sum3, temp, by = c('orf_name','category'), all.x = T)
      col.names <- c(col.names, paste(d, a, c, sep = '_'))
    }
  }
}
colnames(data.fit.sum3) <- col.names
head(data.fit.sum3)

dbWriteTable(conn, 'TRANSLATOME_DEL_FITNESS_MATRIX_DATA', data.fit.sum3, overwrite = T)
data.fit.sum3 <- dbGetQuery(conn, 'select * from TRANSLATOME_DEL_FITNESS_MATRIX_DATA')

##### FITNESS DIFFERENTIAL
data.diff.ref <- NULL
for (d in unique(data.fit.sum2$data)) {
  for (a in unique(data.fit.sum2$attempt[data.fit.sum2$data == d])) {
    for (c in unique(data.fit.sum2$condition[data.fit.sum2$data == d & data.fit.sum2$attempt == a])) {
      temp <- data.fit$fitness[data.fit$data == d & data.fit$attempt == a & data.fit$condition == c &
                                 data.fit$orf_name == 'HO' & !is.na(data.fit$fitness)]
      if (c %in% c('FL','TN','DM')) {
        temp2 <- data.fit$fitness[data.fit$data == d & data.fit$attempt == a & data.fit$condition == 'DM' &
                                    data.fit$orf_name == 'HO' & !is.na(data.fit$fitness)]
      } else {
        temp2 <- data.fit$fitness[data.fit$data == d & data.fit$attempt == a & data.fit$condition == 'YPDA' &
                                    data.fit$orf_name == 'HO' & !is.na(data.fit$fitness)]
      }
      temp.diff.ref <- NULL
      for (i in 1:5000) {
        temp.diff.ref <- c(temp.diff.ref, median(sample(temp, 16) - sample(temp2, 16), na.rm = T))
      }
      data.diff.ref <- rbind(data.diff.ref, data.frame(data = d, attempt = a, condition = c, ref.diff = temp.diff.ref))
    }
  }
}


##### FITNESS DIFFERENTIAL WHEN CYCLING THROUGH P-VALUES
data.diff.fdr <- NULL
for (p_thresh in c(seq(0.1,0.9,0.1),seq(1,20,1))) {
  ul <- 1 - p_thresh/2/100
  ll <- p_thresh/2/100
  
  data.diff.ref.sum <- data.diff.ref %>%
    group_by(data, attempt, condition) %>%
    summarise(ref.diff.ll = quantile(ref.diff, ll, na.rm = T),
              ref.diff.m = median(ref.diff, na.rm = T),
              ref.diff.ul = quantile(ref.diff, ul, na.rm = T),
              .groups = 'keep') %>%
    data.frame()
  
  # colnames(data.fit.sum3)
  data.diff <- cbind(data.fit.sum3[,c("strain_id","orf_name","category")],
                     data.fit.sum3[,c("cs_init_FL","cs_init_TN")] - data.fit.sum3[,c("cs_init_DM")],
                     data.fit.sum3[,c("cs_init_HO","cs_init_HU","cs_init_SA")] - data.fit.sum3[,c("cs_init_YPDA")],
                     data.fit.sum3[,c("cs_valid_FL","cs_valid_TN")] - data.fit.sum3[,c("cs_valid_DM")],
                     data.fit.sum3[,c("cs_valid_HO","cs_valid_HU","cs_valid_SA")] - data.fit.sum3[,c("cs_valid_YPDA")],
                     data.fit.sum3[,c("growth_init_FL","growth_init_TN")] - data.fit.sum3[,c("growth_init_DM")],
                     data.fit.sum3[,c("growth_init_HO","growth_init_HU","growth_init_SA")] - data.fit.sum3[,c("growth_init_YPDA")],
                     data.fit.sum3[,c("growth_valid_FL","growth_valid_TN")] - data.fit.sum3[,c("growth_valid_DM")],
                     data.fit.sum3[,c("growth_valid_HO","growth_valid_HU","growth_valid_SA")] - data.fit.sum3[,c("growth_valid_YPDA")])
  data.diff <- melt(data.diff, id.vars = c("strain_id","orf_name","category"),
                    variable.name = 'sample', value.name = 'diff.fit')
  temp <- data.frame(cbind(sample = data.diff$sample %>% unique(), str_split(data.diff$sample %>% unique(), '_', simplify = T) %>% data.frame()))
  colnames(temp) <- c('sample','data','attempt','condition')
  data.diff <- merge(temp, data.diff, by = 'sample')
  data.diff <- merge(data.diff[,-1], data.diff.ref.sum, by = c('data','attempt','condition'))
  
  data.diff$phenotype[data.diff$diff.fit > data.diff$ref.diff.ul] <- 'Beneficial'
  data.diff$phenotype[data.diff$diff.fit < data.diff$ref.diff.ll] <- 'Deleterious'
  data.diff$phenotype[is.na(data.diff$phenotype) & !is.na(data.diff$diff.fit)] <- 'Neutral'
  
  # data.diff %>%
  #   mutate(delta.fit = diff.fit - ref.diff.m)
  
  ##### PHENOTYPE OVERLAP
  temp.diff.fdr <- merge(merge(merge(data.diff %>%
                                       filter(category == 'Transient', !is.na(phenotype), attempt != 'lin') %>%
                                       group_by(data, category, condition, strain_id, orf_name) %>%
                                       count() %>% filter(n > 1) %>%
                                       group_by(data, condition, category) %>%
                                       count() %>% data.frame(),
                                     data.diff %>%
                                       filter(category == 'Transient', !is.na(phenotype), attempt != 'lin') %>%
                                       group_by(data, category, condition, phenotype, strain_id, orf_name) %>%
                                       count() %>% filter(n > 1) %>%
                                       group_by(data, condition, category) %>%
                                       count() %>% data.frame(),
                                     by = c('data','condition','category'),
                                     suffixes = c('_total','_common_phenotype'),
                                     all = T),
                               merge(data.diff %>%
                                       filter(category == 'Transient', !is.na(phenotype), attempt != 'lin',
                                              diff.fit > ref.diff.m) %>%
                                       group_by(data, category, condition, strain_id, orf_name) %>%
                                       count() %>% filter(n > 1) %>%
                                       group_by(data, condition, category) %>%
                                       count() %>% data.frame(), 
                                     data.diff %>%
                                       filter(category == 'Transient', !is.na(phenotype), attempt != 'lin',
                                              diff.fit <= ref.diff.m) %>%
                                       group_by(data, category, condition, strain_id, orf_name) %>%
                                       count() %>% filter(n > 1) %>%
                                       group_by(data, condition, category) %>%
                                       count() %>% data.frame(),
                                     by = c('data','condition','category'),
                                     suffixes = c('_max_possible_ben', '_max_possible_del'),
                                     all = T),
                               by = c('data','condition','category'),
                               all = T),
                         merge(data.diff %>%
                                 filter(category == 'Transient', !is.na(phenotype), attempt != 'lin',
                                        diff.fit > ref.diff.ul) %>%
                                 group_by(data, category, condition, strain_id, orf_name) %>%
                                 count() %>% filter(n > 1) %>%
                                 group_by(data, condition, category) %>%
                                 count() %>% data.frame(), 
                               data.diff %>%
                                 filter(category == 'Transient', !is.na(phenotype), attempt != 'lin',
                                        diff.fit <= ref.diff.ll) %>%
                                 group_by(data, category, condition, strain_id, orf_name) %>%
                                 count() %>% filter(n > 1) %>%
                                 group_by(data, condition, category) %>%
                                 count() %>% data.frame(),
                               by = c('data','condition','category'),
                               suffixes = c('_beneficial_phenotype', '_deleterious_phenotype'),
                               all = T),
                         by = c('data','condition','category'),
                         all = T)
  temp.diff.fdr$p_value <- p_thresh
  temp.diff.fdr <- temp.diff.fdr %>%
    mutate(beneficial_fdr = n_max_possible_ben * (p_thresh/2)/n_beneficial_phenotype,
           deleterious_fdr = n_max_possible_del * (p_thresh/2)/n_deleterious_phenotype)
  data.diff.fdr <- rbind(data.diff.fdr, temp.diff.fdr)
}

data.diff.fdr %>%
  ggplot(aes(x = p_value, y = beneficial_fdr)) +
  geom_hline(yintercept = 5) +
  geom_line(aes(col = condition)) +
  scale_y_log10() +
  facet_wrap(.~data)

data.diff.fdr %>%
  filter(data == 'growth', p_value == 5)

dbWriteTable(conn, 'TRANSLATOME_DEL_FITNESS_DIFF_PHENOTYPE_COUNTS', data.diff.fdr, overwrite = T)

##### ADDING EFFECT SIZE TO THE DIFF FITNESS EFFECTS
data.diff$es <- NULL
for (d in unique(data.diff.ref$data)) {
  for (a in unique(data.diff.ref$attempt[data.diff.ref$data == d])) {
    for (c in unique(data.diff.ref$condition[data.diff.ref$data == d & data.diff.ref$attempt == a])) {
      temp.ref <- data.diff.ref$ref.diff[data.diff.ref$data == d & data.diff.ref$attempt == a & data.diff.ref$condition == c]
      for (s in unique(data.diff$strain_id[data.diff$data == d & data.diff$attempt == a & data.diff$condition == c &
                                           !is.na(data.diff$diff.fit)])) {
        percentile <- ecdf(temp.ref)
        temp.es <- percentile(data.diff$diff.fit[data.diff$data == d & data.diff$attempt == a & data.diff$condition == c & data.diff$strain_id == s])
        data.diff$es[data.diff$data == d & data.diff$attempt == a & data.diff$condition == c & data.diff$strain_id == s] <- temp.es
      }
    }
  }
}

dbWriteTable(conn, 'TRANSLATOME_DEL_FITNESS_DIFF_SUMMARY_DATA', data.diff, overwrite = T)

# data.diff %>%
#   filter(phenotype == 'Beneficial', category == 'Transient', condition == 'SA', data == 'growth') %>%
#   group_by(data, condition, strain_id, orf_name) %>%
#   count() %>% filter(n > 2) %>% data.frame()
# 
# data.diff %>%
#   filter(category %in% c('Transient','Conserved'), data != 'auc') %>%
#   ggplot(aes(y = diff.fit, x = 1, col = phenotype)) +
#   geom_jitter() +
#   facet_grid(category ~ data * condition)

data.diff2 <- data.fit.sum3[,c('strain_id','orf_name','category')]
col.names <- c('strain_id','orf_name','category')
for (d in unique(data.diff$data)) {
  for (a in unique(data.diff$attempt[data.diff$data == d])) {
    for (c in unique(data.diff$condition[data.diff$data == d & data.diff$attempt == a])) {
      temp <- data.diff[data.diff$data == d & data.diff$attempt == a & data.diff$condition == c,
                        c('strain_id','orf_name','diff.fit','phenotype','es')] %>% filter(strain_id > 0)
      data.diff2 <- merge(data.diff2, temp, by = c('strain_id','orf_name'), all.x = T)
      col.names <- c(col.names, paste(d, a, c, sep = '_'), paste(d, a, c, 'p', sep = '_'), paste(d, a, c, 'es', sep = '_'))
    }
  }
}
colnames(data.diff2) <- col.names
head(data.diff2)

dbWriteTable(conn, 'TRANSLATOME_OE_FITNESS_DIFF_MATRIX_DATA', data.diff2, overwrite = T)


###### COLONY TO COLONY FITNESS CORRELATION
data.fit.hrs

merge(data.fit[data.fit$data == 'cs' & data.fit$attempt == 'init' & data.fit$condition == 'GA' & data.fit$hours == 21,
               c('data','condition','pos','strain_id','orf_name','category','fitness')],
      data.fit[data.fit$data == 'cs' & data.fit$attempt == 'init' & data.fit$condition == 'GA2' & data.fit$hours == 24,
               c('data','condition','pos','strain_id','orf_name','category','fitness')],
      by = c('data','pos','strain_id','orf_name','category'),
      suffix = c('_ga','_ga2')) %>%
  filter(category != 'ConditionControls') %>%
  ggplot(aes(x = fitness_ga, y = fitness_ga2)) +
  geom_point() +
  geom_abline() +
  stat_cor(method = 'pearson') +
  facet_wrap(.~category, nrow = 1)


###### REP TO REP FITNESS CORRELATION
merge(data.fit.sum2[data.fit.sum2$attempt == 'init' & data.fit.sum2$condition == 'GA',
                    c('data','condition','rep','strain_id','orf_name','category','fit.median')],
      data.fit.sum2[data.fit.sum2$attempt == 'init' & data.fit.sum2$condition == 'GA2',
                    c('data','condition','rep','strain_id','orf_name','category','fit.median')],
      by = c('data','rep','strain_id','orf_name','category'),
      suffix = c('_ga','_ga2')) %>%
  filter(category != 'ConditionControls', data == 'cs') %>%
  ggplot(aes(x = fit.median_ga, y = fit.median_ga2)) +
  geom_point() +
  geom_abline() +
  stat_cor(method = 'pearson') +
  facet_wrap(.~category, nrow = 1)


###### STRAIN_ID to STRAIN_ID FITNESS CORRELATION
data.fit.sum3 %>%
  filter(category != 'ConditionControls', category == 'Transient') %>%
  mutate(phenotype = paste(growth_init_SA_p, growth_valid_SA_p)) %>%
  filter(!is.na(category)) %>%
  ggplot(aes(x = growth_init_SA, y = growth_valid_SA)) +
  geom_abline() +
  geom_point() +
  scale_color_discrete(guide = 'none') +
  # geom_smooth(method = lm) +
  stat_cor(method = 'spearman') +
  facet_wrap(.~category, nrow = 1) #+
# coord_cartesian(xlim = c(-0.5,0.5),
#                 ylim = c(-0.5,0.5))
# coord_cartesian(xlim = c(0,1.2),
#                 ylim = c(0,1.2))


###### STRAIN_ID to STRAIN_ID FITNESS DIFF CORRELATION
data.diff2 %>%
  filter(category != 'ConditionControls', category == 'Transient') %>%
  mutate(phenotype = paste(growth_init_SA_p, growth_valid_SA_p)) %>%
  filter(!is.na(category)) %>%
  ggplot(aes(x = growth_init_SA, y = growth_valid_SA)) +
  geom_abline() +
  geom_point() +
  # geom_point(aes(col = phenotype)) +
  scale_color_discrete(guide = 'none') +
  # geom_smooth(method = lm) +
  stat_cor(method = 'spearman') +
  facet_wrap(.~category, nrow = 1) #+
# coord_cartesian(xlim = c(-0.5,0.5),
#                 ylim = c(-0.5,0.5))
# coord_cartesian(xlim = c(0,1.2),
#                 ylim = c(0,1.2))


###### STRAIN_ID to STRAIN_ID FITNESS DIFF ES CORRELATION
data.diff2 %>%
  filter(category != 'ConditionControls', category == 'Transient') %>%
  mutate(phenotype = paste(growth_init_SA_p, growth_valid_SA_p)) %>%
  filter(!is.na(category)) %>%
  ggplot(aes(x = growth_init_SA_es, y = growth_valid_SA_es)) +
  geom_abline() +
  # geom_point() +
  geom_point(aes(col = phenotype)) +
  scale_color_discrete(guide = 'none') +
  stat_cor(method = 'spearman') +
  facet_wrap(.~category, nrow = 1) #+
# coord_cartesian(xlim = c(0.9,1),
#                 ylim = c(0.9,1))
# coord_cartesian(xlim = c(0,1.2),
#                 ylim = c(0,1.2))


##### INTERCHANGING REPLICATES
data.fit.cor <- merge(data.fit, data.fit.hrs, by = c('attempt','condition','hours'))
# use mutate to add a column where fitness values per attempt per condition per rep are randomly exchanged
arrange(data.fit.cor, attempt, condition, rep) %>%
  mutate(fitness_ = sample(fitness))


for (a in unique(data.fit.cor$attempt)) {
  for (c in unique(data.fit.cor$condition[data.fit.cor$attempt == a])) {
    for (r in unique(data.fit.cor$rep[data.fit.cor$attempt == a & data.fit.cor$condition == c])) {
      data.fit.cor$fitness_[data.fit.cor$attempt == a & data.fit.cor$condition == c & data.fit.cor$rep == r] <-
        sample(data.fit.cor$fitness[data.fit.cor$attempt == a & data.fit.cor$condition == c & data.fit.cor$rep == r])
    }
  }
}

hi <- merge(data.fit.cor[data.fit.cor$attempt == 'init', 
                         c('condition','strain_id','orf_name','category','fitness','fitness_','plate_no','plate_col','plate_row')],
            data.fit.cor[data.fit.cor$attempt == 'lin' , 
                         c('condition','strain_id','orf_name','fitness','fitness_','plate_no','plate_col','plate_row')],
            by = c('condition','strain_id','orf_name','plate_no','plate_col','plate_row'), suffixes = c('_init','_lin'))
head(hi)

hi %>%
  ggplot(aes(x = fitness_init, y = fitness_lin)) +
  geom_point() +
  stat_cor(method = 'spearman') +
  facet_grid(condition~category)

hi %>%
  ggplot(aes(x = fitness__init, y = fitness__lin)) +
  geom_point() +
  stat_cor(method = 'spearman') +
  facet_grid(condition~category)

##### ORIGINAL BF VS EXPANSION
data.fit.bf.ori.expan <- merge(data.fit.sum2[,c('attempt','condition','strain_id','orf_name','category','fit.median')] %>%
                                 filter(orf_name %in% c('YNL146W',
                                                        'YJL189W',
                                                        'YHR199C-A',
                                                        'YER014C-A',
                                                        'YER044C-A',
                                                        'YDL061C',
                                                        'YML101C-A',
                                                        'YIL009C-A'), strain_id < 10000),
                               data.fit.sum2[,c('attempt','condition','strain_id','orf_name','fit.median')] %>%
                                 filter(orf_name %in% c('YNL146W',
                                                        'YJL189W',
                                                        'YHR199C-A',
                                                        'YER014C-A',
                                                        'YER044C-A',
                                                        'YDL061C',
                                                        'YML101C-A',
                                                        'YIL009C-A'), strain_id > 10000),
                               by = c('attempt','condition','orf_name'), suffixes = c('_original','_expanded')) 
data.fit.bf.ori.expan$fit.diff <- abs(data.fit.bf.ori.expan$fit.median_expanded - data.fit.bf.ori.expan$fit.median_original)

data.fit.ref.diff <- NULL
for (a in unique(data.fit.sum2$attempt)) {
  for (c in unique(data.fit.sum2$condition[data.fit.sum2$attempt == a])) {
    temp <- NULL
    for (i in 1:1000) {
      bf1 <- sample(data.fit.sum2$fit.median[data.fit.sum2$attempt == a &
                                               data.fit.sum2$condition == c &
                                               data.fit.sum2$orf_name == 'BF_control' &
                                               !is.na(data.fit.sum2$fit.median)], 1)
      bf2 <- sample(data.fit.sum2$fit.median[data.fit.sum2$attempt == a &
                                               data.fit.sum2$condition == c &
                                               data.fit.sum2$orf_name == 'BF_control'&
                                               !is.na(data.fit.sum2$fit.median)], 1)
      temp <- rbind(temp, abs(bf1 - bf2))
    }
    data.fit.ref.diff <- rbind(data.fit.ref.diff, 
                               data.frame(attempt = a, condition = c, ref_diff = temp))
  }
}

data.fit.ref.diff %>%
  filter(attempt != 'valid') %>%
  ggplot(aes(y = ref_diff)) +
  geom_violin(aes(x = 1)) +
  geom_jitter(data = data.fit.bf.ori.expan, 
              aes(y = fit.diff, x = 1, col = category)) +
  scale_y_log10() +
  facet_grid(attempt~condition)

