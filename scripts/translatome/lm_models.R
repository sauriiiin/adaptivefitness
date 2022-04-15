
head(data.fit.sum2)

data.fit.sum4 <- rbind(merge(data.fit.sum2[data.fit.sum2$condition %in% c('FL','TN'),-15],
            data.fit.sum2[data.fit.sum2$condition %in% c('DM'),-15],
            by = c('data','attempt','rep','strain_id','orf_name','category'),
            suffixes = c('_stress','_nostress')),
      merge(data.fit.sum2[data.fit.sum2$condition %in% c('SA','HU','HO'),-15],
            data.fit.sum2[data.fit.sum2$condition %in% c('GA'),-15],
            by = c('data','attempt','rep','strain_id','orf_name','category'),
            suffixes = c('_stress','_nostress')))
head(data.fit.sum4)

fitness.mod1 <- lm(fit.median_stress ~ condition_stress + orf_name,
                   data = data.fit.sum4 %>% filter(data == 'cs', attempt == 'init', category %in% c('Conserved','Transient','Others'), 
                                                   !is.na(fit.median_stress), !is.na(fit.median_nostress)))
fitness.mod2 <- lm(fit.median_stress ~ condition_stress + orf_name + fit.median_nostress,
                   data = data.fit.sum4 %>% filter(data == 'cs', attempt == 'init', category %in% c('Conserved','Transient','Others'),
                                                   !is.na(fit.median_stress), !is.na(fit.median_nostress)))

anova(fitness.mod1)
anova(fitness.mod2)
anova(fitness.mod1, fitness.mod2)

library(lme4)
fitness.mmod1 <- lmer(fit.median_stress ~ condition_stress + orf_name + (condition_stress | orf_name), 
                     data = data.fit.sum4[!is.na(data.fit.sum4$fit.median_stress) & !is.na(data.fit.sum4$fit.median_nostress),] %>%
                       filter(data == 'cs', attempt == 'init', category %in% c('Conserved','Transient','Others'),
                              !is.na(fit.median_stress), !is.na(fit.median_nostress)))
fitness.mmod2 <- lmer(fit.median_stress ~ condition_stress + (condition_stress | orf_name), 
                     data = data.fit.sum4[!is.na(data.fit.sum4$fit.median_stress) & !is.na(data.fit.sum4$fit.median_nostress),] %>%
                       filter(data == 'cs', attempt == 'init', category %in% c('Conserved','Transient','Others'),
                              !is.na(fit.median_stress), !is.na(fit.median_nostress)))
fitness.mmod3 <- lmer(fit.median_stress ~ condition_stress + fit.median_nostress + orf_name + (condition_stress | orf_name), 
                     data = data.fit.sum4[!is.na(data.fit.sum4$fit.median_stress) & !is.na(data.fit.sum4$fit.median_nostress),] %>% 
                       filter(data == 'cs', attempt == 'init', category %in% c('Conserved','Transient','Others'),
                              !is.na(fit.median_stress), !is.na(fit.median_nostress)))

anova(fitness.mmod1)
anova(fitness.mmod2)
anova(fitness.mmod3)


library(nlme)
fitness.mod2 <- lme(fit.median_stress ~ condition_stress + fit.median_nostress + orf_name, random= ~1|category, 
                    data = data.fit.sum4[!is.na(data.fit.sum4$fit.median_stress) & !is.na(data.fit.sum4$fit.median_nostress),] %>% 
                      filter(data == 'cs', attempt == 'init', category %in% c('Conserved','Transient','Others')))
anova(fitness.mod2)

