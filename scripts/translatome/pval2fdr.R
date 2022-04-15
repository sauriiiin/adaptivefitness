# fit.sum <- data.fit.sum2[,c(1:10)]
# ref.data <- data.fit.sum2 %>% filter(orf_name == 'BF_control')
# attempts <- c('init','lin')
# p_values <- c(seq(0.1,0.9,0.1),seq(1,20,1))


pval2fdr <- function(fit.sum, ref.data, attempts, p_values) {
  fit.pheno <- NULL
  fit.fdr <- NULL
  fit.lim <- NULL
  
  fit.sum <- fit.sum %>%
    group_by(data, attempt, condition, strain_id, orf_name, category) %>%
    summarise(fit.median = median(fit.median, na.rm = T), .groups = 'keep') %>%
    data.frame()
  
  for (p_thresh in p_values) {
    ul <- 1 - p_thresh/2/100
    ll <- p_thresh/2/100
    
    temp.lim <- ref.data %>%
      group_by(data, attempt, condition) %>%
      summarise(fit.ll = quantile(fit.median, ll, na.rm = T),
                fit.m = median(fit.median, na.rm = T),
                fit.ul = quantile(fit.median, ul, na.rm = T),
                .groups = 'keep') %>%
      data.frame()
    fit.lim <- rbind(fit.lim, data.frame(cbind(p_value = p_thresh/100, temp.lim)))
    
    temp.fit.sum <- merge(fit.sum, temp.lim, by = c('data','attempt','condition'))
    
    temp.fit.sum$phenotype[temp.fit.sum$fit.median > temp.fit.sum$fit.ul] <- 'Beneficial'
    temp.fit.sum$phenotype[temp.fit.sum$fit.median < temp.fit.sum$fit.ll] <- 'Deleterious'
    temp.fit.sum$phenotype[is.na(temp.fit.sum$phenotype) & !is.na(temp.fit.sum$fit.median)] <- 'Neutral'
    
    temp.fit.pheno <- merge(merge(temp.fit.sum %>%
                                    filter(category %in% c('Transient','Conserved','Others'), !is.na(phenotype)) %>%
                                    group_by(data, attempt, condition) %>%
                                    count(),
                                  temp.fit.sum %>%
                                    filter(category %in% c('Transient','Conserved','Others'), !is.na(phenotype)) %>%
                                    group_by(data, attempt, condition, category) %>%
                                    count(), 
                                  by = c('data','attempt','condition'),
                                  suffixes = c('_total','_per_category'), all = T),
                            merge(merge(temp.fit.sum %>%
                                          filter(category %in% c('Transient','Conserved','Others'), !is.na(phenotype)) %>%
                                          group_by(data, attempt, condition, category) %>%
                                          filter(fit.median > fit.m) %>%
                                          count(), 
                                        temp.fit.sum %>%
                                          filter(category %in% c('Transient','Conserved','Others'), !is.na(phenotype)) %>%
                                          group_by(data, attempt, condition, category) %>%
                                          filter(fit.median < fit.m) %>%
                                          count(),
                                        by = c('data','attempt','condition','category'),
                                        suffixes = c('_max_possible_ben', '_max_possible_del'), all = T),
                                  merge(temp.fit.sum %>%
                                          filter(category %in% c('Transient','Conserved','Others'), !is.na(phenotype)) %>%
                                          group_by(data, attempt, condition, category) %>%
                                          filter(fit.median > fit.ul) %>%
                                          count(),
                                        temp.fit.sum %>%
                                          filter(category %in% c('Transient','Conserved','Others'), !is.na(phenotype)) %>%
                                          group_by(data, attempt, condition, category) %>%
                                          filter(fit.median < fit.ll) %>%
                                          count(),
                                        by = c('data','attempt','condition','category'),
                                        suffixes = c('_beneficial_phenotype', '_deleterious_phenotype'), all = T),
                                  by = c('data','attempt','condition','category'), all = T),
                            by = c('data','attempt','condition','category'), all = T)
    
    temp.fit.fdr <- merge(merge(merge(temp.fit.sum %>%
                                        filter(category %in% c('Transient','Conserved','Others'), !is.na(phenotype), attempt %in% attempts) %>%
                                        group_by(data, category, condition, strain_id, orf_name) %>%
                                        count() %>% filter(n > 1) %>%
                                        group_by(data, condition, category) %>%
                                        count() %>% data.frame(),
                                      temp.fit.sum %>%
                                        filter(category %in% c('Transient','Conserved','Others'), !is.na(phenotype), attempt %in% attempts) %>%
                                        group_by(data, category, condition, phenotype, strain_id, orf_name) %>%
                                        count() %>% filter(n > 1) %>%
                                        group_by(data, condition, category) %>%
                                        count() %>% data.frame(),
                                      by = c('data','condition','category'),
                                      suffixes = c('_total','_common_phenotype'),
                                      all = T),
                                merge(temp.fit.sum %>%
                                        filter(category %in% c('Transient','Conserved','Others'), !is.na(phenotype), attempt %in% attempts,
                                               fit.median > fit.m) %>%
                                        group_by(data, category, condition, strain_id, orf_name) %>%
                                        count() %>% filter(n > 1) %>%
                                        group_by(data, condition, category) %>%
                                        count() %>% data.frame(), 
                                      temp.fit.sum %>%
                                        filter(category %in% c('Transient','Conserved','Others'), !is.na(phenotype), attempt %in% attempts,
                                               fit.median <= fit.m) %>%
                                        group_by(data, category, condition, strain_id, orf_name) %>%
                                        count() %>% filter(n > 1) %>%
                                        group_by(data, condition, category) %>%
                                        count() %>% data.frame(),
                                      by = c('data','condition','category'),
                                      suffixes = c('_max_possible_ben', '_max_possible_del'),
                                      all = T),
                                by = c('data','condition','category'),
                                all = T),
                          merge(merge(temp.fit.sum %>%
                                        filter(category %in% c('Transient','Conserved','Others'), phenotype == 'Beneficial', attempt == attempts[1]) %>%
                                        group_by(data,condition,category) %>%
                                        count() %>% data.frame(),
                                      temp.fit.sum %>%
                                        filter(category %in% c('Transient','Conserved','Others'), phenotype == 'Beneficial', attempt == attempts[2]) %>%
                                        group_by(data,condition,category) %>%
                                        count() %>% data.frame(),
                                      by = c('data','condition','category'),
                                      suffixes = c('_beneficial_attempt1','_beneficial_attempt2')),
                                merge(temp.fit.sum %>%
                                        filter(category %in% c('Transient','Conserved','Others'), !is.na(phenotype), attempt %in% attempts,
                                               fit.median > fit.ul) %>%
                                        group_by(data, category, condition, strain_id, orf_name) %>%
                                        count() %>% filter(n > 1) %>%
                                        group_by(data, condition, category) %>%
                                        count() %>% data.frame(), 
                                      temp.fit.sum %>%
                                        filter(category %in% c('Transient','Conserved','Others'), !is.na(phenotype), attempt %in% attempts,
                                               fit.median <= fit.ll) %>%
                                        group_by(data, category, condition, strain_id, orf_name) %>%
                                        count() %>% filter(n > 1) %>%
                                        group_by(data, condition, category) %>%
                                        count() %>% data.frame(),
                                      by = c('data','condition','category'),
                                      suffixes = c('_beneficial_phenotype', '_deleterious_phenotype'),
                                      all = T),
                                by = c('data','condition','category'), all = T),
                          by = c('data','condition','category'),
                          all = T)
    
    temp.fit.pheno$p_value <- p_thresh/100
    temp.fit.pheno <- temp.fit.pheno %>%
      mutate(beneficial_fdr = n_max_possible_ben * (p_thresh/2)/n_beneficial_phenotype,
             deleterious_fdr = n_max_possible_del * (p_thresh/2)/n_deleterious_phenotype) 
    
    temp.fit.fdr$p_value <- p_thresh/100
    temp.fit.fdr <- temp.fit.fdr %>%
      mutate(beneficial_fdr = n_max_possible_ben * (p_thresh/2)/n_beneficial_phenotype,
             deleterious_fdr = n_max_possible_del * (p_thresh/2)/n_deleterious_phenotype)  
    
    fit.pheno <- rbind(fit.pheno, temp.fit.pheno)
    fit.fdr <- rbind(fit.fdr, temp.fit.fdr)
  }
  output <- list(fit.lim, fit.pheno, fit.fdr)
  return(output)
} 

