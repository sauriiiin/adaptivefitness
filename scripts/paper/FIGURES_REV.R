##### LID RESPONSE TO REVIEWERS
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 11/20/2020

##### INITIALIZE
library(ggplot2)
library(ggridges)
library(gridExtra)
library(grid)
library(tidyverse)
library(ggpubr)
library(stringr)
library(scales)
library(egg)
library(zoo)
library(ggrepel)
library(reshape2)
library(qvalue)
library(RMariaDB)

source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")
out_path = 'figs/paper/';

load(file = sprintf('%sRND2_RESULTS.RData',out_path))

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 7
txt <- 7
lbls <- 9

##### GATHER DATA
rnd2_data <- dbGetQuery(conn, 'select lid.*, mcat.fitness mcat_fitness, mcat.mcat_p, mcat.mcat_stat
                        from
                        (select c.*, d.p lid_p, d.stat lid_stat
                        from 
                        (select b.*, a.fitness lid_fitness
                        from 4C4_FS_RND2_6144_FITNESS a, 4C4_FS_RND2_6144_DATA b
                        where a.orf_name is not NULL and a.orf_name != "BFC100"
                        and a.pos = b.pos and a.hours = b.hours
                        and a.fitness is not NULL
                        and a.hours > 0) c
                        LEFT JOIN
                        (select * from 4C4_FS_RND2_6144_PVALUE) d
                        on
                        c.hours = d.hours and c.orf_name = d.orf_name
                        order by c.hours, c.orf_name) lid
                        JOIN
                        (select c.*, d.p mcat_p, d.stat mcat_stat
                        from 
                        (select b.*, a.fitness
                        from 4C4_FS_RND2_BEAN_6144_FITNESS a, 4C4_FS_RND2_6144_DATA b
                        where a.orf_name is not NULL and a.orf_name != "BFC100"
                        and a.pos = b.pos and a.hours = b.hours
                        and a.fitness is not NULL
                        and a.hours > 0) c
                        LEFT JOIN
                        (select * from 4C4_FS_RND2_BEAN_6144_PVALUE) d
                        on
                        c.hours = d.hours and c.orf_name = d.orf_name
                        order by c.hours, c.orf_name) mcat
                        on lid.pos = mcat.pos and lid.hours = mcat.hours
                        order by hours, plate, col, row')

rnd2_res <- NULL
i <- 1
for (h in sort(unique(rnd2_data$hours))) {
  rnd2_data$lid_pthresh[rnd2_data$hours == h] <- unique(lid.rnd2.res$pthresh[lid.rnd2.res$hours == h & lid.rnd2.res$rep == 16 & lid.rnd2.res$ref == 0.25])
  rnd2_data$mcat_pthresh[rnd2_data$hours == h] <- unique(mcat.rnd2.res$pthresh[mcat.rnd2.res$hours == h])
  for (o in sort(unique(rnd2_data$orf_name[rnd2_data$hours == h & rnd2_data$orf_name != 'BF_control']))){
    rnd2_res$hours[i] <- h
    rnd2_res$orf_name[i] <- o
    rnd2_res$rnd_hrs[i] <- unique(rnd2_data$rnd_hrs[rnd2_data$hours == h & rnd2_data$orf_name == o])
    
    rnd2_res$lid_p[i] <- unique(lid.rnd2.res$p[lid.rnd2.res$hours == h & lid.rnd2.res$rep == 16 & lid.rnd2.res$ref == 0.25 & lid.rnd2.res$orf_name == o])
    rnd2_res$lid_stat[i] <- unique(rnd2_data$lid_stat[rnd2_data$hours == h & rnd2_data$orf_name == o])
    rnd2_res$lid_es[i] <- unique(lid.rnd2.res$es[lid.rnd2.res$hours == h & lid.rnd2.res$rep == 16 & lid.rnd2.res$ref == 0.25 & lid.rnd2.res$orf_name == o])
    rnd2_res$lid_pthresh[i] <- unique(lid.rnd2.res$pthresh[lid.rnd2.res$hours == h & lid.rnd2.res$rep == 16 & lid.rnd2.res$ref == 0.25])
    
    ref_fit <- rnd2_data$lid_fitness[rnd2_data$hours == h & rnd2_data$orf_name == 'BF_control' & !is.na(rnd2_data$lid_fitness)]
    orf_fit <- rnd2_data$lid_fitness[rnd2_data$hours == h & rnd2_data$orf_name == o & !is.na(rnd2_data$lid_fitness)]
    wilcox_stats <- wilcox.test(orf_fit, ref_fit)
    # rnd2_data$lid_wilcox_p[rnd2_data$hours == h & rnd2_data$orf_name == o] <- wilcox_stats$p.value
    rnd2_res$lid_wilcox_p[i] <- wilcox_stats$p.value
    
    rnd2_res$mcat_p[i] <- unique(rnd2_data$mcat_p[rnd2_data$hours == h & rnd2_data$orf_name == o])
    rnd2_res$mcat_stat[i] <- unique(rnd2_data$mcat_stat[rnd2_data$hours == h & rnd2_data$orf_name == o])
    rnd2_res$mcat_es[i] <- unique(mcat.rnd2.res$es[mcat.rnd2.res$hours == h & mcat.rnd2.res$orf_name == o])
    rnd2_res$mcat_pthresh[i] <- unique(mcat.rnd2.res$pthresh[mcat.rnd2.res$hours == h])
    
    ref_fit <- rnd2_data$mcat_fitness[rnd2_data$hours == h & rnd2_data$orf_name == 'BF_control' & !is.na(rnd2_data$mcat_fitness)]
    orf_fit <- rnd2_data$mcat_fitness[rnd2_data$hours == h & rnd2_data$orf_name == o & !is.na(rnd2_data$mcat_fitness)]
    wilcox_stats <- wilcox.test(orf_fit, ref_fit)
    # rnd2_data$mcat_wilcox_p[rnd2_data$hours == h & rnd2_data$orf_name == o] <- wilcox_stats$p.value
    rnd2_res$mcat_wilcox_p[i] <- wilcox_stats$p.value
    
    i <- i + 1
  }
  # rnd2_data$lid_adj_p[rnd2_data$hours == h] <- p.adjust(rnd2_data$lid_wilcox_p[rnd2_data$hours == h], method = "bonferroni")
  # rnd2_data$mcat_adj_p[rnd2_data$hours == h] <- p.adjust(rnd2_data$mcat_wilcox_p[rnd2_data$hours == h], method = "bonferroni")
  rnd2_res$lid_q[rnd2_res$hours == h] <- qvalue(rnd2_res$lid_wilcox_p[rnd2_res$hours == h], fdr.level=0.05, pi0 = 1)$qvalue
  rnd2_res$mcat_q[rnd2_res$hours == h] <- qvalue(rnd2_res$mcat_wilcox_p[rnd2_res$hours == h], fdr.level=0.05, pi0 = 1)$qvalue
}
rnd2_res <- data.frame(rnd2_res)
head(rnd2_res)

rnd2_res$phenotype[rnd2_res$hours < rnd2_res$rnd_hrs] <- 'Beneficial'
rnd2_res$phenotype[rnd2_res$hours > rnd2_res$rnd_hrs] <- 'Deleterious'
rnd2_res$phenotype[rnd2_res$hours == rnd2_res$rnd_hrs] <- 'Neutral'

rnd2_res$lid_phenotype[rnd2_res$lid_p <= rnd2_res$lid_pthresh & rnd2_res$lid_es > 1] <- 'Beneficial'
rnd2_res$lid_phenotype[rnd2_res$lid_p <= rnd2_res$lid_pthresh & rnd2_res$lid_es < 1] <- 'Deleterious'
rnd2_res$lid_phenotype[is.na(rnd2_res$lid_phenotype)] <- 'Neutral'
dim(rnd2_res[rnd2_res$phenotype != rnd2_res$lid_phenotype,])[[1]]

rnd2_res$lid_wilcox_phenotype[rnd2_res$lid_wilcox_p <= 0.05 & rnd2_res$lid_q <= 0.05 & rnd2_res$lid_es > 1] <- 'Beneficial'
rnd2_res$lid_wilcox_phenotype[rnd2_res$lid_wilcox_p <= 0.05 & rnd2_res$lid_q <= 0.05 & rnd2_res$lid_es < 1] <- 'Deleterious'
rnd2_res$lid_wilcox_phenotype[is.na(rnd2_res$lid_wilcox_phenotype)] <- 'Neutral'

rnd2_res$mcat_phenotype[rnd2_res$mcat_p <= rnd2_res$mcat_pthresh & rnd2_res$mcat_es > 1] <- 'Beneficial'
rnd2_res$mcat_phenotype[rnd2_res$mcat_p <= rnd2_res$mcat_pthresh & rnd2_res$mcat_es < 1] <- 'Deleterious'
rnd2_res$mcat_phenotype[is.na(rnd2_res$mcat_phenotype)] <- 'Neutral'
dim(rnd2_res[rnd2_res$phenotype != rnd2_res$mcat_phenotype,])[[1]]

rnd2_res$mcat_wilcox_phenotype[rnd2_res$mcat_wilcox_p <= 0.05 & rnd2_res$mcat_q <= 0.05 & rnd2_res$mcat_es > 1] <- 'Beneficial'
rnd2_res$mcat_wilcox_phenotype[rnd2_res$mcat_wilcox_p <= 0.05 & rnd2_res$mcat_q <= 0.05 & rnd2_res$mcat_es < 1] <- 'Deleterious'
rnd2_res$mcat_wilcox_phenotype[is.na(rnd2_res$mcat_wilcox_phenotype)] <- 'Neutral'


###### QC
plyr::count(lid.rnd2.res[lid.rnd2.res$rep == 16 & lid.rnd2.res$ref == 0.25 & lid.rnd2.res$attempt == 1,],
            vars = c('hours','effect'))
plyr::count(rnd2_res, vars = c('hours','lid_phenotype'))

hello <- merge(lid.rnd2.res[lid.rnd2.res$rep == 16 & lid.rnd2.res$ref == 0.25 & lid.rnd2.res$attempt == 1,],
      rnd2_res, by = c('hours','orf_name'), suffix = c('old','new'))

hi <- hello[hello$effect != hello$lid_phenotype,]

######
head(rnd2_res)

pheno_mismatch <- NULL
pheno_mismatch$hours <- sort(unique(rnd2_res$hours))
pheno_mismatch <- merge(pheno_mismatch,
                        plyr::count(rnd2_res[rnd2_res$phenotype != rnd2_res$lid_phenotype,], vars = 'hours'),
                        by = 'hours', all = T)
pheno_mismatch <- merge(pheno_mismatch,
                        plyr::count(rnd2_res[rnd2_res$phenotype != rnd2_res$lid_wilcox_phenotype,], vars = 'hours'),
                        by = 'hours', all = T)
pheno_mismatch <- merge(pheno_mismatch,
                        plyr::count(rnd2_res[rnd2_res$phenotype != rnd2_res$mcat_phenotype,], vars = 'hours'),
                        by = 'hours', all = T)
pheno_mismatch <- merge(pheno_mismatch,
                        plyr::count(rnd2_res[rnd2_res$phenotype != rnd2_res$mcat_wilcox_phenotype,], vars = 'hours'),
                        by = 'hours', all = T)

colnames(pheno_mismatch) <- c('hours','lid_emp','lid_wilcox','mcat_emp','mcat_wilcox')
pheno_mismatch[is.na(pheno_mismatch)] <- 0

##### PLOT THE PHENOTYPE RESULTS
ggplot(melt(pheno_mismatch, id.vars = 'hours', variable.name = 'method', value.name = 'count'),
       aes(x = hours, y = count + 0.001)) +
  geom_line(aes(col = method)) +
  scale_x_continuous(breaks = seq(1,13,2)) +
  scale_y_continuous(breaks = seq(0,300,50)) +
  scale_color_discrete(name = 'Method') +
  labs(x = 'Time (hours)',
       y = 'Mismatch') +
  coord_cartesian(xlim = c(1,11)) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom')

##### FALSE POSITIVES & NEGATIVES
head(rnd2_res)
rnd2_res$truth[rnd2_res$phenotype != 'Neutral'] <- 'True Positive'
rnd2_res$truth[rnd2_res$phenotype == 'Neutral'] <- 'True Negative'

rnd2_res$lid_results[rnd2_res$phenotype == 'Beneficial' & rnd2_res$lid_phenotype == 'Beneficial'] <- 'True Positive'
rnd2_res$lid_results[rnd2_res$phenotype == 'Deleterious' & rnd2_res$lid_phenotype == 'Deleterious'] <- 'True Positive'
rnd2_res$lid_results[rnd2_res$phenotype == 'Neutral' & rnd2_res$lid_phenotype == 'Neutral'] <- 'True Negative'
rnd2_res$lid_results[rnd2_res$phenotype == 'Neutral' & rnd2_res$lid_phenotype != 'Neutral'] <- 'False Positive'
rnd2_res$lid_results[rnd2_res$phenotype != 'Neutral' & rnd2_res$lid_phenotype == 'Neutral'] <- 'False Negative'
rnd2_res$lid_results[rnd2_res$phenotype != 'Neutral' & rnd2_res$lid_phenotype != 'Neutral' &
                        rnd2_res$phenotype != rnd2_res$lid_phenotype] <- 'Switched Positive'

rnd2_res$lid_wilcox_results[rnd2_res$phenotype == 'Beneficial' & rnd2_res$lid_wilcox_phenotype == 'Beneficial'] <- 'True Positive'
rnd2_res$lid_wilcox_results[rnd2_res$phenotype == 'Deleterious' & rnd2_res$lid_wilcox_phenotype == 'Deleterious'] <- 'True Positive'
rnd2_res$lid_wilcox_results[rnd2_res$phenotype == 'Neutral' & rnd2_res$lid_wilcox_phenotype == 'Neutral'] <- 'True Negative'
rnd2_res$lid_wilcox_results[rnd2_res$phenotype == 'Neutral' & rnd2_res$lid_wilcox_phenotype != 'Neutral'] <- 'False Positive'
rnd2_res$lid_wilcox_results[rnd2_res$phenotype != 'Neutral' & rnd2_res$lid_wilcox_phenotype == 'Neutral'] <- 'False Negative'
rnd2_res$lid_wilcox_results[rnd2_res$phenotype != 'Neutral' & rnd2_res$lid_wilcox_phenotype != 'Neutral' &
                        rnd2_res$phenotype != rnd2_res$lid_wilcox_phenotype] <- 'Switched Positive'

rnd2_res$mcat_results[rnd2_res$phenotype == 'Beneficial' & rnd2_res$mcat_phenotype == 'Beneficial'] <- 'True Positive'
rnd2_res$mcat_results[rnd2_res$phenotype == 'Deleterious' & rnd2_res$mcat_phenotype == 'Deleterious'] <- 'True Positive'
rnd2_res$mcat_results[rnd2_res$phenotype == 'Neutral' & rnd2_res$mcat_phenotype == 'Neutral'] <- 'True Negative'
rnd2_res$mcat_results[rnd2_res$phenotype == 'Neutral' & rnd2_res$mcat_phenotype != 'Neutral'] <- 'False Positive'
rnd2_res$mcat_results[rnd2_res$phenotype != 'Neutral' & rnd2_res$mcat_phenotype == 'Neutral'] <- 'False Negative'
rnd2_res$mcat_results[rnd2_res$phenotype != 'Neutral' & rnd2_res$mcat_phenotype != 'Neutral' &
                        rnd2_res$phenotype != rnd2_res$mcat_phenotype] <- 'Switched Positive'

rnd2_res$mcat_wilcox_results[rnd2_res$phenotype == 'Beneficial' & rnd2_res$mcat_wilcox_phenotype == 'Beneficial'] <- 'True Positive'
rnd2_res$mcat_wilcox_results[rnd2_res$phenotype == 'Deleterious' & rnd2_res$mcat_wilcox_phenotype == 'Deleterious'] <- 'True Positive'
rnd2_res$mcat_wilcox_results[rnd2_res$phenotype == 'Neutral' & rnd2_res$mcat_wilcox_phenotype == 'Neutral'] <- 'True Negative'
rnd2_res$mcat_wilcox_results[rnd2_res$phenotype == 'Neutral' & rnd2_res$mcat_wilcox_phenotype != 'Neutral'] <- 'False Positive'
rnd2_res$mcat_wilcox_results[rnd2_res$phenotype != 'Neutral' & rnd2_res$mcat_wilcox_phenotype == 'Neutral'] <- 'False Negative'
rnd2_res$mcat_wilcox_results[rnd2_res$phenotype != 'Neutral' & rnd2_res$mcat_wilcox_phenotype != 'Neutral' &
                               rnd2_res$phenotype != rnd2_res$mcat_wilcox_phenotype] <- 'Switched Positive'

ggplot(rnd2_res,
       aes(x = hours)) +
  geom_histogram(aes(fill = mcat_results))

truth <- data.frame(plyr::count(rnd2_res, vars = c('hours','orf_name','truth'))[,c(1:3)], method = 'truth')
colnames(truth) <- c('hours','orf_name','results','method') 

lid_results <- data.frame(plyr::count(rnd2_res, vars = c('hours','orf_name','lid_results'))[,c(1:3)], method = 'lid_emp')
colnames(lid_results) <- c('hours','orf_name','results','method')           

lid_wilcox_results <- data.frame(plyr::count(rnd2_res, vars = c('hours','orf_name','lid_wilcox_results'))[,c(1:3)], method = 'lid_wilcox')
colnames(lid_wilcox_results) <- c('hours','orf_name','results','method') 

mcat_results <- data.frame(plyr::count(rnd2_res, vars = c('hours','orf_name','mcat_results'))[,c(1:3)], method = 'mcat_emp')
colnames(mcat_results) <- c('hours','orf_name','results','method') 

mcat_wilcox_results <- data.frame(plyr::count(rnd2_res, vars = c('hours','orf_name','mcat_wilcox_results'))[,c(1:3)], method = 'mcat_wilcox')
colnames(mcat_wilcox_results) <- c('hours','orf_name','results','method') 

all_results <- rbind(truth, lid_results, lid_wilcox_results, mcat_results, mcat_wilcox_results)
head(all_results)

ggplot(all_results,
       aes(x = hours)) +
  geom_bar(aes(fill = results)) +
  facet_wrap(.~method) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom')

ggplot(all_results,
       aes(x = method)) +
  geom_bar(aes(fill = results)) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom')


##### FALSE POSITIVES & NEGATIVES SUMMARY
all_results_sum <- plyr::count(all_results, vars = c('method','results'))
all_results_sum$count <- all_results_sum$freq 
all_results_sum$freq <- all_results_sum$freq/(9.10*11)
all_results_sum$results <- factor(all_results_sum$results,
                              levels = c('True Negative','True Positive','False Negative','False Positive'))
all_results_sum$method <- as.character(all_results_sum$method)
all_results_sum$method[all_results_sum$method == 'truth'] <- 'Truth'
all_results_sum$method[all_results_sum$method == 'lid_emp'] <- 'LID + Empirical Test'
all_results_sum$method[all_results_sum$method == 'lid_wilcox'] <- 'LID + Wilcoxon Ranksum'
all_results_sum$method[all_results_sum$method == 'mcat_emp'] <- 'MCAT + Empirical Test'
all_results_sum$method[all_results_sum$method == 'mcat_wilcox'] <- 'MCAT + Wilcoxon Ranksum'

head(all_results_sum)

truth.pie <- ggplot(all_results_sum[all_results_sum$method == 'Truth',], aes(x = "", y = freq, fill = results)) +
  geom_bar(stat = 'identity', col = 'black') +
  coord_polar("y", start = 0) +
  geom_label_repel(aes(label = sprintf("%0.2f%%",freq)),
                   position = position_stack(vjust = 0.5),
                   colour ='white',
                   label.size = 0.15,
                   show.legend = F) +
  scale_fill_discrete(name = 'Result', drop = F) +
  facet_wrap(.~method) +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(size = titles, hjust = 0.5),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = 'bottom',
        strip.text = element_text(size = titles),
        legend.box.spacing = unit(0.5,"mm"))

test.pie <- ggplot(all_results_sum[all_results_sum$method != 'Truth',], aes(x = "", y = freq, fill = results)) +
  geom_bar(stat = 'identity', col = 'black') +
  coord_polar("y", start = 0) +
  geom_label_repel(aes(label = sprintf("%0.2f%%",freq)),
                   position = position_stack(vjust = 0.5),
                   colour ='white',
                   label.size = 0.15,
                   show.legend = F) +
  scale_fill_discrete(name = 'Result', drop = F) +
  facet_wrap(.~method, ncol = 2) +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(size = titles, hjust = 0.5),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = 'bottom',
        strip.text = element_text(size = titles),
        legend.box.spacing = unit(0.5,"mm"))

perf.pie <- ggpubr::ggarrange(truth.pie, test.pie, nrow = 1,
                  widths = c(1,2), common.legend = T, legend = 'bottom')
ggsave(sprintf("%sr2r/FIGURE_RND_PERFORMACE.jpg",out_path), perf.pie,
       height = 120, width = two.c, units = 'mm',
       dpi = 300)


all_results_sum$sensitivity[all_results_sum$method == 'LID + Empirical Test'] <- all_results_sum$freq[all_results_sum$method == 'LID + Empirical Test' & all_results_sum$results == 'True Positive']/
  sum(all_results_sum$freq[all_results_sum$method == 'LID + Empirical Test' & all_results_sum$results %in% c('True Positive', 'False Negative')]) * 100
all_results_sum$specificity[all_results_sum$method == 'LID + Empirical Test'] <- all_results_sum$freq[all_results_sum$method == 'LID + Empirical Test' & all_results_sum$results == 'True Negative']/
  sum(all_results_sum$freq[all_results_sum$method == 'LID + Empirical Test' & all_results_sum$results %in% c('True Negative', 'False Positive')]) * 100
all_results_sum$ppv[all_results_sum$method == 'LID + Empirical Test'] <- all_results_sum$freq[all_results_sum$method == 'LID + Empirical Test' & all_results_sum$results == 'True Positive']/
  sum(all_results_sum$freq[all_results_sum$method == 'LID + Empirical Test' & all_results_sum$results %in% c('True Positive', 'False Positive')]) * 100
all_results_sum$npv[all_results_sum$method == 'LID + Empirical Test'] <- all_results_sum$freq[all_results_sum$method == 'LID + Empirical Test' & all_results_sum$results == 'True Negative']/
  sum(all_results_sum$freq[all_results_sum$method == 'LID + Empirical Test' & all_results_sum$results %in% c('True Negative', 'False Negative')]) * 100

all_results_sum$sensitivity[all_results_sum$method == 'LID + Wilcoxon Ranksum'] <- all_results_sum$freq[all_results_sum$method == 'LID + Wilcoxon Ranksum' & all_results_sum$results == 'True Positive']/
  sum(all_results_sum$freq[all_results_sum$method == 'LID + Wilcoxon Ranksum' & all_results_sum$results %in% c('True Positive', 'False Negative')]) * 100
all_results_sum$specificity[all_results_sum$method == 'LID + Wilcoxon Ranksum'] <- all_results_sum$freq[all_results_sum$method == 'LID + Wilcoxon Ranksum' & all_results_sum$results == 'True Negative']/
  sum(all_results_sum$freq[all_results_sum$method == 'LID + Wilcoxon Ranksum' & all_results_sum$results %in% c('True Negative', 'False Positive')]) * 100
all_results_sum$ppv[all_results_sum$method == 'LID + Wilcoxon Ranksum'] <- all_results_sum$freq[all_results_sum$method == 'LID + Wilcoxon Ranksum' & all_results_sum$results == 'True Positive']/
  sum(all_results_sum$freq[all_results_sum$method == 'LID + Wilcoxon Ranksum' & all_results_sum$results %in% c('True Positive', 'False Positive')]) * 100
all_results_sum$npv[all_results_sum$method == 'LID + Wilcoxon Ranksum'] <- all_results_sum$freq[all_results_sum$method == 'LID + Wilcoxon Ranksum' & all_results_sum$results == 'True Negative']/
  sum(all_results_sum$freq[all_results_sum$method == 'LID + Wilcoxon Ranksum' & all_results_sum$results %in% c('True Negative', 'False Negative')]) * 100

all_results_sum$sensitivity[all_results_sum$method == 'MCAT + Empirical Test'] <- all_results_sum$freq[all_results_sum$method == 'MCAT + Empirical Test' & all_results_sum$results == 'True Positive']/
  sum(all_results_sum$freq[all_results_sum$method == 'MCAT + Empirical Test' & all_results_sum$results %in% c('True Positive', 'False Negative')]) * 100
all_results_sum$specificity[all_results_sum$method == 'MCAT + Empirical Test']  <- all_results_sum$freq[all_results_sum$method == 'MCAT + Empirical Test' & all_results_sum$results == 'True Negative']/
  sum(all_results_sum$freq[all_results_sum$method == 'MCAT + Empirical Test' & all_results_sum$results %in% c('True Negative', 'False Positive')]) * 100
all_results_sum$ppv[all_results_sum$method == 'MCAT + Empirical Test']  <- all_results_sum$freq[all_results_sum$method == 'MCAT + Empirical Test' & all_results_sum$results == 'True Positive']/
  sum(all_results_sum$freq[all_results_sum$method == 'MCAT + Empirical Test' & all_results_sum$results %in% c('True Positive', 'False Positive')]) * 100
all_results_sum$npv[all_results_sum$method == 'MCAT + Empirical Test']  <- all_results_sum$freq[all_results_sum$method == 'MCAT + Empirical Test' & all_results_sum$results == 'True Negative']/
  sum(all_results_sum$freq[all_results_sum$method == 'MCAT + Empirical Test' & all_results_sum$results %in% c('True Negative', 'False Negative')]) * 100

all_results_sum$sensitivity[all_results_sum$method == 'MCAT + Wilcoxon Ranksum'] <- all_results_sum$freq[all_results_sum$method == 'MCAT + Wilcoxon Ranksum' & all_results_sum$results == 'True Positive']/
  sum(all_results_sum$freq[all_results_sum$method == 'MCAT + Wilcoxon Ranksum' & all_results_sum$results %in% c('True Positive', 'False Negative')]) * 100
all_results_sum$specificity[all_results_sum$method == 'MCAT + Wilcoxon Ranksum'] <- all_results_sum$freq[all_results_sum$method == 'MCAT + Wilcoxon Ranksum' & all_results_sum$results == 'True Negative']/
  sum(all_results_sum$freq[all_results_sum$method == 'MCAT + Wilcoxon Ranksum' & all_results_sum$results %in% c('True Negative', 'False Positive')]) * 100
all_results_sum$ppv[all_results_sum$method == 'MCAT + Wilcoxon Ranksum'] <- all_results_sum$freq[all_results_sum$method == 'MCAT + Wilcoxon Ranksum' & all_results_sum$results == 'True Positive']/
  sum(all_results_sum$freq[all_results_sum$method == 'MCAT + Wilcoxon Ranksum' & all_results_sum$results %in% c('True Positive', 'False Positive')]) * 100
all_results_sum$npv[all_results_sum$method == 'MCAT + Wilcoxon Ranksum'] <- all_results_sum$freq[all_results_sum$method == 'MCAT + Wilcoxon Ranksum' & all_results_sum$results == 'True Negative']/
  sum(all_results_sum$freq[all_results_sum$method == 'MCAT + Wilcoxon Ranksum' & all_results_sum$results %in% c('True Negative', 'False Negative')]) * 100

perf.points <- ggplot(melt(all_results_sum[all_results_sum$method != 'Truth',c(1,4,5)], id.vars = 'method'),
       aes(x = variable, y = value)) +
  geom_line(aes(group = method)) +
  geom_point(aes(col = method)) +
  labs(x = 'Variable', y = 'Value (%)') +
  scale_color_discrete(name = 'Method') +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = 'bottom',
        strip.text = element_text(size = titles),
        legend.box.spacing = unit(0.5,"mm"))
ggsave(sprintf("%sr2r/FIGURE_RND_PERFORMACE2.jpg",out_path), perf.points,
       height = one.5c, width = one.5c, units = 'mm',
       dpi = 300)


##### NUMBER OF REPLICATES
outlier_reps <- dbGetQuery(conn, 'select * from 4C4_FS_CC_6144_FITNESS_STATS')
outlier.histo <- ggplot(outlier_reps[!(outlier_reps$orf_name %in% c('BFC100','BF_control')) &
                      outlier_reps$hours > 0,],
       aes(x = N, y = ..count../910*100)) +
  geom_histogram(binwidth = 0.5) +
  coord_cartesian(ylim = c(0,100)) +
  labs(x = 'Replicates Retained',
       y = 'Proportion of Strains (%)') +
  # scale_y_log10() +
  facet_wrap(.~hours) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%sr2r/FIGURE_OUTLIERS_HISTOGRAM.jpg",out_path), outlier.histo,
       height = one.5c, width = one.5c, units = 'mm',
       dpi = 300)

outlier_data <- dbGetQuery(conn, 'select a.*, b.plate, b.col, b.row
                           from 4C4_FS_CC_6144_FITNESS a, 4C4_pos2coor b
                           where a.pos = b.pos
                           and hours > 0')
outlier_data$orf_name[is.na(outlier_data$orf_name)] <- 'Gap'
outlier_data$fitness[outlier_data$orf_name %in% c('BFC100','BF_control')] <- NA

for (h in sort(unique(outlier_data$hours))) {
  for (o in sort(unique(outlier_data$orf_name[outlier_data$hours == h]))) {
    temp <- outlier_data$fitness[outlier_data$hours == h & outlier_data$orf_name == o]
    m <- median(temp, na.rm = T)
    madev <- mad(temp, na.rm = T)
    ul <- m + 3*madev
    ll <- m - 3*madev
    
    outlier_data$outlier[outlier_data$hours == h & outlier_data$orf_name == o & outlier_data$fitness > ul] <- TRUE
    outlier_data$outlier[outlier_data$hours == h & outlier_data$orf_name == o & outlier_data$fitness < ll] <- TRUE
    
    outlier_data$outlier_count[outlier_data$hours == h & outlier_data$orf_name == o] <- 
      sum(outlier_data$outlier[outlier_data$hours == h & outlier_data$orf_name == o], na.rm = T)
  }
}
outlier_data$outlier[is.na(outlier_data$outlier)] <- FALSE
head(outlier_data)

outlier.violin <- ggplot() +
  geom_violin(data = outlier_data[outlier_data$orf_name %in% unique(outlier_data$orf_name[outlier_data$outlier &
                                                                                            outlier_data$hours == 11.04]) &
                                    !(outlier_data$outlier) &
                                    outlier_data$hours == 11.04,],
              aes(x = orf_name, y = fitness)) +
  geom_jitter(data = outlier_data[outlier_data$orf_name %in% unique(outlier_data$orf_name[outlier_data$outlier &
                                                                                            outlier_data$hours == 11.04]) &
                                    outlier_data$outlier &
                                    outlier_data$hours == 11.04,],
              aes(x = orf_name, y = fitness), col = 'red', size = 0.8,
              width = 0) +
  facet_wrap(.~hours, scales = 'free', ncol = 1) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 90),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%sr2r/FIGURE_OUTLIERS_EXAMPLE.jpg",out_path), outlier.violin,
       height = 150, width = 1000, units = 'mm',
       dpi = 300)

head(outlier_data)
outlier.map <- ggplot() +
  geom_tile(data = outlier_data[outlier_data$hours == 11.04,],
            aes(x = col, y = row), fill = 'white', col = 'black') +
  geom_tile(data = plyr::count(outlier_data[outlier_data$outlier & outlier_data$hours > 7.5,], vars = c('plate','col','row')),
            aes(x = col, y = row, fill = freq), col = 'black') +
  geom_tile(data = outlier_data[outlier_data$orf_name == 'Gap',],
            aes(x = col, y = row), fill = 'black', col = 'black') +
  scale_fill_distiller(name = 'Outlier\nFreq.',
                       palette = "Spectral") +
  scale_y_continuous(trans = 'reverse') +
  facet_wrap(.~plate) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%sr2r/FIGURE_OUTLIERS_MAP.jpg",out_path), outlier.map,
       height = 150, width = two.c, units = 'mm',
       dpi = 300)

##### GROWTH CURVES
outlier_data$colony[outlier_data$orf_name == 'BF_control'] <- 'Reference'
outlier_data$colony[outlier_data$orf_name != 'BF_control'] <- 'Query'
outlier_data$colony[outlier_data$orf_name == 'Gap'] <- 'Gap'

gc_data <- data.frame(outlier_data %>%
  group_by(hours,colony) %>%
  summarise(cs = mean(average, na.rm = T)))

gc.plot <- ggplot(gc_data[gc_data$colony != 'Gap',]) +
  geom_line(aes(x = hours, y = cs, col = colony)) +
  geom_point(aes(x = hours, y = cs, col = colony)) +
  labs(x = 'Time (hours)',
       y = 'log10 Pixel Count') +
  # scale_color_discrete(name = '') +
  scale_y_log10() +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_blank(),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%sr2r/FIGURE_GROWTHCURVES.jpg",out_path), gc.plot,
       height = one.c, width = one.c, units = 'mm',
       dpi = 300)


##### MODIFIED MAP
##### for empirical testing using a different population

map.384 <- dbGetQuery(conn, 'select a.*, b.orf_name
                      from 4C4_pos2coor a, 4C4_pos2orf_name b
                      where a.density = 384 and a.pos = b.pos')
map.384$orf_name[map.384$plate == 1 & !is.na(map.384$orf_name)] <- "REF"

pvalue_4C4 <- dbGetQuery(conn, 'select * from 4C4_FS_CC_6144_PVALUE where hours > 0')
fitness_4C4 <- dbGetQuery(conn, 'select * from 4C4_FS_CC_6144_FITNESS where hours > 0')

fitness_4C4$source <- as.numeric(substr(as.character(fitness_4C4$pos), 3, 6))
fitness_4C4 <- merge(fitness_4C4, map.384, by.x = 'source', by.y = 'pos', suffixes = c('','_new'), all = T)
head(fitness_4C4)

# make ref fitness and bf_control fitness distribution and compare
ggplot(fitness_4C4[fitness_4C4$orf_name_new %in% c('REF','BF_control'),]) +
  geom_density_ridges(aes(x = fitness, y = orf_name_new, group = orf_name_new, fill = orf_name_new),
                      quantile_lines = TRUE,
                      scale = 3, alpha = 0.8, size = 0.3,
                      vline_size = 0.2, vline_color = "black") +
  facet_wrap(.~hours, scale = 'free')
# remove outliers and repeat

# use ref distribution to give pvalues and compare results with pvalue_4C4
