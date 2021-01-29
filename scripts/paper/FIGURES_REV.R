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

load("~/R/Projects/adaptivefitness/output/workspace/r2r.RData")

source("R/functions/isoutlier.R")
source("R/functions/empirical_p.R")
source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")
out_path = 'figs/paper/';

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 7
txt <- 7
lbls <- 9

##### REVIEWER #1

##### l. 101: Expt pipeline reducing source effect
##### FIGURE S1. Prescreens Plate Effect
load(sprintf('%s4C4PSDATA.RData',out_path))

ps.sourceeffect <- ggplot(data.ps,
                         aes(x = source, y = average-median_cs)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0, col = 'red', linetype = 'dashed') +
  labs(x = 'Source Plate',
       y = 'Plate Median Subtracted Pixel Counts') +
  scale_x_discrete(breaks=c("1TL","2TR","3BL","4BR"),
                   labels=c("Top Left","Top Right","Bottom Left","Bottom Right")) +
  # stat_compare_means(label.x = 2.5, label.y = 700, hjust = 0.5, size = 1.5) +
  facet_wrap(.~stage) +
  theme_linedraw()+
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = txt,
                                  margin = margin(2,0,2,0, "mm"))) +
  coord_cartesian(ylim = c(-500,500))
ggsave(sprintf("%sr2r/FIGURE_SOURCE_EFFECT.jpg",out_path), ps.sourceeffect,
       height = one.c, width = one.c, units = 'mm',
       dpi = 300)

##### l. 305: Using some other statistical testing in the virtual plates with random CS distribution
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

##### CONSOLIDATING EMPIRICAL TEST AND WILCOXON
rnd2_res <- NULL
i <- 1
for (h in sort(unique(rnd2_data$hours))) {
  rnd2_data$lid_pthresh[rnd2_data$hours == h] <- unique(lid.rnd2.res$pthresh[lid.rnd2.res$hours == h & lid.rnd2.res$rep == 16 & lid.rnd2.res$ref == 0.25])
  rnd2_data$mcat_pthresh[rnd2_data$hours == h] <- unique(mcat.rnd2.res$pthresh[mcat.rnd2.res$hours == h])
  for (o in sort(unique(rnd2_data$orf_name[rnd2_data$hours == h & rnd2_data$orf_name != 'BF_control']))){
    rnd2_res$hours[i] <- h
    rnd2_res$orf_name[i] <- o
    rnd2_res$rnd_hrs[i] <- unique(rnd2_data$rnd_hrs[rnd2_data$hours == h & rnd2_data$orf_name == o])
    ## LID
    # empirical test results
    rnd2_res$lid_p[i] <- unique(lid.rnd2.res$p[lid.rnd2.res$hours == h & lid.rnd2.res$rep == 16 & lid.rnd2.res$ref == 0.25 & lid.rnd2.res$orf_name == o])
    rnd2_res$lid_stat[i] <- unique(rnd2_data$lid_stat[rnd2_data$hours == h & rnd2_data$orf_name == o])
    rnd2_res$lid_es[i] <- unique(lid.rnd2.res$es[lid.rnd2.res$hours == h & lid.rnd2.res$rep == 16 & lid.rnd2.res$ref == 0.25 & lid.rnd2.res$orf_name == o])
    rnd2_res$lid_pthresh[i] <- unique(lid.rnd2.res$pthresh[lid.rnd2.res$hours == h & lid.rnd2.res$rep == 16 & lid.rnd2.res$ref == 0.25])
    # wilcoxon ranksum test
    ref_fit <- rnd2_data$lid_fitness[rnd2_data$hours == h & rnd2_data$orf_name == 'BF_control' & !is.na(rnd2_data$lid_fitness)]
    ref_fit <- ref_fit[!isoutlier(ref_fit)]
    orf_fit <- rnd2_data$lid_fitness[rnd2_data$hours == h & rnd2_data$orf_name == o & !is.na(rnd2_data$lid_fitness)]
    orf_fit <- orf_fit[!isoutlier(orf_fit)]
    wilcox_stats <- wilcox.test(orf_fit, ref_fit)
    rnd2_res$lid_wilcox_p[i] <- wilcox_stats$p.value
  
    ## MCAT
    # empirical test results
    rnd2_res$mcat_p[i] <- unique(rnd2_data$mcat_p[rnd2_data$hours == h & rnd2_data$orf_name == o])
    rnd2_res$mcat_stat[i] <- unique(rnd2_data$mcat_stat[rnd2_data$hours == h & rnd2_data$orf_name == o])
    rnd2_res$mcat_es[i] <- unique(mcat.rnd2.res$es[mcat.rnd2.res$hours == h & mcat.rnd2.res$orf_name == o])
    rnd2_res$mcat_pthresh[i] <- unique(mcat.rnd2.res$pthresh[mcat.rnd2.res$hours == h])
    # wilcoxon ranksum test
    ref_fit <- rnd2_data$mcat_fitness[rnd2_data$hours == h & rnd2_data$orf_name == 'BF_control' & !is.na(rnd2_data$mcat_fitness)]
    ref_fit <- ref_fit[!isoutlier(ref_fit)]
    orf_fit <- rnd2_data$mcat_fitness[rnd2_data$hours == h & rnd2_data$orf_name == o & !is.na(rnd2_data$mcat_fitness)]
    orf_fit <- orf_fit[!isoutlier(orf_fit)]
    wilcox_stats <- wilcox.test(orf_fit, ref_fit)
    rnd2_res$mcat_wilcox_p[i] <- wilcox_stats$p.value
    
    i <- i + 1
  }
  # multiple hypothesis correction
  rnd2_res$lid_q[rnd2_res$hours == h] <- qvalue(rnd2_res$lid_wilcox_p[rnd2_res$hours == h], fdr.level=0.05, pi0 = 1)$qvalue
  rnd2_res$mcat_q[rnd2_res$hours == h] <- qvalue(rnd2_res$mcat_wilcox_p[rnd2_res$hours == h], fdr.level=0.05, pi0 = 1)$qvalue
}
rnd2_res <- data.frame(rnd2_res)
head(rnd2_res)

# true phenotype
rnd2_res$phenotype[rnd2_res$hours < rnd2_res$rnd_hrs] <- 'Beneficial'
rnd2_res$phenotype[rnd2_res$hours > rnd2_res$rnd_hrs] <- 'Deleterious'
rnd2_res$phenotype[rnd2_res$hours == rnd2_res$rnd_hrs] <- 'Neutral'
# lid empirical test phenotype
rnd2_res$lid_phenotype[rnd2_res$lid_p <= rnd2_res$lid_pthresh & rnd2_res$lid_es > 1] <- 'Beneficial'
rnd2_res$lid_phenotype[rnd2_res$lid_p <= rnd2_res$lid_pthresh & rnd2_res$lid_es < 1] <- 'Deleterious'
rnd2_res$lid_phenotype[is.na(rnd2_res$lid_phenotype)] <- 'Neutral'
# lid wilcoxon test phenotype
rnd2_res$lid_wilcox_phenotype[rnd2_res$lid_wilcox_p <= 0.05 & rnd2_res$lid_q <= 0.05 & rnd2_res$lid_es > 1] <- 'Beneficial'
rnd2_res$lid_wilcox_phenotype[rnd2_res$lid_wilcox_p <= 0.05 & rnd2_res$lid_q <= 0.05 & rnd2_res$lid_es < 1] <- 'Deleterious'
rnd2_res$lid_wilcox_phenotype[is.na(rnd2_res$lid_wilcox_phenotype)] <- 'Neutral'
# mcat empirical test phenotype
rnd2_res$mcat_phenotype[rnd2_res$mcat_p <= rnd2_res$mcat_pthresh & rnd2_res$mcat_es > 1] <- 'Beneficial'
rnd2_res$mcat_phenotype[rnd2_res$mcat_p <= rnd2_res$mcat_pthresh & rnd2_res$mcat_es < 1] <- 'Deleterious'
rnd2_res$mcat_phenotype[is.na(rnd2_res$mcat_phenotype)] <- 'Neutral'
# mcat wilcoxon test phenotype
rnd2_res$mcat_wilcox_phenotype[rnd2_res$mcat_q <= 0.05 & rnd2_res$mcat_es > 1] <- 'Beneficial'
rnd2_res$mcat_wilcox_phenotype[rnd2_res$mcat_q <= 0.05 & rnd2_res$mcat_es < 1] <- 'Deleterious'
rnd2_res$mcat_wilcox_phenotype[is.na(rnd2_res$mcat_wilcox_phenotype)] <- 'Neutral'

# Internal QC
# plyr::count(lid.rnd2.res[lid.rnd2.res$rep == 16 & lid.rnd2.res$ref == 0.25 & lid.rnd2.res$attempt == 1,],
#             vars = c('hours','effect'))
# plyr::count(rnd2_res, vars = c('hours','lid_phenotype'))
# 
# hello <- merge(lid.rnd2.res[lid.rnd2.res$rep == 16 & lid.rnd2.res$ref == 0.25 & lid.rnd2.res$attempt == 1,],
#       rnd2_res, by = c('hours','orf_name'), suffix = c('old','new'))
# 
# hi <- hello[hello$effect != hello$lid_phenotype,]

###### COMPILING ALL PHENOTYPE RESULTS
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

# Plotting the mismatch counts over time
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

# FALSE POSITIVES & NEGATIVES
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

# all results
all_results <- rbind(truth, lid_results, lid_wilcox_results, mcat_results, mcat_wilcox_results)
head(all_results)

# ggplot(all_results,
#        aes(x = hours)) +
#   geom_bar(aes(fill = results)) +
#   facet_wrap(.~method) +
#   theme_linedraw() +
#   theme(plot.title = element_text(size = titles),
#         axis.title = element_text(size = titles),
#         axis.text = element_text(size = txt),
#         legend.title = element_text(size = titles),
#         legend.text = element_text(size = txt),
#         legend.position = 'bottom')
# 
# ggplot(all_results,
#        aes(x = method)) +
#   geom_bar(aes(fill = results)) +
#   theme_linedraw() +
#   theme(plot.title = element_text(size = titles),
#         axis.title = element_text(size = titles),
#         axis.text = element_text(size = txt),
#         legend.title = element_text(size = titles),
#         legend.text = element_text(size = txt),
#         legend.position = 'bottom')

# summary of the performance results
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

# PLOTTING THE PERFORMANCE SUMMARY
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
                  widths = c(1,2), common.legend = T, legend = 'bottom',
                  labels = c('A','B'),
                  font.label = list(face = 'bold', size = lbls, family = "sans"),
                  hjust=-1)
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

perf.points <- ggplot(melt(all_results_sum[all_results_sum$method != 'Truth',c(1,5,6)], id.vars = 'method'),
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

perf.all <- ggpubr::ggarrange(perf.pie, perf.points, ncol = 1,
                              heights = c(2,1),
                              labels = c('','C'),
                              font.label = list(face = 'bold', size = lbls, family = "sans"),
                              hjust=-1)
ggsave(sprintf("%sr2r/FIGURE_RND_PERFORMACE3.jpg",out_path), perf.all,
       height = two.c, width = one.5c, units = 'mm',
       dpi = 300)

##### l. 413 colonies around gaps
gap_data <- dbGetQuery(conn, 'select c.*, d.p, d.stat from 
                       (select a.*, b.plate, b.row, b.col
                       from 4C4_FS_CC_6144_FITNESS a, 4C4_pos2coor b
                       where a.pos = b.pos and a.hours = 11.04) c
                       left join 4C4_FS_CC_6144_PVALUE d
                       on c.orf_name = d.orf_name and c.hours = d.hours
                       order by plate, col, row')
gap_data$phenotype[gap_data$p <= 0.05 & gap_data$stat > 0] <- 'Beneficial'
gap_data$phenotype[gap_data$p <= 0.05 & gap_data$stat < 0] <- 'Deleterious'
gap_data$phenotype[is.na(gap_data$phenotype) &
                     !is.na(gap_data$orf_name) &
                     gap_data$orf_name != 'BF_control' &
                     !is.na(gap_data$fitness)] <- 'Neutral'
gap_data$phenotype[is.na(gap_data$orf_name)] <- 'Gap'
gap_data$phenotype[gap_data$orf_name == 'BF_control'] <- 'BF_control'
# head(gap_data)

ggplot(gap_data[!is.na(gap_data$plate),],
       aes(x = col, y = row)) +
  geom_tile(aes(fill = phenotype)) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Beneficial','Neutral','Deleterious','Gap','BF_control'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107',
                               'Gap'='Red',
                               'BF_control'='grey60')) +
  facet_wrap(.~plate)

##### l. 539: Is outlier removal sound, and what ends up getting removed
##### NUMBER OF REPLICATES AFTER OUTLIER REMOVAL
outlier_reps <- dbGetQuery(conn, 'select * from 4C4_FS_CC_6144_FITNESS_STATS')
outlier.histo <- ggplot(outlier_reps[!(outlier_reps$orf_name %in% c('BFC100','BF_control')) &
                      outlier_reps$hours > 0,],
       aes(x = N, y = ..count../910*100)) +
  geom_histogram(binwidth = 0.5) +
  coord_cartesian(ylim = c(0,100)) +
  labs(x = 'Replicates Retained',
       y = 'Proportion of Mutants (%)') +
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

outlier_reps[outlier_reps$N <= 14 & outlier_reps$hours > 7,] %>%
  # group_by(hours) %>%
  count()/(911*4) * 100

##### WHAT ARE THE OUTLIER VALUES
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

# example
outlier.violin <- ggplot() +
  geom_violin(data = outlier_data[outlier_data$orf_name %in% unique(outlier_data$orf_name[outlier_data$outlier &
                                                                                            outlier_data$hours == 11.04 &
                                                                                            outlier_data$outlier_count > 2]) &
                                    !(outlier_data$outlier) &
                                    outlier_data$hours == 11.04,],
              aes(x = orf_name, y = fitness), draw_quantiles = c(0.25,0.5,0.75)) +
  geom_jitter(data = outlier_data[outlier_data$orf_name %in% unique(outlier_data$orf_name[outlier_data$outlier &
                                                                                            outlier_data$hours == 11.04 &
                                                                                            outlier_data$outlier_count > 2]) &
                                    outlier_data$outlier &
                                    outlier_data$hours == 11.04,],
              aes(x = orf_name, y = fitness), col = 'red', size = 0.8,
              width = 0) +
  facet_wrap(.~hours, scales = 'free', ncol = 1) +
  theme_linedraw() +
  labs(x = 'Mutant',
       y = 'Fitness') +
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
       height = one.5c, width = two.c, units = 'mm',
       dpi = 300)

##### WHERE ON THE PLATE DO THE OUTLIERS LIE
## combined results for around saturation time points
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

##### l. 196: ACTUAL GROWTH CURVES OVER THE 11 TIME POINTS
##### GROWTH CURVES
# pooling all references and mutants together
outlier_data$colony[outlier_data$orf_name == 'BF_control'] <- 'Reference'
outlier_data$colony[outlier_data$orf_name != 'BF_control'] <- 'Mutant'
outlier_data$colony[outlier_data$orf_name == 'Gap'] <- 'Gap'

gc_data <- data.frame(outlier_data %>%
  group_by(hours,colony) %>%
  summarise(cs = mean(average, na.rm = T)))

gc.plot <- ggplot(gc_data[gc_data$colony != 'Gap',]) +
  geom_line(aes(x = hours, y = cs)) +
  geom_point(aes(x = hours, y = cs)) +
  labs(x = 'Time (hours)',
       y = 'log10 mean Pixel Count') +
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


##### MAJOR COMMENT
##### Additional reference source plate suggestion for creating the null distribution

##### Using the condition negative dataset
pvalue_4C4 <- dbGetQuery(conn, 'select * from 4C4_FS_CC_6144_PVALUE where hours > 0')
fitness_4C4 <- dbGetQuery(conn, 'select * from 4C4_FS_CC_6144_FITNESS where hours > 0')
fitness_4C4$source <- as.numeric(substr(as.character(fitness_4C4$pos), 3, 6))
##### Modified map for empirical testing using a different population
map.384 <- dbGetQuery(conn, 'select a.*, b.orf_name
                      from 4C4_pos2coor a, 4C4_pos2orf_name b
                      where a.density = 384 and a.pos = b.pos')
# using plate 2 as the reference plate for empirical testing
map.384$orf_name[map.384$plate == 2 & !is.na(map.384$orf_name)] <- "REF"
fitness_4C4 <- merge(fitness_4C4, map.384, by.x = 'source', by.y = 'pos', suffixes = c('','_new'), all = T)
head(fitness_4C4)

# remove outliers
fitness_4C4$outlier <- FALSE
for (h in sort(unique(fitness_4C4$hours))) {
  for (p in sort(unique(fitness_4C4$source[fitness_4C4$hours == h]))) {
    temp <- fitness_4C4$fitness[fitness_4C4$hours == h & fitness_4C4$source == p]
    m <- median(temp, na.rm = T)
    if (!is.na(m)) {
      madev <- mad(temp, na.rm = T)
      ul <- m + 3*madev
      ll <- m - 3*madev
      
      fitness_4C4$outlier[fitness_4C4$hours == h & fitness_4C4$source == p & fitness_4C4$fitness > ul] <- TRUE
      fitness_4C4$outlier[fitness_4C4$hours == h & fitness_4C4$source == p & fitness_4C4$fitness < ll] <- TRUE
      # fitness_4C4$outlier[fitness_4C4$hours == h & fitness_4C4$source == p & is.na(fitness_4C4$outlier)] <- FALSE
    }
    temp2 <- fitness_4C4$fitness[fitness_4C4$hours == h & fitness_4C4$source == p & !fitness_4C4$outlier]
    fitness_4C4$cs_mean[fitness_4C4$hours == h & fitness_4C4$source == p] <- mean(temp2, na.rm = T)
    fitness_4C4$cs_median[fitness_4C4$hours == h & fitness_4C4$source == p] <- median(temp2, na.rm = T)
    fitness_4C4$cs_std[fitness_4C4$hours == h & fitness_4C4$source == p] <- sd(temp2, na.rm = T)
  }
}

# compare REF and BF_control distribution
# BF_control used for spatial bias correction
# REF will be used for empirical testing
fitness_4C4$title <- sprintf('Time = %0.1f hours', fitness_4C4$hours)
fitness_4C4$title <- factor(fitness_4C4$title, levels = sprintf('Time = %0.1f hours', sort(unique(fitness_4C4$hours))))
new_ref.plot <- ggplot(fitness_4C4[fitness_4C4$orf_name_new %in% c('REF','BF_control') &
                     !fitness_4C4$outlier &
                     !is.na(fitness_4C4$fitness),],
       aes(x = fitness, y = orf_name_new, group = orf_name_new, fill = orf_name_new)) +
  geom_density_ridges(quantile_lines = TRUE,
                      quantiles = c(0.25, 0.5, 0.75),
                      scale = 3, alpha = 0.8, size = 0.3,
                      vline_size = 0.2, vline_color = "black",
                      na.rm = T) +
  labs(x = 'Fitness',
       y = 'Strain') +
  scale_y_discrete(breaks = c('REF','BF_control'),
                   labels = c('REF' = 'Tester', 'BF_control' = 'Reference')) +
  scale_fill_discrete(name = 'Strain',
                      breaks = c('REF','BF_control'),
                      labels = c('REF' = 'Tester', 'BF_control' = 'Reference'),
                      guide = F) +
  facet_wrap(.~title, scale = 'free_x') +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.y = element_text(angle = -45, vjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%sr2r/FIGURE_NEW_REFS.jpg",out_path), new_ref.plot,
       height = two.c, width = two.c, units = 'mm',
       dpi = 300)

# ggplot(fitness_4C4[fitness_4C4$orf_name_new %in% c('REF','BF_control') &
#                      !fitness_4C4$outlier &
#                      !is.na(fitness_4C4$fitness),],
#        aes(x = orf_name_new, y = fitness)) +
#   geom_violin() +
#   facet_wrap(.~title, scales = 'free')

# use REF distribution for empirical testing and compare results with those from using BF_control
fitness_stat_4C4 <- data.frame(fitness_4C4[!is.na(fitness_4C4$fitness),] %>%
  group_by(hours, source, orf_name_new) %>%
  summarise(cs_mean = mean(cs_mean, na.rm = T),
            cs_median = mean(cs_median, na.rm = T),
            cs_std = mean(cs_std, na.rm = T)))

colnames(fitness_stat_4C4) <- c("hours", "source", "orf_name", "cs_mean", "cs_median", "cs_std")

fitness_stat_4C4 <- empirical_p(fitness_stat_4C4, 'REF', 'cs_mean', 'p1', 'es1', c('REF','BF_control'), 0.05)
fitness_stat_4C4 <- empirical_p(fitness_stat_4C4, 'BF_control', 'cs_mean', 'p2', 'es2', c('REF','BF_control'), 0.05)
head(fitness_stat_4C4)

# ggplot(fitness_stat_4C4,
#        aes(x = p1, y = p2)) +
#   geom_point() +
#   stat_cor(method = 'spearman')

##### CONTROLING FOR 5% FPR AND ASSIGNING PHENOTYPES
p_thresh <- data.frame(fitness_stat_4C4 %>%
                         group_by(hours) %>%
                         summarise(p1_thresh = mean(p1_thresh, na.rm = T),
                                   p2_thresh = mean(p2_thresh, na.rm = T)))

fitness_stat_4C4$phenotype1 <- NULL
fitness_stat_4C4$phenotype1[fitness_stat_4C4$p1 <= fitness_stat_4C4$p1_thresh &
                              fitness_stat_4C4$es1 > 0 &
                              !(fitness_stat_4C4$orf_name %in% c('REF','BF_control'))] <- 'Beneficial'
fitness_stat_4C4$phenotype1[fitness_stat_4C4$p1 <= fitness_stat_4C4$p1_thresh &
                              fitness_stat_4C4$es1 < 0 &
                              !(fitness_stat_4C4$orf_name %in% c('REF','BF_control'))] <- 'Deleterious'
fitness_stat_4C4$phenotype1[is.na(fitness_stat_4C4$phenotype1) &
                              !(fitness_stat_4C4$orf_name %in% c('REF','BF_control'))] <- 'Neutral'

fitness_stat_4C4$phenotype2 <- NULL
fitness_stat_4C4$phenotype2[fitness_stat_4C4$p2 <= fitness_stat_4C4$p2_thresh &
                              fitness_stat_4C4$es2 > 0 &
                              !(fitness_stat_4C4$orf_name %in% c('REF','BF_control'))] <- 'Beneficial'
fitness_stat_4C4$phenotype2[fitness_stat_4C4$p2 <= fitness_stat_4C4$p2_thresh &
                              fitness_stat_4C4$es2 < 0 &
                              !(fitness_stat_4C4$orf_name %in% c('REF','BF_control'))] <- 'Deleterious'
fitness_stat_4C4$phenotype2[is.na(fitness_stat_4C4$phenotype2) &
                              !(fitness_stat_4C4$orf_name %in% c('REF','BF_control'))] <- 'Neutral'
head(fitness_stat_4C4)

diff_pheno <- dim(fitness_stat_4C4[fitness_stat_4C4$phenotype1 != fitness_stat_4C4$phenotype2 &
                   !is.na(fitness_stat_4C4$phenotype1),])[[1]]/
  dim(fitness_stat_4C4[!is.na(fitness_stat_4C4$phenotype1),])[[1]] * 100

cn.diff_pheno.plot <- ggplot(plyr::count(fitness_stat_4C4[fitness_stat_4C4$phenotype1 != fitness_stat_4C4$phenotype2 &
                               !is.na(fitness_stat_4C4$phenotype1),],
            vars = c('phenotype1','phenotype2')),
       aes(x = phenotype1, y = phenotype2)) +
  geom_point(aes(size = freq/6688 * 100, fill = freq), col = 'black', shape = 21) +
  scale_x_discrete(breaks = c('Deleterious','Neutral','Beneficial'),
                   limits = c('Deleterious','Neutral','Beneficial')) +
  scale_y_discrete(breaks = c('Deleterious','Neutral','Beneficial'),
                   limits = c('Deleterious','Neutral','Beneficial')) +
  scale_fill_distiller(name = 'Count',
                        palette="Spectral",
                        limits = c(0,65)) +
  scale_size_continuous(name = 'Overall %',
                        limits = c(0,1)) +
  labs(title = 'Condition Negative Dataset (Uniform CS Distribution)',
       subtitle = sprintf('%0.2f%% of Mutants have a Different Phenotype', diff_pheno),
       x = 'Empirical Testing using Tester',
       y = 'Empirical Testing using Reference') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles),
        plot.subtitle = element_text(size = txt-1),
        panel.background = element_rect(fill = 'grey90'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "right",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))

cn.diff_pheno.time <- ggplot(plyr::count(fitness_stat_4C4[fitness_stat_4C4$phenotype1 != fitness_stat_4C4$phenotype2 &
                                                                !is.na(fitness_stat_4C4$phenotype1),], vars = c('hours')),
                             aes(x = hours, y = 0)) +
  geom_point(aes(size = freq/6688 * 100, fill = freq), col = 'black', shape = 21) +
  scale_fill_distiller(name = 'Count',
                        palette="Spectral",
                        limits = c(0,65)) +
  scale_size_continuous(name = 'Overall %',
                        limits = c(0,1)) +
  labs(x = 'Time (hours)') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles),
        plot.subtitle = element_text(size = txt-1),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = 'grey90'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "right",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))

# cn.diff_pheno <- ggpubr::ggarrange(cn.diff_pheno.plot, cn.diff_pheno.time,
#                                    ncol = 1, common.legend = T, legend = 'bottom',
#                                    heights = c(3,1), align = 'v')
# ggsave(sprintf("%sr2r/FIGURE_CN_DIFFPHENO.jpg",out_path), cn.diff_pheno,
#        height = 110, width = 110, units = 'mm',
#        dpi = 300)

##### Sensitivity
head(fitness_stat_4C4)
ref_spe <- NULL
for (p in seq(0,1,0.01)) {
  temp <- plyr::count(fitness_stat_4C4[!is.na(fitness_stat_4C4$p1) & fitness_stat_4C4$p1 <= p,], vars = 'hours')
  temp$pthresh <- p
  ref_spe <- rbind(ref_spe, temp)
}
ref_spe <- data.frame(ref_spe)
head(ref_spe)

# load(sprintf('%sLID_SPECIFICITY_DATA.RData',out_path))
# lid.spe.data <- stats.tmp
ref.specificity <- ggplot() +
  stat_summary(data = ref_spe,
               aes(x = pthresh, y = (1-freq/1139)*100, group = 1),
               fun.data=mean_sdl, fun.args = list(mult=1),
               geom="ribbon", alpha = 0.4) +
  stat_summary(data = ref_spe,
               aes(x = pthresh, y = (1-freq/1139)*100, col = 'Tester'),
               fun=mean, geom="line", lwd =0.7) +
  stat_summary(data = lid.spe.data[lid.spe.data$hours == lid.spe.data$cont_hrs,],
               fun.data=mean_sdl, fun.args = list(mult=1),
               aes(x = p, y = (1-fpr)*100, group = 1),
               geom="ribbon", alpha = 0.4) +
  stat_summary(data = lid.spe.data[lid.spe.data$hours == lid.spe.data$cont_hrs,],
               aes(x = p, y = (1-fpr)*100, col = 'Reference'),
               fun=mean, geom="line", lwd =0.7) +
  labs(x = "p-value",
       y = "Specificity") +
  scale_x_continuous(breaks = seq(-1,1,0.025),
                     minor_breaks = seq(-1,1,0.0125)) +
  scale_y_continuous(breaks = seq(0,200,1),
                     minor_breaks = seq(0,200,0.5),
                     labels = paste(sprintf('%d',seq(0,200,1)),'%', sep = '')) +
  coord_cartesian(xlim = c(0, 0.1), ylim = c(90, 100)) +
  scale_color_manual(name = "Empirical Testing\nUsing",
                     breaks = c('Tester','Reference'),
                     values = c('Tester' = '#1976D2', 'Reference' = '#009688')) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(4, "mm"),
        legend.position = c(0.2,0.2),
        legend.margin = margin(0.5,0.5,0.5,0.5, "mm"))
# ggsave(sprintf("%sr2r/FIGURE_REF_SPECIFICITY.jpg",out_path), ref.specificity,
#        height = one.c, width = one.c, units = 'mm',
#        dpi = 300)

##### Using the condition positive dataset - bimodal distribution
ref_vps <- NULL
for (tr in sort(unique(fitness_4C4$hours))) {
  temp_tr <- fitness_4C4[fitness_4C4$hours == tr,]
  for (tm in sort(unique(fitness_4C4$hours))) {
    temp_tm <- fitness_4C4[fitness_4C4$hours == tm,]
    temp_tm$average[temp_tm$orf_name == 'BF_control' & !is.na(temp_tm$orf_name)] <-
      temp_tr$average[temp_tr$orf_name == 'BF_control' & !is.na(temp_tm$orf_name)]
    temp_tm$average[temp_tm$orf_name_new == 'REF' & !is.na(temp_tm$orf_name_new)] <-
      temp_tr$average[temp_tr$orf_name_new == 'REF' & !is.na(temp_tm$orf_name_new)]
    temp_tm$bg <- temp_tr$bg
    temp_tm$fitness <- temp_tm$average/temp_tm$bg
    temp_tm$tr <- tr
    temp_tm$tm <- tm
    temp_tm$effect <- median(temp_tm$average[temp_tm$orf_name != 'BF_control' & !is.na(temp_tm$orf_name)], na.rm = T)/
      median(temp_tr$average[temp_tr$orf_name == 'BF_control' & !is.na(temp_tm$orf_name)], na.rm = T)
    ref_vps <- rbind(ref_vps, temp_tm)
  }
}

# consolidate stats
ref_stat_vps <- data.frame(ref_vps[!is.na(ref_vps$fitness) & !(ref_vps$outlier),] %>%
                             group_by(tr, tm, effect, source, orf_name_new) %>%
                             summarise(cs_mean = mean(fitness, na.rm = T),
                                       cs_median = mean(fitness, na.rm = T),
                                       cs_std = mean(fitness, na.rm = T)))
colnames(ref_stat_vps) <- c('tr','tm','effect','source','orf_name','cs_mean','cs_median','cs_std')
ref_stat_vps$hours <- ref_stat_vps$tm
ref_stat_vps2 <- NULL
for (h in sort(unique(ref_stat_vps$tr))) {
  temp <- ref_stat_vps[ref_stat_vps$tr == h,]
  temp <- empirical_p(temp, 'REF', 'cs_mean', 'p1', 'es1', c('REF','BF_control'), 0.05)
  ref_stat_vps2 <- rbind(ref_stat_vps2, temp)
}
head(ref_stat_vps2)

ref_stat_vps2$phenotype <- NULL
ref_stat_vps2$phenotype[ref_stat_vps2$p1 <= ref_stat_vps2$p1_thresh &
                          ref_stat_vps2$es1 > 0 &
                          !(ref_stat_vps2$orf_name %in% c('REF','BF_control'))] <- 'Beneficial'
ref_stat_vps2$phenotype[ref_stat_vps2$p1 <= ref_stat_vps2$p1_thresh &
                          ref_stat_vps2$es1 < 0 &
                          !(ref_stat_vps2$orf_name %in% c('REF','BF_control'))] <- 'Deleterious'
ref_stat_vps2$phenotype[is.na(ref_stat_vps2$phenotype) &
                          !(ref_stat_vps2$orf_name %in% c('REF','BF_control'))] <- 'Neutral'

hello2 <- plyr::count(ref_stat_vps2[!(ref_stat_vps2$orf_name %in% c('REF','BF_control')),], vars = c('tr','tm','effect','phenotype'))
hello2 <- data.frame(merge(merge(hello2[hello2$phenotype == 'Beneficial',] %>% group_by(tr, tm, effect) %>% summarise(Beneficial = freq),
                                hello2[hello2$phenotype == 'Neutral',] %>% group_by(tr, tm, effect) %>% summarise(Neutral = freq),
                                by = c('tr','tm','effect'), all = T),
                          hello2[hello2$phenotype == 'Deleterious',] %>% group_by(tr, tm, effect) %>% summarise(Deleterious = freq),
                          by = c('tr','tm','effect'), all = T))
hello2[is.na(hello2)] <- 0

hello2$Beneficial <- hello2$Deleterious + hello2$Beneficial
hello2$Neutral <- hello2$Beneficial + hello2$Neutral

# histogram
t <- 608
ref.sensitivity <- ggplot(hello2[hello2$effect != 1,]) +
  geom_area(aes(x = (effect-1)*100, y = Neutral, fill = 'Neutral'), alpha = 1) +
  geom_area(aes(x = (effect-1)*100, y = Beneficial, fill = 'Beneficial'), alpha = 1) +
  geom_area(aes(x = (effect-1)*100, y = Deleterious, fill = 'Deleterious'), alpha = 1) +
  geom_vline(xintercept = seq(-2,2,0.025)*100, col = '#757575', lwd = 0.5, alpha =0.5) +
  geom_hline(yintercept = c(t * seq(0,1,0.05)), col = '#757575', lwd = 0.5, alpha =0.5) +
  geom_vline(xintercept = c(-5,5), col = 'red', lwd = 1, linetype = 'dashed') +
  geom_text(x=0, y=t*0.9, label="Empirical Testing\nUsing Tester", col = 'white', size = 3) +
  labs(x = 'Fitness Effect',
       y = 'Sensitivity') +
  scale_y_continuous(breaks = c(t * seq(0,1,0.1)),
                     minor_breaks = c(t * seq(0,1,0.05)),
                     labels = c('0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')) +
  scale_x_continuous(breaks = seq(-2,2,0.05)*100,
                     minor_breaks = seq(-2,2,0.025)*100,
                     labels = paste0(round(seq(-2,2,0.05)*100), '%', sep = '')) +
  scale_color_manual(name = 'Phenotype',
                     breaks = c('Beneficial','Neutral','Deleterious'),
                     values = c('Deleterious'='#3F51B5',
                                'Neutral'='#212121',
                                'Beneficial'='#FFC107'),
                     guide = F) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107')) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(4, "mm"),
        legend.position = c(0.87,0.2),
        legend.margin = margin(0.5,0.5,0.5,0.5, "mm")) +
  coord_cartesian(xlim = c(-10,10),
                  ylim = c(0,t))

# lid sensitivity
load(sprintf('%sLID_SENSITIVITY_DATA.RData',out_path))
lid.sen.data <- dat.cnt2

t <- lid.sen.data$Beneficial[1]
lid.sensitivity <- ggplot(lid.sen.data) +
  geom_area(aes(x = (cen-1)*100, y = Neutral, fill = 'Neutral'), alpha = 1) +
  geom_area(aes(x = (cen-1)*100, y = Beneficial, fill = 'Beneficial'), alpha = 1) +
  geom_area(aes(x = (cen-1)*100, y = Deleterious, fill = 'Deleterious'), alpha = 1) +
  geom_vline(xintercept = seq(-2,2,0.025)*100, col = '#757575', lwd = 0.5, alpha =0.5) +
  geom_hline(yintercept = c(t * seq(0,1,0.05)), col = '#757575', lwd = 0.5, alpha =0.5) +
  geom_vline(xintercept = c(-5,5), col = 'red', lwd = 1, linetype = 'dashed') +
  geom_text(x=0, y=t*0.9, label="Empirical Testing\nUsing Reference", col = 'white', size = 3) +
  labs(x = 'Fitness Effect',
       y = 'Sensitivity') +
  scale_y_continuous(breaks = c(t * seq(0,1,0.1)),
                     minor_breaks = c(t * seq(0,1,0.05)),
                     labels = c('0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')) +
  scale_x_continuous(breaks = seq(-2,2,0.05)*100,
                     minor_breaks = seq(-2,2,0.025)*100,
                     labels = paste0(round(seq(-2,2,0.05)*100), '%', sep = '')) +
  scale_color_manual(name = 'Phenotype',
                     breaks = c('Beneficial','Neutral','Deleterious'),
                     values = c('Deleterious'='#3F51B5',
                                'Neutral'='#212121',
                                'Beneficial'='#FFC107'),
                     guide = F) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107')) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(4, "mm"),
        legend.position = c(0.87,0.2),
        legend.margin = margin(0.5,0.5,0.5,0.5, "mm")) +
  coord_cartesian(xlim = c(-10,10),
                  ylim = c(0,t))

ref.perf <- ggpubr::ggarrange(ref.specificity,
                              ggpubr::ggarrange(ref.sensitivity, lid.sensitivity,
                                                nrow = 1,
                                                labels = c('B','C'),
                                                font.label = list(face = 'bold', size = lbls, family = "sans"),
                                                hjust=-1),
                              nrow = 2, heights = c(1.2,2),
                              labels = c('A',''),
                              font.label = list(face = 'bold', size = lbls, family = "sans"),
                              hjust=-1)
ggsave(sprintf("%sr2r/FIGURE_REF_PERFORMANCE.jpg",out_path), ref.perf,
       height = 140, width = two.c, units = 'mm',
       dpi = 300)

##### Using the condition positive dataset - random distribution
fitness_rnd_4C4 <- dbGetQuery(conn, 'select a.*, b.rnd_hrs
                              from 4C4_FS_RND2_6144_FITNESS a, 4C4_FS_RND2_6144_DATA b
                              where a.hours = b.hours and a.pos = b.pos
                              and a.hours > 0 and a.orf_name is not NULL
                              order by hours, orf_name')
fitness_rnd_4C4 <- merge(fitness_rnd_4C4, fitness_4C4[,c(1,3,4,6,12)], by = c('hours','pos'), suffixes = c('','_new'))

fitness_rnd_4C4$average[fitness_rnd_4C4$orf_name != fitness_rnd_4C4$orf_name_new] <-
  fitness_rnd_4C4$average_new[fitness_rnd_4C4$orf_name != fitness_rnd_4C4$orf_name_new]
fitness_rnd_4C4$rnd_hrs[fitness_rnd_4C4$orf_name != fitness_rnd_4C4$orf_name_new] <-
  fitness_rnd_4C4$hours[fitness_rnd_4C4$orf_name != fitness_rnd_4C4$orf_name_new]

fitness_rnd_4C4$fitness <- fitness_rnd_4C4$average/fitness_rnd_4C4$bg

fitness_stat_rnd_4C4 <- data.frame(fitness_rnd_4C4[!is.na(fitness_rnd_4C4$fitness),] %>%
                                 group_by(hours, source, orf_name_new) %>%
                                 summarise(cs_mean = mean(fitness, na.rm = T),
                                           cs_median = mean(fitness, na.rm = T),
                                           cs_std = mean(fitness, na.rm = T)))

colnames(fitness_stat_rnd_4C4) <- c("hours", "source", "orf_name", "cs_mean", "cs_median", "cs_std")

fitness_stat_rnd_4C4 <- empirical_p(fitness_stat_rnd_4C4, 'REF', 'cs_mean', 'p1', 'es1', c('REF','BF_control'), 0.05)
fitness_stat_rnd_4C4 <- empirical_p(fitness_stat_rnd_4C4, 'BF_control', 'cs_mean', 'p2', 'es2', c('REF','BF_control'), 0.05)
head(fitness_stat_rnd_4C4)

fitness_stat_rnd_4C4 <- merge(fitness_stat_rnd_4C4, p_thresh, by = 'hours', suffixes = c('','_org'))

fitness_stat_rnd_4C4$phenotype1 <- NULL
fitness_stat_rnd_4C4$phenotype1[fitness_stat_rnd_4C4$p1 <= fitness_stat_rnd_4C4$p1_thresh_org &
                              fitness_stat_rnd_4C4$es1 > 0 &
                              !(fitness_stat_rnd_4C4$orf_name %in% c('REF','BF_control'))] <- 'Beneficial'
fitness_stat_rnd_4C4$phenotype1[fitness_stat_rnd_4C4$p1 <= fitness_stat_rnd_4C4$p1_thresh_org &
                              fitness_stat_rnd_4C4$es1 < 0 &
                              !(fitness_stat_rnd_4C4$orf_name %in% c('REF','BF_control'))] <- 'Deleterious'
fitness_stat_rnd_4C4$phenotype1[is.na(fitness_stat_rnd_4C4$phenotype1) &
                              !(fitness_stat_rnd_4C4$orf_name %in% c('REF','BF_control'))] <- 'Neutral'

fitness_stat_rnd_4C4$phenotype2 <- NULL
fitness_stat_rnd_4C4$phenotype2[fitness_stat_rnd_4C4$p2 <= fitness_stat_rnd_4C4$p2_thresh_org &
                              fitness_stat_rnd_4C4$es2 > 0 &
                              !(fitness_stat_rnd_4C4$orf_name %in% c('REF','BF_control'))] <- 'Beneficial'
fitness_stat_rnd_4C4$phenotype2[fitness_stat_rnd_4C4$p2 <= fitness_stat_rnd_4C4$p2_thresh_org &
                              fitness_stat_rnd_4C4$es2 < 0 &
                              !(fitness_stat_rnd_4C4$orf_name %in% c('REF','BF_control'))] <- 'Deleterious'
fitness_stat_rnd_4C4$phenotype2[is.na(fitness_stat_rnd_4C4$phenotype2) &
                              !(fitness_stat_rnd_4C4$orf_name %in% c('REF','BF_control'))] <- 'Neutral'
head(fitness_stat_rnd_4C4)

diff_pheno <- dim(fitness_stat_rnd_4C4[fitness_stat_rnd_4C4$phenotype1 != fitness_stat_rnd_4C4$phenotype2 &
                       !is.na(fitness_stat_rnd_4C4$phenotype1),])[[1]]/
  dim(fitness_stat_rnd_4C4[!is.na(fitness_stat_rnd_4C4$phenotype1),])[[1]] * 100

cp.diff_pheno.plot <- ggplot(plyr::count(fitness_stat_rnd_4C4[fitness_stat_rnd_4C4$phenotype1 != fitness_stat_rnd_4C4$phenotype2 &
                                                            !is.na(fitness_stat_rnd_4C4$phenotype1),],
                                         vars = c('phenotype1','phenotype2')),
                             aes(x = phenotype1, y = phenotype2)) +
  geom_point(aes(size = freq/6688 * 100, fill = freq), col = 'black', shape = 21) +
  scale_x_discrete(breaks = c('Deleterious','Neutral','Beneficial'),
                   limits = c('Deleterious','Neutral','Beneficial')) +
  scale_y_discrete(breaks = c('Deleterious','Neutral','Beneficial'),
                   limits = c('Deleterious','Neutral','Beneficial')) +
  scale_fill_distiller(name = 'Count',
                        palette="Spectral",
                        limits = c(0,65)) +
  scale_size_continuous(name = 'Overall %',
                        limits = c(0,1)) +
  labs(title = 'Condition Positive Dataset (Random CS Distribution)',
       subtitle = sprintf('%0.2f%% of Mutants have a Different Phenotype', diff_pheno),
       x = 'Empirical Testing using Tester',
       y = 'Empirical Testing using Reference') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles),
        plot.subtitle = element_text(size = txt-1),
        panel.background = element_rect(fill = 'grey90'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "right",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))

cp.diff_pheno.time <- ggplot(plyr::count(fitness_stat_rnd_4C4[fitness_stat_rnd_4C4$phenotype1 != fitness_stat_rnd_4C4$phenotype2 &
                       !is.na(fitness_stat_rnd_4C4$phenotype1),], vars = c('hours')),
       aes(x = hours, y = 0)) +
  geom_point(aes(size = freq/6688 * 100, fill = freq), col = 'black', shape = 21) +
  scale_fill_distiller(name = 'Count',
                     palette="Spectral",
                     limits = c(0,65)) +
  scale_size_continuous(name = 'Overall %',
                        limits = c(0,1)) +
  labs(x = 'Time (hours)') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles),
        plot.subtitle = element_text(size = txt-1),
        panel.background = element_rect(fill = 'grey90'),
        panel.grid = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "right",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))

# cp.diff_pheno <- ggpubr::ggarrange(cp.diff_pheno.plot, cp.diff_pheno.time,
#                   ncol = 1, common.legend = T, legend = 'bottom',
#                   heights = c(3,1), align = 'v')
# ggsave(sprintf("%sr2r/FIGURE_CP_DIFFPHENO.jpg",out_path), cp.diff_pheno,
#        height = 110, width = 110, units = 'mm',
#        dpi = 300)

ref.diff_pheno <- ggpubr::ggarrange(cn.diff_pheno.plot, cp.diff_pheno.plot,
                                    cn.diff_pheno.time, cp.diff_pheno.time,
                                    labels = c('A','B','C','D'),
                                    font.label = list(face = 'bold', size = lbls, family = "sans"),
                                    hjust=-1,
                                    ncol = 2, nrow = 2,
                                    heights = c(3,1), align = 'v',
                                    common.legend = T, legend = 'bottom')
ggsave(sprintf("%sr2r/FIGURE_REF_DIFFPHENO.jpg",out_path), ref.diff_pheno,
       height = 130, width = 110*2, units = 'mm',
       dpi = 300)

##### l. 476 IMPACT OF BRODER REMOVAL
bor_4C4 <- dbGetQuery(conn, 'select * from 4C4_BOR_6144_FITNESS
                      where hours > 0 and orf_name is not NULL
                      order by hours, orf_name')
bor_4C4$source <- as.numeric(substr(as.character(bor_4C4$pos), 3, 6))
head(bor_4C4)

# remove outliers
bor_4C4$outlier <- FALSE
for (h in sort(unique(bor_4C4$hours))) {
  for (p in sort(unique(bor_4C4$source[bor_4C4$hours == h]))) {
    temp <- bor_4C4$fitness[bor_4C4$hours == h & bor_4C4$source == p]
    m <- median(temp, na.rm = T)
    if (!is.na(m)) {
      madev <- mad(temp, na.rm = T)
      ul <- m + 3*madev
      ll <- m - 3*madev
      
      bor_4C4$outlier[bor_4C4$hours == h & bor_4C4$source == p & bor_4C4$fitness > ul] <- TRUE
      bor_4C4$outlier[bor_4C4$hours == h & bor_4C4$source == p & bor_4C4$fitness < ll] <- TRUE
    }
  }
}

##### CONDITION NEGATIVE ANALYSIS
# consolidate stats
bor_stat_4C4 <- data.frame(bor_4C4[!is.na(bor_4C4$fitness) & !(bor_4C4$outlier),] %>%
                                 group_by(hours, source, orf_name) %>%
                                 summarise(cs_mean = mean(fitness, na.rm = T),
                                           cs_median = mean(fitness, na.rm = T),
                                           cs_std = mean(fitness, na.rm = T)))

bor_stat_4C4 <- empirical_p(bor_stat_4C4, 'BF_control', 'cs_mean', 'p1', 'es1', c('REF','BF_control'), 0.05)
head(bor_stat_4C4)

plyr::count(bor_stat_4C4[!is.na(bor_stat_4C4$p1) & bor_stat_4C4$p1 <= 0.05,], vars = 'hours')/
  plyr::count(bor_stat_4C4[!is.na(bor_stat_4C4$p1),], vars = 'hours')

bor_spe <- NULL
for (p in seq(0,1,0.01)) {
  temp <- plyr::count(bor_stat_4C4[!is.na(bor_stat_4C4$p1) & bor_stat_4C4$p1 <= p,], vars = 'hours')
  temp$pthresh <- p
  bor_spe <- rbind(bor_spe, temp)
}
bor_spe <- data.frame(bor_spe)
head(bor_spe)

load(sprintf('%sLID_SPECIFICITY_DATA.RData',out_path))
lid.spe.data <- stats.tmp

bor.specificity <- ggplot() +
  # stat_summary(data = bor_spe,
  #              aes(x = pthresh, y = (1-freq/1139)*100, group = 1),
  #              fun.data=mean_sdl, fun.args = list(mult=1),
  #              geom="ribbon", alpha = 0.4) +
  stat_summary(data = bor_spe,
               aes(x = pthresh, y = (1-freq/1139)*100, col = 'No'),
               fun=mean, geom="line", lwd =0.7) +
  # stat_summary(data = lid.spe.data[lid.spe.data$hours == lid.spe.data$cont_hrs,],
  #              fun.data=mean_sdl, fun.args = list(mult=1),
  #              aes(x = p, y = (1-fpr)*100, group = 1),
  #              geom="ribbon", alpha = 0.4) +
  stat_summary(data = lid.spe.data[lid.spe.data$hours == lid.spe.data$cont_hrs,],
               aes(x = p, y = (1-fpr)*100, col = 'Yes'),
               fun=mean, geom="line", lwd =0.7) +
  labs(x = "p-value",
       y = "Specificity") +
  scale_x_continuous(breaks = seq(-1,1,0.025),
                     minor_breaks = seq(-1,1,0.0125)) +
  scale_y_continuous(breaks = seq(0,200,1),
                     minor_breaks = seq(0,200,0.5),
                     labels = paste(sprintf('%d',seq(0,200,1)),'%', sep = '')) +
  coord_cartesian(xlim = c(0, 0.1), ylim = c(90, 100)) +
  scale_color_manual(name = "Border Removed",
                     breaks = c('No','Yes'),
                     values = c('No' = '#1976D2', 'Yes' = '#009688')) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(4, "mm"),
        legend.position = c(0.2,0.2),
        legend.margin = margin(0.5,0.5,0.5,0.5, "mm"))
ggsave(sprintf("%sr2r/FIGURE_BORDER_SPECIFICITY.jpg",out_path), bor.specificity,
       height = one.c, width = one.c, units = 'mm',
       dpi = 300)


##### CONDITION POSITIVE ANALYSIS
# building the dataset
bor_vps <- NULL
for (tr in sort(unique(bor_4C4$hours))) {
  temp_tr <- bor_4C4[bor_4C4$hours == tr,]
  for (tm in sort(unique(bor_4C4$hours))) {
    temp_tm <- bor_4C4[bor_4C4$hours == tm,]
    temp_tm$average[temp_tm$orf_name == 'BF_control'] <- temp_tr$average[temp_tr$orf_name == 'BF_control']
    temp_tm$bg <- temp_tr$bg
    temp_tm$fitness <- temp_tm$average/temp_tm$bg
    temp_tm$tr <- tr
    temp_tm$tm <- tm
    temp_tm$effect <- median(temp_tm$average[temp_tm$orf_name != 'BF_control'], na.rm = T)/
      median(temp_tr$average[temp_tr$orf_name == 'BF_control'], na.rm = T)
    bor_vps <- rbind(bor_vps, temp_tm)
  }
}

# consolidate stats
bor_stat_vps <- data.frame(bor_vps[!is.na(bor_vps$fitness) & !(bor_vps$outlier),] %>%
                             group_by(tr, tm, effect, source, orf_name) %>%
                             summarise(cs_mean = mean(fitness, na.rm = T),
                                       cs_median = mean(fitness, na.rm = T),
                                       cs_std = mean(fitness, na.rm = T)))
bor_stat_vps$hours <- bor_stat_vps$tm
bor_stat_vps2 <- NULL
for (h in sort(unique(bor_stat_vps$tr))) {
  temp <- bor_stat_vps[bor_stat_vps$tr == h,]
  temp <- empirical_p(temp, 'BF_control', 'cs_mean', 'p1', 'es1', c('REF','BF_control'), 0.05)
  bor_stat_vps2 <- rbind(bor_stat_vps2, temp)
}
head(bor_stat_vps2)

bor_stat_vps2$phenotype <- NULL
bor_stat_vps2$phenotype[bor_stat_vps2$p1 <= bor_stat_vps2$p1_thresh &
                          bor_stat_vps2$es1 > 0 &
                          !(bor_stat_vps2$orf_name %in% c('REF','BF_control'))] <- 'Beneficial'
bor_stat_vps2$phenotype[bor_stat_vps2$p1 <= bor_stat_vps2$p1_thresh &
                          bor_stat_vps2$es1 < 0 &
                          !(bor_stat_vps2$orf_name %in% c('REF','BF_control'))] <- 'Deleterious'
bor_stat_vps2$phenotype[is.na(bor_stat_vps2$phenotype) &
                          !(bor_stat_vps2$orf_name %in% c('REF','BF_control'))] <- 'Neutral'

hello <- plyr::count(bor_stat_vps2[!(bor_stat_vps2$orf_name %in% c('REF','BF_control')),], vars = c('tr','tm','effect','phenotype'))
hello <- data.frame(merge(merge(hello[hello$phenotype == 'Beneficial',] %>% group_by(tr, tm, effect) %>% summarise(Beneficial = freq),
      hello[hello$phenotype == 'Neutral',] %>% group_by(tr, tm, effect) %>% summarise(Neutral = freq),
      by = c('tr','tm','effect'), all = T),
      hello[hello$phenotype == 'Deleterious',] %>% group_by(tr, tm, effect) %>% summarise(Deleterious = freq),
      by = c('tr','tm','effect'), all = T))
hello[is.na(hello)] <- 0

hello$Beneficial <- hello$Deleterious + hello$Beneficial
hello$Neutral <- hello$Beneficial + hello$Neutral

ggplot(hello,
       aes(x = (effect - 1)* 100, y = freq/1139 * 100)) +
  # geom_point(aes(col = phenotype)) +
  stat_summary(aes(col = phenotype),
               fun=mean, geom="line", lwd =0.7) +
  coord_cartesian(xlim = c(-10,10))

# histogram
t <- 1139
bor.sensitivity <- ggplot(hello[hello$effect != 1,]) +
  geom_area(aes(x = (effect-1)*100, y = Neutral, fill = 'Neutral'), alpha = 1) +
  geom_area(aes(x = (effect-1)*100, y = Beneficial, fill = 'Beneficial'), alpha = 1) +
  geom_area(aes(x = (effect-1)*100, y = Deleterious, fill = 'Deleterious'), alpha = 1) +
  geom_vline(xintercept = seq(-2,2,0.025)*100, col = '#757575', lwd = 0.5, alpha =0.5) +
  geom_hline(yintercept = c(t * seq(0,1,0.05)), col = '#757575', lwd = 0.5, alpha =0.5) +
  geom_vline(xintercept = c(-5,5), col = 'red', lwd = 1, linetype = 'dashed') +
  geom_text(x=0, y=t*0.9, label="LID + Border", col = 'white', size = 3) +
  labs(x = 'Fitness Effect',
       y = 'Sensitivity') +
  scale_y_continuous(breaks = c(t * seq(0,1,0.1)),
                     minor_breaks = c(t * seq(0,1,0.05)),
                     labels = c('0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')) +
  scale_x_continuous(breaks = seq(-2,2,0.05)*100,
                     minor_breaks = seq(-2,2,0.025)*100,
                     labels = paste0(round(seq(-2,2,0.05)*100), '%', sep = '')) +
  scale_color_manual(name = 'Phenotype',
                     breaks = c('Beneficial','Neutral','Deleterious'),
                     values = c('Deleterious'='#3F51B5',
                                'Neutral'='#212121',
                                'Beneficial'='#FFC107'),
                     guide = F) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107')) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(4, "mm"),
        legend.position = c(0.87,0.2),
        legend.margin = margin(0.5,0.5,0.5,0.5, "mm")) +
  coord_cartesian(xlim = c(-10,10),
                  ylim = c(0,t))

bor.perf <- ggpubr::ggarrange(bor.specificity, bor.sensitivity,
                  nrow = 1,
                  labels = c('A','B'),
                  font.label = list(face = 'bold', size = lbls, family = "sans"),
                  hjust=-1)
ggsave(sprintf("%sr2r/FIGURE_BORDER_PERFORMANCE.jpg",out_path), bor.perf,
       height = one.c, width = two.c, units = 'mm',
       dpi = 300)

hours.labs <- sprintf('Time = %0.1f hours', sort(unique(bor_stat_4C4$hours)))
names(hours.labs ) <- sort(unique(bor_stat_4C4$hours))

bor.fitness.density <- ggplot() +
  # geom_density_ridges(data = bor_stat_4C4[bor_stat_4C4$orf_name == 'BF_control',],
  geom_density_ridges(data = bor_stat_4C4,
            aes(x = cs_mean, y = 'Yes'),
            quantile_lines = TRUE,
            quantiles = c(0.25, 0.5, 0.75),
            scale = 0.05, alpha = 0.8, size = 0.2,
            vline_size = 0.2, vline_color = "black",
            na.rm = T) +
  # geom_density_ridges(data = fitness_stat_4C4[fitness_stat_4C4$orf_name == 'BF_control',],
  geom_density_ridges(data = fitness_stat_4C4,
            aes(x = cs_mean, y = 'No'),
            quantile_lines = TRUE,
            quantiles = c(0.25, 0.5, 0.75),
            scale = 0.05, alpha = 0.8, size = 0.2,
            vline_size = 0.2, vline_color = "black",
            na.rm = T) +
  labs(x = 'Fitness',
       y = 'Border') +
  facet_wrap(.~hours, labeller = labeller(hours = hours.labs)) +
  coord_cartesian(xlim = c(0.7,1.3)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(4, "mm"),
        legend.position = c(0.2,0.2),
        legend.margin = margin(0.5,0.5,0.5,0.5, "mm"))
ggsave(sprintf("%sr2r/FIGURE_BORDER_DENSITY.jpg",out_path), bor.fitness.density,
       height = two.c, width = two.c, units = 'mm',
       dpi = 300)


