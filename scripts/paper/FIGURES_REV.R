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
library(RMariaDB)

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
                        on lid.pos = mcat.pos and lid.hours = mcat.hours')


rnd2_data$phenotype[rnd2_data$hours < rnd2_data$rnd_hrs] <- 'Beneficial'
rnd2_data$phenotype[rnd2_data$hours > rnd2_data$rnd_hrs] <- 'Deleterious'
rnd2_data$phenotype[rnd2_data$hours == rnd2_data$rnd_hrs] <- 'Neutral'

rnd2_data$lid_phenotype[rnd2_data$lid_p <= 0.05 & rnd2_data$lid_stat > 0] <- 'Beneficial'
rnd2_data$lid_phenotype[rnd2_data$lid_p <= 0.05 & rnd2_data$lid_stat < 0] <- 'Deleterious'
rnd2_data$lid_phenotype[is.na(rnd2_data$lid_phenotype)] <- 'Neutral'
dim(rnd2_data[rnd2_data$phenotype != rnd2_data$lid_phenotype,])[[1]]

rnd2_data$mcat_phenotype[rnd2_data$mcat_p <= 0.05 & rnd2_data$mcat_stat > 0] <- 'Beneficial'
rnd2_data$mcat_phenotype[rnd2_data$mcat_p <= 0.05 & rnd2_data$mcat_stat < 0] <- 'Deleterious'
rnd2_data$mcat_phenotype[is.na(rnd2_data$mcat_phenotype)] <- 'Neutral'
dim(rnd2_data[rnd2_data$phenotype != rnd2_data$mcat_phenotype,])[[1]]

head(rnd2_data) 

for (h in sort(unique(rnd2_data$hours))) {
  for (o in sort(unique(rnd2_data$orf_name[rnd2_data$hours == h & rnd2_data$orf_name != 'BF_control']))){
    ref_fit <- rnd2_data$lid_fitness[rnd2_data$hours == h & rnd2_data$orf_name == 'BF_control' & !is.na(rnd2_data$lid_fitness)]
    orf_fit <- rnd2_data$lid_fitness[rnd2_data$hours == h & rnd2_data$orf_name == o & !is.na(rnd2_data$lid_fitness)]
    wilcox_stats <- wilcox.test(orf_fit, ref_fit)
    rnd2_data$lid_wilcox_p[rnd2_data$hours == h & rnd2_data$orf_name == o] <- wilcox_stats$p.value
    
    ref_fit <- rnd2_data$mcat_fitness[rnd2_data$hours == h & rnd2_data$orf_name == 'BF_control' & !is.na(rnd2_data$mcat_fitness)]
    orf_fit <- rnd2_data$mcat_fitness[rnd2_data$hours == h & rnd2_data$orf_name == o & !is.na(rnd2_data$mcat_fitness)]
    wilcox_stats <- wilcox.test(orf_fit, ref_fit)
    rnd2_data$mcat_wilcox_p[rnd2_data$hours == h & rnd2_data$orf_name == o] <- wilcox_stats$p.value
  }
  rnd2_data$lid_adj_p[rnd2_data$hours == h] <- p.adjust(rnd2_data$lid_wilcox_p[rnd2_data$hours == h], method = "bonferroni")
  rnd2_data$mcat_adj_p[rnd2_data$hours == h] <- p.adjust(rnd2_data$mcat_wilcox_p[rnd2_data$hours == h], method = "bonferroni")
}

rnd2_data$lid_wilcox_phenotype[rnd2_data$lid_wilcox_p <= 0.05 & rnd2_data$lid_stat > 0] <- 'Beneficial'
rnd2_data$lid_wilcox_phenotype[rnd2_data$lid_wilcox_p <= 0.05 & rnd2_data$lid_stat < 0] <- 'Deleterious'
rnd2_data$lid_wilcox_phenotype[is.na(rnd2_data$lid_wilcox_phenotype)] <- 'Neutral'

rnd2_data$mcat_wilcox_phenotype[rnd2_data$mcat_wilcox_p <= 0.05 & rnd2_data$mcat_stat > 0] <- 'Beneficial'
rnd2_data$mcat_wilcox_phenotype[rnd2_data$mcat_wilcox_p <= 0.05 & rnd2_data$mcat_stat < 0] <- 'Deleterious'
rnd2_data$mcat_wilcox_phenotype[is.na(rnd2_data$mcat_wilcox_phenotype)] <- 'Neutral'

head(rnd2_data)

pheno_mismatch <- NULL
pheno_mismatch$hours <- sort(unique(rnd2_data$hours))
pheno_mismatch <- merge(pheno_mismatch,
                        plyr::count(plyr::count(rnd2_data[rnd2_data$phenotype != rnd2_data$lid_phenotype,],
                                    vars = c('hours', 'orf_name'))[,1:2], vars = 'hours'),
                        by = 'hours', all = T)
pheno_mismatch <- merge(pheno_mismatch,
                        plyr::count(plyr::count(rnd2_data[rnd2_data$phenotype != rnd2_data$lid_wilcox_phenotype,],
                                    vars = c('hours', 'orf_name'))[,1:2], vars = 'hours'),
                        by = 'hours', all = T)
pheno_mismatch <- merge(pheno_mismatch,
                        plyr::count(plyr::count(rnd2_data[rnd2_data$phenotype != rnd2_data$mcat_phenotype,],
                                    vars = c('hours', 'orf_name'))[,1:2], vars = 'hours'),
                        by = 'hours', all = T)
pheno_mismatch <- merge(pheno_mismatch,
                        plyr::count(plyr::count(rnd2_data[rnd2_data$phenotype != rnd2_data$mcat_wilcox_phenotype,],
                                    vars = c('hours', 'orf_name'))[,1:2], vars = 'hours'),
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
head(rnd2_data)
rnd2_data$lid_results[rnd2_data$phenotype == 'Beneficial' & rnd2_data$lid_phenotype == 'Beneficial'] <- 'True Positive'
rnd2_data$lid_results[rnd2_data$phenotype == 'Deleterious' & rnd2_data$lid_phenotype == 'Deleterious'] <- 'True Positive'
rnd2_data$lid_results[rnd2_data$phenotype == 'Neutral' & rnd2_data$lid_phenotype == 'Neutral'] <- 'True Negative'
rnd2_data$lid_results[rnd2_data$phenotype == 'Neutral' & rnd2_data$lid_phenotype != 'Neutral'] <- 'False Positive'
rnd2_data$lid_results[rnd2_data$phenotype != 'Neutral' & rnd2_data$lid_phenotype == 'Neutral'] <- 'False Negative'
rnd2_data$lid_results[rnd2_data$phenotype != 'Neutral' & rnd2_data$lid_phenotype != 'Neutral' &
                        rnd2_data$phenotype != rnd2_data$lid_phenotype] <- 'Switched Positive'

rnd2_data$lid_wilcox_results[rnd2_data$phenotype == 'Beneficial' & rnd2_data$lid_wilcox_phenotype == 'Beneficial'] <- 'True Positive'
rnd2_data$lid_wilcox_results[rnd2_data$phenotype == 'Deleterious' & rnd2_data$lid_wilcox_phenotype == 'Deleterious'] <- 'True Positive'
rnd2_data$lid_wilcox_results[rnd2_data$phenotype == 'Neutral' & rnd2_data$lid_wilcox_phenotype == 'Neutral'] <- 'True Negative'
rnd2_data$lid_wilcox_results[rnd2_data$phenotype == 'Neutral' & rnd2_data$lid_wilcox_phenotype != 'Neutral'] <- 'False Positive'
rnd2_data$lid_wilcox_results[rnd2_data$phenotype != 'Neutral' & rnd2_data$lid_wilcox_phenotype == 'Neutral'] <- 'False Negative'
rnd2_data$lid_wilcox_results[rnd2_data$phenotype != 'Neutral' & rnd2_data$lid_wilcox_phenotype != 'Neutral' &
                        rnd2_data$phenotype != rnd2_data$lid_wilcox_phenotype] <- 'Switched Positive'

rnd2_data$mcat_results[rnd2_data$phenotype == 'Beneficial' & rnd2_data$mcat_phenotype == 'Beneficial'] <- 'True Positive'
rnd2_data$mcat_results[rnd2_data$phenotype == 'Deleterious' & rnd2_data$mcat_phenotype == 'Deleterious'] <- 'True Positive'
rnd2_data$mcat_results[rnd2_data$phenotype == 'Neutral' & rnd2_data$mcat_phenotype == 'Neutral'] <- 'True Negative'
rnd2_data$mcat_results[rnd2_data$phenotype == 'Neutral' & rnd2_data$mcat_phenotype != 'Neutral'] <- 'False Positive'
rnd2_data$mcat_results[rnd2_data$phenotype != 'Neutral' & rnd2_data$mcat_phenotype == 'Neutral'] <- 'False Negative'
rnd2_data$mcat_results[rnd2_data$phenotype != 'Neutral' & rnd2_data$mcat_phenotype != 'Neutral' &
                        rnd2_data$phenotype != rnd2_data$mcat_phenotype] <- 'Switched Positive'

rnd2_data$mcat_wilcox_results[rnd2_data$phenotype == 'Beneficial' & rnd2_data$mcat_wilcox_phenotype == 'Beneficial'] <- 'True Positive'
rnd2_data$mcat_wilcox_results[rnd2_data$phenotype == 'Deleterious' & rnd2_data$mcat_wilcox_phenotype == 'Deleterious'] <- 'True Positive'
rnd2_data$mcat_wilcox_results[rnd2_data$phenotype == 'Neutral' & rnd2_data$mcat_wilcox_phenotype == 'Neutral'] <- 'True Negative'
rnd2_data$mcat_wilcox_results[rnd2_data$phenotype == 'Neutral' & rnd2_data$mcat_wilcox_phenotype != 'Neutral'] <- 'False Positive'
rnd2_data$mcat_wilcox_results[rnd2_data$phenotype != 'Neutral' & rnd2_data$mcat_wilcox_phenotype == 'Neutral'] <- 'False Negative'
rnd2_data$mcat_wilcox_results[rnd2_data$phenotype != 'Neutral' & rnd2_data$mcat_wilcox_phenotype != 'Neutral' &
                               rnd2_data$phenotype != rnd2_data$mcat_wilcox_phenotype] <- 'Switched Positive'

ggplot(rnd2_data,
       aes(x = hours)) +
  geom_histogram(aes(fill = mcat_results))

lid_results <- data.frame(plyr::count(rnd2_data, vars = c('hours','orf_name','lid_results'))[,c(1:3)], method = 'lid_emp')
colnames(lid_results) <- c('hours','orf_name','results','method')           

lid_wilcox_results <- data.frame(plyr::count(rnd2_data, vars = c('hours','orf_name','lid_wilcox_results'))[,c(1:3)], method = 'lid_wilcox')
colnames(lid_wilcox_results) <- c('hours','orf_name','results','method') 

mcat_results <- data.frame(plyr::count(rnd2_data, vars = c('hours','orf_name','mcat_results'))[,c(1:3)], method = 'mcat_emp')
colnames(mcat_results) <- c('hours','orf_name','results','method') 

mcat_wilcox_results <- data.frame(plyr::count(rnd2_data, vars = c('hours','orf_name','mcat_wilcox_results'))[,c(1:3)], method = 'mcat_wilcox')
colnames(mcat_wilcox_results) <- c('hours','orf_name','results','method') 

all_results <- rbind(lid_results, lid_wilcox_results, mcat_results, mcat_wilcox_results)
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
all_results_sum$freq <- all_results_sum$freq/(9.11*11)
head(all_results_sum)

ggplot(all_results_sum, aes(x = "", y = freq, fill = results)) +
  geom_bar(stat = 'identity', col = 'black') +
  coord_polar("y", start = 0) +
  geom_label_repel(aes(label = sprintf("%0.2f%%",freq)),
                   position = position_stack(vjust = 0.5),
                   colour ='white',
                   label.size = 0.15,
                   show.legend = F) +
  scale_fill_discrete(name = 'Method') +
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


all_results_sum$freq[all_results_sum$method == 'lid_emp' & all_results_sum$results == 'True Positive']/
  sum(all_results_sum$freq[all_results_sum$method == 'lid_emp' & all_results_sum$results %in% c('True Positive', 'False Negative')]) * 100
all_results_sum$freq[all_results_sum$method == 'lid_emp' & all_results_sum$results == 'True Negative']/
  sum(all_results_sum$freq[all_results_sum$method == 'lid_emp' & all_results_sum$results %in% c('True Negative', 'False Positive')]) * 100

all_results_sum$freq[all_results_sum$method == 'lid_wilcox' & all_results_sum$results == 'True Positive']/
  sum(all_results_sum$freq[all_results_sum$method == 'lid_wilcox' & all_results_sum$results %in% c('True Positive', 'False Negative')]) * 100
all_results_sum$freq[all_results_sum$method == 'lid_wilcox' & all_results_sum$results == 'True Negative']/
  sum(all_results_sum$freq[all_results_sum$method == 'lid_wilcox' & all_results_sum$results %in% c('True Negative', 'False Positive')]) * 100


all_results_sum$freq[all_results_sum$method == 'mcat_emp' & all_results_sum$results == 'True Positive']/
  sum(all_results_sum$freq[all_results_sum$method == 'mcat_emp' & all_results_sum$results %in% c('True Positive', 'False Negative')]) * 100
all_results_sum$freq[all_results_sum$method == 'mcat_emp' & all_results_sum$results == 'True Negative']/
  sum(all_results_sum$freq[all_results_sum$method == 'mcat_emp' & all_results_sum$results %in% c('True Negative', 'False Positive')]) * 100

all_results_sum$freq[all_results_sum$method == 'mcat_wilcox' & all_results_sum$results == 'True Positive']/
  sum(all_results_sum$freq[all_results_sum$method == 'mcat_wilcox' & all_results_sum$results %in% c('True Positive', 'False Negative')]) * 100
all_results_sum$freq[all_results_sum$method == 'mcat_wilcox' & all_results_sum$results == 'True Negative']/
  sum(all_results_sum$freq[all_results_sum$method == 'mcat_wilcox' & all_results_sum$results %in% c('True Negative', 'False Positive')]) * 100

