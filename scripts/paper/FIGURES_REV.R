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

plyr::count(rnd2_data[(rnd2_data$lid_phenotype != rnd2_data$lid_wilcox_phenotype) &
                        (rnd2_data$phenotype == rnd2_data$lid_wilcox_phenotype),],
            vars = c('hours'))
plyr::count(rnd2_data[rnd2_data$phenotype != rnd2_data$lid_wilcox_phenotype,],
            vars = c('hours'))

plyr::count(rnd2_data[(rnd2_data$mcat_phenotype != rnd2_data$mcat_wilcox_phenotype) &
                        (rnd2_data$phenotype == rnd2_data$mcat_wilcox_phenotype),],
            vars = c('hours'))
plyr::count(rnd2_data[rnd2_data$phenotype != rnd2_data$mcat_wilcox_phenotype,],
            vars = c('hours'))

pheno_mismatch <- NULL
pheno_mismatch$hours <- sort(unique(rnd2_data$hours))
pheno_mismatch <- merge(pheno_mismatch,
                        plyr::count(rnd2_data[rnd2_data$phenotype != rnd2_data$lid_phenotype,],
                                    vars = c('hours')),
                        by = 'hours', all = T)
pheno_mismatch <- merge(pheno_mismatch,
                        plyr::count(rnd2_data[rnd2_data$phenotype != rnd2_data$lid_wilcox_phenotype,],
                                    vars = c('hours')),
                        by = 'hours', all = T)
pheno_mismatch <- merge(pheno_mismatch,
                        plyr::count(rnd2_data[rnd2_data$phenotype != rnd2_data$mcat_phenotype,],
                                    vars = c('hours')),
                        by = 'hours', all = T)
pheno_mismatch <- merge(pheno_mismatch,
                        plyr::count(rnd2_data[rnd2_data$phenotype != rnd2_data$mcat_wilcox_phenotype,],
                                    vars = c('hours')),
                        by = 'hours', all = T)

colnames(pheno_mismatch) <- c('hours','lid_emp','lid_wilcox','mcat_emp','mcat_wilcox')
pheno_mismatch

ggplot(melt(pheno_mismatch, id.vars = 'hours', variable.name = 'method', value.name = 'count'),
       aes(x = hours, y = count + 0.001)) +
  geom_line(aes(col = method)) +
  scale_y_log10()

######
ggplot(rnd2_data[rnd2_data$hours == 11.04,],
       aes(x = orf_name, y = fitness)) +
  geom_boxplot(aes(fill = lid_phenotype),
               outlier.shape = NA) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107'),
                    guide = F) +
  facet_wrap(.~hours)



