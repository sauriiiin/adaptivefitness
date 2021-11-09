##### F28FUR EMPIRICAL PVALUES AND SIG FITNESS CHANGES
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 03/29/2021

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

source("R/functions/empirical_p.R")
source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")
out_path = 'figs/f28fu/repeat/'

load("/home/sbp29/R/Projects/adaptivefitness/figs/f28fu/repeat/f28fu_data.RData")

##### FIGURE SIZE
one.c <- 90 #single column 
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 7
txt <- 7
lbls <- 9

#####
data <- foa.data
ref <- 'BF_control'
value <- 'fitness_clean'
p <- 'emp_p'
es <- 'effsize'
avoid <- c(ref,'BOR','REF','YHR021W-A')
fpr <- c(0.05)

data[,p] <- NULL
data[,es] <- NULL
data[,sprintf('%s_thresh',p)] <- NULL

for (c in unique(data$condition)) {
  for (d in unique(data$density[data$condition == c])) {
    for (h in sort(unique(data$hours[data$condition == c & data$density == d]))) {
      ref_fit <- data[data$condition == c & data$density == d & data$hours == h & data$orf_name == ref,] %>%
        group_by(rep, replicate) %>%
        dplyr::summarise(median = median(fitness_clean, na.rm = T), .groups = "keep") %>%
        data.frame()
      ref_fit <- ref_fit$median[!is.na(ref_fit$median)]
      
      for (o in sort(unique(data$orf_name[!(data$orf_name %in% avoid) & !is.na(data$orf_name)]))) {
        orf_fit <- median(data[data$condition == c & data$density == d & data$hours == h & data$orf_name == o, value], na.rm = T)
        es_val <- (orf_fit - mean(ref_fit, na.rm = T))/mean(ref_fit, na.rm = T) * 100
        
        p_val1 <- sum(orf_fit < ref_fit, na.rm = T)/length(ref_fit)
        p_val2 <- sum(orf_fit > ref_fit, na.rm = T)/length(ref_fit)
        p_val <- min(c(p_val1, p_val2))*2
        data[data$condition == c & data$density == d & data$hours == h & data$orf_name == o & !is.na(data$orf_name), p] <- p_val
        data[data$condition == c & data$density == d & data$hours == h & data$orf_name == o & !is.na(data$orf_name), es] <- es_val
      }
      data[data$condition == c & data$density == d & data$hours == h, sprintf('%s_thresh',p)] <- 
        quantile(data[data$condition == c & data$density == d & data$hours == h, p], probs = fpr, na.rm = T)
    } 
  }
}
head(data)
save(data, file = sprintf('%sf28fur_fs_data.RData',out_path))

data.sum <- data[!(data$orf_name %in% c('BOR','REF','BF_control','YHR021W-A')) & !is.na(data$orf_name),] %>%
  group_by(condition, density, replicate, hours, orf_name) %>%
  summarise(fitness = mean(fitness, na.rm = T),
            emp_p = mean(emp_p, na.rm = T),
            effsize = mean(effsize, na.rm = T),
            p_thresh = mean(emp_p_thresh, na.rm = T),
            .groups = "keep") %>%
  data.frame()

hi <- data.sum[data.sum$emp_p <= 0.05 &
           ((data.sum$density == 6144 & data.sum$hours == 36) |
           (data.sum$density == 1536 & data.sum$hours == 80)),]

length(ref_fit)

##### 5FOA FITNESS (S3)

##### PS2 PHENOTYPE (PS2)

##### PAIRWISE COMPARISIONS
# using median fitness values
# F28FUR S3, PS2, FS & SDPG in all densities

