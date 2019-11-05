##### SEQUENCE PROPERTIES
##### proto-gene and gene sequence properties that might set them apart
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 10/30/2019

##### INITIALIZE
library(seqinr)
library(Peptides)
# Piptide::instaIndex() calculates the instability index proposed by Guruprasad (1990).
# This index predicts the stability of a protein based on its amino acid composition,
# a protein whose instability index is smaller than 40 is predicted as stable, a value
# above 40 predicts that the protein may be unstable.
library(RMariaDB)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggExtra)
library(grid)
library(tidyverse)
library(ggpubr)
library(stringr)
library(reshape2)
source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")
out_path = 'figs/seq_prop/'

data(caitab)
w <- caitab$sc
# codon adaptation index for sc

##### GETTING SEQ DATA
data <- dbGetQuery(conn, "select *, length(seq_nt_atg) len_nt, length(seq_aa) len_aa
              from ORFS_SEQUENCES
              where orf_name in
              (select distinct orf_name from BARFLEX_BACTERIAL_MAP
              where orf_name is not NULL and orf_name != 'BF_control')")

pgs <- dbGetQuery(conn, "select orf_name from PROTOGENES
                  where selected + translated + longer < 3")

data$protogene <- NULL
data$protogene[data$orf_name %in% pgs$orf_name] <- 1
data$protogene[is.na(data$protogene)] <- 0
data$protogene[grepl('[s]',data$orf_name)] <- 2
data$protogene <- as.character(data$protogene)

data$len_nt <- as.numeric(data$len_nt)
data$len_aa <- as.numeric(data$len_aa)

data$gc <- NULL
data$tiny <- NULL
data$small <- NULL
data$aliphatic <- NULL
data$aromatic <- NULL
data$nonpolar <- NULL
data$polar <- NULL
data$charged <- NULL
data$basic <- NULL
data$acidic <- NULL
data$pI <- NULL
data$instability <- NULL
data$cai <- NULL

for (i in 1:length(data$orf_name)) {
  prop_temp <- AAstat(s2c(data$seq_aa[i]), plot = F)
  
  data$gc[i] <- GC(s2c(data$seq_nt[i]))
  data$tiny[i] <- prop_temp$Prop$Tiny
  data$small[i] <- prop_temp$Prop$Small
  data$aliphatic[i] <- prop_temp$Prop$Aliphatic
  data$aromatic[i] <- prop_temp$Prop$Aromatic
  data$nonpolar[i] <- prop_temp$Prop$Non.polar
  data$polar[i] <- prop_temp$Prop$Polar
  data$charged[i] <- prop_temp$Prop$Charged
  data$basic[i] <- prop_temp$Prop$Basic
  data$acidic[i] <- prop_temp$Prop$Acidic
  data$pI[i] <- prop_temp$Pi
  data$instability[i] <- instaIndex(str_remove(data$seq_aa[i], '[*]'))
  data$cai[i] <- cai(s2c(data$seq_nt[i]), w = w)[[1]] 
}

melted <- data[c(1,5:20)]
save(melted,file = sprintf('%spg1.RData',out_path))
melted <- melt(melted, id.vars = c('orf_name','protogene'))

ggplot(melted, aes(x = protogene, y = value)) +
  geom_violin() +
  stat_compare_means(comparisons = list(c("0","1"),c("1","2"),c("0","2"))) +
  facet_wrap(~variable, scales = 'free_y',
             nrow = 5) +
  theme_linedraw()
