##### MODEL FOR BENEFICIENCE
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 12/19/2019

##### INITIALIZE
library(RMariaDB)
library(ggplot2)
library(ggExtra)
library(gridExtra)
library(ggpubr)
library(reshape2)
source("R/functions/initialize.sql.R")

conn <- initialize.sql("saurin_test")

data <- dbGetQuery(conn, 'select id, orf_name, start, end, is_gene, focal_aaseq, focal_seq, orf_class, evidence, len_nt, len_aa, gc, tiny, small, aliphatic, aromatic, nonpolar, polar, charged, basic, acidic, pI, instability, cai, frame_conservation, coding_score, pair_diffs, diverse, colony_size, effect_cs, protogene
from TRANSLATEOME_PROPS a, brian_031918.DATASET_6 b
where a.is_gene = b.orf_name
and b.exp_id = 28')

dat2 <- data[,c(2,8:31)]
dat2 <- dat2[dat2$evidence > 0,]
dat2 <- melt(dat2, id.vars = c('orf_name','orf_class','effect_cs','protogene','colony_size'))

ggplot(dat2) +
  geom_violin(aes(x = effect_cs, y = value, fill = protogene)) +
  facet_wrap(~variable, scale = 'free')

# data[data$protogene == 'proto-gene' & data$effect_cs == 'beneficial',]
