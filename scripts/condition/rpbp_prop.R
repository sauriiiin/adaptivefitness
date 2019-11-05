##### RPBP ANALYSIS
##### Aaron's results from the RiboSeq analysis
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 10/31/2019

##### INITIALIZE
library(seqinr)
# loading CAI
data(caitab)
w <- caitab$sc
#
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

source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

rpbp <- read.table('/home/acwach/Synteny/rpbps_table', header = T)
orfs <- read.table('/home/acwach/Synteny/orfs_comprehensive', header = T, fill = T, sep = ' ')
xtra <- read.table('/home/acwach/Synteny/additional_orf_info',header = T, fill = T, sep = ',')

orfs$id <- orfs$id + 9000000
rpbp$id <- orfs$id
rpbp <- rpbp[,c(86,1:85)]

# dbWriteTable(conn, "RPBP_TABLE", rpbp, overwrite = T)
# dbWriteTable(conn, "RPBP_ORFS_COMPREHENSIVE", orfs, overwrite = T)

##### COMBINING PROPERTIES WITH EXISTING DATABASE
# rpbp[,c(3:86)] <- rpbp[,c(3:86)] > 10
# rpbp[,c(3:86)] <- lapply(rpbp[,c(3:86)], as.numeric)
# rpbp$evidence <- rowSums(rpbp[,c(3:86)])

# pg.can <- orfs[orfs$id %in% rpbp$id[rpbp$evidence > 0],] 
# 
# pg.can$len_nt <- pg.can$end - pg.can$start + 1
# pg.can$len_aa <- pg.can$len_nt/3
# 
# pg.can <- pg.can[pg.can$is_gene == 'X',]
# pg.can <- pg.can[,c("id","focal_aaseq","focal_seq","len_nt","len_aa")]
# 
# colnames(pg.can) <- c('id','seq_aa','seq_nt','len_nt','len_aa')
# pg.can[,c(2:3)] <- lapply(pg.can[,c(2:3)], as.character)
# 
# pg.can$gc <- NULL
# pg.can$tiny <- NULL
# pg.can$small <- NULL
# pg.can$aliphatic <- NULL
# pg.can$aromatic <- NULL
# pg.can$nonpolar <- NULL
# pg.can$polar <- NULL
# pg.can$charged <- NULL
# pg.can$basic <- NULL
# pg.can$acidic <- NULL
# pg.can$pI <- NULL
# pg.can$instability <- NULL
# pg.can$cai <- NULL
# 
# for (i in 1:length(pg.can$id)) {
#   prop_temp <- AAstat(s2c(pg.can$seq_aa[i]), plot = F)
#   
#   pg.can$gc[i] <- GC(s2c(pg.can$seq_nt[i]))
#   pg.can$tiny[i] <- prop_temp$Prop$Tiny
#   pg.can$small[i] <- prop_temp$Prop$Small
#   pg.can$aliphatic[i] <- prop_temp$Prop$Aliphatic
#   pg.can$aromatic[i] <- prop_temp$Prop$Aromatic
#   pg.can$nonpolar[i] <- prop_temp$Prop$Non.polar
#   pg.can$polar[i] <- prop_temp$Prop$Polar
#   pg.can$charged[i] <- prop_temp$Prop$Charged
#   pg.can$basic[i] <- prop_temp$Prop$Basic
#   pg.can$acidic[i] <- prop_temp$Prop$Acidic
#   pg.can$pI[i] <- prop_temp$Pi
#   pg.can$instability[i] <- instaIndex(str_remove(pg.can$seq_aa[i], '[*]'))
#   pg.can$cai[i] <- cai(s2c(pg.can$seq_nt[i]), w = w)[[1]] 
# }
# 
# pg.can <- pg.can[,c(1,4:18)]
# pg.can$protogene <- 3

# load(sprintf('%spg1.RData',out_path))
# props <- rbind(melted[,c(2:17)], pg.can[,c(2:17)])
# save(props,file = sprintf('%spg12.RData',out_path))
load(sprintf('%spg12.RData',out_path))
props <- melt(props, id.vars = c('protogene'))

ggplot(props, aes(x = protogene, y = value)) +
  geom_violin() +
  stat_compare_means(comparisons = list(c("0","1"),c("0","2"),c("0","3"),c("1","2"),c("1","3"),c("2","3"))) +
  facet_wrap(~variable, scales = 'free_y',
             nrow = 5) +
  theme_linedraw()

#####
rpbp[,c(3:86)] <- rpbp[,c(3:86)] > 10
rpbp[,c(3:86)] <- lapply(rpbp[,c(3:86)], as.numeric)
rpbp$evidence <- rowSums(rpbp[,c(3:86)])
data <- cbind(orfs, evidence = rpbp$evidence)
data$protogene[data$is_gene == 'X'] <- 'Y'
data$protogene[is.na(data$protogene)] <- 'N'
data$orf_class <- as.character(data$orf_class)
data$orf_class[data$orf_class == ""] <- 'unannotated'

ggplot(data[data$evidence > 0,]) +
  geom_line(aes(x = evidence, col = orf_class), stat = 'density', trim = T) +
  theme_linedraw()


data$len_nt <- abs(data$end - data$start) + 1
data$len_aa <- data$len_nt/3

data[,c("focal_aaseq","focal_seq")] <- lapply(data[,c("focal_aaseq","focal_seq")], as.character)

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

for (i in 1:length(data$id)) {
  prop_temp <- AAstat(s2c(data$focal_aaseq[i]), plot = F)

  data$gc[i] <- GC(s2c(data$focal_seq[i]))
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
  if (str_length(str_remove(data$focal_aaseq[i], 'X')) > 2) {
    data$instability[i] <- instaIndex(str_remove(data$focal_aaseq[i], 'X'))
  } else {
    data$instability[i] <- NA
  }
  data$cai[i] <- cai(s2c(data$focal_seq[i]), w = w)[[1]]
}

data$grp_len_aa <- as.numeric(cut(data$len_aa, 10))
data$grp_evidence <- as.numeric(cut(data$evidence, 5))
data <- cbind(data, xtra[2:5])
# save(data,file = sprintf('%srpbp_all_props.RData',out_path))
write.csv(data, file = 'TRANSLATEOME_PROPS.csv', row.names = F)
dbWriteTable(conn, "RPBP_ORFS_PROPS", data, overwrite = T)

dat.plt <- data[,c("id","orf_class","evidence","len_nt","len_aa","gc","tiny",
                   "small","aliphatic","aromatic","nonpolar","polar","charged",
                   "basic","acidic","pI","instability","cai","grp_len_aa",
                   "grp_evidence","frame_conservation","coding_score","pair_diffs","diverse")]

# dat.plt <- dat.plt[dat.plt$orf_class != 'pseudogene' & dat.plt$orf_class != 'te',]
# dat.plt <- dat.plt[dat.plt$evidence > 0,]

# dat.plt$grp_len_aa <- as.numeric(cut(dat.plt$len_aa, 10))
# dat.plt$grp_evidence <- as.numeric(cut(dat.plt$evidence, 4))
dat.plt$grp_evidence <- NULL
dat.plt$grp_evidence[dat.plt$evidence == 0] <- 0
dat.plt$grp_evidence[dat.plt$evidence > 0 & dat.plt$evidence < 23] <- 1
dat.plt$grp_evidence[dat.plt$evidence > 22 & dat.plt$evidence < 45] <- 2
dat.plt$grp_evidence[dat.plt$evidence > 44 & dat.plt$evidence < 67] <- 3
dat.plt$grp_evidence[dat.plt$evidence > 66] <- 4

dat.plt$grp_fc <- NULL
dat.plt$grp_fc[dat.plt$frame_conservation <= 0.6] <- 1
dat.plt$grp_fc[dat.plt$frame_conservation > 0.8] <- 3
dat.plt$grp_fc[is.na(dat.plt$grp_fc)] <- 2
# dat.plt$orf_class <- ordered(dat.plt$orf_class, levels = c('unannotated','dubious','uncharacterized','te','pseudogene','verified'))
melted <- melt(dat.plt, id.vars = c("id","orf_class","evidence","grp_len_aa","grp_evidence"))

plt <- ggplot(melted,
       aes(x = grp_evidence, y = value, group = grp_evidence)) +
  geom_jitter(aes(col = orf_class), size = 0.08, alpha = 0.4) +
  geom_violin(fill = NA, size = 1.2,
              draw_quantiles = c(0.05, 0.5, 0.95), trim = T) +
  stat_summary(fun.y = median, geom = 'line', group = 'orf_class', col = 'blue') +
  stat_summary(fun.y = median, geom = 'point', group = 'orf_class', col = 'red') +
  # stat_compare_means(comparisons = list(c("0","1"),c("1","2"),c("0","2"))) +
  facet_wrap(~variable, scales = 'free_y',
             nrow = 5) +
  scale_x_continuous(name = 'ORF RiboSeq Evidence\n(no. of expt. with bayes factor > 10)',
                     breaks = 0:4,
                     labels = c('0','1-22','23-44','45-66','67-88')) +
  scale_color_discrete(name = 'SGD Class') +
  theme_linedraw() +
  theme(legend.position = 'bottom') +
  guides(color = guide_legend(override.aes = list(size=3, alpha = 1)),
         shape = guide_legend(override.aes = list(size=3, alpha = 1)))
ggsave(sprintf("%srpbp_prop.jpg",out_path),
       plt,
       width = 15, height = 15,
       dpi = 500)

evi2sgd <- ggplot(melted) +
  geom_bar(aes(x = evidence, fill = orf_class), stat = 'count', position = 'stack') +
  scale_y_continuous(trans = 'log10') +
  # coord_trans(y="log10") +
  labs(x = 'ORF RiboSeq Evidence\n(no. of expt. with bayes factor > 10)',
       y = 'Log10(count)') +
  scale_fill_discrete(name = 'SGD Class') +
  # annotation_logticks(sides="l") +
  theme_linedraw()
ggsave(sprintf("%sevi2sgd.jpg",out_path),
       evi2sgd,
       width = 10, height = 10,
       dpi = 300)

# hello <- ggplot(melted) +
#   geom_line(aes(x = evidence, col = orf_class), stat = 'density')
# ggsave(sprintf("%shello.jpg",out_path),
#        hello,
#        width = 10, height = 10,
#        dpi = 300)
##### CLUSTERING
# BiocManager::install("M3C")
library(M3C)