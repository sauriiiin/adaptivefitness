##### EXPANDED BARFLEX COLLECTION
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 01/16/2020

##### INITIALIZE
library(RMariaDB)
library(readxl)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(reshape2)
library(grid)
library(tidyverse)
library(stringr)

source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

out_path = 'figs/collection/'

##### GATHER DATA
bf.ori <- dbGetQuery(conn, 'select * from COLLECTIONS.BARFLEX_SPACE_AGAR')
bf.exp22 <- dbGetQuery(conn, 'select * from BARFLEX_SPACE_AGAR where 384plate = 22')
bf.exp <- dbGetQuery(conn, 'select * from COLLECTIONS.PROTOGENE_COLLECTION')
bf.all <- rbind(bf.ori, bf.exp, bf.exp22)
colnames(bf.all) <- c('strain_id','plate','row','col','orf_name')

bf.all$expan[bf.all$plate > 18] <- 1
bf.all$expan[is.na(bf.all$expan)] <- 0

pgs <- dbGetQuery(conn, 'select * from PROTOGENES
                  where selected + longer + translated < 3')
seq <- dbGetQuery(conn, 'select * from ORFS_SEQUENCES')

bf.all$protogene[bf.all$orf_name %in% pgs$orf_name] <- 1
bf.all$protogene[is.na(bf.all$protogene)] <- 0
bf.all$protogene[str_detect(bf.all$orf_name,'smorf')] <- 2
# 6405 ORFs in total - 1840 PGS - 1121 smorfs

for (orf in bf.all$orf_name) {
  bf.all$seq_nt_atg[bf.all$orf_name == orf] = seq$seq_nt_atg[seq$orf_name == orf]
  bf.all$nt_len[bf.all$orf_name == orf] = nchar(seq$seq_nt_atg[seq$orf_name == orf])
  
}

##### VISUALIZING DATA
## ORF SEQ NT LENGHTS
orf.lengths <- ggplot(bf.all) +
  geom_line(aes(x = nt_len, col = as.character(protogene)),
            lwd = 1.2,
            stat = 'density', trim = T) +
  scale_x_continuous(trans = 'log10') +
  labs(title = 'ORF Seq NT Lengths',
       subtitle = 'Expanded BF Collection',
       x = 'log10(NT Length)',
       y = 'Density') +
  theme_linedraw() +
  annotation_logticks(sides = 'b') +
  scale_color_discrete(name = 'ORF',
                       breaks = c(0,1,2),
                       labels = c('Gene','Y-PG','S-PG')) +
  facet_wrap(.~plate, scales = 'free_y',
             nrow = 5, ncol = 5)
ggsave(sprintf('%sEXPANDED_BF_NT_LENS.png',out_path), orf.lengths,
       width = 10, height = 10,
       dpi = 600)

## PG PROPORTIONS / PLATE
pie.dat <- plyr::count(bf.all, vars = c('plate','protogene'))
freq <- plyr::count(bf.all, vars = c('plate'))

for (p in freq$plate) {
 pie.dat$freq[pie.dat$plate == p] <- pie.dat$freq[pie.dat$plate == p]/freq$freq[freq$plate == p] * 100
}

pg.prop <- ggplot(pie.dat,
       aes(x = "", y = freq, fill = as.character(protogene))) +
  geom_bar(stat = 'identity') +
  coord_polar(theta = "y", start = 0) +
  labs(title = 'Proto-gene Proportion',
       subtitle = 'Expanded BF Collection') +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank()) +
  scale_fill_discrete(name = 'ORF',
                       breaks = c(0,1,2),
                       labels = c('Gene','Y-PG','S-PG')) +
  facet_wrap(.~plate,
             nrow = 5, ncol = 5)
ggsave(sprintf('%sEXPANDED_BF_NT_PROP.png',out_path), pg.prop,
       width = 10, height = 10,
       dpi = 600)

## PLATE MAPS
plate.map <- ggplot(bf.all) +
  geom_point(aes(x = col, y = row, col = as.character(protogene)), size = 2) +
  labs(title = 'Plate Maps',
       subtitle = 'Expanded BF Collection',
       x = 'Columns',
       y = 'Rows') +
  scale_x_continuous(breaks = seq(0,24,2),limits = c(1,24)) +
  scale_y_continuous(breaks = seq(0,26,2),limits = c(16,1),trans = 'reverse') +
  theme_linedraw() +
  scale_color_discrete(name = 'ORF',
                       breaks = c(0,1,2),
                       labels = c('Gene','Y-PG','S-PG')) +
  facet_wrap(.~plate,
             nrow = 5, ncol = 5)
ggsave(sprintf('%sEXPANDED_BF_MAPS.png',out_path), plate.map,
       width = 13, height = 10,
       dpi = 600)

plate.map2 <- ggplot(bf.all) +
  geom_point(aes(x = col, y = row, col = as.character(protogene)), size = 2) +
  geom_text(aes(x = col, y = row, label = orf_name, angle = 30), size = 0.6) +
  labs(title = 'Plate Maps',
       subtitle = 'Expanded BF Collection',
       x = 'Columns',
       y = 'Rows') +
  scale_x_continuous(breaks = seq(0,24,2),limits = c(1,24)) +
  scale_y_continuous(breaks = seq(0,26,2),limits = c(16,1),trans = 'reverse') +
  theme_linedraw() +
  scale_color_discrete(name = 'ORF',
                       breaks = c(0,1,2),
                       labels = c('Gene','Y-PG','S-PG')) +
  facet_wrap(.~plate,
             nrow = 5, ncol = 5)
ggsave(sprintf('%sEXPANDED_BF_MAPS_ORFS.png',out_path), plate.map2,
       width = 13, height = 10,
       dpi = 1000)
