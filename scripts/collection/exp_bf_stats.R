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

source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

out_path = 'figs/collection/'

##### GATHER DATA
bf.ori <- dbGetQuery(conn, 'select * from COLLECTIONS.BARFLEX_SPACE_AGAR')
bf.exp <- dbGetQuery(conn, 'select * from COLLECTIONS.PROTOGENE_COLLECTION')
bf.all <- rbind(bf.ori, bf.exp)
colnames(bf.all) <- c('strain_id','plate','row','col','orf_name')

bf.all$expan[bf.all$plate > 18] <- 1
bf.all$expan[is.na(bf.all$expan)] <- 0

pgs <- dbGetQuery(conn, 'select * from PROTOGENES
                  where selected + longer + translated < 3')
seq <- dbGetQuery(conn, 'select * from ORFS_SEQUENCES')

bf.all$protogene[bf.all$orf_name %in% pgs$orf_name] <- 1
bf.all$protogene[is.na(bf.all$protogene)] <- 0

for (orf in bf.all$orf_name) {
  bf.all$seq_nt_atg[bf.all$orf_name == orf] = seq$seq_nt_atg[seq$orf_name == orf]
  bf.all$nt_len[bf.all$orf_name == orf] = nchar(seq$seq_nt_atg[seq$orf_name == orf])
}

ggplot(bf.all) +
  geom_line(aes(x = nt_len, col = as.character(protogene)),
            lwd = 1.2,
            stat = 'density', trim = T) +
  labs(title = 'ORF Seq NT Lengths',
       subtitle = 'Expanded BF Collection',
       x = 'NT Length',
       y = 'Density') +
  theme_linedraw() +
  scale_color_discrete(name = 'PG') +
  facet_wrap(.~plate, scales = 'free')
ggsave(sprintf('%sEXPANDED_BF_NT_LENS.png',out_path),
       width = 10, height = 10,
       dpi = 300)

ggplot(bf.all,
       aes(x = protogene, y = "", fill = as.character(protogene))) +
  geom_bar(stat = 'identity') +
  coord_polar(theta = "y", start = 0) +
  labs(title = 'Proto-gene Proportion',
       subtitle = 'Expanded BF Collection') +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank()) +
  scale_fill_discrete(name = 'PG') +
  facet_wrap(.~plate)
ggsave(sprintf('%sEXPANDED_BF_NT_PROP.png',out_path),
       width = 10, height = 10,
       dpi = 300)

