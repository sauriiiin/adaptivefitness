library(RMariaDB)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggExtra)
library(grid)
library(tidyverse)
library(ggpubr)
library(stringr)
source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

out_path = 'figs/collection/';

seq.dat <- dbGetQuery(conn, 'select * from COLLECTION')
seq.dat$protogene <- as.character(seq.dat$protogene)

ggplot(seq.dat) +
  geom_line(aes(x = length_aa, col = protogene), stat = 'density', lwd = 1.3, trim = T) +
  scale_x_continuous(minor = seq(0,10000,50), trans = 'log10') + 
  scale_color_discrete(name = 'ORF',
                       breaks = c('0','1','2'),
                       labels = c('0' = 'Annotated\nGene\n','1' = 'Annotated\nProto-gene\n','2' = 'Unannotated\nProto-gene\n')) +
  annotation_logticks(sides = 'b') +
  labs(title = 'Amino-acid seq length distribution',
       subtitle = 'of expanded BarFLEX collection',
       x = 'Log10(Length)',
       y = 'Density') +
  theme_linedraw() +
  guides(color = guide_legend(override.aes = list(size=4)),
         shape = guide_legend(override.aes = list(size=4)))
ggsave(sprintf("%slength_aa.jpg",out_path),
       width = 8, height = 7,
       dpi = 300)



prop.dat <- dbGetQuery(conn, 'select a.orf_name, b.aromaticity, b.stability, b.pI, a.length_nt, a.length_aa, a.protogene
      from COLLECTION a, NT_PROPERTIES b
      where a.orf_name = b.orf_name')
prop.dat$protogene <- as.character(prop.dat$protogene)

ggplot(prop.dat) +
  geom_line(aes(x = pI, col = protogene), stat = 'density', lwd = 1.3, trim = T) +
  scale_x_continuous(minor = seq(0,20,0.5)) +
  scale_color_discrete(name = 'ORF',
                       breaks = c('0','1','2'),
                       labels = c('0' = 'Annotated\nGene\n','1' = 'Annotated\nProto-gene\n','2' = 'Unannotated\nProto-gene\n')) +
  # annotation_logticks(sides = 'b') +
  labs(title = 'Isoelectric Point distribution',
       subtitle = 'of expanded BarFLEX collection',
       x = 'pH',
       y = 'Density') +
  theme_linedraw() +
  guides(color = guide_legend(override.aes = list(size=4)),
         shape = guide_legend(override.aes = list(size=4)))
ggsave(sprintf("%spI_aa.jpg",out_path),
       width = 8, height = 7,
       dpi = 300)


hel.dat <- dbGetQuery(conn, 'select * from COLLECTION a, NT_TMHMM b
      where a.orf_name = b.orf_name')
hel.dat$protogene <- as.character(hel.dat$protogene)

ggplot(hel.dat[hel.dat$predhel < 6,]) +
  geom_violin(aes(x = protogene, y = predhel, fill = protogene)) +
  scale_x_discrete(labels = c('0' = 'Annotated\nGene\n','1' = 'Annotated\nProto-gene\n','2' = 'Unannotated\nProto-gene\n')) +
  scale_fill_discrete(name = 'ORF',
                       breaks = c('0','1','2'),
                       labels = c('0' = 'Annotated\nGene\n','1' = 'Annotated\nProto-gene\n','2' = 'Unannotated\nProto-gene\n'),
                      guide = F) +
  labs(title = 'TMHMM Predicted TM Domains',
       subtitle = 'in expanded BarFLEX collection',
       x = 'ORF',
       y = 'Density') +
  theme_linedraw() +
  coord_cartesian(ylim = c(0,5))
ggsave(sprintf("%shel_aa.jpg",out_path),
       width = 7, height = 7,
       dpi = 300)
