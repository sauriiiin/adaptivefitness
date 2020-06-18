##### OESP1 Analysis
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 02/11/2020
##### Analyzing the first mini-screen using OE collection and LID

##### INITIALIZE
library(RMariaDB)
library(ggplot2)
library(ggExtra)
library(gridExtra)
library(ggpubr)
library(stringr)
library(reshape2)
library(ggrepel)
source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

out_path <- 'figs/OESP1/'
expt <- 'OESP1'

pgs <- dbGetQuery(conn, 'select * from PROTOGENES
                  where selected + longer + translated < 3')
seq <- dbGetQuery(conn, 'select * from ORFS_SEQUENCES')

# bf.exp1 <- dbGetQuery(conn, 'select * from BARFLEX_SPACE_AGAR')
# bf.exp2 <- dbGetQuery(conn, 'select * from COLLECTIONS.PROTOGENE_COLLECTION')
# bf.all <- rbind(bf.exp1, bf.exp2)
# colnames(bf.all) <- c('strain_id','plate','row','col','orf_name')
# 
# for (orf in bf.all$orf_name[!is.na(bf.all$orf_name) & bf.all$orf_name != 'BF_control']) {
#   bf.all$seq_nt_atg[bf.all$orf_name == orf] = seq$seq_nt_atg[seq$orf_name == orf]
#   bf.all$nt_len[bf.all$orf_name == orf] = nchar(seq$seq_nt_atg[seq$orf_name == orf])
#   
# }
# save(bf.all, file = sprintf('%sbfall.RData',out_path))

load(sprintf('%sbfall.RData',out_path))

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 7
txt <- 5
lbls <- 9

##### PLATEMAPS
# platemap <- dbGetQuery(conn, sprintf('select b.strain_id, a.*
#                       from %s_pos2coor a, %s_pos2strainid b
#                       where a.pos = b.pos',expt,expt))
# platemap$colony[platemap$strain_id == 0]  <- 'NULL'
# platemap$colony[platemap$strain_id == -1] <- 'Reference'
# platemap$colony[platemap$strain_id == -2] <- 'Border'
# platemap$colony[is.na(platemap$colony)]  <- 'Query'
# 
# pl.map <- ggplot(platemap, aes(x = col, y = row, col = colony)) +
#   geom_point(size = 0.5) +
#   scale_y_continuous(trans = 'reverse') +
#   scale_color_discrete(name = 'Colony',
#                        breaks = c('Border','Reference','Query','NULL')) +
#   theme_linedraw() +
#   theme(axis.title = element_text(size = titles),
#         axis.text = element_text(size = txt),
#         axis.text.x = element_text(angle = 45, hjust = 1),
#         strip.text = element_text(size = txt,
#                                   margin = margin(0.5,0,0.5,0, "mm")),
#         legend.position = 'bottom') +
#   guides(color = guide_legend(override.aes = list(size=4)),
#          shape = guide_legend(override.aes = list(size=4))) +
#   facet_wrap(.~density*plate,
#              ncol = 8,
#              scales = 'free')
# ggsave(sprintf("%s%s_platemaps.jpg",out_path,expt), pl.map,
#        width = two.c, height = one.c, unit = "mm",
#        dpi = 600)


##### FINAL SCREEN
stage <- 'FS'
density <- 6144

fitstats <- dbGetQuery(conn, sprintf('select a.*, b.p, b.stat, b.es
                                   from %s_%s_%d_FITNESS_STATS a, %s_%s_%d_PVALUE b
                                   where a.hours = b.hours and a.orf_name = b.orf_name',
                                   expt,stage,density,
                                   expt,stage,density))

fitstats$protogene[fitstats$orf_name %in% pgs$orf_name] <- 1
fitstats$protogene[is.na(fitstats$protogene)] <- 0
fitstats$protogene[str_detect(fitstats$orf_name,'smorf')] <- 2
fitstats$protogene <- as.factor(fitstats$protogene)

fitstats$effect[fitstats$p <= 0.05 & fitstats$stat > 0] <- 'Beneficial'
fitstats$effect[fitstats$p <= 0.05 & fitstats$stat < 0] <- 'Deleterious'
fitstats$effect[is.na(fitstats$effect)] <- 'Neutral'

fitstats$nt_len <- NULL
fitstats$plate <- NULL
for (sid in fitstats$strain_id) {
  fitstats$nt_len[fitstats$strain_id == sid] <- mean(bf.all$nt_len[bf.all$strain_id == sid])
  fitstats$plate[fitstats$strain_id == sid] <- max(bf.all$plate[bf.all$strain_id == sid])
}

fit.plt <- ggplot(fitstats[fitstats$hours == 22,]) +
  geom_line(aes(x = cs_mean, col = protogene), stat = 'density', lwd = 1.2,
            trim = T) +
  scale_color_discrete(name = 'ORF Class',
                       breaks = c(0,1,2),
                       labels = c('Gene','Y-PG','S-PG')) +
  labs(title = 'Time = 22 hours',
       x = 'Fitness',
       y = 'Density') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        strip.text = element_text(size = txt,
                                  margin = margin(0.5,0,0.5,0, "mm")),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'right') +
  guides(color = guide_legend(override.aes = list(size=4)),
         shape = guide_legend(override.aes = list(size=4)))
# ggsave(sprintf("%s%s_fitness.jpg",out_path,expt), fit.plt,
#        width = one.c, height = one.c, unit = "mm",
#        dpi = 600)

fitpie <- plyr::count(fitstats, vars = c('hours','protogene','effect'))

for (h in fitpie$hours) {
  for (p in fitpie$protogene[fitpie$hours == h]) {
    fitpie$total[fitpie$hours == h & fitpie$protogene == p] <-
      sum(fitpie$freq[fitpie$hours == h & fitpie$protogene == p])
  }
}

fitpie$name[fitpie$protogene == 0] <- 'Gene'
fitpie$name[fitpie$protogene == 1] <- 'Y-PG'
fitpie$name[fitpie$protogene == 2] <- 'S-PG'

fitpie$name <- factor(fitpie$name, levels=c('Gene','Y-PG','S-PG'))

effect.pie <- ggplot(data = fitpie[fitpie$hours > 0,], aes(x = "", y = freq/total, fill = effect)) +
  geom_bar(stat = 'identity') +
  coord_polar("y", start = 0) +
  geom_label_repel(aes(label = sprintf("%0.2f %%",freq/total * 100)),
                   col = 'white',
                   position = position_stack(vjust = 0.5),
                   label.size = 0.15,
                   show.legend = F,
                   fontface = 'bold',
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.5, "lines"),
                   segment.color = 'grey50',
                   seed = 10) +
  scale_fill_manual(name = 'Effects',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107')) +
  facet_wrap(.~hours*name,
             ncol = 3) +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.text =  element_blank(),
        strip.text = element_text(size = txt,
                                  margin = margin(0.5,0,0.5,0, "mm")),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'right')
ggsave(sprintf("%s%s_effectpie.jpg",out_path,expt), effect.pie,
       width = two.c, height = two.c, unit = "mm",
       dpi = 600)
  
  
# plyr::count(fitstats[fitstats$hours == 17,], vars = c('protogene'))

##### WHERE ARE THE COLONIES COMING FROM?
effect.pg.plt <- plyr::count(fitstats, vars = c('hours','plate','protogene','effect'))
for (h in effect.pg.plt$hours) {
  for (pl in effect.pg.plt$plate[effect.pg.plt$hours == h])
    # for (p in effect.pg.plt$protogene[effect.pg.plt$hours == h & effect.pg.plt$plate == pl]) {
      effect.pg.plt$total[effect.pg.plt$hours == h & effect.pg.plt$plate == pl] <-
        sum(effect.pg.plt$freq[effect.pg.plt$hours == h & effect.pg.plt$plate == pl])
    # }
}

effect.pg.plt$name[effect.pg.plt$protogene == 0] <- 'Gene'
effect.pg.plt$name[effect.pg.plt$protogene == 1] <- 'Y-PG'
effect.pg.plt$name[effect.pg.plt$protogene == 2] <- 'S-PG'

effect.pg.plt$name <- factor(effect.pg.plt$name, levels=c('Gene','Y-PG','S-PG'))

efc.pp <- ggplot(effect.pg.plt[effect.pg.plt$hours > 0,]) +
  geom_bar(aes(x = as.factor(plate), y = freq/total, fill = effect, alpha = name),
           stat="identity") +
  scale_fill_manual(name = 'Effects',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107')) +
  scale_alpha_manual(name = 'ORF',
                    breaks = c('Gene','Y-PG','S-PG'),
                    values = c('Gene'=1,
                               'Y-PG'=0.7,
                               'S-PG'=0.4)) +
  facet_wrap(.~hours, ncol = 3) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text =  element_text(size = txt),
        strip.text = element_text(size = txt,
                                  margin = margin(0.5,0,0.5,0, "mm")),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'right')
ggsave(sprintf("%s%s_effectpgplt.jpg",out_path,expt), efc.pp,
       width = two.c, height = two.c, unit = "mm",
       dpi = 300)

effect.plt <- plyr::count(fitstats, vars = c('hours','plate','effect'))
for (h in effect.plt$hours) {
  for (pl in effect.plt$plate[effect.plt$hours == h])
    effect.plt$total[effect.plt$hours == h & effect.plt$plate == pl] <-
      sum(effect.plt$freq[effect.plt$hours == h & effect.plt$plate == pl])
}

efc.p <- ggplot(effect.plt[effect.plt$hours > 0,]) +
  geom_bar(aes(x = as.factor(plate), y = freq/total*100, fill = effect),
           stat="identity") +
  scale_fill_manual(name = 'Effects',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107')) +
  facet_wrap(.~hours, ncol = 3) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text =  element_text(size = txt),
        strip.text = element_text(size = txt,
                                  margin = margin(0.5,0,0.5,0, "mm")),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'right')
ggsave(sprintf("%s%s_effectplt.jpg",out_path,expt), efc.p,
       width = two.c, height = two.c, unit = "mm",
       dpi = 300)

plyr::count(fitstats, vars = c('hours','nt_len','effect'))

ggplot(fitstats[fitstats$hours > 0,]) +
  geom_histogram(aes(x = nt_len, y = (..count..)/sum(..count..) * 100, fill = effect), binwidth = 100) +
  scale_fill_manual(name = 'Effects',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107')) +
  # scale_y_continuous(trans = 'log') +
  facet_grid(.~hours) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text =  element_text(size = txt),
        strip.text = element_text(size = txt,
                                  margin = margin(0.5,0,0.5,0, "mm")),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'right')
ggsave(sprintf("%s%s_effectlen.jpg",out_path,expt),
       width = two.c, height = one.c, unit = "mm",
       dpi = 300)
