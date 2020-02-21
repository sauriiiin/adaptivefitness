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
library(reshape2)
source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

out_path <- 'figs/OESP1/'
expt <- 'OESP1'

pgs <- dbGetQuery(conn, 'select * from PROTOGENES
                  where selected + longer + translated < 3')
seq <- dbGetQuery(conn, 'select * from ORFS_SEQUENCES')

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 7
txt <- 5
lbls <- 9

##### PLATEMAPS
platemap <- dbGetQuery(conn, sprintf('select b.strain_id, a.*
                      from %s_pos2coor a, %s_pos2strainid b
                      where a.pos = b.pos',expt,expt))
platemap$colony[platemap$strain_id == 0]  <- 'NULL'
platemap$colony[platemap$strain_id == -1] <- 'Reference'
platemap$colony[platemap$strain_id == -2] <- 'Border'
platemap$colony[is.na(platemap$colony)]  <- 'Query'

pl.map <- ggplot(platemap, aes(x = col, y = row, col = colony)) +
  geom_point(size = 0.5) +
  scale_y_continuous(trans = 'reverse') +
  scale_color_discrete(name = 'Colony',
                       breaks = c('Border','Reference','Query','NULL')) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = txt,
                                  margin = margin(0.5,0,0.5,0, "mm")),
        legend.position = 'bottom') +
  guides(color = guide_legend(override.aes = list(size=4)),
         shape = guide_legend(override.aes = list(size=4))) +
  facet_wrap(.~density*plate,
             ncol = 8,
             scales = 'free')
ggsave(sprintf("%s%s_platemaps.jpg",out_path,expt), pl.map,
       width = two.c, height = one.c, unit = "mm",
       dpi = 600)


##### FINAL SCREEN
stage <- 'FS'
density <- 6144

fitdat <- dbGetQuery(conn, sprintf('select a.*, b.p, b.stat, b.es
                                   from %s_%s_%d_FITNESS_STATS a, %s_%s_%d_PVALUE b
                                   where a.hours = b.hours and a.orf_name = b.orf_name',
                                   expt,stage,density,
                                   expt,stage,density))

fitdat$protogene[fitdat$orf_name %in% pgs$orf_name] <- 1
fitdat$protogene[is.na(fitdat$protogene)] <- 0
fitdat$protogene[str_detect(fitdat$orf_name,'smorf')] <- 2
fitdat$protogene <- as.factor(fitdat$protogene)

fitdat$effect[fitdat$p <= 0.05 & fitdat$stat > 0] <- 'Beneficial'
fitdat$effect[fitdat$p <= 0.05 & fitdat$stat < 0] <- 'Deleterious'
fitdat$effect[is.na(fitdat$effect)] <- 'Neutral'

fit.plt <- ggplot(fitdat[fitdat$hours == 22,]) +
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
ggsave(sprintf("%s%s_fitness.jpg",out_path,expt), fit.plt,
       width = one.c, height = one.c, unit = "mm",
       dpi = 600)

fitpie <- plyr::count(fitdat, vars = c('hours','protogene','effect'))

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
  
  
plyr::count(fitdat[fitdat$hours == 17,], vars = c('protogene'))

