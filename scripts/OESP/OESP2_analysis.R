##### OESP2 Analysis
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 03/29/2020
##### Analyzing the second mini-screen using OE collection and LID

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
expt <- 'OESP2_OLD'

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 7
txt <- 5
lbls <- 9

##### PLATEMAPS
plate384 <- dbGetQuery(conn, 'select b.*, a.plate 384plate from OESP1_pos2coor a, OESP1_pos2coor b
                        where a.density = 384 and b.density = 384 and a.pos = b.pos
                        order by b.plate, b.col, b.row')

plate1536 <- dbGetQuery(conn, 'select b.*, a.plate 384plate from OESP1_pos2coor a, OESP1_pos2coor b
                    where a.density = 384 and b.density = 1536 and
                    (a.pos = b.pos - 100000
                    or a.pos = b.pos - 200000
                    or a.pos = b.pos - 300000
                    or a.pos = b.pos - 400000
                    or a.pos = b.pos - 500000
                    or a.pos = b.pos - 600000
                    or a.pos = b.pos - 700000
                    or a.pos = b.pos - 800000)
                    order by b.plate, b.col, b.row')

plate6144 <- dbGetQuery(conn, 'select b.*, a.plate 384plate from OESP1_pos2coor a, OESP1_pos2coor b
                        where a.density = 384 and b.density = 6144 and
                        (
                        a.pos = b.pos - 1100000
                        or a.pos = b.pos - 1200000
                        or a.pos = b.pos - 1300000
                        or a.pos = b.pos - 1400000
                        or a.pos = b.pos - 1500000
                        or a.pos = b.pos - 1600000
                        or a.pos = b.pos - 1700000
                        or a.pos = b.pos - 1800000
                        or a.pos = b.pos - 2100000
                        or a.pos = b.pos - 2200000
                        or a.pos = b.pos - 2300000
                        or a.pos = b.pos - 2400000
                        or a.pos = b.pos - 2500000
                        or a.pos = b.pos - 2600000
                        or a.pos = b.pos - 2700000
                        or a.pos = b.pos - 2800000
                        or a.pos = b.pos - 3100000
                        or a.pos = b.pos - 3200000
                        or a.pos = b.pos - 3300000
                        or a.pos = b.pos - 3400000
                        or a.pos = b.pos - 3500000
                        or a.pos = b.pos - 3600000
                        or a.pos = b.pos - 3700000
                        or a.pos = b.pos - 3800000
                        or a.pos = b.pos - 4100000
                        or a.pos = b.pos - 4200000
                        or a.pos = b.pos - 4300000
                        or a.pos = b.pos - 4400000
                        or a.pos = b.pos - 4500000
                        or a.pos = b.pos - 4600000
                        or a.pos = b.pos - 4700000
                        or a.pos = b.pos - 4800000
                        or a.pos = b.pos - 5100000
                        or a.pos = b.pos - 5200000
                        or a.pos = b.pos - 5300000
                        or a.pos = b.pos - 5400000
                        or a.pos = b.pos - 5500000
                        or a.pos = b.pos - 5600000
                        or a.pos = b.pos - 5700000
                        or a.pos = b.pos - 5800000
                        or a.pos = b.pos - 6100000
                        or a.pos = b.pos - 6200000
                        or a.pos = b.pos - 6300000
                        or a.pos = b.pos - 6400000
                        or a.pos = b.pos - 6500000
                        or a.pos = b.pos - 6600000
                        or a.pos = b.pos - 6700000
                        or a.pos = b.pos - 6800000
                        or a.pos = b.pos - 7100000
                        or a.pos = b.pos - 7200000
                        or a.pos = b.pos - 7300000
                        or a.pos = b.pos - 7400000
                        or a.pos = b.pos - 7500000
                        or a.pos = b.pos - 7600000
                        or a.pos = b.pos - 7700000
                        or a.pos = b.pos - 7800000
                        or a.pos = b.pos - 8100000
                        or a.pos = b.pos - 8200000
                        or a.pos = b.pos - 8300000
                        or a.pos = b.pos - 8400000
                        or a.pos = b.pos - 8500000
                        or a.pos = b.pos - 8600000
                        or a.pos = b.pos - 8700000
                        or a.pos = b.pos - 8800000
                        )
                        order by b.plate, b.col, b.row')

plates <- rbind(plate384, plate1536, plate6144)

ggplot(plates) +
  geom_tile(aes(x = col, y = row, fill = as.factor(`384plate`))) +
  scale_y_reverse() +
  scale_fill_manual(name = '384-Density\nPlate',
                    labels = c("1"="REF-1","2"="REF-2","3"="BF-11","4"="BF-22",
                               "5"="BF-23","6"="BF-25","7"="BF-27","8"="BF-28"),
                    values = c("1"="black","2"="black","3"="#388E3C","4"="#FFC107",
                               "5"="grey90","6"="#E040FB","7"="#FF5722","8"="#448AFF")) +
  facet_wrap(.~density * plate,
             scales = 'free',
             ncol = 4) +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = 'bottom',
        strip.text = element_text(size = txt,
                                 margin = margin(0.15,0,0.15,0, "mm")),
        legend.box.spacing = unit(0.5,"mm"))

ggsave(sprintf("%sOESP_PM.jpg",out_path),
       height = two.c, width = 160, units = 'mm',
       dpi = 600)

for (pl in unique(plates$plate[plates$density == 6144])) {
  ggplot(plates[plates$density == 6144 & plates$plate == pl,]) +
    geom_tile(aes(x = col, y = row, fill = as.factor(`384plate`)), col = 'black') +
    # scale_y_reverse() +
    scale_fill_manual(name = '384-Density\nPlate',
                      labels = c("1"="REF-1","2"="REF-2","3"="BF-11","4"="BF-22",
                                 "5"="BF-23","6"="BF-25","7"="BF-27","8"="BF-28"),
                      values = c("1"="black","2"="black","3"="#388E3C","4"="#FFC107",
                                 "5"="grey90","6"="#E040FB","7"="#FF5722","8"="#448AFF"),
                      guide = F) +
    # facet_wrap(.~density * plate,
    #            scales = 'free',
    #            ncol = 4) +
    theme_linedraw() +
    coord_cartesian(xlim = c(4.5,8.5),
                    ylim = c(8.5,4.5)) +
    theme(axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.title = element_text(size = titles),
          legend.text = element_text(size = txt),
          legend.key.size = unit(3, "mm"),
          legend.position = 'bottom',
          strip.text = element_text(size = txt,
                                    margin = margin(0.15,0,0.15,0, "mm")),
          legend.box.spacing = unit(0.5,"mm"))
  
  ggsave(sprintf("%sOESP_6144_%d_PM.jpg",out_path, pl),
         height = one.c, width = one.c, units = 'mm',
         dpi = 600)
}


##### LOOKING AT RAW DATA
data <- dbGetQuery(conn, 'select b.*, a.average
                   from OESP2_OLD_FS_6144_JPEG a, OESP1_pos2coor b
                   where a.hours = 24 and a.pos = b.pos
                   order by hours, plate, col, row')

ggplot(data) +
  geom_line(aes(x  = average), stat = 'density', trim = T) +
  facet_wrap(.~plate, nrow = 2)

ggplot(data[data$average >= 600,]) +
  geom_tile(aes(x = col, y = row, col = average)) +
  facet_wrap(.~plate, nrow = 2) 


##### COMPARING RESULTS BETWEEN OESP1 and OESP2
# OESP1 = 1419 (91, 66 del, 25 neut)
# OESP2 = 1336 (8, 6 del, 2 neut)

fdata <- dbGetQuery(conn, 'select a.strain_id, a.orf_name,
                    a.cs_mean oesp1_cs_mean, a.cs_median oesp1_cs_median,
                    b.cs_mean oesp2_cs_mean, b.cs_median oesp2_cs_median
                    from OESP1_FS_6144_FITNESS_STATS a, OESP2_FS_6144_FITNESS_STATS b
                    where a.strain_id = b.strain_id and a.hours = 22 and b.hours = 24')

# fdata$oesp1_cs_mean[is.na(fdata$oesp1_cs_mean) & !is.na(fdata$oesp2_cs_mean)] <- 0
# fdata$oesp2_cs_mean[is.na(fdata$oesp2_cs_mean) & !is.na(fdata$oesp1_cs_mean)] <- 0

ggplot(fdata,
       aes(x = oesp1_cs_mean, y = oesp2_cs_mean)) +
  geom_point() +
  geom_abline(col = 'red', linetype = 'dashed') +
  geom_point(col = '#303F9F', alpha = 0.7) +
  geom_density_2d(col = 'red', lwd = 0.1) +
  annotate("text", x = 0.25, y = 1,
           label = lm_eqn(data.frame(average = fdata$oesp1_cs_mean, bg = fdata$oesp2_cs_mean)),
           size = 2, parse = T,
           hjust = 0) +
  labs(x = 'Pilot #1',
       y = 'Pilot #2')
  coord_cartesian(xlim = c(0, 1.1),
                  ylim = c(0,1.1)) +
  theme_linedraw()



pdata <- dbGetQuery(conn, 'select a.strain_id, a.orf_name,
                    a.p oesp1_p, a.stat oesp1_stat, a.es oesp1_es,
                    b.p oesp2_p, b.stat oesp2_stat, b.es oesp2_es
                    from OESP1_FS_6144_PVALUE a, OESP2_FS_6144_PVALUE b 
                    where a.hours = 22 and b.hours = 24 
                    and a.strain_id = b.strain_id')

ggplot(pdata,
       aes(x = oesp1_es, y = oesp2_es)) +
  geom_point() +
  geom_abline(col = 'red', linetype = 'dashed') +
  geom_point(col = '#303F9F', alpha = 0.7) +
  geom_density_2d(col = 'red', lwd = 0.1) +
  annotate("text", x = -0.5, y = 0,
           label = lm_eqn(data.frame(average = pdata$oesp1_es, bg = pdata$oesp2_es)),
           size = 2, parse = T,
           hjust = 0)

ggplot(pdata) +
  geom_vline(xintercept = 0, col = 'red', linetype = 'dashed') +
  geom_hline(yintercept = 0, col = 'red', linetype = 'dashed') +
  geom_point(aes(x = oesp1_es, y = oesp2_es)) +
  coord_cartesian(xlim = c(-0.1,0.1),
                  ylim = c(-0.1,0.1)) +
  theme_linedraw()

pdata$oesp1_effect[pdata$oesp1_p <= 0.05 & pdata$oesp1_stat < 0] <- 'Deleterious'
pdata$oesp1_effect[pdata$oesp1_p <= 0.05 & pdata$oesp1_stat > 0] <- 'Beneficial'
pdata$oesp1_effect[is.na(pdata$oesp1_effect)] <- 'Neutral'
pdata$oesp2_effect[pdata$oesp2_p <= 0.05 & pdata$oesp2_stat < 0] <- 'Deleterious'
pdata$oesp2_effect[pdata$oesp2_p <= 0.05 & pdata$oesp2_stat > 0] <- 'Beneficial'
pdata$oesp2_effect[is.na(pdata$oesp2_effect)] <- 'Neutral'

pdata$effect <- paste(strtrim(pdata$oesp1_effect,3), strtrim(pdata$oesp2_effect,3), sep = '/')
pdata$es <- rowMeans(cbind(pdata$oesp1_es, pdata$oesp2_es), na.rm = T)

pie.per <- plyr::count(pdata, vars = c('oesp1_effect'))
pie.per$per <- pie.per$freq/sum(pie.per$freq) * 100


ggplot(data = pie.per, aes(x = "", y = per, fill = oesp1_effect)) +
  geom_bar(stat = 'identity') +
  coord_polar("y", start = 0) +
  geom_label_repel(aes(label = sprintf("%0.2f%%\n(%d)",per, freq)),
                   position = position_stack(vjust = 0.5),
                   label.size = 0.15,
                   box.padding = 0.4,
                   seed = 100,
                   show.legend = F) +
  scale_fill_discrete(name = 'Effect (P1/P2)') +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = 'bottom',
        strip.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm"))

ggsave(sprintf("%sOESP_EFFECTS.jpg",out_path),
       height = one.c, width = one.c, units = 'mm',
       dpi = 300)


ggplot(pdata,
       aes(x = effect, y = es)) +
  geom_violin()

##### OESP2 vs OESP2_OLD
pdata <- dbGetQuery(conn, 'select a.strain_id, a.orf_name,
                    a.p oesp1_p, a.stat oesp1_stat, a.es oesp1_es,
                    b.p oesp2_p, b.stat oesp2_stat, b.es oesp2_es
                    from OESP2_FS_6144_PVALUE a, OESP2_OLD_FS_6144_PVALUE b 
                    where a.hours = 24 and b.hours = 24 
                    and a.strain_id = b.strain_id')

# ggplot(pdata) +
#   geom_point(aes(x = oesp1_es, y = oesp2_es))
# 
# ggplot(pdata) +
#   geom_vline(xintercept = 0, col = 'red', linetype = 'dashed') +
#   geom_hline(yintercept = 0, col = 'red', linetype = 'dashed') +
#   geom_point(aes(x = oesp1_es, y = oesp2_es)) +
#   coord_cartesian(xlim = c(-0.1,0.1),
#                   ylim = c(-0.1,0.1)) +
#   theme_linedraw()

pdata$oesp1_effect[pdata$oesp1_p <= 0.05 & pdata$oesp1_stat < 0] <- 'Deleterious'
pdata$oesp1_effect[pdata$oesp1_p <= 0.05 & pdata$oesp1_stat > 0] <- 'Beneficial'
pdata$oesp1_effect[is.na(pdata$oesp1_effect)] <- 'Neutral'
pdata$oesp2_effect[pdata$oesp2_p <= 0.05 & pdata$oesp2_stat < 0] <- 'Deleterious'
pdata$oesp2_effect[pdata$oesp2_p <= 0.05 & pdata$oesp2_stat > 0] <- 'Beneficial'
pdata$oesp2_effect[is.na(pdata$oesp2_effect)] <- 'Neutral'

pdata$effect <- paste(strtrim(pdata$oesp1_effect,3), strtrim(pdata$oesp2_effect,3), sep = '/')

pie.per <- plyr::count(pdata, vars = c('effect'))
pie.per$per <- pie.per$freq/sum(pie.per$freq) * 100


ggplot(data = pie.per, aes(x = "", y = per, fill = effect)) +
  geom_bar(stat = 'identity') +
  coord_polar("y", start = 0) +
  geom_label_repel(aes(label = sprintf("%0.2f%%\n(%d)",per, freq)),
                   position = position_stack(vjust = 0.5),
                   label.size = 0.15,
                   box.padding = 0.4,
                   seed = 100,
                   show.legend = F) +
  scale_fill_discrete(name = 'Effect (P2/P1-O)') +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = 'bottom',
        strip.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm"))

ggsave(sprintf("%sOESP2_1O_EFFECTS.jpg",out_path),
       height = one.c, width = one.c, units = 'mm',
       dpi = 300)


##### OESP1 vs OESP2_OLD
pdata <- dbGetQuery(conn, 'select a.strain_id, a.orf_name,
                    a.p oesp1_p, a.stat oesp1_stat, a.es oesp1_es,
                    b.p oesp2_p, b.stat oesp2_stat, b.es oesp2_es
                    from OESP1_FS_6144_PVALUE a, OESP2_OLD_FS_6144_PVALUE b 
                    where a.hours = 22 and b.hours = 24 
                    and a.strain_id = b.strain_id')

# ggplot(pdata) +
#   geom_point(aes(x = oesp1_es, y = oesp2_es))
# 
# ggplot(pdata) +
#   geom_vline(xintercept = 0, col = 'red', linetype = 'dashed') +
#   geom_hline(yintercept = 0, col = 'red', linetype = 'dashed') +
#   geom_point(aes(x = oesp1_es, y = oesp2_es)) +
#   coord_cartesian(xlim = c(-0.1,0.1),
#                   ylim = c(-0.1,0.1)) +
#   theme_linedraw()

pdata$oesp1_effect[pdata$oesp1_p <= 0.05 & pdata$oesp1_stat < 0] <- 'Deleterious'
pdata$oesp1_effect[pdata$oesp1_p <= 0.05 & pdata$oesp1_stat > 0] <- 'Beneficial'
pdata$oesp1_effect[is.na(pdata$oesp1_effect)] <- 'Neutral'
pdata$oesp2_effect[pdata$oesp2_p <= 0.05 & pdata$oesp2_stat < 0] <- 'Deleterious'
pdata$oesp2_effect[pdata$oesp2_p <= 0.05 & pdata$oesp2_stat > 0] <- 'Beneficial'
pdata$oesp2_effect[is.na(pdata$oesp2_effect)] <- 'Neutral'

pdata$effect <- paste(strtrim(pdata$oesp1_effect,3), strtrim(pdata$oesp2_effect,3), sep = '/')

pie.per <- plyr::count(pdata, vars = c('effect'))
pie.per$per <- pie.per$freq/sum(pie.per$freq) * 100


ggplot(data = pie.per, aes(x = "", y = per, fill = effect)) +
  geom_bar(stat = 'identity') +
  coord_polar("y", start = 0) +
  geom_label_repel(aes(label = sprintf("%0.2f%%\n(%d)",per, freq)),
                   position = position_stack(vjust = 0.5),
                   label.size = 0.15,
                   box.padding = 0.4,
                   seed = 100,
                   show.legend = F) +
  scale_fill_discrete(name = 'Effect (P1/P1-O)') +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = 'bottom',
        strip.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm"))

ggsave(sprintf("%sOESP1_1O_EFFECTS.jpg",out_path),
       height = one.c, width = one.c, units = 'mm',
       dpi = 300)

##### OESP1, OESP2 & OESP2_OLD
pdata <- dbGetQuery(conn, 'select a.strain_id, a.orf_name,
                    a.p oesp1_p, a.stat oesp1_stat, a.es oesp1_es,
                    b.p oesp2_p, b.stat oesp2_stat, b.es oesp2_es,
                    c.p oesp2o_p, c.stat oesp2o_stat, c.es oesp2o_es
                    from OESP1_FS_6144_PVALUE a, OESP2_FS_6144_PVALUE b, OESP2_OLD_FS_6144_PVALUE c 
                    where a.hours = 22 and b.hours = 24 and c.hours = 24
                    and a.strain_id = b.strain_id and b.strain_id = c.strain_id')

# ggplot(pdata) +
#   geom_point(aes(x = oesp1_es, y = oesp2_es))
# 
# ggplot(pdata) +
#   geom_vline(xintercept = 0, col = 'red', linetype = 'dashed') +
#   geom_hline(yintercept = 0, col = 'red', linetype = 'dashed') +
#   geom_point(aes(x = oesp1_es, y = oesp2_es)) +
#   coord_cartesian(xlim = c(-0.1,0.1),
#                   ylim = c(-0.1,0.1)) +
#   theme_linedraw()

pdata$oesp1_effect[pdata$oesp1_p <= 0.05 & pdata$oesp1_stat < 0] <- 'Deleterious'
pdata$oesp1_effect[pdata$oesp1_p <= 0.05 & pdata$oesp1_stat > 0] <- 'Beneficial'
pdata$oesp1_effect[is.na(pdata$oesp1_effect)] <- 'Neutral'
pdata$oesp2_effect[pdata$oesp2_p <= 0.05 & pdata$oesp2_stat < 0] <- 'Deleterious'
pdata$oesp2_effect[pdata$oesp2_p <= 0.05 & pdata$oesp2_stat > 0] <- 'Beneficial'
pdata$oesp2_effect[is.na(pdata$oesp2_effect)] <- 'Neutral'
pdata$oesp2o_effect[pdata$oesp2o_p <= 0.05 & pdata$oesp2o_stat < 0] <- 'Deleterious'
pdata$oesp2o_effect[pdata$oesp2o_p <= 0.05 & pdata$oesp2o_stat > 0] <- 'Beneficial'
pdata$oesp2o_effect[is.na(pdata$oesp2o_effect)] <- 'Neutral'

pdata$effect <- paste(strtrim(pdata$oesp1_effect,3), strtrim(pdata$oesp2_effect,3), strtrim(pdata$oesp2o_effect,3), sep = '/')

pie.per <- plyr::count(pdata, vars = c('effect'))
pie.per$per <- pie.per$freq/sum(pie.per$freq) * 100
write.csv(pie.per, file = 'oesp_pie.csv')

ggplot(data = pie.per, aes(x = "", y = per, fill = effect)) +
  geom_bar(stat = 'identity') +
  coord_polar("y", start = 0) +
  geom_label_repel(aes(label = sprintf("%0.2f%%\n(%d)",per, freq)),
                   position = position_stack(vjust = 0.5),
                   label.size = 0.15,
                   box.padding = 0.4,
                   seed = 100,
                   show.legend = F) +
  scale_fill_manual(name = 'Effect (P1/P2/P3)',
                      values = c('Del/Del/Del' = '#FF5252',
                                 'Del/Neu/Del' = '#7C4DFF',
                                 'Neu/Del/Del' = '#448AFF',
                                 'Neu/Neu/Del' = '#388E3C',
                                 'Del/Del/Neu' = '#9C27B0',
                                 'Del/Neu/Neu' = '#303F9F',
                                 'Neu/Del/Neu' = '#00BCD4',
                                 'Neu/Neu/Neu' = '#00796B')) +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = 'bottom',
        strip.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm"))

ggsave(sprintf("%sOESP123_EFFECTS.jpg",out_path),
       height = one.c, width = one.c, units = 'mm',
       dpi = 300)

## 3D Scatter Plot of all the Effect Sizes
library(plotly)
plot_ly(x=pdata$oesp1_es, y=pdata$oesp2_es, z=pdata$oesp2o_es, type="scatter3d", mode="markers", color=pdata$effect)

## ES and Effect
pdata$es <- rowMeans(cbind(pdata$oesp1_es, pdata$oesp2_es, pdata$oesp2o_es), na.rm = T)

ggplot(pdata) +
  geom_point(aes(x = 1:dim(pdata)[1], y = es, col = effect))

pdata$effect2[pdata$oesp1_effect == 'Deleterious' &
                pdata$oesp2_effect == 'Deleterious' &
                pdata$oesp2o_effect == 'Deleterious'] <- 'Deleterious'
pdata$effect2[pdata$oesp1_effect == 'Neutral' &
                pdata$oesp2_effect == 'Neutral' &
                pdata$oesp2o_effect == 'Neutral'] <- 'Neutral'
pdata$effect2[is.na(pdata$effect2)] <- 'Variable'

ggplot(pdata) +
  geom_hline(yintercept = quantile(pdata$es[pdata$effect2 == 'Variable'], c(0.05,0.95)),
             col = 'red', linetype = 'dashed') +
  geom_violin(aes(x = effect2, y = es)) +
  scale_x_discrete(breaks = c('Deleterious','Variable','Neutral'),
                   limits = c('Deleterious','Variable','Neutral')) +
  scale_y_continuous(breaks = seq(-2,2,0.1),
                     minor_breaks = seq(-2,2,0.05)) +
  labs(x = 'Fitness Effect',
       y = 'Mean Effect Size') +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = 'bottom',
        strip.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm"))

ggsave(sprintf("%sOESP_ESvsEFFECTS.jpg",out_path),
       height = one.c, width = one.c, units = 'mm',
       dpi = 300)


pie.per <- plyr::count(pdata, vars = c('effect2'))
pie.per$per <- pie.per$freq/sum(pie.per$freq) * 100

pie.per$effect2 <- factor(pie.per$effect2, levels = c('Deleterious','Variable','Neutral'))

ggplot(data = pie.per, aes(x = "", y = per, fill = effect2)) +
  geom_bar(stat = 'identity') +
  coord_polar("y", start = 0) +
  geom_label_repel(aes(label = sprintf("%0.2f%%\n(%d)",per, freq)),
                   position = position_stack(vjust = 0.5),
                   label.size = 0.15,
                   box.padding = 0.4,
                   seed = 100,
                   show.legend = F) +
  scale_fill_manual(name = 'Effect',
                    values = c('Deleterious' = '#FF5252',
                               'Variable' = 'grey50',
                               'Neutral' = '#00796B')) +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = 'bottom',
        strip.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm"))

ggsave(sprintf("%sOESP123_EFFECTS2.jpg",out_path),
       height = one.c, width = one.c, units = 'mm',
       dpi = 300)
