##### SAN DIEGO ANALYSIS
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 09/10/2019

##### INITIALIZE
library(RMariaDB)
library(ggplot2)
library(ggExtra)
library(gridExtra)
library(ggpubr)
source("R/functions/initialize.sql.R")

##### GET DATA
conn <- initialize.sql("saurin_test")
bf_met <- dbGetQuery(conn, "select * from brian_031918.`FITNESS_v2_BF1-MET_SC_screen`
                      where exp_id = 30 and hours = 68 and fitness is not NULL
                      order by orf_name")
bf_metmet <- dbGetQuery(conn, "select * from brian_031918.`FITNESS_v2_BF1-MET_SD_screen`
                      where exp_id = 39 and hours = 68 and fitness is not NULL
                      order by orf_name")

ko_met <- dbGetQuery(conn, "select * from brian_031918.`FITNESS_v2_KO1-MET_SC_screen`
                      where exp_id = 8 and hours = 84 and fitness is not NULL
                      order by orf_name")
ko_metmet <- dbGetQuery(conn, "select * from brian_031918.`FITNESS_v2_KO1-MET_SD_screen`
                      where exp_id = 20 and hours = 68 and fitness is not NULL
                      order by orf_name")

##### GROWTH IN DIFFERENT MEDIA
ggplot() +
  geom_line(data = bf_met[bf_met$orf_name != 'BF_control',],
            aes(x = fitness, col = 'Cases'), stat = 'density',
            lwd = 1.2) +
  geom_line(data = bf_met[bf_met$orf_name == 'BF_control',],
            aes(x = fitness, col = 'Control'), stat = 'density',
            lwd = 1.2)

ggplot() +
  geom_line(data = bf_metmet[bf_metmet$orf_name != 'BF_control',],
            aes(x = fitness, col = 'Cases'), stat = 'density',
            lwd = 1.2) +
  geom_line(data = bf_metmet[bf_metmet$orf_name == 'BF_control',],
            aes(x = fitness, col = 'Control'), stat = 'density',
            lwd = 1.2)


ggplot() +
  geom_line(data = ko_met[ko_met$orf_name != 'KO_control',],
            aes(x = fitness, col = 'Cases'), stat = 'density',
            lwd = 1.2) +
  geom_line(data = ko_met[ko_met$orf_name == 'KO_control',],
            aes(x = fitness, col = 'Control'), stat = 'density',
            lwd = 1.2)

ggplot() +
  geom_line(data = ko_metmet[ko_metmet$orf_name != 'KO_control',],
            aes(x = fitness, col = 'Cases'), stat = 'density',
            lwd = 1.2) +
  geom_line(data = ko_metmet[ko_metmet$orf_name == 'KO_control',],
            aes(x = fitness, col = 'Control'), stat = 'density',
            lwd = 1.2)

##### CHANGE IN REF FITNESS WITH PINS
pin1 <- dbGetQuery(conn, 'select * from brian_031918.FITNESS_v2_BF3_pin1
              where exp_id in (50,52,54,56)
              order by exp_id, orf_name')
pin1$colony[pin1$orf_name == 'BF_control' & !is.na(pin1$orf_name) & pin1$orf_name != 'null'] <- 'BFC'
pin1$colony[pin1$orf_name != 'BF_control' & !is.na(pin1$orf_name) & pin1$orf_name != 'null'] <- 'MUT'

pin2 <- dbGetQuery(conn, 'select * from brian_031918.FITNESS_v2_BF3_pin2
              where exp_id in (50,52,54,56)
              order by exp_id, orf_name')
pin2$colony[pin2$orf_name == 'BF_control' & !is.na(pin2$orf_name) & pin2$orf_name != 'null'] <- 'BFC'
pin2$colony[pin2$orf_name != 'BF_control' & !is.na(pin2$orf_name) & pin2$orf_name != 'null'] <- 'MUT'

pin3 <- dbGetQuery(conn, 'select * from brian_031918.FITNESS_v2_BF3_pin3
              where exp_id in (50,52,54,56)
              order by exp_id, orf_name')
pin3$colony[pin3$orf_name == 'BF_control' & !is.na(pin3$orf_name) & pin3$orf_name != 'null'] <- 'BFC'
pin3$colony[pin3$orf_name != 'BF_control' & !is.na(pin3$orf_name) & pin3$orf_name != 'null'] <- 'MUT'

pin4 <- dbGetQuery(conn, 'select * from brian_031918.FITNESS_v2_BF3_pin4
              where exp_id in (50,52,54,56)
              order by exp_id, orf_name')
pin4$colony[pin4$orf_name == 'BF_control' & !is.na(pin4$orf_name) & pin4$orf_name != 'null'] <- 'BFC'
pin4$colony[pin4$orf_name != 'BF_control' & !is.na(pin4$orf_name) & pin4$orf_name != 'null'] <- 'MUT'

ggplot(pin2[!is.na(pin2$colony),]) +
  geom_line(aes(x = average, col = colony), stat = 'density')

strain.fit <- ggplot() +
  geom_line(data = pin1[pin1$colony == 'MUT',],
            aes(x = fitness, col = 'Pin1'), stat = 'density', lwd = 1.2) +
  geom_line(data = pin2[pin2$colony == 'MUT',],
            aes(x = fitness, col = 'Pin2'), stat = 'density', lwd = 1.2) +
  geom_line(data = pin3[pin3$colony == 'MUT',],
            aes(x = fitness, col = 'Pin3'), stat = 'density', lwd = 1.2) +
  geom_line(data = pin4[pin4$colony == 'MUT',],
            aes(x = fitness, col = 'Pin4'), stat = 'density', lwd = 1.2) +
  labs(x = 'Fitness', y = 'Density',
       title = '', subtitle = 'Fitness Distribution') +
  theme_linedraw() +
  scale_color_discrete(name = 'Pin Level') +
  coord_cartesian(xlim = c(0,2.0))

strain.raw <- ggplot() +
  geom_line(data = pin1[pin1$colony == 'MUT',],
            aes(x = average, col = 'Pin1'), stat = 'density', lwd = 1.2) +
  geom_line(data = pin2[pin2$colony == 'MUT',],
            aes(x = average, col = 'Pin2'), stat = 'density', lwd = 1.2) +
  geom_line(data = pin3[pin3$colony == 'MUT',],
            aes(x = average, col = 'Pin3'), stat = 'density', lwd = 1.2) +
  geom_line(data = pin4[pin4$colony == 'MUT',],
            aes(x = average, col = 'Pin4'), stat = 'density', lwd = 1.2) +
  labs(x = 'Pix Count', y = 'Density',
       title = 'Strains', subtitle = 'Raw Distribution') +
  theme_linedraw() +
  scale_color_discrete(name = 'Pin Level') +
  coord_cartesian(xlim = c(0,1000))


ref.fit <- ggplot() +
  geom_line(data = pin1[pin1$colony == 'BFC',],
            aes(x = fitness, col = 'Pin1'), stat = 'density', lwd = 1.2) +
  geom_line(data = pin2[pin2$colony == 'BFC',],
            aes(x = fitness, col = 'Pin2'), stat = 'density', lwd = 1.2) +
  geom_line(data = pin3[pin3$colony == 'BFC',],
            aes(x = fitness, col = 'Pin3'), stat = 'density', lwd = 1.2) +
  geom_line(data = pin4[pin4$colony == 'BFC',],
            aes(x = fitness, col = 'Pin4'), stat = 'density', lwd = 1.2) +
  labs(x = 'Fitness', y = 'Density',
       title = '', subtitle = 'Fitness Distribution') +
  theme_linedraw() +
  scale_color_discrete(name = 'Pin Level') +
  coord_cartesian(xlim = c(0,2.0))

ref.raw <- ggplot() +
  geom_line(data = pin1[pin1$colony == 'BFC',],
            aes(x = average, col = 'Pin1'), stat = 'density', lwd = 1.2) +
  geom_line(data = pin2[pin2$colony == 'BFC',],
            aes(x = average, col = 'Pin2'), stat = 'density', lwd = 1.2) +
  geom_line(data = pin3[pin3$colony == 'BFC',],
            aes(x = average, col = 'Pin3'), stat = 'density', lwd = 1.2) +
  geom_line(data = pin4[pin4$colony == 'BFC',],
            aes(x = average, col = 'Pin4'), stat = 'density', lwd = 1.2) +
  labs(x = 'Pix Count', y = 'Density',
       title = 'Reference', subtitle = 'Raw Distribution') +
  theme_linedraw() +
  scale_color_discrete(name = 'Pin Level') +
  coord_cartesian(xlim = c(0,1000))

ggarrange(strain.raw, strain.fit, ref.raw, ref.fit,
          nrow = 2, ncol = 2,
          common.legend = T, legend = 'bottom')


ref.raw2 <- ggplot() +
  geom_line(data = pin1[pin1$colony == 'BFC',],
            aes(x = average - median(pin1$average[pin1$colony == 'BFC'], na.rm = T), col = 'Pin1'), stat = 'density', lwd = 1.2) +
  geom_line(data = pin2[pin2$colony == 'BFC',],
            aes(x = average - median(pin2$average[pin2$colony == 'BFC'], na.rm = T), col = 'Pin2'), stat = 'density', lwd = 1.2) +
  geom_line(data = pin3[pin3$colony == 'BFC',],
            aes(x = average - median(pin3$average[pin3$colony == 'BFC'], na.rm = T), col = 'Pin3'), stat = 'density', lwd = 1.2) +
  geom_line(data = pin4[pin4$colony == 'BFC',],
            aes(x = average - median(pin4$average[pin4$colony == 'BFC'], na.rm = T), col = 'Pin4'), stat = 'density', lwd = 1.2) +
  labs(x = 'Pix Count', y = 'Density',
       title = 'Reference', subtitle = 'Raw Distribution') +
  theme_linedraw() +
  scale_color_discrete(name = 'Pin Level') +
  coord_cartesian(xlim = c(-600,600))

strain.raw2 <- ggplot() +
  geom_line(data = pin1[pin1$colony == 'MUT',],
            aes(x = average - median(pin1$average[pin1$colony == 'BFC'], na.rm = T), col = 'Pin1'), stat = 'density', lwd = 1.2) +
  geom_line(data = pin2[pin2$colony == 'MUT',],
            aes(x = average - median(pin2$average[pin2$colony == 'BFC'], na.rm = T), col = 'Pin2'), stat = 'density', lwd = 1.2) +
  geom_line(data = pin3[pin3$colony == 'MUT',],
            aes(x = average - median(pin3$average[pin3$colony == 'BFC'], na.rm = T), col = 'Pin3'), stat = 'density', lwd = 1.2) +
  geom_line(data = pin4[pin4$colony == 'MUT',],
            aes(x = average - median(pin4$average[pin4$colony == 'BFC'], na.rm = T), col = 'Pin4'), stat = 'density', lwd = 1.2) +
  labs(x = 'Pix Count', y = 'Density',
       title = 'Strains', subtitle = 'Raw Distribution') +
  theme_linedraw() +
  scale_color_discrete(name = 'Pin Level') +
  coord_cartesian(xlim = c(-600,600))

ggarrange(strain.raw2, ref.raw2,
          common.legend = T, legend = 'bottom')


sd(pin1$average[pin1$colony == 'BFC'], na.rm = T)
sd(pin2$average[pin2$colony == 'BFC'], na.rm = T)
sd(pin3$average[pin3$colony == 'BFC'], na.rm = T)
sd(pin4$average[pin4$colony == 'BFC'], na.rm = T)

##### QC FOR ANOVA USING 4C3 DATA
fourC3 <- dbGetQuery(conn, 'select b.*, a.average
          from 4C3_GA_1536_JPEG a, 4C3_pos2coor1536 b
          where a.pos = b.pos and average is not NULL')

colnames(fourC3) <- c('pos','plate','row','col','average')
fourC3$plate <- as.character(fourC3$plate)

ggplot(fourC3) +
  geom_line(aes(x = average, col = plate), stat = 'density')

anova(lm(average ~ plate, fourC3))
lm(plate ~ average, fourC3)

wilcox.test(average ~ plate, data = fourC3,
            subset = plate %in% c('1','2'))

