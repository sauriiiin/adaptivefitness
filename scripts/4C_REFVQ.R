##### INITIALIZE
library(RMariaDB)
library(ggplot2)
library(ggExtra)
library(gridExtra)
library(ggpubr)
source("R/functions/initialize.sql.R")

##### GET/SET DATA
expt_name = '4C3_GA3'
expt = 'FS1-GA3'
out_path = 'figs/';
density = 6144;

##### CHECK POSITION WISE VARIABILITY
conn <- initialize.sql("saurin_test")
tablename_fit = sprintf('%s_%d_FITNESS',expt_name,density);

alldat <- dbGetQuery(conn, sprintf('select * from %s', tablename_fit))

alldat$colony[alldat$orf_name == 'BF_control'] = 'Reference'
alldat$colony[alldat$orf_name != 'BF_control'] = 'Query'

hr = 18
med <- quantile(alldat$average[alldat$hours == hr], 0.5, na.rm = T)[[1]]
lim.low <- med - 4*sd(alldat$average[alldat$hours == hr], na.rm = T)
lim.hig <- med + 4*sd(alldat$average[alldat$hours == hr], na.rm = T)

es.fit <- mean(alldat$fitness[alldat$hours == hr & alldat$colony == 'Query'], na.rm = T)/
  mean(alldat$fitness[alldat$hours == hr & alldat$colony == 'Reference'], na.rm = T)
es.pix <- mean(alldat$average[alldat$hours == hr & alldat$colony == 'Query'], na.rm = T)/
  mean(alldat$average[alldat$hours == hr & alldat$colony == 'Reference'], na.rm = T)

f <- ggplot(alldat[alldat$hours == hr,]) +
  # geom_histogram(aes(x = fitness, fill = colony), bins = 40, alpha = 0.7) +
  geom_line(aes(x = fitness, col = colony), stat = 'density', lwd = 1.2) +
  coord_cartesian(xlim = c(0.7,1.3)) +
  labs(subtitle = sprintf('ES = %0.3f', es.fit),
       x = 'Fitness',
       y = 'Density') +
  scale_color_manual(name = 'Strain',
                     breaks = c('Reference', 'Query'),
                     values = c('Reference' = '#303F9F',
                                'Query' = '#FF5252'),
                     labels = c('Ref','Query'),
                     drop = T) +
  theme_linedraw()

cs <- ggplot(alldat[alldat$hours == hr,]) +
  # geom_histogram(aes(x = average, fill = colony), bins = 30, alpha = 0.7) +
  geom_line(aes(x = average, col = colony), stat = 'density', lwd = 1.2) +
  coord_cartesian(xlim = c(lim.low,lim.hig)) +
  labs(subtitle = sprintf('ES = %0.3f', es.pix),
       x = 'Colony Size (pix)',
       y = 'Density') +
  scale_color_manual(name = 'Strain',
                     breaks = c('Reference', 'Query'),
                     values = c('Reference' = '#303F9F',
                                'Query' = '#FF5252'),
                     labels = c('Ref','Query'),
                     drop = T) +
  theme_linedraw()

annotate_figure(ggarrange(cs, f,
                          align = 'v',
                          common.legend = T, legend = 'right'),
                top = expt)
ggsave(sprintf('%s%s_REFVQ_DIS.png',out_path,expt_name),
       height = 7, width = 14)

#####
dbDisconnect(conn)
