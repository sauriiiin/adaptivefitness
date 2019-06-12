##### LID PAPER FIGURES
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 05/18/2019

##### INITIALIZE
library(RMariaDB)
library(readxl)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(grid)
library(tidyverse)
library(egg)
source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

##### INITIALIZE
expt_name = '4C3_GA1'
expt = 'FS1-1'
out_path = 'figs/lid_paper/';
density = 6144;

tablename_pval = sprintf('%s_%d_PVALUE',expt_name,density)
hours = dbGetQuery(conn, sprintf('select distinct hours from %s order by hours asc', tablename_pval))
pvals = seq(0,1,0.01)

tablename_fit = sprintf('%s_%d_FITNESS',expt_name,density);
tablename_rfit = sprintf('%s_RAW_%d_FITNESS',expt_name,density);
tablename_nfit = sprintf('%s_NIL_%d_FITNESS',substr(expt_name,1,6),density);
tablename_p2o = '4C3_pos2orf_name1';
tablename_bpos = '4C3_borderpos';

p2c_info = NULL
p2c_info[1] = '4C3_pos2coor6144'
p2c_info[2] = '6144plate'
p2c_info[3] = '6144col'
p2c_info[4] = '6144row'

p2c = dbGetQuery(conn, sprintf('select * from %s a order by a.%s, a.%s, a.%s',
                               p2c_info[1],
                               p2c_info[2],
                               p2c_info[3],
                               p2c_info[4]))

n_plates = dbGetQuery(conn, sprintf('select distinct %s from %s a order by %s asc',
                                    p2c_info[2],
                                    p2c_info[1],
                                    p2c_info[2]))
hr = hours[[1]][9]
pl = n_plates[[1]][1]

fitdat = dbGetQuery(conn, sprintf('select c.*, a.orf_name, a.hours, a.bg, a.average, a.fitness,
                                  b.bg nbg, b.average naverage, b.fitness nfitness
                                  from %s a, %s b, %s c
                                  where a.hours = %d and a.hours = b.hours
                                  and a.pos = b.pos and b.pos = c.pos
                                  and c.%s = %d order by c.%s, c.%s',
                                  tablename_fit,tablename_nfit,
                                  p2c_info[1],hr,p2c_info[2],
                                  pl,p2c_info[3],p2c_info[4]))
fitdat$bg[is.na(fitdat$average)] = NA

fitdat$source[fitdat$`6144row`%%2==1 & fitdat$`6144col`%%2==1] = 'TL'
fitdat$source[fitdat$`6144row`%%2==0 & fitdat$`6144col`%%2==1] = 'BL'
fitdat$source[fitdat$`6144row`%%2==1 & fitdat$`6144col`%%2==0] = 'TR'
fitdat$source[fitdat$`6144row`%%2==0 & fitdat$`6144col`%%2==0] = 'BR'

fitdat$colony[fitdat$orf_name == 'BF_control'] = 'Reference'
fitdat$colony[fitdat$orf_name != 'BF_control'] = 'Query'
fitdat$colony[is.na(fitdat$orf_name)] = 'Gap'

rfitdat = dbGetQuery(conn, sprintf('select c.*, a.orf_name, a.hours, a.bg, a.average, a.fitness,
                                  b.bg nbg, b.average naverage, b.fitness nfitness
                                  from %s a, %s b, %s c
                                  where a.hours = %d and a.hours = b.hours
                                  and a.pos = b.pos and b.pos = c.pos
                                  and c.%s = %d order by c.%s, c.%s',
                                  tablename_rfit,tablename_nfit,
                                  p2c_info[1],hr,p2c_info[2],
                                  pl,p2c_info[3],p2c_info[4]))
rfitdat$bg[is.na(rfitdat$average)] = NA

rfitdat$source[rfitdat$`6144row`%%2==1 & rfitdat$`6144col`%%2==1] = 'TL'
rfitdat$source[rfitdat$`6144row`%%2==0 & rfitdat$`6144col`%%2==1] = 'BL'
rfitdat$source[rfitdat$`6144row`%%2==1 & rfitdat$`6144col`%%2==0] = 'TR'
rfitdat$source[rfitdat$`6144row`%%2==0 & rfitdat$`6144col`%%2==0] = 'BR'

rfitdat$colony[rfitdat$orf_name == 'BF_control'] = 'Reference'
rfitdat$colony[rfitdat$orf_name != 'BF_control'] = 'Query'
rfitdat$colony[is.na(rfitdat$orf_name)] = 'Gap'

p2c_info_384 = NULL
p2c_info_384[1] = '4C3_pos2coor384'
p2c_info_384[2] = '384plate'
p2c_info_384[3] = '384col'
p2c_info_384[4] = '384row'

lay.384 <- dbGetQuery(conn, sprintf('select * from %s a, %s b
                                    where a.pos = b.pos
                                    order by b.%s, b.%s, b.%s',
                                    tablename_p2o, p2c_info_384[1],
                                    p2c_info_384[2], p2c_info_384[3], p2c_info_384[4]))

lay.384$colony[lay.384$orf_name == 'BF_control'] = 'Reference'
lay.384$colony[lay.384$orf_name != 'BF_control'] = 'Query'
lay.384$colony[is.na(lay.384$orf_name)] = 'Gap'

stk.lay.384 <- lay.384
stk.lay.384$colony[stk.lay.384$colony != 'Gap'] = 'Stock'

p2c_info_1536 = NULL
p2c_info_1536[1] = '4C3_pos2coor1536'
p2c_info_1536[2] = '1536plate'
p2c_info_1536[3] = '1536col'
p2c_info_1536[4] = '1536row'

lay.1536 <- dbGetQuery(conn, sprintf('select * from %s a, %s b
                                    where a.pos = b.pos
                                    order by b.%s, b.%s, b.%s',
                                     tablename_p2o, p2c_info_1536[1],
                                     p2c_info_1536[2], p2c_info_1536[3], p2c_info_1536[4]))

lay.1536$colony[lay.1536$orf_name == 'BF_control'] = 'Reference'
lay.1536$colony[lay.1536$orf_name != 'BF_control'] = 'Query'
lay.1536$colony[is.na(lay.1536$orf_name)] = 'Gap'

lay.6144 <- dbGetQuery(conn, sprintf('select * from %s a, %s b
                                    where a.pos = b.pos
                                    order by b.%s, b.%s, b.%s',
                                     tablename_p2o, p2c_info[1],
                                     p2c_info[2], p2c_info[3], p2c_info[4]))

lay.6144$colony[lay.6144$orf_name == 'BF_control'] = 'Reference'
lay.6144$colony[lay.6144$orf_name != 'BF_control'] = 'Query'
lay.6144$colony[is.na(lay.6144$orf_name)] = 'Gap'

##### FIGURE SIZE
max.w <- 20
max.h <- 25
ttle <- 1
stle <- 10/12

##### TEXT SIZES
plt.ttle <- 12
plt.stle <- 10
lax.ttle <- 10
lax.txt <- 9

##### SHAPE SIZE
sz.384 <- 5
sz.1536 <- 3
sz.6144 <- 1.5

##### FIGURE 1
col.lim <- 12
row.lim <- 8

gly1 <- ggplot(lay.384[lay.384$`384plate` == 1,]) +
  geom_point(aes(x = `384col`, y = `384row`, col = colony),
             shape = 15,
             size = sz.384) +
  scale_color_manual(name = 'Colony Type',
                     breaks = c('Stock','Reference', 'Query', 'Gap'),
                     values = c("Stock" = "#689F38",
                                "Reference" = "#FFC107",
                                "Query" = "#303F9F",
                                "Gap" = "#FF5252"),
                     guide = F) +
  labs(title = "",
       subtitle = "") +
  scale_x_continuous(breaks = seq(1,col.lim,1),limits = c(1,col.lim)) +
  scale_y_continuous(breaks = seq(1,row.lim,1),limits = c(row.lim,1),trans = 'reverse') +
  theme_linedraw() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_blank(),
        legend.title = element_blank()) +
  coord_cartesian(xlim = c(1,col.lim),
                  ylim = c(row.lim,1))

gly2 <- ggplot(lay.384[lay.384$`384plate` == 2,]) +
  geom_point(aes(x = `384col`, y = `384row`, col = colony),
             shape = 15,
             size = sz.384) +
  scale_color_manual(name = 'Colony Type',
                     breaks = c('Stock','Reference', 'Query', 'Gap'),
                     values = c("Stock" = "#689F38",
                                "Reference" = "#FFC107",
                                "Query" = "#303F9F",
                                "Gap" = "#FF5252"),
                     guide = F) +
  labs(title = "",
       subtitle = "") +
  scale_x_continuous(breaks = seq(1,col.lim,1),limits = c(1,col.lim)) +
  scale_y_continuous(breaks = seq(1,row.lim,1),limits = c(row.lim,1),trans = 'reverse') +
  theme_linedraw() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_blank(),
        legend.title = element_blank()) +
  coord_cartesian(xlim = c(1,col.lim),
                  ylim = c(row.lim,1))

gly3 <- ggplot(lay.384[lay.384$`384plate` == 3,]) +
  geom_point(aes(x = `384col`, y = `384row`, col = colony),
             shape = 15,
             size = sz.384) +
  scale_color_manual(name = 'Colony Type',
                     breaks = c('Stock','Reference', 'Query', 'Gap'),
                     values = c("Stock" = "#689F38",
                                "Reference" = "#FFC107",
                                "Query" = "#303F9F",
                                "Gap" = "#FF5252"),
                     guide = F) +
  labs(title = "",
       subtitle = "") +
  scale_x_continuous(breaks = seq(1,col.lim,1),limits = c(1,col.lim)) +
  scale_y_continuous(breaks = seq(1,row.lim,1),limits = c(row.lim,1),trans = 'reverse') +
  theme_linedraw() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_blank(),
        legend.title = element_blank()) +
  coord_cartesian(xlim = c(1,col.lim),
                  ylim = c(row.lim,1))

gly4 <- ggplot(lay.384[lay.384$`384plate` == 4,]) +
  geom_point(aes(x = `384col`, y = `384row`, col = colony),
             shape = 15,
             size = sz.384) +
  scale_color_manual(name = 'Colony Type',
                     breaks = c('Stock','Reference', 'Query', 'Gap'),
                     values = c("Stock" = "#689F38",
                                "Reference" = "#FFC107",
                                "Query" = "#303F9F",
                                "Gap" = "#FF5252"),
                     guide = F) +
  labs(title = "",
       subtitle = "") +
  scale_x_continuous(breaks = seq(1,col.lim,1),limits = c(1,col.lim)) +
  scale_y_continuous(breaks = seq(1,row.lim,1),limits = c(row.lim,1),trans = 'reverse') +
  theme_linedraw() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_blank(),
        legend.title = element_blank()) +
  coord_cartesian(xlim = c(1,col.lim),
                  ylim = c(row.lim,1))

# fig1.s1 <- ggarrange(gly1, gly2, gly3, gly4,
#                      nrow = 1)
# ggsave(sprintf("%sfigure1_s1.jpg",out_path),
#        fig1.s1,
#        width = max.w, height = max.w/4 * 3/4,
#        dpi = 300)

wc1 <- ggplot(lay.384[lay.384$`384plate` == 1,]) +
  geom_point(aes(x = `384col`, y = `384row`, col = colony),
             shape = 19,
             size = sz.384) +
  scale_color_manual(name = 'Colony Type',
                     breaks = c('Reference', 'Query', 'Gap'),
                     values = c("Reference" = "#FFC107",
                                "Query" = "#303F9F",
                                "Gap" = "#FF5252"),
                     guide = F) +
  labs(title = "",
       subtitle = "") +
  scale_x_continuous(breaks = seq(1,col.lim,1),limits = c(1,col.lim)) +
  scale_y_continuous(breaks = seq(1,row.lim,1),limits = c(row.lim,1),trans = 'reverse') +
  theme_linedraw() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_blank(),
        legend.title = element_blank()) +
  coord_cartesian(xlim = c(1,col.lim),
                  ylim = c(row.lim,1))

wc2 <- ggplot(lay.384[lay.384$`384plate` == 2,]) +
  geom_point(aes(x = `384col`, y = `384row`, col = colony),
             shape = 19,
             size = sz.384) +
  scale_color_manual(name = 'Colony Type',
                     breaks = c('Reference', 'Query', 'Gap'),
                     values = c("Reference" = "#FFC107",
                                "Query" = "#303F9F",
                                "Gap" = "#FF5252"),
                     guide = F) +
  labs(title = "",
       subtitle = "") +
  scale_x_continuous(breaks = seq(1,col.lim,1),limits = c(1,col.lim)) +
  scale_y_continuous(breaks = seq(1,row.lim,1),limits = c(row.lim,1),trans = 'reverse') +
  theme_linedraw() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_blank(),
        legend.title = element_blank()) +
  coord_cartesian(xlim = c(1,col.lim),
                  ylim = c(row.lim,1))

wc3 <- ggplot(lay.384[lay.384$`384plate` == 3,]) +
  geom_point(aes(x = `384col`, y = `384row`, col = colony),
             shape = 19,
             size = sz.384) +
  scale_color_manual(name = 'Colony Type',
                     breaks = c('Reference', 'Query', 'Gap'),
                     values = c("Reference" = "#FFC107",
                                "Query" = "#303F9F",
                                "Gap" = "#FF5252"),
                     guide = F) +
  labs(title = "",
       subtitle = "") +
  scale_x_continuous(breaks = seq(1,col.lim,1),limits = c(1,col.lim)) +
  scale_y_continuous(breaks = seq(1,row.lim,1),limits = c(row.lim,1),trans = 'reverse') +
  theme_linedraw() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_blank(),
        legend.title = element_blank()) +
  coord_cartesian(xlim = c(1,col.lim),
                  ylim = c(row.lim,1))

wc4 <- ggplot(lay.384[lay.384$`384plate` == 4,]) +
  geom_point(aes(x = `384col`, y = `384row`, col = colony),
             shape = 19,
             size = sz.384) +
  scale_color_manual(name = 'Colony Type',
                     breaks = c('Reference', 'Query', 'Gap'),
                     values = c("Reference" = "#FFC107",
                                "Query" = "#303F9F",
                                "Gap" = "#FF5252"),
                     guide = F) +
  labs(title = "",
       subtitle = "") +
  scale_x_continuous(breaks = seq(1,col.lim,1),limits = c(1,col.lim)) +
  scale_y_continuous(breaks = seq(1,row.lim,1),limits = c(row.lim,1),trans = 'reverse') +
  theme_linedraw() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_blank(),
        legend.title = element_blank()) +
  coord_cartesian(xlim = c(1,col.lim),
                  ylim = c(row.lim,1))

# fig1.s2 <- ggarrange(wc1, wc2, wc3, wc4,
#                      nrow = 1)
# ggsave(sprintf("%sfigure1_s2.jpg",out_path),
#        fig1.s2,
#        width = max.w, height = max.w/4 * 3/4,
#        dpi = 300)

ps1 <- ggplot(lay.1536[lay.1536$`1536plate` == 1,]) +
  geom_point(aes(x = `1536col`, y = `1536row`, col = colony),
             shape = 19,
             size = sz.1536) +
  scale_color_manual(name = 'Colony Type',
                     breaks = c('Reference', 'Query', 'Gap'),
                     values = c("Reference" = "#FFC107",
                                "Query" = "#303F9F",
                                "Gap" = "#FF5252"),
                     guide = F) +
  labs(title = "",
       subtitle = "") +
  scale_x_continuous(breaks = seq(1,col.lim*2,1),limits = c(1,col.lim*2)) +
  scale_y_continuous(breaks = seq(1,row.lim*2,1),limits = c(row.lim*2,1),trans = 'reverse') +
  theme_linedraw() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_blank(),
        legend.title = element_blank()) +
  coord_cartesian(xlim = c(1,col.lim*2),
                  ylim = c(row.lim*2,1))

ps2 <- ggplot(lay.1536[lay.1536$`1536plate` == 2,]) +
  geom_point(aes(x = `1536col`, y = `1536row`, col = colony),
             shape = 19,
             size = sz.1536) +
  scale_color_manual(name = 'Colony Type',
                     breaks = c('Reference', 'Query', 'Gap'),
                     values = c("Reference" = "#FFC107",
                                "Query" = "#303F9F",
                                "Gap" = "#FF5252"),
                     guide = F) +
  labs(title = "",
       subtitle = "") +
  scale_x_continuous(breaks = seq(1,col.lim*2,1),limits = c(1,col.lim*2)) +
  scale_y_continuous(breaks = seq(1,row.lim*2,1),limits = c(row.lim*2,1),trans = 'reverse') +
  theme_linedraw() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_blank(),
        legend.title = element_blank()) +
  coord_cartesian(xlim = c(1,col.lim*2),
                  ylim = c(row.lim*2,1))

ps3 <- ggplot(lay.1536[lay.1536$`1536plate` == 3,]) +
  geom_point(aes(x = `1536col`, y = `1536row`, col = colony),
             shape = 19,
             size = sz.1536) +
  scale_color_manual(name = 'Colony Type',
                     breaks = c('Reference', 'Query', 'Gap'),
                     values = c("Reference" = "#FFC107",
                                "Query" = "#303F9F",
                                "Gap" = "#FF5252"),
                     guide = F) +
  labs(title = "",
       subtitle = "") +
  scale_x_continuous(breaks = seq(1,col.lim*2,1),limits = c(1,col.lim*2)) +
  scale_y_continuous(breaks = seq(1,row.lim*2,1),limits = c(row.lim*2,1),trans = 'reverse') +
  theme_linedraw() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_blank(),
        legend.title = element_blank()) +
  coord_cartesian(xlim = c(1,col.lim*2),
                  ylim = c(row.lim*2,1))

ps4 <- ggplot(lay.1536[lay.1536$`1536plate` == 4,]) +
  geom_point(aes(x = `1536col`, y = `1536row`, col = colony),
             shape = 19,
             size = sz.1536) +
  scale_color_manual(name = 'Colony Type',
                     breaks = c('Reference', 'Query', 'Gap'),
                     values = c("Reference" = "#FFC107",
                                "Query" = "#303F9F",
                                "Gap" = "#FF5252"),
                     guide = F) +
  labs(title = "",
       subtitle = "") +
  scale_x_continuous(breaks = seq(1,col.lim*2,1),limits = c(1,col.lim*2)) +
  scale_y_continuous(breaks = seq(1,row.lim*2,1),limits = c(row.lim*2,1),trans = 'reverse') +
  theme_linedraw() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_blank(),
        legend.title = element_blank()) +
  coord_cartesian(xlim = c(1,col.lim*2),
                  ylim = c(row.lim*2,1))

# fig1.s3 <- ggarrange(ps1, ps2, ps3, ps4,
#                      nrow = 1)
# ggsave(sprintf("%sfigure1_s3.jpg",out_path),
#        fig1.s3,
#        width = max.w, height = max.w/4 * 3/4,
#        dpi = 300)

fs1 <- ggplot(lay.6144[lay.6144$`6144plate` == 1,]) +
  geom_point(aes(x = `6144col`, y = `6144row`, col = colony),
             shape = 19,
             size = sz.6144) +
  scale_color_manual(name = 'Colony Type',
                     breaks = c('Stock','Reference', 'Query', 'Gap'),
                     values = c("Stock" = "#689F38",
                                "Reference" = "#FFC107",
                                "Query" = "#303F9F",
                                "Gap" = "#FF5252"),
                     drop = F,
                     guide = F) +
  labs(title = "",
       subtitle = "") +
  scale_x_continuous(breaks = seq(0,col.lim*4,2),limits = c(1,col.lim*4)) +
  scale_y_continuous(breaks = seq(0,row.lim*4,2),limits = c(row.lim*4,1),trans = 'reverse') +
  theme_linedraw() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_blank(),
        legend.title = element_blank()) +
  coord_cartesian(xlim = c(1,col.lim*4),
                  ylim = c(row.lim*4,1))

lay.6144[6144*2+1,4] = 2
lay.6144[6144*2+1,7] = 'Stock' # done to add 'Stock' colony type to legend

fs2 <- ggplot(lay.6144[lay.6144$`6144plate` == 2,]) +
  geom_point(aes(x = `6144col`, y = `6144row`, col = colony,
                 shape = colony),
             size = sz.6144) +
  scale_color_manual(name = 'Colony Type',
                     breaks = c('Reference', 'Query', 'Gap'),
                     values = c("Stock" = "#689F38",
                                "Reference" = "#FFC107",
                                "Query" = "#303F9F",
                                "Gap" = "#FF5252"),
                     labels = c('Reference', 'Query', 'Gap'),
                     drop = FALSE,
                     guide = F) +
  scale_shape_manual(name = 'Plate Type',
                     breaks = c('Stock','Reference'),
                     values = c("Stock" = 15,
                                "Reference" = 19,
                                "Query" = 19,
                                "Gap" = 19),
                     labels = c('Deep Well','Agar'),
                     drop = FALSE,
                     guide = F) +
  labs(title = "",
       subtitle = "") +
  scale_x_continuous(breaks = seq(0,col.lim*4,2),limits = c(1,col.lim*4)) +
  scale_y_continuous(breaks = seq(0,row.lim*4,2),limits = c(row.lim*4,1),trans = 'reverse') +
  theme_linedraw() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_blank(),
        legend.title = element_blank()) +
  coord_cartesian(xlim = c(1,col.lim*4),
                  ylim = c(row.lim*4,1))

# fig1.s4 <- ggarrange(fs1, fs2,
#                      nrow = 1)
# ggsave(sprintf("%sfigure1_s4.jpg",out_path),
#        fig1.s4,
#        width = max.w/2, height = max.w/4 * 3/4,
#        dpi = 300)

fig1 <- ggarrange(gly1,gly2,gly3,gly4,
                  wc1,wc2,wc3,wc4,
                  ps1,ps2,ps3,ps4,
                  fs1,fs2,
                  nrow = 4,
                  ncol = 4)
ggsave(sprintf("%sfigure1.jpg",out_path),
       fig1,
       width = max.w * 4/6, height = max.w * .80 * 4/6,
       dpi = 300)

##### FIGURE 2
f2.min = 200
f2.max = 700

fig2.s1 <- ggplot(data = rfitdat, aes(x = `6144col`, y = `6144row`)) +
  geom_point(aes(x = `6144col`, y = `6144row`,col = average,
                 shape = colony,
                 alpha = colony),
             size = sz.6144/3, na.rm = T) +
  labs(title = "",
       subtitle = "",
       x = "",
       y = "") +
  scale_x_continuous(breaks = seq(1,96,1)) +
  scale_y_continuous(breaks = seq(1,64,1),trans = 'reverse') +
  scale_color_distiller(name = "PIX",
                        limits = c(f2.min,f2.max),
                        palette = "Set1",
                        guide = F) +
  scale_shape_manual(name="Colony Kind",
                     values=c("Gap"=15,"Query"=15,"Reference"=15),
                     breaks=c("Reference","Query","Gap"),guide=F) +
  scale_alpha_manual(values=c("Gap"=1,"Query"=1,"Reference"=1),
                     guide=F) +
  theme_linedraw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size=lax.txt),
        legend.title = element_text(size=lax.ttle),
        legend.position = "right") +
  coord_cartesian(xlim = c(1,96),
                  ylim = c(64,1))

ggsave(sprintf("%sfigure2_s1.jpg",out_path),
       fig2.s1,
       width = max.w * 4/6 * 1/4, height = max.w * .80 * 4/6 * 1/4,
       dpi = 300)

step2.tl <- ggplot(data = rfitdat[rfitdat$source == "TL",], aes(x = `6144col`, y = `6144row`)) +
  geom_point(aes(x = `6144col`, y = `6144row`,col = average,
                 shape = colony,
                 alpha = colony),
             size = sz.1536/2, na.rm = T) +
  labs(title = "",
       subtitle = "",
       x = "",
       y = "") +
  scale_x_continuous(breaks = seq(1,96,1)) +
  scale_y_continuous(breaks = seq(1,64,1),trans = 'reverse') +
  scale_color_distiller(name = "PIX",
                        limits = c(f2.min,f2.max),
                        palette = "Set1") +
  scale_shape_manual(name="Colony Kind",
                     values=c("Gap"=15,"Query"=15,"Reference"=15),
                     breaks=c("Reference","Query","Gap"),
                     guide = F) +
  scale_alpha_manual(values=c("Gap"=1,"Query"=1,"Reference"=1),
                     guide=F) +
  theme_linedraw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size=lax.txt),
        legend.title = element_text(size=lax.ttle),
        legend.position = "right") +
  coord_cartesian(xlim = c(1,96),
                  ylim = c(64,1))

step2.tr <- ggplot(data = rfitdat[rfitdat$source == "TR",], aes(x = `6144col`, y = `6144row`)) +
  geom_point(aes(x = `6144col`, y = `6144row`,col = average,
                 shape = colony,
                 alpha = colony),
             size = sz.1536/2, na.rm = T) +
  labs(title = "",
       subtitle = "",
       x = "",
       y = "") +
  scale_x_continuous(breaks = seq(1,96,1)) +
  scale_y_continuous(breaks = seq(1,64,1),trans = 'reverse') +
  scale_color_distiller(name = "PIX",
                        limits = c(f2.min,f2.max),
                        palette = "Set1") +
  scale_shape_manual(name="Colony Kind",
                     values=c("Gap"=15,"Query"=15,"Reference"=15),
                     breaks=c("Reference","Query","Gap"),guide=F) +
  scale_alpha_manual(values=c("Gap"=1,"Query"=1,"Reference"=1),
                     guide=F) +
  theme_linedraw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size=lax.txt),
        legend.title = element_text(size=lax.ttle),
        legend.position = "right") +
  coord_cartesian(xlim = c(1,96),
                  ylim = c(64,1))

step2.bl <- ggplot(data = rfitdat[rfitdat$source == "BL",], aes(x = `6144col`, y = `6144row`)) +
  geom_point(aes(x = `6144col`, y = `6144row`,col = average,
                 shape = colony,
                 alpha = colony),
             size = sz.1536/2, na.rm = T) +
  labs(title = "",
       subtitle = "",
       x = "",
       y = "") +
  scale_x_continuous(breaks = seq(1,96,1)) +
  scale_y_continuous(breaks = seq(1,64,1),trans = 'reverse') +
  scale_color_distiller(name = "PIX",
                        limits = c(f2.min,f2.max),
                        palette = "Set1") +
  scale_shape_manual(name="Colony Kind",
                     values=c("Gap"=15,"Query"=15,"Reference"=15),
                     breaks=c("Reference","Query","Gap"),guide=F) +
  scale_alpha_manual(values=c("Gap"=1,"Query"=1,"Reference"=1),
                     guide=F) +
  theme_linedraw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size=lax.txt),
        legend.title = element_text(size=lax.ttle),
        legend.position = "right") +
  coord_cartesian(xlim = c(1,96),
                  ylim = c(64,1))

step2.br <- ggplot(data = rfitdat[rfitdat$source == "BR",], aes(x = `6144col`, y = `6144row`)) +
  geom_point(aes(x = `6144col`, y = `6144row`,col = average,
                 shape = colony,
                 alpha = colony),
             size = sz.1536/2, na.rm = T) +
  labs(title = "",
       subtitle = "",
       x = "",
       y = "") +
  scale_x_continuous(breaks = seq(1,96,1)) +
  scale_y_continuous(breaks = seq(1,64,1),trans = 'reverse') +
  scale_color_distiller(name = "PIX",
                        limits = c(f2.min,f2.max),
                        palette = "Set1") +
  scale_shape_manual(name="Colony Kind",
                     values=c("Gap"=15,"Query"=15,"Reference"=15),
                     breaks=c("Reference","Query","Gap"),guide=F) +
  scale_alpha_manual(values=c("Gap"=1,"Query"=1,"Reference"=1),
                     guide=F) +
  theme_linedraw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size=lax.txt),
        legend.title = element_text(size=lax.ttle),
        legend.position = "right") +
  coord_cartesian(xlim = c(1,96),
                  ylim = c(64,1))

fig2.s2 <- ggarrange(step2.tl,step2.tr,step2.bl,step2.br,
          nrow = 2, ncol = 2,
          common.legend = T, legend = F)

ggsave(sprintf("%sfigure2_s2.jpg",out_path),
       fig2.s2,
       width = max.w * 4/6 * 1/2, height = max.w * .80 * 4/6 * 1/2,
       dpi = 300)

step3.tl <- ggplot(data = rfitdat[rfitdat$source == "TL",], aes(x = `6144col`, y = `6144row`)) +
  geom_point(aes(x = `6144col`, y = `6144row`,col = bg,
                 shape = colony,
                 alpha = colony),
             size = sz.1536/2, na.rm = T) +
  labs(title = "",
       subtitle = "",
       x = "",
       y = "") +
  scale_x_continuous(breaks = seq(1,96,1)) +
  scale_y_continuous(breaks = seq(1,64,1),trans = 'reverse') +
  scale_color_distiller(name = "PIX",
                        limits = c(f2.min,f2.max),
                        palette = "Set1",
                        guide = F) +
  scale_shape_manual(name="Colony Kind",
                     values=c("Gap"=15,"Query"=15,"Reference"=15),
                     breaks=c("Reference","Query","Gap"),guide=F) +
  scale_alpha_manual(values=c("Gap"=1,"Query"=1,"Reference"=1),
                     guide=F) +
  theme_linedraw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size=lax.txt),
        legend.title = element_text(size=lax.ttle),
        legend.position = "right") +
  coord_cartesian(xlim = c(1,96),
                  ylim = c(64,1))

step3.tr <- ggplot(data = rfitdat[rfitdat$source == "TR",], aes(x = `6144col`, y = `6144row`)) +
  geom_point(aes(x = `6144col`, y = `6144row`,col = bg,
                 shape = colony,
                 alpha = colony),
             size = sz.1536/2, na.rm = T) +
  labs(title = "",
       subtitle = "",
       x = "",
       y = "") +
  scale_x_continuous(breaks = seq(1,96,1)) +
  scale_y_continuous(breaks = seq(1,64,1),trans = 'reverse') +
  scale_color_distiller(name = "PIX",
                        limits = c(f2.min,f2.max),
                        palette = "Set1",
                        guide = F) +
  scale_shape_manual(name="Colony Kind",
                     values=c("Gap"=15,"Query"=15,"Reference"=15),
                     breaks=c("Reference","Query","Gap"),guide=F) +
  scale_alpha_manual(values=c("Gap"=1,"Query"=1,"Reference"=1),
                     guide=F) +
  theme_linedraw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size=lax.txt),
        legend.title = element_text(size=lax.ttle),
        legend.position = "right") +
  coord_cartesian(xlim = c(1,96),
                  ylim = c(64,1))

step3.bl <- ggplot(data = rfitdat[rfitdat$source == "BL",], aes(x = `6144col`, y = `6144row`)) +
  geom_point(aes(x = `6144col`, y = `6144row`,col = bg,
                 shape = colony,
                 alpha = colony),
             size = sz.1536/2, na.rm = T) +
  labs(title = "",
       subtitle = "",
       x = "",
       y = "") +
  scale_x_continuous(breaks = seq(1,96,1)) +
  scale_y_continuous(breaks = seq(1,64,1),trans = 'reverse') +
  scale_color_distiller(name = "PIX",
                        limits = c(f2.min,f2.max),
                        palette = "Set1",
                        guide = F) +
  scale_shape_manual(name="Colony Kind",
                     values=c("Gap"=15,"Query"=15,"Reference"=15),
                     breaks=c("Reference","Query","Gap"),guide=F) +
  scale_alpha_manual(values=c("Gap"=1,"Query"=1,"Reference"=1),
                     guide=F) +
  theme_linedraw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size=lax.txt),
        legend.title = element_text(size=lax.ttle),
        legend.position = "right") +
  coord_cartesian(xlim = c(1,96),
                  ylim = c(64,1))

step3.br <- ggplot(data = rfitdat[rfitdat$source == "BR",], aes(x = `6144col`, y = `6144row`)) +
  geom_point(aes(x = `6144col`, y = `6144row`,col = bg,
                 shape = colony,
                 alpha = colony),
             size = sz.1536/2, na.rm = T) +
  labs(title = "",
       subtitle = "",
       x = "",
       y = "") +
  scale_x_continuous(breaks = seq(1,96,1)) +
  scale_y_continuous(breaks = seq(1,64,1),trans = 'reverse') +
  scale_color_distiller(name = "PIX",
                        limits = c(f2.min,f2.max),
                        palette = "Set1",
                        guide = F) +
  scale_shape_manual(name="Colony Kind",
                     values=c("Gap"=15,"Query"=15,"Reference"=15),
                     breaks=c("Reference","Query","Gap"),guide=F) +
  scale_alpha_manual(values=c("Gap"=1,"Query"=1,"Reference"=1),
                     guide=F) +
  theme_linedraw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size=lax.txt),
        legend.title = element_text(size=lax.ttle),
        legend.position = "right") +
  coord_cartesian(xlim = c(1,96),
                  ylim = c(64,1))

fig2.s3 <- ggarrange(step3.tl,step3.tr,step3.bl,step3.br,
                     nrow = 2, ncol = 2,
                     common.legend = T, legend = F)

ggsave(sprintf("%sfigure2_s3.jpg",out_path),
       fig2.s3,
       width = max.w * 4/6 * 1/2, height = max.w * .80 * 4/6 * 1/2,
       units = "cm",
       dpi = 300)

fig2.s4 <- ggplot(data = rfitdat, aes(x = `6144col`, y = `6144row`)) +
  geom_point(aes(x = `6144col`, y = `6144row`,col = bg,
                 shape = colony,
                 alpha = colony),
             size = sz.6144/3, na.rm = T) +
  labs(title = "",
       subtitle = "",
       x = "",
       y = "") +
  scale_x_continuous(breaks = seq(1,96,1)) +
  scale_y_continuous(breaks = seq(1,64,1),trans = 'reverse') +
  scale_color_distiller(name = "PIX",
                        limits = c(f2.min,f2.max),
                        palette = "Set1",
                        guide = F) +
  scale_shape_manual(name="Colony Kind",
                     values=c("Gap"=15,"Query"=15,"Reference"=15),
                     breaks=c("Reference","Query","Gap"),guide=F) +
  scale_alpha_manual(values=c("Gap"=1,"Query"=1,"Reference"=1),
                     guide=F) +
  theme_linedraw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size=lax.txt),
        legend.title = element_text(size=lax.ttle),
        legend.position = "right") +
  coord_cartesian(xlim = c(1,96),
                  ylim = c(64,1))

ggsave(sprintf("%sfigure2_s4.jpg",out_path),
       fig2.s4,
       width = max.w * 4/6 * 1/4, height = max.w * .80 * 4/6 * 1/4,
       dpi = 300)


# legend <- cowplot::get_legend(fig2.s4)
# grid.newpage()
# grid.draw(legend)
# ggsave(sprintf("%sfigure2_lgnd.jpg",out_path),
#        legend,
#        width = 2,height = 4,
#        limitsize = F)

##### FIGURE 3
obs <- ggplot(data = fitdat, aes(x = `6144col`, y = `6144row`)) +
  geom_point(aes(x = `6144col`, y = `6144row`,col = average,
                 shape = colony,
                 alpha = colony),
             size = sz.6144/3, na.rm = T) +
  labs(title = "Observed Colony Size (Pixel Count)",
       x = "",
       y = "") +
  scale_x_continuous(breaks = seq(1,96,1)) +
  scale_y_continuous(breaks = seq(1,64,1),trans = 'reverse') +
  scale_color_distiller(name = "PIX",
                        limits = c(min,max),
                        palette = "Set1") +
  scale_shape_manual(name="Colony Kind",
                     values=c("Gap"=15,"Query"=15,"Reference"=15),
                     breaks=c("Reference","Query","Gap"),guide=F) +
  scale_alpha_manual(values=c("Gap"=1,"Query"=1,"Reference"=1),
                     guide=F) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.text = element_text(size=15),
        legend.title = element_text(size=20,face="bold"),
        legend.position = "right",
        plot.title = element_text(size=25,hjust = 0),
        plot.subtitle = element_text(size=20,hjust = 0)) +
  coord_cartesian(xlim = c(1,96),
                  ylim = c(64,1))

##### 3B
pre <- ggplot(data = fitdat, aes(x = `6144col`, y = `6144row`)) +
  geom_point(aes(x = `6144col`, y = `6144row`,col = bg,
                 shape = colony,
                 alpha = colony),size = sz,na.rm = T) +
  labs(title = "Predicted Colony Size (Pixel Count)",
       x = "",
       y = "") +
  scale_x_continuous(breaks = seq(1,96,1)) +
  scale_y_continuous(breaks = seq(1,64,1),trans = 'reverse') +
  scale_color_distiller(name = "PIX",
                        limits = c(min,max),
                        palette = "Set1") +
  scale_shape_manual(name="Colony Kind",
                     values=c("Gap"=15,"Query"=15,"Reference"=15),
                     breaks=c("Reference","Query","Gap"),guide=F) +
  scale_alpha_manual(values=c("Gap"=1,"Query"=1,"Reference"=1),
                     guide=F) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.text = element_text(size=15),
        legend.title = element_text(size=20,face="bold"),
        legend.position = "right",
        plot.title = element_text(size=25,hjust = 0),
        plot.subtitle = element_text(size=20,hjust = 0)) +
  coord_cartesian(xlim = c(1,96),
                  ylim = c(64,1))

##### 3C
fit <- ggplot(data = fitdat, aes(x = `6144col`, y = `6144row`)) +
  geom_point(aes(x = `6144col`, y = `6144row`,col = fitness,
                 shape = colony,
                 alpha = colony),size = sz,na.rm = T) +
  labs(title = "Fitness Distribution",
       x = "",
       y = "") +
  scale_x_continuous(breaks = seq(1,96,1)) +
  scale_y_continuous(breaks = seq(1,64,1),trans = 'reverse') +
  scale_color_distiller(name = "FIT",
                        limits = c(0.7,1.3),
                        breaks = c(0.70,0.85,1.00,1.15,1.30),
                        palette = "Set1") +
  scale_shape_manual(name="Colony Kind",
                     values=c("Gap"=15,"Query"=15,"Reference"=15),
                     breaks=c("Reference","Query","Gap"),guide=F) +
  scale_alpha_manual(values=c("Gap"=1,"Query"=1,"Reference"=1),
                     guide=F) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.text = element_text(size=15),
        legend.title = element_text(size=20,face="bold"),
        legend.position = "right",
        plot.title = element_text(size=25,hjust = 0),
        plot.subtitle = element_text(size=20,hjust = 0)) +
  coord_cartesian(xlim = c(1,96),
                  ylim = c(64,1))

##### FINAL FIGURE 3
fig3 <- ggarrange(obs, pre, fit,
                  nrow = 1)
ggsave(sprintf("%sfigure3.jpg",out_path),
       fig3,
       width = 30,height = 6.5)

##### FIGURE 4
##### 4A
raw <- ggplot(data = fitdat, aes(x=average, col = source)) +
  geom_density(lwd = 1.2) + 
  scale_colour_manual(name="Source",
                      values=c("TL"="#D32F2F","TR"="#536DFE","BL"="#388E3C","BR"="#795548"),
                      breaks=c("TL","TR","BL","BR"),
                      labels=c("Top Left","Top Right","Bottom Left","Bottom Right"),
                      guide = F) +
  scale_x_continuous(breaks = seq(0,1000,50),
                     minor_breaks = seq(0,1000,25),
                     limits = c(min,max)) +
  scale_y_continuous(breaks = seq(0,1,0.002),
                     minor_breaks = seq(0,1,0.001),
                     limits = c(0,0.02)) +
  labs(title = 'A. Raw colony sizes',
       x = 'Observed Pixel Count', y = 'Density') +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=20),
        legend.text = element_text(size=15),
        legend.title = element_text(size=20),
        legend.position = 'bottom',
        legend.background = element_rect(fill="lightblue", 
                                         size=0.5, linetype="solid"),
        plot.title = element_text(size=25,hjust = -0.15)) +
  coord_cartesian(ylim = c(0,0.018))

##### 4B
src.nrm <- ggplot(data = fitdat, aes(x=fitness, col = source)) +
  geom_density(lwd = 1.2) + 
  scale_colour_manual(name="Source",
                      values=c("TL"="#D32F2F","TR"="#536DFE","BL"="#388E3C","BR"="#795548"),
                      breaks=c("TL","TR","BL","BR"),
                      labels=c("Top Left","Top Right","Bottom Left","Bottom Right")) +
  scale_x_continuous(breaks = seq(0,2,0.05),
                     minor_breaks = seq(0,2,0.025),
                     limits = c(0.7,1.3)) +
  scale_y_continuous(breaks = seq(0,15,1),
                     minor_breaks = seq(0,15,0.5),
                     limits = c(0,12)) +
  labs(title= 'B. With Source Normalization',
       x = 'Fitness', y = 'Density') +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=20),
        legend.text = element_text(size=15),
        legend.title = element_text(size=20),
        legend.position = 'bottom',
        legend.background = element_rect(fill="lightblue", 
                                         size=0.5, linetype="solid"),
        plot.title = element_text(size=25,hjust = -0.13)) +
  coord_cartesian(ylim = c(0,11))

##### 4C
no.src.nrm <- ggplot(data = fitdat, aes(x=nfitness, col = source)) +
  geom_density(lwd = 1.2) + 
  scale_colour_manual(name="Source",
                      values=c("TL"="#D32F2F","TR"="#536DFE","BL"="#388E3C","BR"="#795548"),
                      breaks=c("TL","TR","BL","BR"),
                      labels=c("Top Left","Top Right","Bottom Left","Bottom Right"),
                      guide = F) +
  scale_x_continuous(breaks = seq(0,2,0.05),
                     minor_breaks = seq(0,2,0.025),
                     limits = c(0.7,1.3)) +
  scale_y_continuous(breaks = seq(0,15,1),
                     minor_breaks = seq(0,15,0.5),
                     limits = c(0,12)) +
  labs(title= 'C. Without Source Normalization',
       x = 'Fitness', y = 'Density') +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=20),
        legend.text = element_text(size=15),
        legend.title = element_text(size=20),
        legend.position = 'bottom',
        legend.background = element_rect(fill="lightblue", 
                                         size=0.5, linetype="solid"),
        plot.title = element_text(size=25,hjust = -0.13)) +
  coord_cartesian(ylim = c(0,11))

##### FINAL FIGURE 4
fig4 <- ggarrange(raw, src.nrm, no.src.nrm,
                      nrow = 1)
ggsave(sprintf("%sfigure4.jpg",out_path),
       fig4,
       width = 30,height = 11)


##### FIGURE S2
pl = 1
vpdat = dbGetQuery(conn, sprintf('select *
                                  from %s a, %s b
                                  where a.pos = b.pos
                                  and a.hours in (10,14,18)
                                  and b.%s = %d
                                  order by a.hours, b.%s, b.%s',
                                  tablename_fit,
                                  p2c_info[1],p2c_info[2],
                                  pl,p2c_info[3],p2c_info[4]))
vpdat$average[is.na(vpdat$average)] = 0

vpdat$colony[vpdat$orf_name == 'BF_control'] = 'Reference'
vpdat$colony[vpdat$orf_name != 'BF_control'] = 'Query'
vpdat$colony[is.na(vpdat$orf_name)] = 'Gap'

vp10 <- vpdat[vpdat$hours == 10,]
vp14 <- vpdat[vpdat$hours == 14,]
vp18 <- vpdat[vpdat$hours == 18,]

vp10$average[vp10$colony == 'Reference'] <- vp14$average[vp14$colony == 'Reference']
vp10$bg <- vp14$bg
vp10$fitness <- vp10$average/vp10$bg

vp18$average[vp18$colony == 'Reference'] <- vp14$average[vp14$colony == 'Reference']
vp18$bg <- vp14$bg
vp18$fitness <- vp18$average/vp18$bg

fit.range <- c(0.5,1.5)

v10 <- ggplot(vp10) +
  geom_point(aes(x = `6144col`, y = `6144row`,
                 col = colony,
                 size = fitness),
             shape = 19) +
  scale_color_manual(name = 'Colony Type',
                     breaks = c('Reference', 'Query', 'Gap'),
                     values = c("Reference" = "#00796B",
                                "Query" = "#673AB7",
                                "Gap" = "#FF5252"),
                     guide = F) +
  scale_size(limits = fit.range, range = c(0, 11), guide = F) +
  labs(title = "A. t(R) = 14 hours, t(Q) = 10 hours",
       x = 'Columns',
       y = 'Rows') +
  scale_x_continuous(breaks = seq(0,96,2),limits = c(1,96)) +
  scale_y_continuous(breaks = seq(0,64,2),limits = c(64,1),trans = 'reverse') +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.title.y =  element_blank(),
        legend.text = element_text(size=15),
        legend.title = element_text(size=18),
        legend.position = 'right',
        plot.title = element_text(size=20,hjust = 0)) +
  coord_cartesian(xlim = c(32,56),
                  ylim = c(38,22))

v14 <- ggplot(vp14) +
  geom_point(aes(x = `6144col`, y = `6144row`,
                 col = colony,
                 size = fitness),
             shape = 19) +
  scale_color_manual(name = 'Colony Type',
                     breaks = c('Reference', 'Query'),
                     values = c("Reference" = "#00796B",
                                "Query" = "#673AB7",
                                "Gap" = "#FF5252"),
                     guide = F) +
  scale_size(limits = fit.range, range = c(0, 11), guide = F) +
  labs(title = "C. t(R) = 14 hours, t(Q) = 14 hours",
       x = 'Columns',
       y = 'Rows') +
  scale_x_continuous(breaks = seq(0,96,2),limits = c(1,96)) +
  scale_y_continuous(breaks = seq(0,64,2),limits = c(64,1),trans = 'reverse') +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.title.y =  element_blank(),
        legend.text = element_text(size=15),
        legend.title = element_text(size=18),
        legend.position = 'right',
        plot.title = element_text(size=20,hjust = 0)) +
  coord_cartesian(xlim = c(32,56),
                  ylim = c(38,22))

v18 <- ggplot(vp18) +
  geom_point(aes(x = `6144col`, y = `6144row`,
                 col = colony,
                 size = fitness),
             shape = 19) +
  scale_color_manual(name = 'Colony Type',
                     breaks = c('Reference', 'Query'),
                     values = c("Reference" = "#00796B",
                                "Query" = "#673AB7",
                                "Gap" = "#FF5252"),
                     guide = F) +
  scale_size(limits = fit.range, range = c(0, 11), guide = F) +
  labs(title = "B. t(R) = 14 hours, t(Q) = 18 hours",
       x = 'Columns',
       y = 'Rows') +
  scale_x_continuous(breaks = seq(0,96,2),limits = c(1,96)) +
  scale_y_continuous(breaks = seq(0,64,2),limits = c(64,1),trans = 'reverse') +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.title.y =  element_blank(),
        legend.text = element_text(size=15),
        legend.title = element_text(size=18),
        legend.position = 'bottom',
        plot.title = element_text(size=20,hjust = 0)) +
  guides(color = guide_legend(override.aes = list(size=6))) +
  coord_cartesian(xlim = c(32,56),
                  ylim = c(38,22))

fig.s2 <- ggarrange(v10, v18, v14,
                     nrow = 1)
ggsave(sprintf("%sfigure_s2.jpg",out_path),
       fig.s2,
       width = 30,height = 7.2)

##### FIGURE S3
pl = 1
sca.dat = dbGetQuery(conn, sprintf('select *
                                  from %s a, %s b
                                  where a.pos = b.pos
                                  and b.%s = %d
                                  order by a.hours, b.%s, b.%s',
                                  tablename_fit,
                                  p2c_info[1],p2c_info[2],
                                  pl,p2c_info[3],p2c_info[4]))

sca.dat$bg[is.na(sca.dat$average)] = NA
min = min(sca.dat$average, na.rm=T)
max = max(sca.dat$average, na.rm=T)

sca.dat$se <- sca.dat$average - sca.dat$bg

sca.dat$source[sca.dat$`6144row`%%2==1 & sca.dat$`6144col`%%2==1] = 'TL'
sca.dat$source[sca.dat$`6144row`%%2==0 & sca.dat$`6144col`%%2==1] = 'BL'
sca.dat$source[sca.dat$`6144row`%%2==1 & sca.dat$`6144col`%%2==0] = 'TR'
sca.dat$source[sca.dat$`6144row`%%2==0 & sca.dat$`6144col`%%2==0] = 'BR'

sca.dat$colony[sca.dat$orf_name == 'BF_control'] = 'Reference'
sca.dat$colony[sca.dat$orf_name != 'BF_control'] = 'Query'
sca.dat$colony[is.na(sca.dat$orf_name)] = 'Gap'

sca.dat$outlier = NA
sca.dat$neigh = NA
sca.dat$diff = NA

for (sr in unique(sca.dat$source)) {
  temp <- sca.dat[sca.dat$source == sr,]
  for (hr in unique(temp$hours)) {
    se.m <- mean(temp$se[temp$hours == hr],na.rm=T)
    se.s <- sd(temp$se[temp$hours == hr],na.rm=T)
    temp$outlier[temp$se > se.m + 3*se.s & temp$hours == hr] = 'Bigger'
    temp$outlier[temp$se < se.m - 3*se.s & temp$hours == hr] = 'Smaller'
    temp$outlier[is.na(temp$outlier) & temp$hours == hr] = 'Normal'

    for (i in seq(2,length(unique(temp$`6144col`)))) {
      col <- unique(temp$`6144col`)[i]
      lf <- tail(temp$`6144col`[temp$`6144col` < col & temp$orf_name == 'BF_control' & !is.na(temp$orf_name)],1)
      rt <- temp$`6144col`[temp$`6144col` > col & temp$orf_name == 'BF_control' & !is.na(temp$orf_name)][1]
      for (ii in seq(2,length(unique(temp$`6144row`[temp$`6144col` == col])))) {
        row <- unique(temp$`6144row`[temp$`6144col` == col])[ii]
        up <- tail(temp$`6144row`[temp$`6144row` < row & temp$orf_name == 'BF_control' & !is.na(temp$orf_name)],1)
        dw <- temp$`6144row`[temp$`6144row` > row & temp$orf_name == 'BF_control' & !is.na(temp$orf_name)][1]
        if (!is.na(sca.dat$average[sca.dat$`6144col` == col & sca.dat$`6144row` == row & sca.dat$hours == hr])) {
          a <- sca.dat$average[sca.dat$`6144col` == col & sca.dat$`6144row` == row & sca.dat$hours == hr]
          u <- sca.dat$average[sca.dat$`6144col` == col & sca.dat$`6144row` == up & sca.dat$hours == hr]
          d <- sca.dat$average[sca.dat$`6144col` == col & sca.dat$`6144row` == dw & sca.dat$hours == hr]
          l <- sca.dat$average[sca.dat$`6144col` == lf & sca.dat$`6144row` == row & sca.dat$hours == hr]
          r <- sca.dat$average[sca.dat$`6144col` == rt & sca.dat$`6144row` == row & sca.dat$hours == hr]
          sca.dat$var[sca.dat$`6144col` == col & sca.dat$`6144row` == row & sca.dat$hours == hr] <-
            sd(c(a,u,d,l,r),na.rm = T)/mean(c(a,u,d,l,r),na.rm = T)
          sca.dat$neigh[sca.dat$`6144col` == col & sca.dat$`6144row` == row & sca.dat$hours == hr] <-  mean(c(u,d,l,r),na.rm = T)
          sca.dat$diff[sca.dat$`6144col` == col & sca.dat$`6144row` == row & sca.dat$hours == hr] <- a - mean(c(u,d,l,r),na.rm = T)
        }
      }
    }
  }
  sca.dat$outlier[sca.dat$pos %in% temp$pos] = temp$outlier
}

# save(sca.dat, file = "output/sca.dat.CS.RData")
# load('output/sca.dat.RData')

# sca.dat <- sca.dat[sca.dat$hours > 0,]
  
sca.tl <- ggplot(sca.dat[sca.dat$source == 'TL',]) +
  geom_abline(intercept = c(0)) +
  geom_point(aes(x=average, y=bg, col = source, shape = colony),
             alpha = 0.5, size = 3) +
  geom_point(aes(x=average, y=bg, col = outlier, shape = colony,alpha = outlier),
             size = 3) +
  scale_colour_manual(name="Outlier Type",
                      values=c("TL"="#D32F2F","TR"="#536DFE","BL"="#388E3C","BR"="#795548",
                               "Bigger"="#FFC107","Smaller"="#7B1FA2","Normal"="transparent"),
                      breaks=c("Bigger","Smaller"),
                      labels=c("Underpredicted","Overpredicted")) +
  scale_shape_manual(name="Colony Type",
                     values=c("Reference"=18,"Query"=15,"Gap"=1),
                     breaks=c("Reference","Query")) +
  scale_alpha_manual(values=c("Smaller"=0.9,"Bigger"=0.9,"Normal"=0),
                     guide = F) +
  labs(title = "S3. Accuracy of background prediction using LID",
       subtitle = sprintf("Top Left | RMSE = %0.3f",
                          sqrt(sum(abs(sca.dat$bg[sca.dat$source == 'TL'] - sca.dat$average[sca.dat$source == 'TL'])^2,na.rm = T)/
                                 length(sca.dat$bg[!is.na(sca.dat$average) & sca.dat$source == 'TL']))),
       x = "",
       y = "Predicted Colony Size (Pixel Count)") +
  scale_x_continuous(breaks = seq(0,1000,100),
                     minor_breaks = seq(0,1000,25)) +
  scale_y_continuous(breaks = seq(0,1000,100),
                     minor_breaks = seq(0,1000,25)) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=20),
        legend.position = c(0.8,0.2),
        legend.background = element_blank(),
        legend.text = element_text(size=15),
        legend.title =  element_text(size=15),
        plot.title = element_text(size=25,hjust = 0),
        plot.subtitle = element_text(size=20,hjust = 0.5)) +
  guides(color = guide_legend(override.aes = list(size=3, alpha = 1)),
         shape = guide_legend(override.aes = list(size=3))) +
  coord_cartesian(xlim = c(0,600),
                  ylim = c(0,600))

sca.tr <- ggplot(sca.dat[sca.dat$source == 'TR',]) +
  geom_abline(intercept = c(0)) +
  geom_point(aes(x=average, y=bg, col = source, shape = colony),
             alpha = 0.5, size = 3) +
  geom_point(aes(x=average, y=bg, col = outlier, shape = colony,alpha = outlier),
             size = 3) +
  scale_colour_manual(name="Outlier Type",
                      values=c("TL"="#D32F2F","TR"="#536DFE","BL"="#388E3C","BR"="#795548",
                               "Bigger"="#FFC107","Smaller"="#7B1FA2","Normal"="transparent"),
                      breaks=c("Bigger","Smaller"),
                      labels=c("Underpredicted","Overpredicted")) +
  scale_shape_manual(name="Colony Type",
                     values=c("Reference"=18,"Query"=15,"Gap"=1),
                     breaks=c("Reference","Query")) +
  scale_alpha_manual(values=c("Smaller"=0.9,"Bigger"=0.9,"Normal"=0),
                     guide = F) +
  labs(title = "",
       subtitle = sprintf("Top Right | RMSE = %0.3f",
                          sqrt(sum(abs(sca.dat$bg[sca.dat$source == 'TR'] - sca.dat$average[sca.dat$source == 'TR'])^2,na.rm = T)/
                                 length(sca.dat$bg[!is.na(sca.dat$average) & sca.dat$source == 'TR']))),
       x = "",
       y = "") +
  scale_x_continuous(breaks = seq(0,1000,100),
                     minor_breaks = seq(0,1000,25)) +
  scale_y_continuous(breaks = seq(0,1000,100),
                     minor_breaks = seq(0,1000,25)) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=20),
        legend.position = c(0.8,0.2),
        legend.background = element_blank(),
        legend.text = element_text(size=15),
        legend.title =  element_text(size=15),
        plot.title = element_text(size=25,hjust = -1),
        plot.subtitle = element_text(size=20,hjust = 0.5)) +
  guides(color = guide_legend(override.aes = list(size=3, alpha = 1)),
         shape = guide_legend(override.aes = list(size=3))) +
  coord_cartesian(xlim = c(0,600),
                  ylim = c(0,600))

sca.bl <- ggplot(sca.dat[sca.dat$source == 'BL',]) +
  geom_abline(intercept = c(0)) +
  geom_point(aes(x=average, y=bg, col = source, shape = colony),
             alpha = 0.5, size = 3) +
  geom_point(aes(x=average, y=bg, col = outlier, shape = colony,alpha = outlier),
             size = 3) +
  scale_colour_manual(name="Outlier Type",
                      values=c("TL"="#D32F2F","TR"="#536DFE","BL"="#388E3C","BR"="#795548",
                               "Bigger"="#FFC107","Smaller"="#7B1FA2","Normal"="transparent"),
                      breaks=c("Bigger","Smaller"),
                      labels=c("Underpredicted","Overpredicted")) +
  scale_shape_manual(name="Colony Type",
                     values=c("Reference"=18,"Query"=15,"Gap"=1),
                     breaks=c("Reference","Query")) +
  scale_alpha_manual(values=c("Smaller"=0.9,"Bigger"=0.9,"Normal"=0),
                     guide = F) +
  labs(title = "",
       subtitle = sprintf("Bottom Left | RMSE = %0.3f",
                          sqrt(sum(abs(sca.dat$bg[sca.dat$source == 'BL'] - sca.dat$average[sca.dat$source == 'BL'])^2,na.rm = T)/
                                 length(sca.dat$bg[!is.na(sca.dat$average) & sca.dat$source == 'BL']))),
       x = "Observed Colony Size (Pixel Count)",
       y = "Predicted Colony Size (Pixel Count)") +
  scale_x_continuous(breaks = seq(0,1000,100),
                     minor_breaks = seq(0,1000,25)) +
  scale_y_continuous(breaks = seq(0,1000,100),
                     minor_breaks = seq(0,1000,25)) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=20),
        legend.position = c(0.8,0.2),
        legend.background = element_blank(),
        legend.text = element_text(size=15),
        legend.title =  element_text(size=15),
        plot.title = element_text(size=25,hjust = -1),
        plot.subtitle = element_text(size=20,hjust = 0.5)) +
  guides(color = guide_legend(override.aes = list(size=3, alpha = 1)),
         shape = guide_legend(override.aes = list(size=3))) +
  coord_cartesian(xlim = c(0,600),
                  ylim = c(0,600))

sca.br <- ggplot(sca.dat[sca.dat$source == 'BR',]) +
  geom_abline(intercept = c(0)) +
  geom_point(aes(x=average, y=bg, col = source, shape = colony),
             alpha = 0.5, size = 3) +
  geom_point(aes(x=average, y=bg, col = outlier, shape = colony,alpha = outlier),
             size = 3) +
  scale_colour_manual(name="Outlier Type",
                      values=c("TL"="#D32F2F","TR"="#536DFE","BL"="#388E3C","BR"="#795548",
                               "Bigger"="#FFC107","Smaller"="#7B1FA2","Normal"="transparent"),
                      breaks=c("Bigger","Smaller"),
                      labels=c("Underpredicted","Overpredicted")) +
  scale_shape_manual(name="Colony Type",
                     values=c("Reference"=18,"Query"=15,"Gap"=1),
                     breaks=c("Reference","Query")) +
  scale_alpha_manual(values=c("Smaller"=0.9,"Bigger"=0.9,"Normal"=0),
                     guide = F) +
  labs(title = "",
       subtitle = sprintf("Bottom Right | RMSE = %0.3f",
                          sqrt(sum(abs(sca.dat$bg[sca.dat$source == 'BR'] - sca.dat$average[sca.dat$source == 'BR'])^2,na.rm = T)/
                                 length(sca.dat$bg[!is.na(sca.dat$average) & sca.dat$source == 'BR']))),
       x = "Observed Colony Size (Pixel Count)",
       y = "") +
  scale_x_continuous(breaks = seq(0,1000,100),
                     minor_breaks = seq(0,1000,25)) +
  scale_y_continuous(breaks = seq(0,1000,100),
                     minor_breaks = seq(0,1000,25)) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=20),
        legend.position = c(0.8,0.2),
        legend.background = element_blank(),
        legend.text = element_text(size=15),
        legend.title =  element_text(size=15),
        plot.title = element_text(size=25,hjust = -1),
        plot.subtitle = element_text(size=20,hjust = 0.5)) +
  guides(color = guide_legend(override.aes = list(size=3, alpha = 1)),
         shape = guide_legend(override.aes = list(size=3))) +
  coord_cartesian(xlim = c(0,600),
                  ylim = c(0,600))

fig.s3 <- ggarrange(sca.tl, sca.tr,
                    sca.bl, sca.br,
                    nrow = 2)
ggsave(sprintf("%sfigureS3.jpg",out_path),
       fig.s3,
       width = 20,height = 20)

##### FIGURE 5
load("/home/sbp29/R/Projects/proto_plots/rawdata/4C3_GA1_FDR.RData")
data = dbGetQuery(conn, sprintf('select * from %s',tablename_pval))
data = data[data$orf_name != 'BFC100',]
pdata = data.frame()
i = 1

for (hr in unique(data$hours)) {
  temp = data[data$hours == hr,]
  n = dim(temp)[1]
  for (p in pvals) {
    pdata = rbind(pdata, c(hr, p, sum(temp$p <= p)/n))
  }
  if (i == 1) {
    colnames(pdata) = c('hours','p','fpr')
    g = ggplot() + 
      geom_line(data = pdata[pdata$hours == hr,], aes(x = p, y = fpr, col = as.character(hours)), lwd = 1.5)
  } else {
    g = g + geom_line(data = pdata[pdata$hours == hr,], aes(x = p, y = fpr, col = as.character(hours)), lwd = 1.5)
  }
  i = i + 1
}

##### A
fpr <- ggplot(stats.tmp[stats.tmp$hours == stats.tmp$cont_hrs,]) +
  geom_abline(col = 'red', linetype = 'dashed', lwd = 1) +
  geom_line(aes(x = p, y = fpr, col = hours), lwd = 1.2) +
  labs(title = "False positive rate",
       subtitle = "",
       x = "p-value cut-off",
       y = "False Positive Rate") +
  scale_x_continuous(breaks = seq(0,1,0.1),
                     minor_breaks = seq(0,1,0.05),
                     limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0,1,0.1),
                     minor_breaks = seq(0,1,0.05),
                     limits = c(0,1)) +
  scale_color_manual(name = "Hours",
                     breaks=c("13","14","16","17","18"),
                     values=c("13"="#D32F2F","14"="#536DFE","16"="#388E3C","17"="#795548","18"="#00BCD4",
                              "0"="transparent","8"="transparent","9"="transparent","10"="transparent","11"="transparent")) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=8),
        axis.title.x = element_text(size=10),
        axis.text.y = element_text(size=8),
        axis.title.y = element_text(size=10),
        # axis.title = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10),
        legend.position = "bottom",
        plot.title = element_text(size=12),
        plot.subtitle = element_text(size=10)) +
  guides(color = guide_legend(override.aes = list(size=6)),
         shape = guide_legend(override.aes = list(size=6))) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

##### B
fpr.zoom <- ggplot(stats.tmp[stats.tmp$hours == stats.tmp$cont_hrs,]) +
  geom_abline(col = 'red', linetype = 'dashed', lwd = 1) +
  geom_line(aes(x = p, y = fpr, col = hours), lwd = 1.2) +
  labs(title = "False positive rate",
       subtitle = "for p < 0.02",
       x = "p-value cut-off",
       y = "False Positive Rate") +
  scale_x_continuous(breaks = seq(0,1,0.05),
                     minor_breaks = seq(0,1,0.01),
                     limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0,1,0.05),
                     minor_breaks = seq(0,1,0.01),
                     limits = c(0,1)) +
  scale_color_manual(name = "Hours",
                     breaks=c("13","14","16","17","18"),
                     values=c("13"="#D32F2F","14"="#536DFE","16"="#388E3C","17"="#795548","18"="#00BCD4",
                              "0"="transparent","8"="transparent","9"="transparent","10"="transparent","11"="transparent")) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=8),
        axis.title.x = element_text(size=10),
        axis.text.y = element_text(size=8),
        axis.title.y = element_text(size=10),
        # axis.title = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10),
        legend.position = "bottom",
        plot.title = element_text(size=12),
        plot.subtitle = element_text(size=10)) +
  guides(color = guide_legend(override.aes = list(size=6)),
         shape = guide_legend(override.aes = list(size=6))) +
  coord_cartesian(xlim = c(0, 0.2), ylim = c(0, 0.2))

##### C
d <- read_excel("rawdata/cdata_ga1_clean.xlsx",col_types = "numeric")
d.smooth <- read_excel("rawdata/cdata_ga1_clean_smooth.xlsx",col_types = "numeric")

pow <- ggplot() +
  geom_line(data=d.smooth,aes(x=es, y=power, col='Smoothed Fit'),
            lwd = 1.2) +
  geom_point(data=d,aes(x=es, y=power, col='Data Point'),
             size = 4) +
  labs(title = "C. Change in power with effect size",
       x = "Effect Size",
       y = "Power") +
  scale_x_continuous(breaks = seq(0,10,0.05),
                     minor_breaks = seq(0,10,0.025)) +
  scale_y_continuous(breaks = seq(-20,140,10),
                     minor_breaks = seq(-20,140,5)) +
  scale_color_manual(name = "Result",
                     breaks=c("Data Point", "Smoothed Fit"),
                     values=c("Data Point"="#303F9F",
                              "Smoothed Fit"="#FF4081")) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=20),
        # axis.title = element_blank(),
        legend.text = element_text(size=15),
        legend.title = element_text(size=20),
        legend.position = "bottom",
        plot.title = element_text(size=25,hjust = -0.17),
        plot.subtitle = element_text(size=20,hjust = -0.1)) +
  coord_cartesian(xlim = c(0.8, 1.2), ylim = c(0, 100))

##### D
data <- read.csv(file="rawdata/4C3_GA1_TEMP_P_ES_17.csv", header=TRUE, sep=",")
hr = 17
alldat <- data.frame()

data <- data[data$orf_name != 'BFC100',]
data$es <- abs(data$es)

for (es in c(0.05,0.1,0.15,0.2,0.5)) {
  data.temp <- data[data$es <= es,]
  fpdat <- data.temp[data.temp$hours == 17,]
  tpdat <- data.temp[data.temp$hours != 17,]
  fp <- NULL
  tp <- NULL
  for (i in 1:length(unique(data$p))) {
    p <- unique(data$p)[i]
    fp[i] <- dim(fpdat[fpdat$p <= p,])[1]/dim(fpdat)[1] * 100
    tp[i] <- dim(tpdat[tpdat$p <= p,])[1]/dim(tpdat)[1] * 100
  }
  plotdat <- data.frame(unique(data$p),fp,tp)
  colnames(plotdat) <- c('p','FalsePositive','TruePositive')
  plotdat$ES <- es*100
  
  alldat <- rbind(alldat,plotdat)
}

lw = 1.2

roc <- ggplot() +
  geom_line(data = alldat[alldat$ES==50,],
            aes(x = FalsePositive, y = p*100, col = 'p'),
            linetype = "dotdash", lwd = lw) +
  geom_line(data = alldat[alldat$ES==5,],
            aes(x = FalsePositive, y = TruePositive, col = '5%'),
            lwd = lw) +
  geom_line(data = alldat[alldat$ES==10,],
            aes(x = FalsePositive, y = TruePositive, col = '10%'),
            lwd = lw) +
  geom_line(data = alldat[alldat$ES==15,],
            aes(x = FalsePositive, y = TruePositive, col = '15%'),
            lwd = lw) +
  geom_line(data = alldat[alldat$ES==20,],
            aes(x = FalsePositive, y = TruePositive, col = '20%'),
            lwd = lw) +
  geom_line(data = alldat[alldat$ES==50,],
            aes(x = FalsePositive, y = TruePositive, col = '50%'),
            lwd = lw) +
  labs(title = "D. ROC curve",
       x = "False Positive Rate",
       y = "True Positive Rate") +
  scale_x_continuous(breaks = seq(0,100,10),
                     minor_breaks = seq(0,100,5),
                     limits = c(0,100)) +
  scale_y_continuous(breaks = seq(0,100,10),
                     minor_breaks = seq(0,100,5),
                     limits = c(0,100),
                     sec.axis = sec_axis(~./100,
                                         breaks = seq(0,1,0.1),
                                         name='p-value')) +
  scale_colour_manual(name="Effect Size Threshold",
                      breaks=c("5%","10%","15%","20%","50%"),
                      values=c("5%"="#D32F2F","10%"="#536DFE","15%"="#388E3C",
                               "20%"="#795548","50%"="#00BCD4","p"="#FFA000")) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=20),
        axis.line.y.right = element_line(color = "#FFA000"),
        axis.text.y.right = element_text(size=15, color = "#FFA000"),
        axis.title.y.right = element_text(size=20, color = "#FFA000"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=20),
        legend.position = "bottom",
        plot.title = element_text(size=25,hjust = -0.1))

##### FINAL FIG5
fig5 <- ggarrange(fpr, fpr.zoom, pow, roc,
                  nrow = 2)
ggsave(sprintf("%sfigure5.jpg",out_path),
       fig5,
       width = 20,height = 20)

#####

ggplot(sca.dat[sca.dat$source == 'TL' & sca.dat$hours == 10,]) +
  geom_point(aes(x=se, y=var, col = source, shape = colony),
             alpha = 0.5, size = 3) +
  geom_point(aes(x=se, y=var, col = outlier, shape = colony),
             alpha = 0.5, size = 3)
  


  
