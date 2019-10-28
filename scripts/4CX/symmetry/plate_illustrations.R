##### PLATE ILLUSTRATIONS
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 05/07/2019

##### INITIALIZE
#install.packages("RMariaDB")
#install.packages("ggplot2")
library(RMariaDB)
library(ggplot2)
source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

##### GET/SET DATA
expt_name = '4C3_GA1'
expt = 'FS1-1'
out_path = 'figs/illustrations/';
density = 6144;

tablename_fit = sprintf('%s_%d_FITNESS',expt_name,density);
tablename_p2o = '4C3_pos2orf_name1';
tablename_bpos = '4C3_borderpos';

p2c_info = NULL
p2c_info[1] = '4C3_pos2coor6144'
p2c_info[2] = '6144plate'
p2c_info[3] = '6144col'
p2c_info[4] = '6144row'

hours = dbGetQuery(conn, sprintf('select distinct hours from %s order by hours asc', tablename_fit))

p2c = dbGetQuery(conn, sprintf('select * from %s a order by a.%s, a.%s, a.%s',
                               p2c_info[1],
                               p2c_info[2],
                               p2c_info[3],
                               p2c_info[4]))
p2c384 = dbGetQuery(conn, 'select 384plate, 384row, 384col
                     from 3P2_pos2coor384
                     where 384plate = 1
                     order by 384col asc, 384row asc')

n_plates = dbGetQuery(conn, sprintf('select distinct %s from %s a order by %s asc',
                                    p2c_info[2],
                                    p2c_info[1],
                                    p2c_info[2]))

hr = hours[[1]][9]
pl = n_plates[[1]][1]

#####

fitdat = dbGetQuery(conn, sprintf('select * from %s a, %s b where a.hours = %d and a.pos = b.pos and b.%s = %d order by b.%s, b.%s',
                                  tablename_fit,p2c_info[1],hr,p2c_info[2],
                                  pl,p2c_info[3],p2c_info[4]))
fitdat$bg[is.na(fitdat$average)] = NA

fitdat$source[fitdat$`6144row`%%2==1 & fitdat$`6144col`%%2==1] = 'TL'
fitdat$source[fitdat$`6144row`%%2==0 & fitdat$`6144col`%%2==1] = 'BL'
fitdat$source[fitdat$`6144row`%%2==1 & fitdat$`6144col`%%2==0] = 'TR'
fitdat$source[fitdat$`6144row`%%2==0 & fitdat$`6144col`%%2==0] = 'BR'

# fitdat <- fitdat[fitdat$`6144row` %in% 23:(23+15) & fitdat$`6144col` %in% 41:(41+23),]
# 
# fitdat <- cbind(fitdat,p2c384)
# 
# fitdat$source[fitdat$`384row`%%2==1 & fitdat$`384col`%%2==1] = 'TL'
# fitdat$source[fitdat$`384row`%%2==0 & fitdat$`384col`%%2==1] = 'BL'
# fitdat$source[fitdat$`384row`%%2==1 & fitdat$`384col`%%2==0] = 'TR'
# fitdat$source[fitdat$`384row`%%2==0 & fitdat$`384col`%%2==0] = 'BR'

fitdat$colony[fitdat$orf_name == 'BF_control'] = 'Reference'
fitdat$colony[fitdat$orf_name != 'BF_control'] = 'Query'
fitdat$colony[is.na(fitdat$orf_name)] = 'Gap'

##### FIGURE 2. BIOINFORMATIC WORK FLOW
ggplot(data = fitdat, aes(x = `6144col`, y = `6144row`)) +
  # geom_point(aes(x = `6144col`, y = `6144row`),shape = 20,size=0.001) +
  geom_point(aes(x = `6144col`, y = `6144row`,
                 col = source,
                 size = bg,
                 shape = colony,
                 alpha = colony),na.rm = T) +
  scale_size_continuous(guide=F) +
  # labs(title = "",
  #      x = "",
  #      y = "") +
  scale_x_continuous(breaks = seq(1,96,1),limits = c(1,96)) +
  scale_y_continuous(breaks = seq(1,64,1),limits = c(64,1),trans = 'reverse') +
  scale_colour_manual(name="Source",
                     values=c("TL"="#D32F2F","TR"="#536DFE","BL"="#388E3C","BR"="#795548"),
                     breaks=c("TL","TR","BL","BR"),
                     labels=c("Top Left","Top Right","Bottom Left","Bottom Right"),guide=F) +
  scale_shape_manual(name="Colony Kind",
                     values=c("Gap"=1,"Query"=15,"Reference"=18),
                     # values=c("Gap"=1,"Query"=18,"Reference"=18),
                     breaks=c("Reference","Query","Gap"),guide=F) +
  scale_size_continuous(range = c(2, 3),guide=F) +
  scale_alpha_manual(values=c("Gap"=1,"Query"=1,"Reference"=1),
                     guide=F) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.text = element_text(size=13),
        legend.title = element_text(size=15,face="bold"),
        legend.position = "right",
        plot.title = element_text(size=20,hjust = 0.5),
        plot.subtitle = element_text(size=13,hjust = 0.5))
ggsave(sprintf("%s%s_step1.png",
               out_path,expt_name,hr,pl),
       width = 3,height = 2)


##### FIGURE 3 OVERVIEW OF RESULTS
min = min(fitdat$average, na.rm=T)
max = max(fitdat$average, na.rm=T)

fitdat$se = (fitdat$average - fitdat$bg)/fitdat$average*100

ggplot(data = fitdat, aes(x = `6144col`, y = `6144row`)) +
  geom_point(aes(x = `6144col`, y = `6144row`,col = average,
                 shape = colony,
                 alpha = colony),size = 7,na.rm = T) +
  # geom_point(aes(x = `384col`, y = `384row`,
  #                col = average,
  #                shape = colony,
  #                alpha = colony),size = 7,na.rm = T) +
  labs(title = "A. Observed Colony Size (Pixel Count)",
       x = "",
       y = "") +
  scale_x_continuous(breaks = seq(1,96,1),limits = c(5,92)) +
  scale_y_continuous(breaks = seq(1,64,1),limits = c(60,5),trans = 'reverse') +
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
        plot.subtitle = element_text(size=20,hjust = 0))
ggsave(sprintf("%s%s_FIG2A.png",
               out_path,expt_name,hr,pl),
       width = 10,height = 6)

ggplot(data = fitdat, aes(x = `6144col`, y = `6144row`)) +
  geom_point(aes(x = `6144col`, y = `6144row`,col = bg,
                 shape = colony,
                 alpha = colony),size = 7,na.rm = T) +
  # geom_point(aes(x = `384col`, y = `384row`,
  #                col = average,
  #                shape = colony,
  #                alpha = colony),size = 7,na.rm = T) +
  labs(title = "B. Predicted Colony Size (Pixel Count)",
       x = "",
       y = "") +
  scale_x_continuous(breaks = seq(1,96,1),limits = c(5,92)) +
  scale_y_continuous(breaks = seq(1,64,1),limits = c(60,5),trans = 'reverse') +
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
        plot.subtitle = element_text(size=20,hjust = 0))
ggsave(sprintf("%s%s_FIG2B.png",
               out_path,expt_name,hr,pl),
       width = 10,height = 6)

ggplot(data = fitdat, aes(x = `6144col`, y = `6144row`)) +
  geom_point(aes(x = `6144col`, y = `6144row`,col = fitness,
                 shape = colony,
                 alpha = colony),size = 7,na.rm = T) +
  # geom_point(aes(x = `384col`, y = `384row`,
  #                col = average,
  #                shape = colony,
  #                alpha = colony),size = 7,na.rm = T) +
  labs(title = "C. Fitness Landscape",
       x = "",
       y = "") +
  scale_x_continuous(breaks = seq(1,96,1),limits = c(5,92)) +
  scale_y_continuous(breaks = seq(1,64,1),limits = c(60,5),trans = 'reverse') +
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
        plot.subtitle = element_text(size=20,hjust = 0))
ggsave(sprintf("%s%s_FIG2C.png",
               out_path,expt_name,hr,pl),
       width = 10,height = 6)

