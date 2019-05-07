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

fitdat$source[fitdat$`6144row`%%2==1 & fitdat$`6144col`%%2==1] = 'TL'
fitdat$source[fitdat$`6144row`%%2==0 & fitdat$`6144col`%%2==1] = 'BL'
fitdat$source[fitdat$`6144row`%%2==1 & fitdat$`6144col`%%2==0] = 'TR'
fitdat$source[fitdat$`6144row`%%2==0 & fitdat$`6144col`%%2==0] = 'BR'

fitdat$colony[fitdat$orf_name == 'BF_control'] = 'Reference'
fitdat$colony[fitdat$orf_name != 'BF_control'] = 'Query'
fitdat$colony[is.na(fitdat$orf_name)] = 'Gap'

ggplot(data = fitdat, aes(x = `6144col`, y = `6144row`)) +
  geom_point(aes(x = `6144col`, y = `6144row`),shape = 20,size=0.5) +
  geom_point(aes(x = `6144col`, y = `6144row`,
                 col = source,
                 size = fitness,
                 shape = colony),alpha=0.9,na.rm = T) +
  scale_size_continuous(guide=F) +
  labs(x = "Column",
       y = "Row") +
  scale_x_continuous(breaks = seq(1,96,1), minor_breaks = seq(-5,100,1)) +
  scale_y_continuous(breaks = seq(1,64,1), minor_breaks = seq(-5,70,1), trans = 'reverse') +
  scale_colour_manual(name="Source",
                     values=c("TL"="#D32F2F","TR"="#536DFE","BL"="#388E3C","BR"="#00BCD4"),
                     breaks=c("TL","TR","BL","BR"),
                     labels=c("Top Left","Top Right","Bottom Left","Bottom Right")) +
  scale_shape_manual(name="Colony Kind",
                     values=c(15,17,19),
                     breaks=c("Reference","Query","Gap")) +
  theme_light() +
  theme(axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7))
ggsave(sprintf("%s%s_6144_%d%d.png",
               out_path,expt_name,hr,pl),
       width = 21,height = 14)








