##### PROTOGENE COLLECTION QC
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 05/01/2019

##### INITIALIZE
#install.packages("RMariaDB")
#install.packages("ggplot2")
library(RMariaDB)
library(ggplot2)
source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

##### GET/SET DATA
pc <- dbGetQuery(conn,
                 'select *,
                  ((b.384plate-1)*384 + ((b.384col-1)*16) + b.384row) as pos
                  from PROTOGENE_COLLECTION b')
pc$pos <- as.integer(pc$pos)

pos2coor384 <- dbGetQuery(conn,
                         'select * from pos2coor384
                         where 384plate in (1,2,3,4,5,6)
                         order by 384plate, 384col, 384row')
pos2coorPC <- pos2coor384
pos2coorPC$`384plate` <- pos2coorPC$`384plate` + 22
pos2coorPC$pos <- (pos2coorPC$`384plate`-1)*384 + ((pos2coorPC$`384col`-1)*16) + pos2coorPC$`384row`

qc.res <- dbGetQuery(conn,
                     'select * from QC_SEQ_ANALYSIS
                     where qc in (6,7)')
qc12345 <- dbGetQuery(conn,
                      'select b.*,
                      ((b.384plate-1)*384 + ((b.384col-1)*16) + b.384row) as pos, a.score,
                      a.qc, a.flag, a.comments
                      from QC_SEQ_ANALYSIS a, PROTOGENE_COLLECTION b
                      where a.qc < 6 and a.strain_id = b.strain_id')
qc12345$pos <- as.integer(qc12345$pos)

pos2coorPC$seq[pos2coorPC$pos %in% qc.res$pos | pos2coorPC$pos %in% qc12345$pos] = 1

pos2coorPC$flag[pos2coorPC$pos %in% qc.res$pos[qc.res$flag == 2] | pos2coorPC$pos %in% qc12345$pos[qc12345$flag == 2]] = 2

pos2coorPC$mismatch[pos2coorPC$pos %in% qc.res$pos[qc.res$comments == 'switched to other orf']] = 1

pos2coorPC$qc[pos2coorPC$pos %in% qc.res$pos] = 'QC67'
pos2coorPC$qc[pos2coorPC$pos %in% qc12345$pos] = 'QC12345'

pos2coorPC$result[pos2coorPC$pos %in% pc$pos] = 'ORF Present'
pos2coorPC$result[pos2coorPC$pos %in% qc.res$pos | pos2coorPC$pos %in% qc12345$pos] = 'Sequenced'
pos2coorPC$result[pos2coorPC$pos %in% qc.res$pos[qc.res$flag == 2] | pos2coorPC$pos %in% qc12345$pos[qc12345$flag == 2]] = 'Flagged'
pos2coorPC$result[pos2coorPC$pos %in% qc.res$pos[qc.res$comments == 'switched to other orf']] = 'Mismatch'

##### VISUALIZE DATA

pl <- 28

ggplot(data = pos2coorPC[pos2coorPC$`384plate` == pl,], aes(x = `384col`, y = `384row`)) +
  geom_point(aes(x = `384col`, y = `384row`),shape = 20,size=1,na.rm = T) +
  geom_point(aes(x = `384col`, y = `384row`,col = result, shape = result),size=3,na.rm = T) +
  #geom_point(aes(size = average, col = prediction, shape = colsize),alpha=0.9) +
  #scale_size_continuous(name="Results") +
  labs(title = "PROTOGENE COLLECTION",
       subtitle = sprintf("Plate %d",pl),
       x = " Column",
       y = "Row") +
  scale_x_continuous(breaks = seq(1,24,1), minor_breaks = seq(-5,30,1)) +
  scale_y_continuous(breaks = seq(1,16,1), minor_breaks = seq(-5,22,1), trans = 'reverse') +
  scale_colour_manual(name="Result",
                      values=c("Flagged"="red","ORF Present"="black","Sequenced"="green","Mismatch" = "blue"),
                      breaks=c("ORF Present","Sequenced","Mismatch","Flagged")) +
  scale_shape_manual(name="Colony Kind",
                     values=c(8,19,19,19),
                     breaks=c("ORF Present","Sequenced","Mismatch","Flagged"),
                     guide=F) +
  theme_light() +
  theme(axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7))
ggsave(sprintf("figs/pc_plates/Plate%d.png",pl),
       width = 21,height = 14)















