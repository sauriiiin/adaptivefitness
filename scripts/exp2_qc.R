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
pos2coorPC$result[pos2coorPC$pos %in% qc12345$pos[qc12345$comments == 'switched to other orf']] = 'Mismatch'

##### VISUALIZE DATA

for (pl in 23:28) {
  ggplot(data = pos2coorPC[pos2coorPC$`384plate` == pl,], aes(x = `384col`, y = `384row`)) +
    geom_point(aes(x = `384col`, y = `384row`),shape = 20,size=1,na.rm = T) +
    geom_point(aes(x = `384col`, y = `384row`,col = result, shape = result),size=5,na.rm = T) +
    #geom_point(aes(size = average, col = prediction, shape = colsize),alpha=0.9) +
    #scale_size_continuous(name="Results") +
    labs(title = "PROTOGENE COLLECTION",
         subtitle = sprintf("Plate %d",pl),
         x = " Column",
         y = "Row") +
    scale_x_continuous(breaks = seq(1,24,1), minor_breaks = seq(-5,30,1)) +
    scale_y_continuous(breaks = seq(1,16,1), minor_breaks = seq(-5,22,1), trans = 'reverse') +
    scale_colour_manual(name="Result",
                        values=c("Flagged"="red","ORF Present"="grey90","Sequenced"="green","Mismatch" = "blue"),
                        breaks=c("ORF Present","Sequenced","Mismatch","Flagged")) +
    scale_shape_manual(name="Colony Kind",
                       values=c("ORF Present"=19,"Sequenced"=19,"Mismatch"=19,"Flagged"=8),
                       breaks=c("ORF Present","Sequenced","Mismatch","Flagged"),
                       guide=F) +
    theme_dark() +
    theme(axis.text.x = element_text(size=7),
          axis.text.y = element_text(size=7))
  ggsave(sprintf("figs/pc_plates/Plate%d.png",pl),
         width = 15,height = 10) 
}


##### GET/SET 96 DATA
bc <- dbGetQuery(conn,
                 'select ((b.96plate-1)*96 + ((b.96col-1)*8) + b.96row) as pos
                 from BACTERIAL_SPACE_MAP b
                 where strain_id > 0')
bc$pos <- as.integer(bc$pos)

pc <- dbGetQuery(conn,
                  'select ((b.96plate-1)*96 + ((b.96col-1)*8) + b.96row) as pos
                  from PROTOGENE_COLLECTION a, BACTERIAL_SPACE_MAP b
                  where a.strain_id = b.strain_id')
pc$pos <- as.integer(pc$pos)

pos2coor96 <- dbGetQuery(conn,
                          'select * from PT_pos2coor96
                          order by 96col, 96row')
pos2coorPC <- NULL
for (i in 1:29) {
  pos2coor96$`96plate` <- i
  pos2coorPC <- rbind(pos2coorPC,pos2coor96)
}
pos2coorPC$pos <- (pos2coorPC$`96plate`-1)*96 + ((pos2coorPC$`96col`-1)*8) + pos2coorPC$`96row`

qc.res <- dbGetQuery(conn,
                      'select ((b.96plate-1)*96 + ((b.96col-1)*8) + b.96row) as pos,
                      a.score, a.qc, a.flag, a.comments
                      from QC_SEQ_ANALYSIS a, BACTERIAL_SPACE_MAP b
                      where a.strain_id = b.strain_id
                      order by b.96plate, b.96col, b.96row')
qc.res$pos <- as.integer(qc.res$pos)

pos2coorPC$seq[pos2coorPC$pos %in% qc.res$pos] = 1
pos2coorPC$flag[pos2coorPC$pos %in% qc.res$pos[qc.res$flag == 2]] = 2
pos2coorPC$mismatch[pos2coorPC$pos %in% qc.res$pos[qc.res$comments == 'switched to other orf']] = 1
pos2coorPC$result[pos2coorPC$pos %in% bc$pos] = 'ORF Present'
pos2coorPC$result[pos2coorPC$pos %in% pc$pos] = 'EXP'
pos2coorPC$result[pos2coorPC$pos %in% qc.res$pos] = 'Sequenced'
pos2coorPC$result[pos2coorPC$pos %in% qc.res$pos[qc.res$flag == 2]] = 'Flagged'
pos2coorPC$result[pos2coorPC$pos %in% qc.res$pos[qc.res$comments == 'switched to other orf']] = 'Mismatch'

##### VISUALIZE DATA

for (pl in 1:29) {
  ggplot(data = pos2coorPC[pos2coorPC$`96plate` == pl,], aes(x = `96col`, y = `96row`)) +
    geom_point(aes(x = `96col`, y = `96row`),shape = 20,size=1,na.rm = T) +
    geom_point(aes(x = `96col`, y = `96row`,col = result, shape = result),size=5,na.rm = T) +
    #geom_point(aes(size = average, col = prediction, shape = colsize),alpha=0.9) +
    #scale_size_continuous(name="Results") +
    labs(title = "Entry clone collection",
         subtitle = sprintf("Plate %d",pl),
         x = " Column",
         y = "Row") +
    scale_x_continuous(breaks = seq(1,12,1), minor_breaks = seq(-2,14,1)) +
    scale_y_continuous(breaks = seq(1,8,1), minor_breaks = seq(-2,10,1), trans = 'reverse') +
    scale_colour_manual(name="Result",
                        values=c("Flagged"="#D32F2F","ORF Present"="#9E9E9E",
                                 "EXP"="#9C27B0","Sequenced"="#388E3C","Mismatch" = "#536DFE"),
                        breaks=c("ORF Present","EXP","Sequenced","Mismatch","Flagged")) +
    scale_shape_manual(name="Colony Kind",
                       values=c("ORF Present"=19,"EXP"=19,"Sequenced"=19,"Mismatch"=19,"Flagged"=8),
                       breaks=c("ORF Present","EXP","Sequenced","Mismatch","Flagged"),
                       guide=F) +
    theme_linedraw() +
    theme(axis.text.x = element_text(size=7),
          axis.text.y = element_text(size=7))
  ggsave(sprintf("figs/pc_plates/ec_plates/ECPlate%d.png",pl),
         width = 15,height = 10) 
}












