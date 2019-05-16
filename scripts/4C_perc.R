
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 05/16/2019

##### INITIALIZE
library(RMariaDB)
library(ggplot2)
library(gridExtra)
source("R/functions/initialize.sql.R")

##### GET/SET DATA
expt_name = '4C3_GA1'
expt = 'FS1-1'
out_path = 'figs/perc/';
density = 6144;

data <- read.csv(file="rawdata/perc.csv", header=TRUE, sep=",")
data$pl[data$i < 5] = "One"; data$pl[data$i > 4] = "Two"
data$source[data$i == 1 | data$i == 5] = 'TL'
data$source[data$i == 2 | data$i == 6] = 'TR'
data$source[data$i == 3 | data$i == 7] = 'BL'
data$source[data$i == 4 | data$i == 8] = 'BR'

data$a_diff <- data$avg_ul - data$avg_ll
data$f_diff <- data$f_ul - data$f_ll
data$m_diff <- data$mck_ul - data$mck_ll
data$mf_diff <- data$m_ul - data$m_ll

##### PERC PLOT
ggplot(data, aes(x = data$avg_sd, y = data$a_diff, col = pl, fill = source),size = 4) +
  geom_point(shape = 21, lwd = 2) +
  labs(title = 'Pixel Count',
       x = 'Standard Deviation',
       y = 'Upper Range - Lower Range') +
  scale_colour_manual(name="Plate",
                      values=c("One"="#212121","Two"="#FFA000")) +
  scale_fill_manual(name="Source",
                      values=c("TL"="#D32F2F","TR"="#536DFE","BL"="#388E3C","BR"="#795548"),
                      breaks=c("TL","TR","BL","BR"),
                      labels=c("Top Left","Top Right","Bottom Left","Bottom Right")) +
  scale_x_continuous(breaks = seq(0,100,10),
                     minor_breaks = seq(0,100,5),
                     limits = c(10,90)) +
  scale_y_continuous(breaks = seq(-100,100,25),
                     minor_breaks = seq(-100,100,12.5)) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=10),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=10),
        axis.title.y = element_text(size=15),
        legend.text = element_text(size=10),
        legend.title = element_text(size=15),
        legend.position = 'bottom',
        plot.title = element_text(size=15,hjust = 0.5))
ggsave(sprintf("%s%s_PERC_PIX.png",
               out_path,expt_name),
       width = 10,height = 10.5)


ggplot(data, aes(x = source, y = mf_diff, col = pl, fill = source, size = mock_sd)) +
  geom_point(shape = 21) +
  labs(title = 'Mock Fitness',
       x = 'Source Plate',
       y = 'Upper Range - Lower Range') +
  scale_size(name = 'Pix Range Diff',range = c(2, 6)) +
  scale_colour_manual(name="Plate",
                      values=c("One"="#212121","Two"="#FFA000")) +
  scale_fill_manual(name="Source",
                    values=alpha(c("TL"="#D32F2F","TR"="#536DFE","BL"="#388E3C","BR"="#795548"),0.2),
                    breaks=c("TL","TR","BL","BR"),
                    labels=c("Top Left","Top Right","Bottom Left","Bottom Right")) +
  # scale_x_continuous(breaks = seq(0,100,10),
  #                    minor_breaks = seq(0,100,5),
  #                    limits = c(10,90)) +
  # scale_y_continuous(breaks = seq(-100,100,25),
  #                    minor_breaks = seq(-100,100,12.5)) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=10),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=10),
        axis.title.y = element_text(size=15),
        legend.text = element_text(size=10),
        legend.title = element_text(size=15),
        legend.position = 'right',
        plot.title = element_text(size=15,hjust = 0.5))
ggsave(sprintf("%s%s_PERC_sourceMF.png",
               out_path,expt_name),
       width = 11,height = 10)

##### CHECK POSITION WISE VARIABILITY
conn <- initialize.sql("saurin_test")

expt = 'FS1-1'
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

for (hr in hours[[1]]) {
  for (pl in n_plates[[1]]) {
    fitdat = dbGetQuery(conn, sprintf('select a.*, b.*
                                      from %s a, %s b
                                      where a.hours = %d
                                      and a.pos = b.pos
                                      and b.%s = %d order by b.%s, b.%s',
                                      tablename_fit,
                                      p2c_info[1],hr,p2c_info[2],
                                      pl,p2c_info[3],p2c_info[4]))
    
    fitdat$source[fitdat$`6144row`%%2==1 & fitdat$`6144col`%%2==1] = 'TL'
    fitdat$source[fitdat$`6144row`%%2==0 & fitdat$`6144col`%%2==1] = 'BL'
    fitdat$source[fitdat$`6144row`%%2==1 & fitdat$`6144col`%%2==0] = 'TR'
    fitdat$source[fitdat$`6144row`%%2==0 & fitdat$`6144col`%%2==0] = 'BR'
    
    var_avg <- NULL
    source <- NULL
    for (col in unique(fitdat$`6144col`)[7:90]) {
      for (row in unique(fitdat$`6144row`)[7:58]) {
        if (!is.na(fitdat$orf_name[fitdat$`6144col` == col & fitdat$`6144row` == row])) {
          a <- fitdat$average[fitdat$`6144col` == col & fitdat$`6144row` == row]
          u <- fitdat$average[fitdat$`6144col` == col & fitdat$`6144row` == row - 2]
          d <- fitdat$average[fitdat$`6144col` == col & fitdat$`6144row` == row + 2]
          l <- fitdat$average[fitdat$`6144col` == col - 2 & fitdat$`6144row` == row]
          r <- fitdat$average[fitdat$`6144col` == col + 2 & fitdat$`6144row` == row]
          fitdat$var[fitdat$`6144col` == col & fitdat$`6144row` == row] = sd(c(a,u,d,l,r),na.rm = T)/mean(c(a,u,d,l,r),na.rm = T)
          var_avg <- c(var_avg, sd(c(a,u,d,l,r),na.rm = T)/mean(c(a,u,d,l,r),na.rm = T))
          source <- c(source, fitdat$source[fitdat$`6144col` == col & fitdat$`6144row` == row])
        }
      }
    }
    var_avg <- data.frame(var_avg, source)
    colnames(var_avg) <- c("var", "source")
    ggplot(var_avg, aes(x=var, col=source)) +
      geom_density(lwd = 2) +
      labs(title = "Variability b/w Neighbors",
           subtitle = sprintf("%s | %d hours | Plate %d", expt, hr, pl),
           x = "Neighbor Variance",
           y = "Density") +
      scale_colour_manual(name="Source",
                          values=c("TL"="#D32F2F","TR"="#536DFE","BL"="#388E3C","BR"="#795548"),
                          breaks=c("TL","TR","BL","BR"),
                          labels=c("Top Left","Top Right","Bottom Left","Bottom Right")) +
      scale_x_continuous(breaks = seq(0,0.5,0.05),
                         minor_breaks = seq(0,0.5,0.025)) +
      scale_y_continuous(breaks = seq(0,50,10),
                         minor_breaks = seq(0,50,5)) +
      theme_linedraw() +
      coord_cartesian(xlim = c(0, 0.20))
    ggsave(sprintf("%s%s_VAR_NEIGH_AVG_%d_%d.png",
                   out_path,expt_name,
                   hr,pl),
           width = 11,height = 10)
  }
}



