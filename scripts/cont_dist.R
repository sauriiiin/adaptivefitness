##### CONTROL DISTRIBUTION (4C)
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 04/26/2019

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
out_path = 'figs/';
density = 6144;

tablename_fit = sprintf('%s_%d_FITNESS',expt_name,density);
tablename_p2o = '4C3_pos2orf_name1';
tablename_bpos = '4C3_borderpos';

p2c_info = NULL
p2c_info[1] = '4C3_pos2coor6144'
p2c_info[2] = '6144plate'
p2c_info[3] = '6144col'
p2c_info[4] = '6144row'

cont.name  = 'BF_control'
dbFetch(query, n=-1)

query = dbSendQuery(conn, sprintf('select distinct hours from %s order by hours asc', tablename_fit))
hours = dbFetch(query,n=-1)

query = dbSendQuery(conn, sprintf('select * from %s a order by a.%s, a.%s, a.%s',
                                  p2c_info[1],
                                  p2c_info[2],
                                  p2c_info[3],
                                  p2c_info[4]))
p2c = dbFetch(query,n=-1)

query = dbSendQuery(conn, sprintf('select distinct %s from %s a order by %s asc',
                                  p2c_info[2],
                                  p2c_info[1],
                                  p2c_info[2]))
n_plates = dbFetch(query, n=-1) 

##### LOCATE THE BIG SMALL COLONIES

for (hr in 7:length(hours$hours)) {
  for (pl in n_plates$`6144plate`) {
    query = dbSendQuery(conn, sprintf('select * from %s a, %s b where a.hours = %d and a.pos = b.pos and b.%s = %d order by b.%s, b.%s',
                                      tablename_fit,p2c_info[1],hours$hours[hr],p2c_info[2],
                                      pl,p2c_info[3],p2c_info[4]))
    fitdat = dbFetch(query,n=-1)
    fitdat$fitness[fitdat$orf_name != 'BF_control'] = NA
    lims = quantile(fitdat$fitness, c(.025, .975),na.rm = TRUE) 
    
    fitdat$colsize[fitdat$fitness < lims[[1]]] = "Very Small"
    fitdat$colsize[fitdat$fitness > lims[[2]]] = "Very Big"
    fitdat$colsize[is.na(fitdat$colsize)] = "Normal"
    
    ggplot(data = fitdat, aes(x = `6144col`, y = `6144row`)) +
      geom_point(aes(size = fitness, col = colsize),alpha=0.7) +
      scale_size_continuous(guide=F) +
      scale_color_discrete(name="Colony Size") +
      theme_light() +
      labs(title = sprintf("%s",expt_name),
           subtitle = sprintf("Plate %d @ %d hrs",pl,hours$hours[hr]),
           x = " Column",
           y = "Row") +
      scale_x_continuous(breaks = seq(1,96,1), minor_breaks = seq(-5,100,1)) +
      scale_y_continuous(breaks = seq(1,64,1), minor_breaks = seq(-5,70,1)) +
      theme(axis.text.x = element_text(size=5),
            axis.text.y = element_text(size=5))
    
    ggsave(sprintf("%scolsize/colsize%s_%d%d.png",
                   out_path,expt_name,hours$hours[hr],pl),
           width = 12,height = 8)
    
    # query = dbSendQuery(conn, sprintf("select a.* from %s a, %s b where a.hours = %d and a.pos = b.pos and b.%s = %d and orf_name = '%s' order by b.%s, b.%s",
    #                                   tablename_fit,p2c_info[1],hours$hours[hr],p2c_info[2],
    #                                   pl,cont.name,p2c_info[3],p2c_info[4]))
    # contpos = dbFetch(query,n=-1)
  }
}

