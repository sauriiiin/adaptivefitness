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

##### LOCATE THE BIG SMALL COLONIES

for (hr in 9:length(hours$hours)) {
  for (pl in n_plates$`6144plate`) {
    fitdat = dbGetQuery(conn, sprintf('select * from %s a, %s b where a.hours = %d and a.pos = b.pos and b.%s = %d order by b.%s, b.%s',
                                      tablename_fit,p2c_info[1],hours$hours[hr],p2c_info[2],
                                      pl,p2c_info[3],p2c_info[4]))
    fitdat$fitness[fitdat$orf_name != 'BF_control'] = NA
    lims = quantile(fitdat$fitness, c(.025, .975),na.rm = TRUE) 
    
    fitdat$source[fitdat$`6144row`%%2==1 & fitdat$`6144col`%%2==1] = 'TL'
    fitdat$source[fitdat$`6144row`%%2==0 & fitdat$`6144col`%%2==1] = 'BL'
    fitdat$source[fitdat$`6144row`%%2==1 & fitdat$`6144col`%%2==0] = 'TR'
    fitdat$source[fitdat$`6144row`%%2==0 & fitdat$`6144col`%%2==0] = 'BR'
    
    fitdat$colsize[fitdat$fitness < lims[[1]]] = "Very Small"
    fitdat$colsize[fitdat$fitness > lims[[2]]] = "Very Big"
    fitdat$colsize[is.na(fitdat$colsize)] = "Normal"
    
    fitdat$se = fitdat$average - fitdat$bg
    fitdat$se[fitdat$orf_name != 'BF_control'] = NA
    fitdat$prediction[fitdat$se < 0] = 'Overestimate'
    fitdat$prediction[fitdat$se > 0] = 'Underestimate'
    fitdat$prediction[fitdat$se == 0] = 'Exact'
    
    fitdat$Gap[is.na(fitdat$orf_name)] = 'Gap'
    
    # ggplot(data = fitdat, aes(x = `6144col`, y = `6144row`)) +
    #   geom_point(aes(x = `6144col`, y = `6144row`),shape = 20,size=0.5) +
    #   geom_point(aes(x = `6144col`, y = `6144row`, shape = Gap),na.rm = T) +
    #   geom_point(aes(size = average, col = prediction, shape = colsize),alpha=0.9) +
    #   scale_size_continuous(name="Pixel Count") +
    #   labs(title = sprintf("%s Reference Colonies",expt),
    #        subtitle = sprintf("Plate %d @ %d hrs",pl,hours$hours[hr]),
    #        x = " Column",
    #        y = "Row") +
    #   scale_x_continuous(breaks = seq(1,96,1), minor_breaks = seq(-5,100,1)) +
    #   scale_y_continuous(breaks = seq(1,64,1), minor_breaks = seq(-5,70,1), trans = 'reverse') +
    #   scale_colour_manual(name="Prediction",
    #                      values=c("Overestimate"="#F44336","Underestimate"="#536DFE","Exact"="#4CAF50"),
    #                      breaks=c("Underestimate","Exact","Overestimate")) +
    #   scale_shape_manual(name="Colony Kind",
    #                      values=c(8,1,19,19)) +
    #   theme_light() +
    #   theme(axis.text.x = element_text(size=7),
    #         axis.text.y = element_text(size=7))
    # ggsave(sprintf("%scolsize/posse%s_%d%d.png",
    #                out_path,expt_name,hours$hours[hr],pl),
    #        width = 21,height = 14)
    
    for (sr in unique(fitdat$source)) {
      ggplot(data = fitdat[fitdat$source==sr,], aes(x = `6144col`, y = `6144row`)) +
        geom_point(aes(x = `6144col`, y = `6144row`),shape = 20,size=0.5) +
        geom_point(aes(x = `6144col`, y = `6144row`, shape = Gap),na.rm = T) +
        geom_point(aes(size = average, col = prediction, shape = colsize),alpha=0.9) +
        scale_size_continuous(name="Pixel Count") +
        labs(title = sprintf("%s Reference Colonies",expt),
             subtitle = sprintf("Plate %d %s @ %d hrs",pl,sr,hours$hours[hr]),
             x = " Column",
             y = "Row") +
        scale_x_continuous(breaks = seq(1,96,1), minor_breaks = seq(-5,100,1)) +
        scale_y_continuous(breaks = seq(1,64,1), minor_breaks = seq(-5,70,1), trans = 'reverse') +
        scale_colour_manual(name="Prediction",
                            values=c("Overestimate"="#F44336","Underestimate"="#536DFE","Exact"="#4CAF50"),
                            breaks=c("Underestimate","Exact","Overestimate")) +
        scale_shape_manual(name="Colony Kind",
                           values=c(8,1,19,19)) +
        theme_light() +
        theme(axis.text.x = element_text(size=7),
              axis.text.y = element_text(size=7))
      ggsave(sprintf("%scolsize/posse%s_%d%d%s.png",
                     out_path,expt_name,hours$hours[hr],pl,sr),
             width = 21,height = 14)
    }
  }
}

